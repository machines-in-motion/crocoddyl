///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2021, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifdef CROCODDYL_WITH_MULTITHREADING
#include <omp.h>
#endif  // CROCODDYL_WITH_MULTITHREADING

#include <iostream>
#include <iomanip>

#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/solvers/fadmm.hpp"

namespace crocoddyl {

// const std::vector<boost::shared_ptr<ConstraintModelAbstract> >& constraint_models
SolverFADMM::SolverFADMM(boost::shared_ptr<ShootingProblem> problem, 
                          const std::vector<boost::shared_ptr<ConstraintModelAbstract>>& constraint_models)
    : SolverDDP(problem){
      
      // std::cout << tmp << std::endl;
      const std::size_t T = this->problem_->get_T();
      const std::size_t ndx = problem_->get_ndx();
      // std::cout << "ndx" << ndx << std::endl;
      fs_try_.resize(T + 1);
      fs_flat_.resize(ndx*(T + 1));
      fs_flat_.setZero();
      dx_.resize(T+1);
      lag_mul_.resize(T+1);
      du_.resize(T);
      KKT_ = 0.;

      const std::vector<boost::shared_ptr<ActionModelAbstract> >& models = problem_->get_runningModels();
      for (std::size_t t = 0; t < T; ++t) {
        const boost::shared_ptr<ActionModelAbstract>& model = models[t];
        const std::size_t nu = model->get_nu();
        dx_[t].resize(ndx); du_[t].resize(nu);
        fs_try_[t].resize(ndx);
        lag_mul_[t].resize(ndx); 
        lag_mul_[t].setZero();
        dx_[t].setZero();
        du_[t] = Eigen::VectorXd::Zero(nu);
        fs_try_[t] = Eigen::VectorXd::Zero(ndx);
      }
      lag_mul_.back().resize(ndx);
      lag_mul_.back().setZero();
      dx_.back().resize(ndx);
      dx_.back().setZero();
      fs_try_.back().resize(ndx);
      fs_try_.back() = Eigen::VectorXd::Zero(ndx);
      

      const std::size_t n_alphas = 10;
      alphas_.resize(n_alphas);
      for (std::size_t n = 0; n < n_alphas; ++n) {
        alphas_[n] = 1. / pow(2., static_cast<double>(n));
      }
      if (th_stepinc_ < alphas_[n_alphas - 1]) {
        th_stepinc_ = alphas_[n_alphas - 1];
        std::cerr << "Warning: th_stepinc has higher value than lowest alpha value, set to "
                  << std::to_string(alphas_[n_alphas - 1]) << std::endl;
      }
    }

SolverFADMM::~SolverFADMM() {}

bool SolverFADMM::solve(const std::vector<Eigen::VectorXd>& init_xs, const std::vector<Eigen::VectorXd>& init_us,
                       const std::size_t maxiter, const bool is_feasible, const double reginit) {

  START_PROFILER("SolverFADMM::solve");
  if (problem_->is_updated()) {
    resizeData();
  }
  xs_try_[0] = problem_->get_x0();  // it is needed in case that init_xs[0] is infeasible
  setCandidate(init_xs, init_us, false);

  if (std::isnan(reginit)) {
    xreg_ = reg_min_;
    ureg_ = reg_min_;
  } else {
    xreg_ = reginit;
    ureg_ = reginit;
  }

  

  for (iter_ = 0; iter_ < maxiter; ++iter_) {


    was_feasible_ = false;
    bool recalcDiff = true;

    while (true) {
      try {
        computeDirection(recalcDiff);
      } 
      catch (std::exception& e) {
        recalcDiff = false;
        increaseRegularization();
        if (xreg_ == reg_max_) {
          return false;
        } else {
          continue;
        }
      }
      break;
    }

    // We need to recalculate the derivatives when the step length passes
    for (std::vector<double>::const_iterator it = alphas_.begin(); it != alphas_.end(); ++it) {
      steplength_ = *it;
      try {
        merit_try_ = tryStep(steplength_);
      } catch (std::exception& e) {
        continue;
      }
      // Heuristic line search criteria
      if(use_heuristic_line_search_){
        if (cost_ > cost_try_ || gap_norm_ > gap_norm_try_ ) {
          setCandidate(xs_try_, us_try_, false);
          recalcDiff = true;
          break;
        }
      }
      // Line-search criteria using merit function
      else{
        if (merit_ > merit_try_) {
          setCandidate(xs_try_, us_try_, false);
          recalcDiff = true;
          break;
        }
      }
    }

    if (steplength_ > th_stepdec_) {
      decreaseRegularization();
    }
    if (steplength_ <= th_stepinc_) {
      increaseRegularization();
      if (xreg_ == reg_max_) {
        STOP_PROFILER("SolverFADMM::solve");
        return false;
      }
    }
    stoppingCriteria();

    // const std::size_t n_callbacks = callbacks_.size();
    // for (std::size_t c = 0; c < n_callbacks; ++c) {
    //   CallbackAbstract& callback = *callbacks_[c];
    //   callback(*this);
    // }
    if(with_callbacks_){
      printCallbacks();
    }

    // KKT termination criteria
    if(use_kkt_criteria_){
      KKT_ = 0.;
      checkKKTConditions();
      if (KKT_  <= termination_tol_) {
        STOP_PROFILER("SolverFADMM::solve");
        return true;
      }
    }  
    // Old criteria
    else {
      if (x_grad_norm_  +  u_grad_norm_ < termination_tol_ ){
        STOP_PROFILER("SolverFADMM::solve");
        return true;
      }
    }
  }
  STOP_PROFILER("SolverFADMM::solve");
  return false;
}


void SolverFADMM::computeDirection(const bool recalcDiff){
  START_PROFILER("SolverFADMM::computeDirection");
  if (recalcDiff) {
    cost_ = calcDiff();
  }
  gap_norm_ = 0;
  const std::size_t T = problem_->get_T();
  for (std::size_t t = 0; t < T; ++t) {
    gap_norm_ += fs_[t].lpNorm<1>();   
  }
  gap_norm_ += fs_.back().lpNorm<1>();   

  merit_ = cost_ + mu_*gap_norm_;

  backwardPass();
  forwardPass();

  STOP_PROFILER("SolverFADMM::computeDirection");

}

void SolverFADMM::checkKKTConditions(){
  const std::size_t T = problem_->get_T();
  const std::size_t ndx = problem_->get_ndx();
  const std::vector<boost::shared_ptr<ActionDataAbstract> >& datas = problem_->get_runningDatas();
  for (std::size_t t = 0; t < T; ++t) {
    const boost::shared_ptr<ActionDataAbstract>& d = datas[t];
    KKT_ = std::max(KKT_, (d->Lx + d->Fx.transpose() * lag_mul_[t+1] - lag_mul_[t]).lpNorm<Eigen::Infinity>());
    KKT_ = std::max(KKT_, (d->Lu + d->Fu.transpose() * lag_mul_[t+1]).lpNorm<Eigen::Infinity>());
    fs_flat_.segment(t*ndx, ndx) = fs_[t];
  }
  fs_flat_.tail(ndx) = fs_.back();
  const boost::shared_ptr<ActionDataAbstract>& d_ter = problem_->get_terminalData();
  KKT_ = std::max(KKT_, (d_ter->Lx - lag_mul_.back()).lpNorm<Eigen::Infinity>());
  KKT_ = std::max(KKT_, fs_flat_.lpNorm<Eigen::Infinity>());
}


void SolverFADMM::forwardPass(){
    START_PROFILER("SolverFADMM::forwardPass");
    x_grad_norm_ = 0; u_grad_norm_ = 0;

    const std::size_t T = problem_->get_T();
    const std::vector<boost::shared_ptr<ActionDataAbstract> >& datas = problem_->get_runningDatas();
    for (std::size_t t = 0; t < T; ++t) {
      const boost::shared_ptr<ActionDataAbstract>& d = datas[t];
      lag_mul_[t].noalias() = Vxx_[t] * dx_[t] + Vx_[t];
      du_[t].noalias() = -K_[t]*(dx_[t]) - k_[t];
      dx_[t+1].noalias() = (d->Fx - (d->Fu * K_[t]))*(dx_[t]) - (d->Fu * (k_[t])) + fs_[t+1];
      x_grad_norm_ += dx_[t].lpNorm<1>(); // assuming that there is no gap in the initial state
      u_grad_norm_ += du_[t].lpNorm<1>();
    }
    lag_mul_.back() = Vxx_.back() * dx_.back() + Vx_.back();
    x_grad_norm_ += dx_.back().lpNorm<1>(); // assuming that there is no gap in the initial state
    x_grad_norm_ = x_grad_norm_/(T+1);
    u_grad_norm_ = u_grad_norm_/T; 
    STOP_PROFILER("SolverFADMM::forwardPass");

}


double SolverFADMM::tryStep(const double steplength) {
    if (steplength > 1. || steplength < 0.) {
        throw_pretty("Invalid argument: "
                    << "invalid step length, value is between 0. to 1.");
    }
    START_PROFILER("SolverFADMM::tryStep");
    cost_try_ = 0.;
    merit_try_ = 0;
    gap_norm_try_ = 0;

    const std::size_t T = problem_->get_T();
    const std::vector<boost::shared_ptr<ActionModelAbstract> >& models = problem_->get_runningModels();
    const std::vector<boost::shared_ptr<ActionDataAbstract> >& datas = problem_->get_runningDatas();

    for (std::size_t t = 0; t < T; ++t) {
        const boost::shared_ptr<ActionModelAbstract>& m = models[t];
        const boost::shared_ptr<ActionDataAbstract>& d = datas[t];
        const std::size_t nu = m->get_nu();

        // error = x + dx - f(x + dx, u + du)
        // std::cout << dx_.size() << std::endl;
        m->get_state()->integrate(xs_[t], steplength * dx_[t], xs_try_[t]); 
        if (nu != 0) {
            us_try_[t].noalias() = us_[t] + steplength * du_[t];
        }        
        m->calc(d, xs_try_[t], us_try_[t]);        
        cost_try_ += d->cost;

        if (t > 0){
          const boost::shared_ptr<ActionDataAbstract>& d_prev = datas[t-1];
          m->get_state()->diff(xs_try_[t], d_prev->xnext, fs_try_[t-1]);
          gap_norm_try_ += fs_try_[t-1].lpNorm<1>(); 
        } 

        if (raiseIfNaN(cost_try_)) {
          STOP_PROFILER("SolverFADMM::tryStep");
          throw_pretty("step_error");
        }   
    }

    // Terminal state update
    const boost::shared_ptr<ActionModelAbstract>& m_ter = problem_->get_terminalModel();
    const boost::shared_ptr<ActionDataAbstract>& d_ter = problem_->get_terminalData();
    m_ter->get_state()->integrate(xs_.back(), steplength * dx_.back(), xs_try_.back()); 
    m_ter->calc(d_ter, xs_try_.back());
    cost_try_ += d_ter->cost;

    const boost::shared_ptr<ActionModelAbstract>& m = models[T-1];
    const boost::shared_ptr<ActionDataAbstract>& d = datas[T-1];
    
    m->get_state()->diff(xs_try_.back(), d->xnext, fs_try_[T-1]);
    gap_norm_try_ += fs_try_[T-1].lpNorm<1>(); 

    merit_try_ = cost_try_ + mu_*gap_norm_try_;

    if (raiseIfNaN(cost_try_)) {
        STOP_PROFILER("SolverFADMM::tryStep");
        throw_pretty("step_error");
    }

    STOP_PROFILER("SolverFADMM::tryStep");

    return merit_try_;
}

void SolverFADMM::printCallbacks(){
  if (this->get_iter() % 10 == 0) {
    std::cout << "iter     merit         cost         grad      step    ||gaps||        KKT";
    std::cout << std::endl;
  }
  std::cout << std::setw(4) << this->get_iter() << "  ";
  std::cout << std::scientific << std::setprecision(5) << this->get_merit() << "  ";
  std::cout << std::scientific << std::setprecision(5) << this->get_cost() << "  ";
  std::cout << this->get_xgrad_norm() + this->get_ugrad_norm() << "  ";
  std::cout << std::fixed << std::setprecision(4) << this->get_steplength() << "  ";
  std::cout << std::scientific << std::setprecision(5) << this->get_gap_norm() << "  ";
  std::cout << std::scientific << std::setprecision(5) << KKT_;
  std::cout << std::endl;
  std::cout << std::flush;
}

void SolverFADMM::setCallbacks(bool inCallbacks){
  with_callbacks_ = inCallbacks;
}

const bool SolverFADMM::getCallbacks(){
  return with_callbacks_;
}
// double SolverFADMM::get_th_acceptnegstep() const { return th_acceptnegstep_; }

// void SolverFADMM::set_th_acceptnegstep(const double th_acceptnegstep) {
//   if (0. > th_acceptnegstep) {
//     throw_pretty("Invalid argument: "
//                  << "th_acceptnegstep value has to be positive.");
//   }
//   th_acceptnegstep_ = th_acceptnegstep;
// }

}  // namespace crocoddyl