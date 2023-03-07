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

#include "crocoddyl/core/utils/exception.hpp"
#include "crocoddyl/core/solvers/gnms.hpp"

namespace crocoddyl {

SolverGNMS::SolverGNMS(boost::shared_ptr<ShootingProblem> problem)
    : SolverDDP(problem){
      
      const std::size_t T = this->problem_->get_T();
      const std::size_t ndx = problem_->get_ndx();
      // std::cout << "ndx" << ndx << std::endl;
      fs_try_.resize(T + 1);
      dx_.resize(T+1);
      du_.resize(T);
      const std::vector<boost::shared_ptr<ActionModelAbstract> >& models = problem_->get_runningModels();
      for (std::size_t t = 0; t < T; ++t) {
        const boost::shared_ptr<ActionModelAbstract>& model = models[t];
        const std::size_t nu = model->get_nu();
        dx_[t].resize(ndx); du_[t].resize(nu);
        fs_try_[t].resize(ndx);
        dx_[t].setZero();
        du_[t] = Eigen::VectorXd::Zero(nu);
        fs_try_[t] = Eigen::VectorXd::Zero(ndx);
      }
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

SolverGNMS::~SolverGNMS() {}

bool SolverGNMS::solve(const std::vector<Eigen::VectorXd>& init_xs, const std::vector<Eigen::VectorXd>& init_us,
                       const std::size_t maxiter, const bool is_feasible, const double reginit) {

  START_PROFILER("SolverGNMS::solve");
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
        // std::cout << "Merit Try : " << merit_try_ << "   Merit   " << merit_ <<  "grad norm " <<  x_grad_norm_  +  u_grad_norm_ << std::endl;

      } catch (std::exception& e) {
        // std::cout << "blows up" << std::endl;
        continue;
      }
      // std::cout << "Merit Try : " << merit_try_ << "   Merit   " << merit_ <<  "grad norm " <<  x_grad_norm_  +  u_grad_norm_ << std::endl;

      if (merit_ > merit_try_) {
        setCandidate(xs_try_, us_try_, false);
        recalcDiff = true;
        break;
      }

    }
    // std::cout << "iter "<< iter_ << " Merit : " << merit_ << "   cost   " << cost_ <<  "  gap norm " <<  gap_norm_  << "  step length "<< steplength_ << std::endl;
    // std::cout << "iter "<< iter_ << " Merit_try : " << merit_try_ << "   cost_try   " << cost_try_ <<  "  gap norm try " <<  gap_norm_try_
    //                                                                                            << std::endl;

    if (steplength_ > th_stepdec_) {
      decreaseRegularization();
    }
    if (steplength_ <= th_stepinc_) {
      increaseRegularization();
      if (xreg_ == reg_max_) {
        STOP_PROFILER("SolverGNMS::solve");
        return false;
      }
    }
    stoppingCriteria();

    const std::size_t n_callbacks = callbacks_.size();
    for (std::size_t c = 0; c < n_callbacks; ++c) {
      CallbackAbstract& callback = *callbacks_[c];
      callback(*this);
    }

    if (x_grad_norm_  +  u_grad_norm_ < termination_tol_ ) {
      STOP_PROFILER("SolverGNMS::solve");
      return true;
    }
  }
  STOP_PROFILER("SolverGNMS::solve");
  return false;
}


void SolverGNMS::computeDirection(const bool recalcDiff){
  START_PROFILER("SolverGNMS::computeDirection");
  // if (recalcDiff) {
  cost_ = calcDiff();
  // }
  gap_norm_ = 0;
  const std::size_t T = problem_->get_T();

  for (std::size_t t = 0; t < T; ++t) {
    gap_norm_ += fs_[t].lpNorm<1>();   
  }
  gap_norm_ += fs_.back().lpNorm<1>();   
  
  merit_ = cost_ + mu_*gap_norm_;

  backwardPass();
  forwardPass();

  STOP_PROFILER("SolverGNMS::computeDirection");

}

void SolverGNMS::forwardPass(){
    START_PROFILER("SolverGNMS::forwardPass");
    STOP_PROFILER("SolverGNMS::forwardPass");
    x_grad_norm_ = 0; u_grad_norm_ = 0;

    const std::size_t T = problem_->get_T();
    const std::vector<boost::shared_ptr<ActionDataAbstract> >& datas = problem_->get_runningDatas();
    for (std::size_t t = 0; t < T; ++t) {
      const boost::shared_ptr<ActionDataAbstract>& d = datas[t];
      du_[t].noalias() = -K_[t]*(dx_[t]) - k_[t];
      dx_[t+1].noalias() = (d->Fx - (d->Fu * K_[t]))*(dx_[t]) - (d->Fu * (k_[t])) + fs_[t+1];
      
      x_grad_norm_ += dx_[t].lpNorm<1>(); // assuming that there is no gap in the initial state
      u_grad_norm_ += du_[t].lpNorm<1>();

    }

    x_grad_norm_ += dx_.back().lpNorm<1>(); // assuming that there is no gap in the initial state
    x_grad_norm_ = x_grad_norm_/T;
    u_grad_norm_ = u_grad_norm_/T; 
}


double SolverGNMS::tryStep(const double steplength) {
    if (steplength > 1. || steplength < 0.) {
        throw_pretty("Invalid argument: "
                    << "invalid step length, value is between 0. to 1.");
    }
    START_PROFILER("SolverGNMS::tryStep");
    cost_try_ = 0.;
    merit_try_ = 0;
    gap_norm_try_ = 0;

    const std::size_t T = problem_->get_T();
    const std::vector<boost::shared_ptr<ActionModelAbstract> >& models = problem_->get_runningModels();
    const std::vector<boost::shared_ptr<ActionDataAbstract> >& datas = problem_->get_runningDatas();
    for (std::size_t t = 0; t < T; ++t) {
        const boost::shared_ptr<ActionModelAbstract>& m = models[t];
        const std::size_t nu = m->get_nu();

        // error = x + dx - f(x + dx, u + du)
        // std::cout << dx_.size() << std::endl;
        m->get_state()->integrate(xs_[t], steplength * dx_[t], xs_try_[t]); 
        if (nu != 0) {
            us_try_[t].noalias() = us_[t] + steplength * du_[t];
        }        
    }

    // Terminal state update
    const boost::shared_ptr<ActionModelAbstract>& m = problem_->get_terminalModel();
    m->get_state()->integrate(xs_.back(), steplength * dx_.back(), xs_try_.back()); 

    for (std::size_t t = 0; t < T; ++t) {
        const boost::shared_ptr<ActionModelAbstract>& m = models[t];
        const boost::shared_ptr<ActionDataAbstract>& d = datas[t];
        
        m->calc(d, xs_try_[t], us_try_[t]);        
        // error = x + dx - f(x + dx, u + du)
        m->get_state()->diff(xs_try_[t+1], d->xnext, fs_try_[t]);

        cost_try_ += d->cost;
        gap_norm_try_ += fs_try_[t].lpNorm<1>(); 
        
        if (raiseIfNaN(cost_try_)) {
          STOP_PROFILER("SolverGNMS::tryStep");
          throw_pretty("step_error");
        }   
    }
    const boost::shared_ptr<ActionModelAbstract>& mter = problem_->get_terminalModel();
    const boost::shared_ptr<ActionDataAbstract>& d = problem_->get_terminalData();
    mter->calc(d, xs_try_.back());
    cost_try_ += d->cost;

    merit_try_ = cost_try_ + mu_*gap_norm_try_;

    // std::cout << "try step cost_try " << cost_try_ << " gap_norm_try_" << gap_norm_try_ << std::endl;


    if (raiseIfNaN(cost_try_)) {
        STOP_PROFILER("SolverGNMS::tryStep");
        throw_pretty("step_error");
    }

    STOP_PROFILER("SolverGNMS::tryStep");

    return merit_try_;
}

// double SolverGNMS::get_th_acceptnegstep() const { return th_acceptnegstep_; }

// void SolverGNMS::set_th_acceptnegstep(const double th_acceptnegstep) {
//   if (0. > th_acceptnegstep) {
//     throw_pretty("Invalid argument: "
//                  << "th_acceptnegstep value has to be positive.");
//   }
//   th_acceptnegstep_ = th_acceptnegstep;
// }

}  // namespace crocoddyl
