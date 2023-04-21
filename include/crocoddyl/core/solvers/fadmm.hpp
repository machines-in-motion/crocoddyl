///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2021, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_SOLVERS_FADMM_HPP_
#define CROCODDYL_CORE_SOLVERS_FADMM_HPP_

#include <Eigen/Cholesky>
#include <vector>

#include "crocoddyl/core/solvers/ddp.hpp"
#include "crocoddyl/core/constraint-base.hpp"


namespace crocoddyl {

/**
 * @brief fadmm solver
 *
 * The fadmm solver computes an optimal trajectory and control commands by iterates running `backwardPass()` and
 * `forwardPass()`. The backward pass accepts infeasible guess as described in the `SolverDDP::backwardPass()`.
 * Additionally, the forward pass handles infeasibility simulations that resembles the numerical behaviour of
 * a multiple-shooting formulation, i.e.:
 * \f{eqnarray}
 *   \mathbf{\hat{x}}_0 &=& \mathbf{\tilde{x}}_0 - (1 - \alpha)\mathbf{\bar{f}}_0,\\
 *   \mathbf{\hat{u}}_k &=& \mathbf{u}_k + \alpha\mathbf{k}_k + \mathbf{K}_k(\mathbf{\hat{x}}_k-\mathbf{x}_k),\\
 *   \mathbf{\hat{x}}_{k+1} &=& \mathbf{f}_k(\mathbf{\hat{x}}_k,\mathbf{\hat{u}}_k) - (1 -
 * \alpha)\mathbf{\bar{f}}_{k+1}.
 * \f}
 * Note that the forward pass keeps the gaps \f$\mathbf{\bar{f}}_s\f$ open according to the step length \f$\alpha\f$
 * that has been accepted. This solver has shown empirically greater globalization strategy. Additionally, the
 * expected improvement computation considers the gaps in the dynamics:
 * \f{equation}
 *   \Delta J(\alpha) = \Delta_1\alpha + \frac{1}{2}\Delta_2\alpha^2,
 * \f}
 * with
 * \f{eqnarray}
 *   \Delta_1 = \sum_{k=0}^{N-1} \mathbf{k}_k^\top\mathbf{Q}_{\mathbf{u}_k} +\mathbf{\bar{f}}_k^\top(V_{\mathbf{x}_k} -
 *   V_{\mathbf{xx}_k}\mathbf{x}_k),\nonumber\\ \Delta_2 = \sum_{k=0}^{N-1}
 *   \mathbf{k}_k^\top\mathbf{Q}_{\mathbf{uu}_k}\mathbf{k}_k + \mathbf{\bar{f}}_k^\top(2 V_{\mathbf{xx}_k}\mathbf{x}_k
 * - V_{\mathbf{xx}_k}\mathbf{\bar{f}}_k). \f}
 *
 * For more details about the feasibility-driven differential dynamic programming algorithm see:
 * \include mastalli-icra20.bib
 *
 * \sa `SolverDDP()`, `backwardPass()`, `forwardPass()`, `expectedImprovement()` and `updateExpectedImprovement()`
 */
class SolverFADMM : public SolverDDP {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * @brief Initialize the fadmm solver
   *
   * @param[in] problem  shooting problem
   */
  explicit SolverFADMM(boost::shared_ptr<ShootingProblem> problem, const std::vector<boost::shared_ptr<ConstraintModelAbstract>>& constraint_models);
  virtual ~SolverFADMM();

  virtual bool solve(const std::vector<Eigen::VectorXd>& init_xs = DEFAULT_VECTOR,
                     const std::vector<Eigen::VectorXd>& init_us = DEFAULT_VECTOR, const std::size_t maxiter = 100,
                     const bool is_feasible = false, const double regInit = 1e-9);

  /**
   * @copybrief SolverAbstract::expectedImprovement
   *
   * This function requires to first run `updateExpectedImprovement()`. The expected improvement computation considers
   * the gaps in the dynamics: \f{equation} \Delta J(\alpha) = \Delta_1\alpha + \frac{1}{2}\Delta_2\alpha^2, \f} with
   * \f{eqnarray}
   *   \Delta_1 = \sum_{k=0}^{N-1} \mathbf{k}_k^\top\mathbf{Q}_{\mathbf{u}_k} +\mathbf{\bar{f}}_k^\top(V_{\mathbf{x}_k}
   * - V_{\mathbf{xx}_k}\mathbf{x}_k),\nonumber\\ \Delta_2 = \sum_{k=0}^{N-1}
   *   \mathbf{k}_k^\top\mathbf{Q}_{\mathbf{uu}_k}\mathbf{k}_k + \mathbf{\bar{f}}_k^\top(2
   * V_{\mathbf{xx}_k}\mathbf{x}_k
   * - V_{\mathbf{xx}_k}\mathbf{\bar{f}}_k). \f}
   */
  // virtual const Eigen::Vector2d& expectedImprovement();

  /**
   * @brief Update internal values for computing the expected improvement
   */
  void updateExpectedImprovement();

  virtual void forwardPass();
  virtual void backwardPass();

  /**
   * @brief Computes the merit function, gaps at the given xs, us along with delta x and delta u
   */
  virtual void computeDirection(const bool recalcDiff);

  virtual double tryStep(const double stepLength);

  virtual void calculate(const bool recalc = true);

  // virtual void set_constraints(const std::vector<boost::shared_ptr<ConstraintModelAbstract>>& constraint_models){
  //   constraint_models_ = constraint_models;
  // };

  /**
   * @brief Compute the KKT conditions residual
   */
  virtual void checkKKTConditions();

  const std::vector<Eigen::VectorXd>& get_xs_try() const { return xs_try_; };
  const std::vector<Eigen::VectorXd>& get_us_try() const { return us_try_; };
  
  const std::vector<Eigen::VectorXd>& get_xs() const { return xs_; };
  const std::vector<Eigen::VectorXd>& get_us() const { return us_; };
  
  const std::vector<Eigen::VectorXd>& get_xs_tilde() const { return dxtilde_; };
  const std::vector<Eigen::VectorXd>& get_us_tilde() const { return dutilde_; };

  const std::vector<Eigen::VectorXd>& get_y() const { return y_; };
  const std::vector<Eigen::VectorXd>& get_z() const { return z_; };

  const double get_gap_norm() const { return gap_norm_; };
  const double get_xgrad_norm() const { return x_grad_norm_; };
  const double get_ugrad_norm() const { return u_grad_norm_; };
  const double get_merit() const { return merit_; };
  const bool get_use_kkt_criteria() const { return use_kkt_criteria_; };
  const bool get_use_heuristic_line_search() const { return use_heuristic_line_search_; };
  const double get_mu() const { return mu_; };
  const double get_termination_tolerance() const { return termination_tol_; };
  const int get_max_qp_iters(){ return max_qp_iters_; };


  void printCallbacks();
  void setCallbacks(bool inCallbacks);
  const bool getCallbacks();


  void set_mu(double mu) { mu_ = mu; };
  void set_termination_tolerance(double tol) { termination_tol_ = tol; };
  void set_use_kkt_criteria(bool inBool) { use_kkt_criteria_ = inBool; };
  void set_use_heuristic_line_search(bool inBool) { use_heuristic_line_search_ = inBool; };
  
  void update_lagrangian_parameters();
  void update_rho_sparse(int iter);

  void set_max_qp_iters(int iters){ max_qp_iters_ = iters; };

 public:
  using SolverDDP::xs_try_;
  using SolverDDP::us_try_;
  using SolverDDP::cost_try_;
  std::vector<Eigen::VectorXd> fs_try_;                                //!< Gaps/defects between shooting nodes
  std::vector<Eigen::VectorXd> dx_;                                    //!< the descent direction for x
  std::vector<Eigen::VectorXd> du_;                                    //!< the descent direction for u
  std::vector<Eigen::VectorXd> lag_mul_;                               //!< the Lagrange multiplier of the dynamics constraint
  Eigen::VectorXd fs_flat_;                                            //!< Gaps/defects between shooting nodes (1D array)
  double KKT_ = std::numeric_limits<double>::infinity();               //!< KKT conditions residual
  bool use_heuristic_line_search_ = false;                              //!< Use heuristic line search

  std::vector<Eigen::VectorXd> dxtilde_;                                    //!< the descent direction for x
  std::vector<Eigen::VectorXd> dutilde_;                                    //!< the descent direction for u

  // ADMM parameters
  std::vector<Eigen::VectorXd> y_;                                    //!< lagrangian dual variable
  std::vector<Eigen::VectorXd> z_;                                    //!< second admm variable
  std::vector<Eigen::VectorXd> z_prev_;                               //!< second admm variable previous
  std::vector<Eigen::VectorXd> z_relaxed_;                           //!< relaxed step of z
  std::vector<Eigen::VectorXd> rho_vec_;                              //!< rho vector

 protected:
  double merit_ = 0;                                           //!< merit function at nominal traj
  double merit_try_ = 0;                                       //!< merit function for the step length tried
  double x_grad_norm_ = 0;                                     //!< 1 norm of the delta x
  double u_grad_norm_ = 0;                                     //!< 1 norm of the delta u
  double gap_norm_ = 0;                                        //!< 1 norm of the gaps
  double constraint_norm_ = 0;                                 //!< 1 norm of constraint violation
  double gap_norm_try_ = 0;                                    //!< 1 norm of the gaps
  double cost_ = 0;                                            //!< cost function
  double mu_ = 1e0;                                            //!< penalty no constraint violation
  double termination_tol_ = 1e-8;                              //!< Termination tolerance
  bool with_callbacks_ = false;                                //!< With callbacks
  bool use_kkt_criteria_ = true;                               //!< Use KKT conditions as termination criteria 
  double sigma_ = 1e-6; // proximal term
  double alpha_ = 1.6; // relaxed step size
  int max_qp_iters_ = 10; // max qp iters

  double rho_estimate_sparse_ = 0.0; // rho estimate
  double rho_sparse_ = 1e-1; // rho
  double rho_min_ = 1e-6; // rho min
  double rho_max_ = 1e6; // rho max
  int rho_update_interval_ = 25; // frequency of update of rho
  double adaptive_rho_tolerance_ = 5; 
  double eps_abs_ = 1e-3; // absolute termination criteria
  double eps_rel_ = 1e-3; // relative termination criteria

  double norm_primal_ = 0.0; // norm primal residual
  double norm_dual_ = 0.0; // norm dual residual
  double norm_primal_rel_ = 0.0; // norm primal relative residual
  double norm_dual_rel_ = 0.0; // norm dual relative residual



 private:
  double th_acceptnegstep_;  //!< Threshold used for accepting step along ascent direction
  const std::vector<boost::shared_ptr<ConstraintModelAbstract>> cmodels_;
  std::vector<boost::shared_ptr<ConstraintDataAbstract>> cdatas_;
  Eigen::VectorXd dual_vecx;
  Eigen::VectorXd dual_vecu;

};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_SOLVERS_FADMM_HPP_
