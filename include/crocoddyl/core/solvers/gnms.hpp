///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2021, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_SOLVERS_GNMS_HPP_
#define CROCODDYL_CORE_SOLVERS_GNMS_HPP_

#include <Eigen/Cholesky>
#include <vector>

#include "crocoddyl/core/solvers/ddp.hpp"


namespace crocoddyl {

/**
 * @brief GNMS solver
 *
 * The GNMS solver computes an optimal trajectory and control commands by iterates running `backwardPass()` and
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
class SolverGNMS : public SolverDDP {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * @brief Initialize the GNMS solver
   *
   * @param[in] problem  shooting problem
   */
  explicit SolverGNMS(boost::shared_ptr<ShootingProblem> problem);
  virtual ~SolverGNMS();

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
  /**
   * @brief Computes the merit function, gaps at the given xs, us along with delta x and delta u
   */
  virtual void computeDirection(const bool recalcDiff);

  virtual double tryStep(const double stepLength);
  
  virtual void checkKKTConditions();

  const std::vector<Eigen::VectorXd>& get_xs_try() const { return xs_try_; };
  const std::vector<Eigen::VectorXd>& get_us_try() const { return us_try_; };
  
  const double get_gap_norm() const { return gap_norm_; };
  const double get_xgrad_norm() const { return x_grad_norm_; };
  const double get_ugrad_norm() const { return u_grad_norm_; };
  const double get_merit() const { return merit_; };
  const bool get_use_kkt_criteria() const { return use_kkt_criteria_; }

  void printCallbacks();
  void setCallbacks(bool inCallbacks);
  const bool getCallbacks();


  void set_mu(double mu) { mu_ = mu; };
  void set_termination_tolerance(double tol) { termination_tol_ = tol; };
  void set_use_kkt_criteria(bool inBool) { use_kkt_criteria_ = inBool; };
 public:
  using SolverDDP::xs_try_;
  using SolverDDP::us_try_;
  using SolverDDP::cost_try_;
  std::vector<Eigen::VectorXd> fs_try_;                               //!< Gaps/defects between shooting nodes
  std::vector<Eigen::VectorXd> dx_;                                    //!< the descent direction for x
  std::vector<Eigen::VectorXd> du_;                                    //!< the descent direction for u
  std::vector<Eigen::VectorXd> lag_mul_;                               //!< the Lagrange multiplier of the dynamics constraint

 protected:
  double merit_ = 0; // merit function at nominal traj
  double merit_try_ = 0; // merit function for the step length tried
  double x_grad_norm_ = 0; // 1 norm of the delta x
  double u_grad_norm_ = 0; // 1 norm of the delta u
  double gap_norm_ = 0; // 1 norm of the gaps
  double gap_norm_try_ = 0; // 1 norm of the gaps
  double cost_ = 0; // cost function
  double mu_ = 1e0; // penalty no constraint violation
  double termination_tol_ = 1e-3;
  bool with_callbacks_ = false;
  double KKT_ = std::numeric_limits<double>::infinity(); // KKT conditions residual
  bool use_kkt_criteria_ = true;

 private:
  double th_acceptnegstep_;  //!< Threshold used for accepting step along ascent direction
};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_SOLVERS_GNMS_HPP_
