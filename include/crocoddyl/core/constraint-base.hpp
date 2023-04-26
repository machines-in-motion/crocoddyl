///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2020-2022, University of Edinburgh, Heriot-Watt University
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_CONSTRAINT_BASE_HPP_
#define CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/state-base.hpp"
#include "crocoddyl/core/data-collector-base.hpp"
#include "crocoddyl/core/residual-base.hpp"

namespace crocoddyl {

/**
 * @brief Abstract class for constraint models
 *
 * A constraint model defines both: inequality \f$\mathbf{g}(\mathbf{x}, \mathbf{u})\in\mathbb{R}^{ng}\f$
 * and equality \f$\mathbf{h}(\mathbf{x}, \mathbf{u})\in\mathbb{R}^{nh}\f$ constraints. The constraint
 * function depends on the state point \f$\mathbf{x}\in\mathcal{X}\f$, which lies in the state manifold
 * described with a `nx`-tuple, its velocity \f$\dot{\mathbf{x}}\in T_{\mathbf{x}}\mathcal{X}\f$ that belongs to
 * the tangent space with `ndx` dimension, and the control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$.
 *
 * The main computations are carried out in `calc()` and `calcDiff()` routines. `calc()` computes the
 * constraint residual and `calcDiff()` computes the Jacobians of the constraint function. Concretely speaking,
 * `calcDiff()` builds a linear approximation of the constraint function with the form:
 * \f$\mathbf{g_x}\in\mathbb{R}^{ng\times ndx}\f$, \f$\mathbf{g_u}\in\mathbb{R}^{ng\times nu}\f$,
 * \f$\mathbf{h_x}\in\mathbb{R}^{nh\times ndx}\f$ \f$\mathbf{h_u}\in\mathbb{R}^{nh\times nu}\f$. Additionally, it is
 * important to note that `calcDiff()` computes the derivatives using the latest stored values by `calc()`. Thus, we
 * need to first run `calc()`.
 *
 * \sa `calc()`, `calcDiff()`, `createData()`
 */
template <typename _Scalar>
class ConstraintModelAbstractTpl {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ConstraintDataAbstractTpl<Scalar> ConstraintDataAbstract;
  typedef StateAbstractTpl<Scalar> StateAbstract;
  typedef typename MathBase::VectorXs VectorXs;

  /**
   * @copybrief Initialize the constraint model
   *
   * @param[in] state  State of the multibody system
   * @param[in] nu     Dimension of control vector
   * @param[in] ng     Number of inequality constraints
   * @param[in] nh     Number of equality constraints
   */
  ConstraintModelAbstractTpl(boost::shared_ptr<StateAbstract> state, const std::size_t nc, const std::size_t nu, 
                                                                    const VectorXs& lb, const VectorXs& ub);

  virtual ~ConstraintModelAbstractTpl();

  /**
   * @brief Compute the constraint value
   *
   * @param[in] data  Constraint data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
  virtual void calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u) = 0;

  /**
   * @brief Compute the Jacobian of the constraint
   *
   * It computes the Jacobian of the constraint function. It assumes that `calc()` has
   * been run first.
   *
   * @param[in] data  Constraint data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
   virtual void calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                        const Eigen::Ref<const VectorXs>& u) = 0;


  /**
   * @brief Create the constraint data
   *
   * The default data contains objects to store the values of the constraint, residual vector and their first
   * derivatives. However, it is possible to specialize this function is we need to create additional data, for
   * instance, to avoid dynamic memory allocation.
   *
   * @param data  Data collector
   * @return the constraint data
   */
  virtual boost::shared_ptr<ConstraintDataAbstract> createData();

  /**
   * @brief Update the lower and upper bounds the upper bound of constraint
   */
  void update_bounds(const VectorXs& lower, const VectorXs& upper);

  /**
   * @brief Remove the bounds of the constraint
   */
  void remove_bounds();

  /**
   * @brief Return the state
   */
  const boost::shared_ptr<StateAbstract>& get_state() const;

  /**
   * @brief Return the lower bound of the constraint
   */
  const VectorXs& get_lb() const;

  /**
   * @brief Return the upper bound of the constraint
   */
  const VectorXs& get_ub() const;

  /**
   * @brief Return the dimension of the control input
   */
  std::size_t get_nu() const;

  /**
   * @brief Return the number of inequality constraints
   */
  std::size_t get_nc() const;

  // /**
  //  * @brief Return the number of equality constraints
  //  */
  // std::size_t get_nh() const;

  /**
   * @brief Print information on the constraint model
   */
  // template <class Scalar>
  // friend std::ostream& operator<<(std::ostream& os, const CostModelAbstractTpl<Scalar>& model);

  /**
   * @brief Print relevant information of the constraint model
   *
   * @param[out] os  Output stream object
   */
  virtual void print(std::ostream& os) const;

 private:

 protected:
  boost::shared_ptr<StateAbstract> state_;             //!< State description
  VectorXs lb_;                                        //!< Lower bound of the constraint
  VectorXs ub_;                                        //!< Upper bound of the constraint
  std::size_t nc_;                                     //!< number of constraints
  std::size_t nu_;                                     //!< number of constraints
  VectorXs unone_;                                     //!< No control vector
};

template <typename _Scalar>
struct ConstraintDataAbstractTpl {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef DataCollectorAbstractTpl<Scalar> DataCollectorAbstract;
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;

  template <template <typename Scalar> class Model>
  ConstraintDataAbstractTpl(Model<Scalar>* const model)
      : c(model->get_nc()),
        Cx(model->get_nc(), model->get_state()->get_ndx()),
        Cu(model->get_nc(), model->get_nu()) {

    c.setZero();
    Cx.setZero();
    Cu.setZero();
  }
  virtual ~ConstraintDataAbstractTpl() {}

  boost::shared_ptr<ResidualDataAbstract> residual;  //!< Residual data
  VectorXs c;                                        //!< Inequality constraint values
  MatrixXs Cx;                                       //!< Jacobian of the inequality constraint
  MatrixXs Cu;                                       //!< Jacobian of the inequality constraint
};

}  // namespace crocoddyl

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "crocoddyl/core/constraint-base.hxx"

#endif  // CROCODDYL_CORE_CONSTRAINT_BASE_HPP_
