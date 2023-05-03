// state constraints

#ifndef CROCODDYL_CORE_NO_CONSTRAINT_HPP_
#define CROCODDYL_CORE_NO_CONSTRAINT_HPP_

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/constraint-base.hpp"

namespace crocoddyl{
/**
 * @brief Residual-based constraint
 *
 * This constraint function uses a residual model to define equality / inequality constraint as \f[
 * \mathbf{\underline{r}} \leq \mathbf{r}(\mathbf{x}, \mathbf{u}) \leq \mathbf{\bar{r}} \f] where
 * \f$\mathbf{r}(\cdot)\f$ describes the residual function, and \f$\mathbf{\underline{r}}\f$, \f$\mathbf{\bar{r}}\f$
 * are the lower and upper bounds, respectively. We can define element-wise equality constraints by defining the same
 * value for both: lower and upper values. Additionally, if we do not define the bounds, then it is assumed that
 * \f$\mathbf{\underline{r}}=\mathbf{\bar{r}}=\mathbf{0}\f$.
 *
 * The main computations are carring out in `calc` and `calcDiff` routines. `calc` computes the constraint residual and
 * `calcDiff` computes the Jacobians of the constraint function. Concretely speaking, `calcDiff` builds
 * a linear approximation of the constraint function with the form: \f$\mathbf{g_x}\in\mathbb{R}^{ng\times ndx}\f$,
 * \f$\mathbf{g_u}\in\mathbb{R}^{ng\times nu}\f$, \f$\mathbf{h_x}\in\mathbb{R}^{nh\times ndx}\f$
 * \f$\mathbf{h_u}\in\mathbb{R}^{nh\times nu}\f$.
 * Additionally, it is important to note that `calcDiff()` computes the derivatives using the latest stored values by
 * `calc()`. Thus, we need to run first `calc()`.
 *
 * \sa `ConstraintModelAbstractTpl`, `calc()`, `calcDiff()`, `createData()`
 */

template <typename _Scalar>
class NoConstraintModelTpl : public ConstraintModelAbstractTpl<_Scalar>{

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef _Scalar Scalar;
        typedef MathBaseTpl<Scalar> MathBase;
        typedef ConstraintModelAbstractTpl<Scalar> Base;
        typedef NoConstraintDataTpl<Scalar> Data;
        typedef ConstraintDataAbstractTpl<Scalar> ConstraintDataAbstract;
        typedef typename MathBase::VectorXs VectorXs;
        typedef typename MathBase::MatrixXs MatrixXs;

    NoConstraintModelTpl(boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nu, const std::string name);

    virtual ~NoConstraintModelTpl();
    
    virtual void calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u);

//   /**
//    * @brief Compute the Jacobian of the constraint
//    *
//    * It computes the Jacobian of the constraint function. It assumes that `calc()` has
//    * been run first.
//    *
//    * @param[in] data  Constraint data
//    * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
//    * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
//    */
    virtual void calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                        const Eigen::Ref<const VectorXs>& u);
  
    virtual boost::shared_ptr<ConstraintDataAbstract> createData();
    virtual void print(std::ostream& os) const;

};

template <typename _Scalar>
struct NoConstraintDataTpl : public ConstraintDataAbstractTpl<_Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ConstraintDataAbstractTpl<Scalar> Base;

  template <template <typename Scalar> class Model>
  NoConstraintDataTpl(Model<Scalar>* const model) : Base(model) {}

//   using Base::VectorXs c;                                        //!< Inequality constraint values
//   using Base::MatrixXs Cx;                                       //!< Jacobian of the inequality constraint
//   using Base::MatrixXs Cu;                                       //!< Jacobian of the inequality constraint
};

}
#include "crocoddyl/core/constraints/no_constraint.hxx"
#endif  // CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

