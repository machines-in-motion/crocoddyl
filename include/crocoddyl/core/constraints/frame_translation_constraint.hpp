// state constraints

#ifndef CROCODDYL_CORE_FRAME_CONSTRAINT_HPP_
#define CROCODDYL_CORE_FRAME_CONSTRAINT_HPP_

#include <pinocchio/algorithm/frames.hpp>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/constraint-base.hpp"
#include "crocoddyl/multibody/states/multibody.hpp"
#include "crocoddyl/core/integ-action-base.hpp"
#include "crocoddyl/core/integrator/euler.hpp"
#include "crocoddyl/multibody/actions/contact-fwddyn.hpp"

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
class FrameTranslationConstraintModelTpl : public ConstraintModelAbstractTpl<_Scalar>{

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef _Scalar Scalar;
        typedef MathBaseTpl<Scalar> MathBase;
        typedef ConstraintModelAbstractTpl<Scalar> Base;
        typedef StateMultibodyTpl<Scalar> StateMultibody;
        typedef FrameTranslationConstraintDataTpl<Scalar> Data;
        typedef IntegratedActionDataEulerTpl<Scalar> IADEuler;
        typedef DifferentialActionDataContactFwdDynamicsTpl<Scalar> DADContact;
        typedef ConstraintDataAbstractTpl<Scalar> ConstraintDataAbstract;
        typedef typename MathBase::VectorXs VectorXs;
        typedef typename MathBase::MatrixXs MatrixXs;

    FrameTranslationConstraintModelTpl(boost::shared_ptr<StateMultibody> state, std::size_t nu, std::size_t fid,  
                                                                    const VectorXs& lb, const VectorXs& ub);

    virtual ~FrameTranslationConstraintModelTpl();
    
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


    private:
        MatrixXs Ix_;                
        MatrixXs Iu_;               
        const std::size_t fid_;
        boost::shared_ptr<typename StateMultibody::PinocchioModel> pin_model_;
};

template <typename _Scalar>
struct FrameTranslationConstraintDataTpl : public ConstraintDataAbstractTpl<_Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ConstraintDataAbstractTpl<Scalar> Base;
  typedef typename MathBase::Matrix6xs Matrix6xs;

  template <template <typename Scalar> class Model>
  FrameTranslationConstraintDataTpl(Model<Scalar>* const model) : Base(model), fJf(6, model->get_state()->get_nv()) {
    fJf.setZero();
  }
  
  pinocchio::DataTpl<Scalar>* pinocchio;
  Matrix6xs fJf; 
};

}
#include "crocoddyl/core/constraints/frame_translation_constraint.hxx"
#endif  // CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

