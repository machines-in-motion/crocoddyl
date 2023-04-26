// state constraints

#ifndef CROCODDYL_CORE_STATE_CONSTRAINT_HPP_
#define CROCODDYL_CORE_STATE_CONSTRAINT_HPP_

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/constraint-base.hpp"


namespace crocoddyl{

template <typename _Scalar>

class StateConstraintModelTpl : public ConstraintModelAbstractTpl<_Scalar>{

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef _Scalar Scalar;
        typedef MathBaseTpl<Scalar> MathBase;
        typedef ConstraintModelAbstractTpl<Scalar> Base;
        typedef ConstraintDataAbstractTpl<Scalar> ConstraintDataAbstract;
        typedef typename MathBase::VectorXs VectorXs;
        typedef typename MathBase::MatrixXs MatrixXs;
        typedef StateAbstractTpl<Scalar> StateAbstract;


     StateConstraintModelTpl(boost::shared_ptr<StateAbstract> state, const std::size_t nc, const std::size_t nu, 
                                                                    const VectorXs& lb, const VectorXs& ub);

    virtual ~StateConstraintModelTpl();
    
    void calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
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

  
    boost::shared_ptr<ConstraintDataAbstract> createData();

    private:
        MatrixXs Ix;                
        MatrixXs Iu;               


};

}
#include "crocoddyl/core/constraints/state_constraint.hxx"
#endif  // CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

