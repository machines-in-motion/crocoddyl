// state constraints

#ifndef CROCODDYL_CORE_CONSTRAINT_STACK_HPP_
#define CROCODDYL_CORE_CONSTRAINT_STACK_HPP_

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/constraint-base.hpp"

namespace crocoddyl{

template <typename _Scalar>
class ConstraintStackTpl :  public ConstraintModelAbstractTpl<_Scalar>{

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef _Scalar Scalar;
        typedef MathBaseTpl<Scalar> MathBase;
        typedef ConstraintModelAbstractTpl<Scalar> Base;
        typedef ConstraintDataStackTpl<Scalar> Data;
        typedef ConstraintDataAbstractTpl<Scalar> ConstraintDataAbstract;
        typedef typename MathBase::VectorXs VectorXs;
        typedef typename MathBase::MatrixXs MatrixXs;

    ConstraintStackTpl(std::vector<boost::shared_ptr<Base>> constraint_models, 
                                boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nc, std::size_t nu, const std::string name);

    virtual ~ConstraintStackTpl();
    
    virtual void calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u);

    virtual void calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                        const Eigen::Ref<const VectorXs>& u);
  
    virtual boost::shared_ptr<ConstraintDataAbstract> createData();

    virtual const std::vector<boost::shared_ptr<Base>> get_constraints();

    private:
        VectorXs tmp_lb_ = Eigen::VectorXd::Zero(4);                                        //!< Lower bound of the constraint
        VectorXs tmp_ub_ = Eigen::VectorXd::Zero(4);                                        //!< Upper bound of the constraint
        int nb_c;
        std::vector<boost::shared_ptr<Base>> constraint_models_;
        std::vector<boost::shared_ptr<ConstraintDataAbstract>> constraint_datas_;

};

template <typename _Scalar>
struct ConstraintDataStackTpl : public ConstraintDataAbstractTpl<_Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ConstraintDataAbstractTpl<Scalar> Base;

  template <template <typename Scalar> class Model>
  ConstraintDataStackTpl(Model<Scalar>* const model) : Base(model) {}

//   using Base::VectorXs c;                                        //!< Inequality constraint values
//   using Base::MatrixXs Cx;                                       //!< Jacobian of the inequality constraint
//   using Base::MatrixXs Cu;                                       //!< Jacobian of the inequality constraint
};


}
#include "crocoddyl/core/constraints/constraint_stack.hxx"
#endif  // CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

