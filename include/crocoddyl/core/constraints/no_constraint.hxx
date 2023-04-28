#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
NoConstraintModelTpl<Scalar>::NoConstraintModelTpl(boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nu)

                                                        : Base(state, 0, nu){

    }

template <typename Scalar>
NoConstraintModelTpl<Scalar>::~NoConstraintModelTpl() {}

template <typename Scalar>
void NoConstraintModelTpl<Scalar>::calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
                                        
    }

template <typename Scalar>
void NoConstraintModelTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
    }

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar>> NoConstraintModelTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
void NoConstraintModelTpl<Scalar>::print(std::ostream& os) const {
  os << "ConstraintModelResidual {""}";
}

}