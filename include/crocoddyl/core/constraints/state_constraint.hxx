#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
StateConstraintModelTpl<Scalar>::StateConstraintModelTpl(boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nu, 
                                                                    const VectorXs& lb, const VectorXs& ub)

                                                        : Base(state, state->get_nx(), nu, lb, ub){

        Iu.resize(this->nc_, this->nu_); Ix.resize(this->nc_, this->state_->get_nx());
        Iu.setZero();
        Ix = MatrixXs::Identity(this->nc_, this->state_->get_nx());
    }

template <typename Scalar>
StateConstraintModelTpl<Scalar>::~StateConstraintModelTpl() {}

template <typename Scalar>
void StateConstraintModelTpl<Scalar>::calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
                    
        data->c = x;
                    
    }

template <typename Scalar>
void StateConstraintModelTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
        data->Cx = Ix;
        data->Cu = Iu;
                    
    }

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar>> StateConstraintModelTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
void StateConstraintModelTpl<Scalar>::print(std::ostream& os) const {
  os << "ConstraintModelResidual {""}";
}

}