#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
ControlConstraintModelTpl<Scalar>::ControlConstraintModelTpl(boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nu, 
                                                                    const VectorXs& lb, const VectorXs& ub)

                                                        : Base(state, nu, nu, lb, ub){

        Iu.resize(this->nc_, this->nu_); Ix.resize(this->nc_, this->state_->get_nx());
        Ix.setZero();
        Iu = MatrixXs::Identity(this->nc_, this->nu_);
    }

template <typename Scalar>
ControlConstraintModelTpl<Scalar>::~ControlConstraintModelTpl() {}

template <typename Scalar>
void ControlConstraintModelTpl<Scalar>::calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
                    
        data->c = u;
                    
    }

template <typename Scalar>
void ControlConstraintModelTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
        data->Cx = Ix;
        data->Cu = Iu;
                    
    }

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar>> ControlConstraintModelTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
void ControlConstraintModelTpl<Scalar>::print(std::ostream& os) const {
  os << "ConstraintModelResidual {""}";
}

}