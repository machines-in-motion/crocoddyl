#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
FrameTranslationConstraintModelTpl<Scalar>::FrameTranslationConstraintModelTpl(boost::shared_ptr<StateMultibody> state, std::size_t nu, std::size_t fid, 
                                                                    const VectorXs& lb, const VectorXs& ub)

                                                      : Base(state, 3, nu, lb, ub), fid_(fid), pin_model_(*state->get_pinocchio().get())
{
        // pinocchio::DataTpl<Scalar> pin_data_tmp(pin_model_);
        // pin_data_ = pin_data_tmp;
        Iu_.resize(this->nc_, this->nu_); Ix_.resize(this->nc_, this->state_->get_nx());
        Iu_.setZero(); Ix_.setZero();
    }

template <typename Scalar>
FrameTranslationConstraintModelTpl<Scalar>::~FrameTranslationConstraintModelTpl() {}

template <typename Scalar>
void FrameTranslationConstraintModelTpl<Scalar>::calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
        
        Data_croc* d = dynamic_cast<Data_croc*>(croc_data.get());

        // data->c = d->differential->pinocchio->oMf[fid_].translation;
                    
    }

template <typename Scalar>
void FrameTranslationConstraintModelTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){

        // pinocchio::getFrameJacobian(*pin_model_.get(), croc_data->differential->pinocchio, fid_, pinocchio::LOCAL, d->fJf);
        data->Cx = Ix_;
        data->Cu = Iu_;
                    
    }

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar>> FrameTranslationConstraintModelTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
void FrameTranslationConstraintModelTpl<Scalar>::print(std::ostream& os) const {
  os << "ConstraintModelResidual {""}";
}

}