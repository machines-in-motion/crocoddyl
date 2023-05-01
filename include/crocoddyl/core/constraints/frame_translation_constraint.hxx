#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
FrameTranslationConstraintModelTpl<Scalar>::FrameTranslationConstraintModelTpl(boost::shared_ptr<StateMultibody> state, std::size_t nu, std::size_t fid, 
                                                                    const VectorXs& lb, const VectorXs& ub)

                                                      : Base(state, 3, nu, lb, ub), fid_(fid), pin_model_(state->get_pinocchio())
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
        
        IADEuler* iad_euler = static_cast<IADEuler*>(croc_data.get());
        DADContact* dad_contact = static_cast<DADContact*>(iad_euler->differential.get());
        data->c = dad_contact->pinocchio.oMf[fid_].translation();
                    
    }

template <typename Scalar>
void FrameTranslationConstraintModelTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){


        IADEuler* iad_euler = static_cast<IADEuler*>(croc_data.get());
        DADContact* dad_contact = static_cast<DADContact*>(iad_euler->differential.get());
        Data* derived_data = static_cast<Data*>(data.get());
        pinocchio::getFrameJacobian(*pin_model_.get(), *derived_data->pinocchio, fid_, pinocchio::LOCAL, derived_data->fJf);
        data->Cx.leftCols(pin_model_->nv) = derived_data->pinocchio->oMf[fid_].rotation() * derived_data->fJf.template topRows<3>();
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