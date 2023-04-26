#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
StateConstraintModelTpl<Scalar>::StateConstraintModelTpl(boost::shared_ptr<StateAbstract> state, const std::size_t nc, const std::size_t nu, 
                                                                    const VectorXs& lb, const VectorXs& ub)

                                                        : Base(state, nc, nu, lb, ub), Ix(nc, nc), Iu(nc, nu){

        Iu.setZero();
        Ix = MatrixXs::Identity(nc, nc);
    }

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
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar> > ConstraintModelAbstractTpl<Scalar>::createData() {
  return boost::allocate_shared<ConstraintDataAbstract>(Eigen::aligned_allocator<ConstraintDataAbstract>(), this);
}
}