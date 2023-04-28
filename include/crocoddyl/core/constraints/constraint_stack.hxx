#include "crocoddyl/core/utils/exception.hpp"


namespace crocoddyl{

template <typename Scalar>
ConstraintStackTpl<Scalar>::ConstraintStackTpl(std::vector<boost::shared_ptr<Base>> constraint_models,  
                                                  boost::shared_ptr<typename Base::StateAbstract> state, std::size_t nc, std::size_t nu) :
              constraint_models_(constraint_models),
              Base(state, nc, nu)
        {

          nb_c = constraint_models_.size(),

          constraint_datas_.resize(nb_c);
          int count = 0;
          for (std::size_t t = 0; t < nb_c; ++t){
            const boost::shared_ptr<ConstraintModelAbstract>& cmodel = constraint_models_[t]; 
            int  nc = cmodel->get_nc();

            const auto ub = cmodel->get_ub(); const auto lb = cmodel->get_lb();
            this->lb_.segment(count,nc) =  lb;
            this->ub_.segment(count,nc) =  ub;
            count += nc;
            constraint_datas_[t] = this->constraint_models_[t]->createData();
          }
    }

template <typename Scalar>
ConstraintStackTpl<Scalar>::~ConstraintStackTpl() {}

template <typename Scalar>
void ConstraintStackTpl<Scalar>::calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){
                    
        int count = 0;
        for (std::size_t t = 0; t < nb_c; ++t){
            const boost::shared_ptr<ConstraintModelAbstract>& cmodel = constraint_models_[t]; 
            int  nc = cmodel->get_nc();
            constraint_models_[t]->calc(constraint_datas_[t], croc_data, x, u);
            data->c.segment(count, nc) = constraint_datas_[t]->c;
            count += nc;
          }
            
    }

template <typename Scalar>
void ConstraintStackTpl<Scalar>::calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u){

        int count = 0;
        for (std::size_t t = 0; t < nb_c; ++t){
            const boost::shared_ptr<ConstraintModelAbstract>& cmodel = constraint_models_[t]; 
            int  nc = cmodel->get_nc();

            constraint_models_[t]->calcDiff(constraint_datas_[t], croc_data, x, u);
            data->Cx.middleRows(count, nc) = constraint_datas_[t]->Cx;
            data->Cu.middleRows(count, nc) = constraint_datas_[t]->Cu;
            count += nc;
          }
                    
    }

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar>> ConstraintStackTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

}