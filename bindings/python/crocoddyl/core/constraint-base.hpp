///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2020-2021, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef BINDINGS_PYTHON_CROCODDYL_CORE_CONSTRAINT_BASE_HPP_
#define BINDINGS_PYTHON_CROCODDYL_CORE_CONSTRAINT_BASE_HPP_

#include "crocoddyl/core/constraint-base.hpp"
#include "crocoddyl/core/utils/exception.hpp"

namespace crocoddyl {
namespace python {

class ConstraintModelAbstract_wrap : public ConstraintModelAbstract, public bp::wrapper<ConstraintModelAbstract> {
 public:

  ConstraintModelAbstract_wrap(boost::shared_ptr<StateAbstract> state, const std::size_t nc, const std::size_t nu, 
                              const VectorXs& lb,const VectorXs& ub)
      : ConstraintModelAbstract(state, nc, nu, lb, ub), bp::wrapper<ConstraintModelAbstract>() {}

  void calc(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const Eigen::VectorXd>& x,
            const Eigen::Ref<const Eigen::VectorXd>& u) {
    if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
      throw_pretty("Invalid argument: "
                   << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
    }
    // if (static_cast<std::size_t>(u.size()) != state_->get_nu()) {
    //   throw_pretty("Invalid argument: "
    //                << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
    // }
    return bp::call<void>(this->get_override("calc").ptr(), data, (Eigen::VectorXd)x, (Eigen::VectorXd)u);
  }

  void calcDiff(const boost::shared_ptr<ConstraintDataAbstract>& data, const boost::shared_ptr<ActionDataAbstract>& croc_data, const Eigen::Ref<const Eigen::VectorXd>& x,
                const Eigen::Ref<const Eigen::VectorXd>& u) {
    if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
      throw_pretty("Invalid argument: "
                   << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
    }
    // if (static_cast<std::size_t>(u.size()) != nu_) {
    //   throw_pretty("Invalid argument: "
    //                << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
    // }
    return bp::call<void>(this->get_override("calcDiff").ptr(), data, (Eigen::VectorXd)x, (Eigen::VectorXd)u);
  }

  boost::shared_ptr<ConstraintDataAbstract> default_createData() {
    return this->ConstraintModelAbstract::createData();
  }
};

}  // namespace python
}  // namespace crocoddyl

#endif  // BINDINGS_PYTHON_CROCODDYL_CORE_CONSTRAINT_BASE_HPP_
