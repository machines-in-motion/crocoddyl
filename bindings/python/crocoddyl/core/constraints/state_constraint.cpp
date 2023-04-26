///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2020-2023, University of Edinburgh, Heriot-Watt University
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "python/crocoddyl/core/core.hpp"
#include "crocoddyl/core/constraints/state_constraint.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "python/crocoddyl/utils/printable.hpp"
#include "python/crocoddyl/utils/vector-converter.hpp"

namespace crocoddyl {
namespace python {

void exposeStateConstraint() {
  
  bp::register_ptr_to_python<boost::shared_ptr<StateConstraintModel> >();

  bp::class_<StateConstraintModel, bp::bases<ConstraintModelAbstract> >(
      "StateConstraintModel",
      "Abstract multibody constraint models.\n\n"
      "A constraint model defines both: inequality g(x,u) and equality h(x, u) constraints.\n"
      "The constraint function depends on the state point x, which lies in the state manifold\n"
      "described with a nx-tuple, its velocity xd that belongs to the tangent space with ndx dimension,\n"
      "and the control input u.",
      bp::init<boost::shared_ptr<StateAbstract>, std::size_t, std::size_t, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&>(
          bp::args("self", "state", "nc", "nu", "lb", "ub"),
          "Initialize the constraint model.\n\n"
          ":param state: state description\n"
          ":param nu: dimension of control vector (default state.nv)\n"
          ":param ng: number of inequality constraints\n"
          ":param nh: number of equality constraints"))
        .def<void (StateConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&, 
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &StateConstraintModel::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
        .def<void (StateConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&, 
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &StateConstraintModel::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)");
}

}  // namespace python
}  // namespace crocoddyl
