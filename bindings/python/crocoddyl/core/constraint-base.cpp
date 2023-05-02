///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2020-2023, University of Edinburgh, Heriot-Watt University
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/core/constraint-base.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "python/crocoddyl/utils/printable.hpp"
#include "python/crocoddyl/utils/vector-converter.hpp"

namespace crocoddyl {
namespace python {

void exposeConstraintAbstract() {
  
  typedef boost::shared_ptr<ConstraintModelAbstract> ConstraintModelPtr;
  typedef boost::shared_ptr<ConstraintDataAbstract> ConstraintDataPtr;
  StdVectorPythonVisitor<ConstraintModelPtr, std::allocator<ConstraintModelPtr>, true>::expose("StdVec_ConstraintModel");
  StdVectorPythonVisitor<ConstraintDataPtr, std::allocator<ConstraintDataPtr>, true>::expose("StdVec_ConstraintData");

  bp::register_ptr_to_python<boost::shared_ptr<ConstraintModelAbstract> >();
  bp::class_<ConstraintModelAbstract_wrap, boost::noncopyable>(
      "ConstraintModelAbstract",
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
      .def("calc", pure_virtual(&ConstraintModelAbstract_wrap::calc), bp::args("self", "data", "croc_data", "x", "u"),
           "Compute the constraint value.\n\n"
           ":param data: constraint data\n"
           ":param x: state point (dim. state.nx)\n"
           ":param u: control input (dim. nu)")

      .def("calcDiff", pure_virtual(&ConstraintModelAbstract_wrap::calcDiff), bp::args("self", "data", "croc_data", "x", "u"),
           "Compute the Jacobians of the constraint function.\n\n"
           "It computes the Jacobians of the constraint function.\n"
           "It assumes that calc has been run first.\n"
           ":param data: constraint data\n"
           ":param x: state point (dim. state.nx)\n"
           ":param u: control input (dim. nu)\n")

      .def("createData", &ConstraintModelAbstract_wrap::default_createData,
           bp::args("self"))
      .def("updateBounds", &ConstraintModelAbstract_wrap::update_bounds, bp::args("self", "lower", "upper"),
           "Update the lower and upper bounds.\n\n"
           ":param lower: lower bound\n"
           ":param upper: upper bound")
      .def("removeBounds", &ConstraintModelAbstract_wrap::remove_bounds, bp::args("self"), "Remove the bounds.")
      .add_property(
          "state",
          bp::make_function(&ConstraintModelAbstract_wrap::get_state, bp::return_value_policy<bp::return_by_value>()),
          "state description")
      .add_property("lb", bp::make_function(&ConstraintModelAbstract_wrap::get_lb, bp::return_internal_reference<>()),
                          bp::make_function(&ConstraintModelAbstract_wrap::set_lb),
                    "lower bound of constraint")

      .add_property("ub", bp::make_function(&ConstraintModelAbstract_wrap::get_ub, bp::return_internal_reference<>()),
                          bp::make_function(&ConstraintModelAbstract_wrap::set_ub),
                    "upper bound of constraint")
      .add_property("nc", bp::make_function(&ConstraintModelAbstract_wrap::get_nc), "dimension of constraint vector")
      .add_property("nu", bp::make_function(&ConstraintModelAbstract_wrap::get_nu), "dimension of control vector")
      .def(CopyableVisitor<ConstraintModelAbstract_wrap>())
      .def(PrintableVisitor<ConstraintModelAbstract>());

  bp::register_ptr_to_python<boost::shared_ptr<ConstraintDataAbstract> >();

  bp::class_<ConstraintDataAbstract, boost::noncopyable>(
      "ConstraintDataAbstract", "Abstract class for constraint data.\n\n",
      bp::init<ConstraintModelAbstract*>(
          bp::args("self", "model"),
          "Create common data shared between constraint models.\n\n"
          ":param model: constraint model\n"
          ":param data: shared data")[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()])
      .add_property("c", bp::make_getter(&ConstraintDataAbstract::c, bp::return_internal_reference<>()),
                    bp::make_setter(&ConstraintDataAbstract::c), "inequality constraint residual")
      .add_property("Cx", bp::make_getter(&ConstraintDataAbstract::Cx, bp::return_internal_reference<>()),
                    bp::make_setter(&ConstraintDataAbstract::Cx), "Jacobian of the inequality constraint")
      .add_property("Cu", bp::make_getter(&ConstraintDataAbstract::Cu, bp::return_internal_reference<>()),
                    bp::make_setter(&ConstraintDataAbstract::Cu), "Jacobian of the inequality constraint")
      .def(CopyableVisitor<ConstraintDataAbstract>());
}

}  // namespace python
}  // namespace crocoddyl
