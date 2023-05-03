#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "crocoddyl/core/constraints/constraint_stack.hpp"
#include "python/crocoddyl/utils/deprecate.hpp"

namespace crocoddyl {
namespace python {

void exposeConstraintStack() {
  bp::register_ptr_to_python<boost::shared_ptr<ConstraintStack> >();

  bp::class_<ConstraintStack, bp::bases<ConstraintModelAbstract>>(
      "ConstraintStack",
      "This defines equality / inequality constraints based on a residual vector and its bounds.",
    bp::init<std::vector<boost::shared_ptr<ConstraintModelAbstract>>, boost::shared_ptr<StateAbstract>, std::size_t, std::size_t, const std::string>
                                    (bp::args("self","constraint_models","state","nc","nu", "name"),
                                "Initialize the residual constraint model as an inequality constraint.\n\n"
                                ":param constraint_models: list of constraint models for time step\n"
                                 ":param constraint_models: list of constraint models for time step\n"
                                  ":param constraint_models: list of constraint models for time step\n."))
      .def<void (ConstraintStack::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &ConstraintStack::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
      .def<void (ConstraintStack::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &ConstraintStack::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the derivatives of the residual constraint.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)\n")
        .def("get_constraints", &ConstraintStack::get_constraints,
           bp::args("self"))
      .def("createData", &ConstraintStack::createData,
           bp::args("self"),
           "Create the residual constraint data.\n\n"
           "Each constraint model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined constraint.\n"
           ":param data: shared data\n"
           ":return constraint data.");

    bp::register_ptr_to_python<boost::shared_ptr<ConstraintDataStack> >();

  bp::class_<ConstraintDataStack, bp::bases<ConstraintDataAbstract> >(
      "StateConstraintData", "Data for residual constraint.\n\n",
      bp::init<ConstraintStack*>(
          bp::args("self", "model"),
          "Create residual constraint data.\n\n"
          ":param model: residual constraint model\n"
          ":param data: shared data")[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 3> >()])
      .def(CopyableVisitor<ConstraintDataStack>());

}


}  // namespace python
}  // namespace crocoddyl