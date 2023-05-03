#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "crocoddyl/core/constraints/control_constraint.hpp"
#include "python/crocoddyl/utils/deprecate.hpp"

namespace crocoddyl {
namespace python {

void exposeControlConstraint() {
  bp::register_ptr_to_python<boost::shared_ptr<ControlConstraintModel> >();

  bp::class_<ControlConstraintModel, bp::bases<ConstraintModelAbstract> >(
      "ControlConstraintModel",
      "This defines equality / inequality constraints based on a residual vector and its bounds.",
      bp::init<boost::shared_ptr<StateAbstract>, std::size_t, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&, const std::string>
                                    (bp::args("self", "state", "nu", "lb", "ub", "name"),
                                "Initialize the residual constraint model as an inequality constraint.\n\n"
                                ":param state: state description\n"
                                ":param residual: residual model\n"
                                ":param lower: lower bound\n"
                                ":param upper: upper bound\n"
                                ":param name: name of the constraint"))
      .def<void (ControlConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &ControlConstraintModel::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
      .def<void (ControlConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &ControlConstraintModel::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the derivatives of the residual constraint.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)\n")
      .def("createData", &ControlConstraintModel::createData,
           bp::args("self"),
           "Create the residual constraint data.\n\n"
           "Each constraint model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined constraint.\n"
           ":param data: shared data\n"
           ":return constraint data.");

  bp::register_ptr_to_python<boost::shared_ptr<ControlConstraintData> >();

  bp::class_<ControlConstraintData, bp::bases<ConstraintDataAbstract> >(
      "ControlConstraintData", "Data for residual constraint.\n\n",
      bp::init<ControlConstraintModel*>(
          bp::args("self", "model"),
          "Create residual constraint data.\n\n"
          ":param model: residual constraint model\n"
          ":param data: shared data")[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 2> >()])
      .def(CopyableVisitor<ControlConstraintData>());


}


}  // namespace python
}  // namespace crocoddyl