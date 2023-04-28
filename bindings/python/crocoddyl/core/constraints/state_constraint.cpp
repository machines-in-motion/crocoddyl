#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "crocoddyl/core/constraints/state_constraint.hpp"
#include "python/crocoddyl/utils/deprecate.hpp"

namespace crocoddyl {
namespace python {

void exposeStateConstraint() {
  bp::register_ptr_to_python<boost::shared_ptr<StateConstraintModel> >();

  bp::class_<StateConstraintModel, bp::bases<ConstraintModelAbstract> >(
      "StateConstraintModel",
      "This defines equality / inequality constraints based on a residual vector and its bounds.",
      bp::init<boost::shared_ptr<StateAbstract>, std::size_t, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&>
                                    (bp::args("self", "state", "nu", "lb", "ub"),
                                "Initialize the residual constraint model as an inequality constraint.\n\n"
                                ":param state: state description\n"
                                ":param lower: lower bound\n"
                                ":param upper: upper bound"))
      .def<void (StateConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &StateConstraintModel::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
      .def<void (StateConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &StateConstraintModel::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the derivatives of the residual constraint.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)\n")
      .def("createData", &StateConstraintModel::createData,
           bp::args("self"),
           "Create the residual constraint data.\n\n"
           "Each constraint model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined constraint.\n"
           ":param data: shared data\n"
           ":return constraint data.");

  bp::register_ptr_to_python<boost::shared_ptr<StateConstraintData> >();

  bp::class_<StateConstraintData, bp::bases<ConstraintDataAbstract> >(
      "StateConstraintData", "Data for residual constraint.\n\n",
      bp::init<StateConstraintModel*>(
          bp::args("self", "model"),
          "Create residual constraint data.\n\n"
          ":param model: residual constraint model\n"
          ":param data: shared data")[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 2> >()])
      .def(CopyableVisitor<StateConstraintData>());


}


}  // namespace python
}  // namespace crocoddyl