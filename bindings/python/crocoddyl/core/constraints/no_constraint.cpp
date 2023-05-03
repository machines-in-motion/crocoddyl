#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "crocoddyl/core/constraints/no_constraint.hpp"
#include "python/crocoddyl/utils/deprecate.hpp"

namespace crocoddyl {
namespace python {

void exposeNoConstraint() {
  bp::register_ptr_to_python<boost::shared_ptr<NoConstraintModel> >();

  bp::class_<NoConstraintModel, bp::bases<ConstraintModelAbstract> >(
      "NoConstraintModel",
      "This defines equality / inequality constraints based on a residual vector and its bounds.",
      bp::init<boost::shared_ptr<StateAbstract>, std::size_t, std::string>
                                    (bp::args("self", "state", "nu"),
                                "Initialize the residual constraint model as an inequality constraint.\n\n"
                                ":param state: state description\n"
                                ":param lower: lower bound\n"
                                ":param upper: upper bound"))
      .def<void (NoConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &NoConstraintModel::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
      .def<void (NoConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &NoConstraintModel::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the derivatives of the residual constraint.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)\n")
      .def("createData", &NoConstraintModel::createData,
           bp::args("self"),
           "Create the residual constraint data.\n\n"
           "Each constraint model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined constraint.\n"
           ":param data: shared data\n"
           ":return constraint data.");

  bp::register_ptr_to_python<boost::shared_ptr<NoConstraintData> >();

  bp::class_<NoConstraintData, bp::bases<ConstraintDataAbstract> >(
      "NoConstraintData", "Data for residual constraint.\n\n",
      bp::init<NoConstraintModel*>(
          bp::args("self", "model"),
          "Create residual constraint data.\n\n"
          ":param model: residual constraint model\n"
          ":param data: shared data")[bp::with_custodian_and_ward<1, 2, bp::with_custodian_and_ward<1, 2> >()])
      .def(CopyableVisitor<NoConstraintData>());


}


}  // namespace python
}  // namespace crocoddyl