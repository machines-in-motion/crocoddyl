#include "pinocchio/container/boost-container-limits.hpp"
#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/utils/copyable.hpp"
#include "crocoddyl/core/constraints/frame_translation_constraint.hpp"
#include "python/crocoddyl/utils/deprecate.hpp"

namespace crocoddyl {
namespace python {

void exposeFrameTranslationConstraint() {
  bp::register_ptr_to_python<boost::shared_ptr<FrameTranslationConstraintModel> >();

  bp::class_<FrameTranslationConstraintModel, bp::bases<ConstraintModelAbstract> >(
      "FrameTranslationConstraintModel",
      "This defines equality / inequality constraints based on a residual vector and its bounds.",
      bp::init<boost::shared_ptr<StateMultibody>, std::size_t, std::size_t, const Eigen::Ref<const Eigen::VectorXd>&, const Eigen::Ref<const Eigen::VectorXd>&, std::string>
                                    (bp::args("self", "state", "nu", "fid", "lb", "ub", "name"),
                                "Initialize the residual constraint model as an inequality constraint.\n\n"
                                ":param state: state description\n"
                                ":param lower: lower bound\n"
                                ":param upper: upper bound\n"
                                ":param name: name of the constraint"))
      .def<void (FrameTranslationConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &FrameTranslationConstraintModel::calc, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the residual constraint.\n\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)")
      .def<void (FrameTranslationConstraintModel::*)(const boost::shared_ptr<ConstraintDataAbstract>&,
                                            const boost::shared_ptr<ActionDataAbstract>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&,
                                             const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &FrameTranslationConstraintModel::calcDiff, bp::args("self", "data", "croc_data", "x", "u"),
          "Compute the derivatives of the residual constraint.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: constraint data\n"
          ":param croc_data: croco data\n"
          ":param x: state point (dim. state.nx)\n"
          ":param u: control input (dim. nu)\n")
      .def("createData", &FrameTranslationConstraintModel::createData,
           bp::args("self"),
           "Create the residual constraint data.\n\n"
           "Each constraint model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined constraint.\n"
           ":param data: shared data\n"
           ":return constraint data.");

  bp::register_ptr_to_python<boost::shared_ptr<FrameTranslationConstraintData> >();

  bp::class_<FrameTranslationConstraintData, bp::bases<ConstraintDataAbstract> >(
      "FrameTranslationConstraintData", "Data for residual constraint.\n\n",
      bp::init<FrameTranslationConstraintModel*>(
          bp::args("self", "model"),
          "Create residual constraint data.\n\n"
          ":param model: frame translation residual model"))
      .add_property("pinocchio",
                    bp::make_getter(&FrameTranslationConstraintData::pinocchio, bp::return_internal_reference<>()),
                    "pinocchio data")
      .add_property("fJf", bp::make_getter(&FrameTranslationConstraintData::fJf, bp::return_internal_reference<>()),
                    "local Jacobian of the frame")
      .def(CopyableVisitor<FrameTranslationConstraintData>());
}


}  // namespace python
}  // namespace crocoddyl