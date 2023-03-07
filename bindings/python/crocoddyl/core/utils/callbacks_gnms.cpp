///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "python/crocoddyl/core/core.hpp"
// #include "python/crocoddyl/core/core.hpp"
#include "crocoddyl/core/utils/callbacks_gnms.hpp"

namespace crocoddyl {
namespace python {

void exposeCallbacksGNMS() {
  bp::register_ptr_to_python<boost::shared_ptr<CallbackVerboseGNMS> >();

  bp::enum_<VerboseLevel>("VerboseLevel").value("_1", _1).value("_2", _2);

  bp::class_<CallbackVerboseGNMS, bp::bases<CallbackAbstract> >(
      "CallbackVerboseGNMS", "Callback function for printing the solver values.",
      bp::init<bp::optional<VerboseLevel> >(bp::args("self", "level"),
                                            "Initialize the differential verbose callback.\n\n"
                                            ":param level: verbose level (default _1)"))
      .def("__call__", &CallbackVerboseGNMS::operator(), bp::args("self", "solver"),
           "Run the callback function given a solver.\n\n"
           ":param solver: solver to be diagnostic");
}

}  // namespace python
}  // namespace crocoddyl
