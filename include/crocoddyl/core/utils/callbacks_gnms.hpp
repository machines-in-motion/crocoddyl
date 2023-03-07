///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019, LAAS-CNRS
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_CORE_UTILS_CALLBACKS_GNMS_HPP_
#define CROCODDYL_CORE_UTILS_CALLBACKS_GNMS_HPP_

#include <iostream>
#include <iomanip>

#include "crocoddyl/core/utils/callbacks.hpp"

namespace crocoddyl {

// enum VerboseLevel { _1 = 0, _2 };
class CallbackVerboseGNMS : public CallbackVerbose {
 public:
  explicit CallbackVerboseGNMS(VerboseLevel level = _1);
  ~CallbackVerboseGNMS();

  virtual void operator()(SolverAbstract& solver);

 private:
  VerboseLevel level;
};

}  // namespace crocoddyl

#endif  // CROCODDYL_CORE_UTILS_CALLBACKS_GNMS_HPP_
