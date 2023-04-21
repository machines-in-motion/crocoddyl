///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, The University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "python/crocoddyl/core/core.hpp"
#include "crocoddyl/core/solvers/fadmm.hpp"

namespace crocoddyl {
namespace python {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SolverFADMM_solves, SolverFADMM::solve, 0, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(SolverFADMM_computeDirections, SolverDDP::computeDirection, 0, 1)

void exposeSolverFADMM() {
  bp::register_ptr_to_python<boost::shared_ptr<SolverFADMM> >();

  bp::class_<SolverFADMM, bp::bases<SolverDDP> >(
      "SolverFADMM",
      "GNMS solver.\n\n"
      "The GNMS solver computes an optimal trajectory and control commands by iterates\n"
      "running backward and forward passes. The backward-pass updates locally the\n"
      "quadratic approximation of the problem and computes descent direction,\n"
      "and the forward-pass rollouts this new policy by integrating the system dynamics\n"
      "along a tuple of optimized control commands U*.\n"
      ":param shootingProblem: shooting problem (list of action models along trajectory.)",
      bp::init<boost::shared_ptr<ShootingProblem>, std::vector<boost::shared_ptr<ConstraintModelAbstract>> >(bp::args("self", "problem", "constraint_models"),
                                                    "Initialize the vector dimension.\n\n"
                                                    ":param problem: shooting problem."))
      .def("solve", &SolverFADMM::solve,
           SolverFADMM_solves(
               bp::args("self", "init_xs", "init_us", "maxiter", "isFeasible", "regInit"),
               "Compute the optimal trajectory xopt, uopt as lists of T+1 and T terms.\n\n"
               "From an initial guess init_xs,init_us (feasible or not), iterate\n"
               "over computeDirection and tryStep until stoppingCriteria is below\n"
               "threshold. It also describes the globalization strategy used\n"
               "during the numerical optimization.\n"
               ":param init_xs: initial guess for state trajectory with T+1 elements (default [])\n"
               ":param init_us: initial guess for control trajectory with T elements (default []).\n"
               ":param maxiter: maximum allowed number of iterations (default 100).\n"
               ":param isFeasible: true if the init_xs are obtained from integrating the init_us (rollout) (default "
               "False).\n"
               ":param regInit: initial guess for the regularization value. Very low values are typical\n"
               "                used with very good guess points (init_xs, init_us) (default None).\n"
               ":returns the optimal trajectory xopt, uopt and a boolean that describes if convergence was reached."))
     //  .def("updateExpectedImprovement", &SolverFADMM::updateExpectedImprovement,
     //       bp::return_value_policy<bp::copy_const_reference>(), bp::args("self"),
     //       "Update the expected improvement model\n\n")
     // .def("increaseRegularization", &solverFDDP::increaseRegularization, bp::args("self"),
     //       "Increase regularization")

      .def("calc", &SolverFADMM::calculate, bp::args("self", "recalc"),
           "")
      .def("update_lagrangian_parameters", &SolverFADMM::update_lagrangian_parameters, bp::args("self"),
           "")
      .def("forwardPass", &SolverFADMM::forwardPass, bp::args("self"),
           "")
      .def("backwardPass", &SolverFADMM::backwardPass, bp::args("self"),
           "")
      .def("update_rho_sparse", &SolverFADMM::update_rho_sparse, bp::args("self", "iter"),
           "")
      .def("computeDirection", &SolverFADMM::computeDirection, bp::args("self", "recalcDiff"),
           "")

      .def_readwrite("xs_try", &SolverFADMM::xs_try_, "xs try")
      .def_readwrite("us_try", &SolverFADMM::us_try_, "us try")  
      .def_readwrite("cost_try", &SolverFADMM::cost_try_, "cost try")
      .def_readwrite("fs_try", &SolverFADMM::fs_try_, "fs_try")
      .def_readwrite("lag_mul", &SolverFADMM::lag_mul_, "lagrange multipliers")
      .def_readwrite("KKT", &SolverFADMM::KKT_, "KKT residual")

      .add_property("with_callbacks", bp::make_function(&SolverFADMM::getCallbacks), bp::make_function(&SolverFADMM::setCallbacks),
                    "Activates the callbacks when true (default: False)")
      .add_property("use_kkt_criteria", bp::make_function(&SolverFADMM::set_use_kkt_criteria), bp::make_function(&SolverFADMM::get_use_kkt_criteria),
                    "Use the KKT residual condition as a termination criteria (default: True)")
      .add_property("mu", bp::make_function(&SolverFADMM::get_mu), bp::make_function(&SolverFADMM::set_mu),
                    "Penalty term for dynamic violation in the merit function (default: 1.)")
      .add_property("xs", make_function(&SolverFADMM::get_xs, bp::return_value_policy<bp::copy_const_reference>()), "xs")
      .add_property("us", make_function(&SolverFADMM::get_us, bp::return_value_policy<bp::copy_const_reference>()), "us")
      .add_property("dx_tilde", make_function(&SolverFADMM::get_xs_tilde, bp::return_value_policy<bp::copy_const_reference>()), "dx_tilde")
      .add_property("du_tilde", make_function(&SolverFADMM::get_us_tilde, bp::return_value_policy<bp::copy_const_reference>()), "du_tilde")
      
      .add_property("y", make_function(&SolverFADMM::get_y, bp::return_value_policy<bp::copy_const_reference>()), "y")
      .add_property("z", make_function(&SolverFADMM::get_z, bp::return_value_policy<bp::copy_const_reference>()), "z")

      .add_property("mu", bp::make_function(&SolverFADMM::get_mu), bp::make_function(&SolverFADMM::set_mu),
                    "Penalty term for dynamic violation in the merit function (default: 1.)")
      .add_property("use_heuristic_line_search", bp::make_function(&SolverFADMM::get_use_heuristic_line_search), bp::make_function(&SolverFADMM::set_use_heuristic_line_search),
                    "Use the heuristic line search criteria (default: False)")
      .add_property("termination_tolerance", bp::make_function(&SolverFADMM::get_termination_tolerance), bp::make_function(&SolverFADMM::set_termination_tolerance),
                    "Termination criteria to exit the iteration (default: 1e-8)")
      .add_property("max_qp_iters", bp::make_function(&SolverFADMM::get_max_qp_iters), bp::make_function(&SolverFADMM::set_max_qp_iters),
                    "get and set max qp iters");
     //  .add_property("th_acceptNegStep", bp::make_function(&SolverFADMM::get_th_acceptnegstep),
     //                bp::make_function(&SolverFADMM::set_th_acceptnegstep),
     //                "threshold for step acceptance in ascent direction");
     
}

}  // namespace python
}  // namespace crocoddyl
