///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2020-2022, University of Edinburgh, Heriot-Watt University
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <boost/core/demangle.hpp>

namespace crocoddyl {

template <typename Scalar>
ConstraintModelAbstractTpl<Scalar>::ConstraintModelAbstractTpl(boost::shared_ptr<StateAbstract> state,
                                                               const std::size_t nc, const std::size_t nu,
                                                               const VectorXs& lb, 
                                                               const VectorXs& ub)
    : state_(state),
      lb_(lb),
      ub_(ub),
      nc_(nc),
      nu_(nu) {}

template <typename Scalar>
ConstraintModelAbstractTpl<Scalar>::~ConstraintModelAbstractTpl() {}

template <typename Scalar>
boost::shared_ptr<ConstraintDataAbstractTpl<Scalar> > ConstraintModelAbstractTpl<Scalar>::createData() {
  return boost::allocate_shared<ConstraintDataAbstract>(Eigen::aligned_allocator<ConstraintDataAbstract>(), this);
}

template <typename Scalar>
void ConstraintModelAbstractTpl<Scalar>::update_bounds(const VectorXs& lower, const VectorXs& upper) {
  if (static_cast<std::size_t>(upper.size()) != nc_ ||
      static_cast<std::size_t>(lower.size()) != nc_) {
    throw_pretty("Invalid argument: the dimension of the lower/upper bound is not the same to ng.")
  }
  if (((upper - lower).array() <= 0.).any()) {
    throw_pretty("Invalid argument: the upper bound is not higher than the lower bound.")
  }
  if ((lb_.array() == std::numeric_limits<Scalar>::infinity()).any() ||
      (lb_.array() == std::numeric_limits<Scalar>::max()).any()) {
    throw_pretty("Invalid argument: the lower bound cannot contain a positive infinity/max value");
  }
  if ((ub_.array() == -std::numeric_limits<Scalar>::infinity()).any() ||
      (ub_.array() == -std::numeric_limits<Scalar>::infinity()).any()) {
    throw_pretty("Invalid argument: the lower bound cannot contain a negative infinity/min value");
  }
  lb_ = lower;
  ub_ = upper;
}

template <typename Scalar>
void ConstraintModelAbstractTpl<Scalar>::remove_bounds() {
  lb_.setConstant(-std::numeric_limits<Scalar>::infinity());
  ub_.setConstant(std::numeric_limits<Scalar>::infinity());
}

template <typename Scalar>
void ConstraintModelAbstractTpl<Scalar>::print(std::ostream& os) const {
  os << boost::core::demangle(typeid(*this).name());
}

template <typename Scalar>
const boost::shared_ptr<StateAbstractTpl<Scalar> >& ConstraintModelAbstractTpl<Scalar>::get_state() const {
  return state_;
}

template <typename Scalar>
const typename MathBaseTpl<Scalar>::VectorXs& ConstraintModelAbstractTpl<Scalar>::get_lb() const {
  return lb_;
}

template <typename Scalar>
const typename MathBaseTpl<Scalar>::VectorXs& ConstraintModelAbstractTpl<Scalar>::get_ub() const {
  return ub_;
}

template <typename Scalar>
std::size_t ConstraintModelAbstractTpl<Scalar>::get_nc() const {
  return nc_;
}

// template <typename Scalar>
// std::size_t ConstraintModelAbstractTpl<Scalar>::get_ng() const {
//   return ng_;
// }

template <typename Scalar>
std::size_t ConstraintModelAbstractTpl<Scalar>::get_nu() const {
  return nu_;
}

template <class Scalar>
std::ostream& operator<<(std::ostream& os, const ConstraintModelAbstractTpl<Scalar>& model) {
  model.print(os);
  return os;
}

}  // namespace crocoddyl
