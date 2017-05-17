/*
 * Copyright (c) 2016, Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Helen Oleynikova, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Rik Bähnemann, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Marija Popovic, ASL, ETH Zurich, Switzerland
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "mav_trajectory_generation_ros/feasibility_base.h"

#include <cmath>
#include <limits>

#include <Eigen/Geometry>

#include <mav_msgs/default_values.h>

namespace mav_trajectory_generation {

std::string getInputFeasibilityResultName(InputFeasibilityResult fr) {
  switch (fr) {
    case InputFeasibilityResult::kInputFeasible:
      return "Feasible";
    case InputFeasibilityResult::kInputIndeterminable:
      return "Indeterminable";
    case InputFeasibilityResult::kInputInfeasibleThrustHigh:
      return "InfeasibleThrustHigh";
    case InputFeasibilityResult::kInputInfeasibleThrustLow:
      return "InfeasibleThrustLow";
    case InputFeasibilityResult::kInputInfeasibleVelocity:
      return "InfeasibleVelocity";
    case InputFeasibilityResult::kInputInfeasibleRollPitchRates:
      return "InfeasibleRollPitchRates";
    case InputFeasibilityResult::kInputInfeasibleYawRates:
      return "InfeasibleYawRates";
    case InputFeasibilityResult::kInputInfeasibleYawAcc:
      return "InfeasibleYawAcc";
  }
  return "Unknown!";
}

InputConstraints::InputConstraints()
    : f_min_(0.0),
      f_max_(std::numeric_limits<double>::max()),
      v_max_(std::numeric_limits<double>::max()),
      omega_xy_max_(std::numeric_limits<double>::max()),
      omega_z_max_(std::numeric_limits<double>::max()),
      omega_z_dot_max_(std::numeric_limits<double>::max()) {}

InputConstraints::InputConstraints(double f_min, double f_max, double v_max,
                                   double omega_xy_max, double omega_z_max,
                                   double omega_z_dot_max) {
setFMin(f_min);
setFMax(f_max);
setVMax(v_max);
setOmegaXYMax(omega_xy_max);
setOmegaZMax(omega_z_max);
setOmegaZDotMax(omega_z_dot_max);
}

void InputConstraints::setFMin(double f_min) {
  f_min = std::abs(f_min);
  f_max_ = f_min > f_max_ ? f_min : f_max_;
  f_min_ = f_min;
}

void InputConstraints::setFMax(double f_max) {
  f_max = std::abs(f_max);
  f_min_ = f_max < f_min_ ? f_max : f_min_;
  f_max_ = f_max;
}

void InputConstraints::setVMax(double v_max) { v_max_ = std::abs(v_max); }
void InputConstraints::setOmegaXYMax(double omega_xy_max) {
  omega_xy_max_ = std::abs(omega_xy_max);
}
void InputConstraints::setOmegaZMax(double omega_z_max) {
  omega_z_max_ = std::abs(omega_z_max);
}
void InputConstraints::setOmegaZDotMax(double omega_z_dot_max) {
  omega_z_dot_max_ = std::abs(omega_z_dot_max);
}
void InputConstraints::setDefaultValues() {
  setFMin(0.5 * mav_msgs::kGravity);
  setFMax(1.5 * mav_msgs::kGravity);
  setVMax(3.0);
  setOmegaXYMax(M_PI / 2.0);
  setOmegaZMax(M_PI / 2.0);
  setOmegaZDotMax(2.0 * M_PI);
}

HalfPlane::HalfPlane(const Eigen::Vector3d& point,
                     const Eigen::Vector3d& normal)
    : point_(point), normal_(normal) {
  CHECK_GT(normal.norm(), 0.0) << "Invalid normal.";
  normal_.normalize();
}

HalfPlane::HalfPlane(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                     const Eigen::Vector3d& c) {
  point_ = a;
  normal_ = (b - a).cross(c - a).normalized();
}

FeasibilityBase::FeasibilityBase()
    : gravity_((Eigen::Vector3d() << 0.0, 0.0, mav_msgs::kGravity).finished()) {
}

FeasibilityBase::FeasibilityBase(const InputConstraints& input_constraints)
    : input_constraints_(input_constraints),
      gravity_((Eigen::Vector3d() << 0.0, 0.0, mav_msgs::kGravity).finished()) {
}

InputFeasibilityResult FeasibilityBase::checkInputFeasibility(
    const Trajectory& trajectory) {
  InputFeasibilityResult result = InputFeasibilityResult::kInputIndeterminable;
  for (const Segment segment : trajectory.segments()) {
    result = checkInputFeasibility(segment);
    if (result != InputFeasibilityResult::kInputFeasible) {
      return result;
    }
  }
  return result;
}

bool FeasibilityBase::checkHalfPlaneFeasibility(const Trajectory& trajectory) {
  for (const Segment segment : trajectory.segments()) {
    if (!checkHalfPlaneFeasibility(segment)) {
      return false;
    }
  }
  return true;
}

}  // namespace mav_trajectory_generation