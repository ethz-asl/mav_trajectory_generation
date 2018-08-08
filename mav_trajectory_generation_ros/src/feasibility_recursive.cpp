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

#include "mav_trajectory_generation_ros/feasibility_recursive.h"

#include <algorithm>
#include <cmath>

#include <Eigen/Core>

#include <mav_trajectory_generation/motion_defines.h>

namespace mav_trajectory_generation {
FeasibilityRecursive::Settings::Settings() : min_section_time_s_(0.05) {}

FeasibilityRecursive::FeasibilityRecursive(const Settings& settings)
    : FeasibilityBase(), settings_(settings) {}
FeasibilityRecursive::FeasibilityRecursive(
    const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints) {}
FeasibilityRecursive::FeasibilityRecursive(
    const Settings& settings, const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints), settings_(settings) {}

InputFeasibilityResult FeasibilityRecursive::checkInputFeasibility(
    const Segment& segment) const {
  // Check user input.
  if (!(segment.D() == 3 || segment.D() == 4)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }

  // Find roots to determine single axis minima / maxima:
  Roots roots_acc, roots_jerk, roots_snap;
  if (input_constraints_.hasConstraint(InputConstraintType::kVMax)) {
    roots_acc.resize(3);
    for (size_t i = 0; i < 3; i++) {
      if (!segment[i].getRoots(derivative_order::ACCELERATION,
              &roots_acc[i])) {
        return InputFeasibilityResult::kInputIndeterminable;
      }
    }
  }

  if (input_constraints_.hasConstraint(InputConstraintType::kFMin) ||
      input_constraints_.hasConstraint(InputConstraintType::kFMax) ||
      input_constraints_.hasConstraint(InputConstraintType::kOmegaXYMax)) {
    roots_jerk.resize(3);
    for (size_t i = 0; i < 3; i++) {
      if (!segment[i].getRoots(derivative_order::JERK,
              &roots_jerk[i])) {
        return InputFeasibilityResult::kInputIndeterminable;
      }
    }
  }

  if (input_constraints_.hasConstraint(InputConstraintType::kOmegaXYMax)) {
    roots_snap.resize(3);
    for (size_t i = 0; i < 3; i++) {
      if (!segment[i].getRoots(derivative_order::SNAP,
              &roots_snap[i])) {
        return InputFeasibilityResult::kInputIndeterminable;
      }
    }
  }

  // Recursive test for velocity, acceleration, and roll and pitch rate
  // feasibility.
  double t_1 = 0.0;
  double t_2 = segment.getTime();
  InputFeasibilityResult result = recursiveFeasibility(
      segment, roots_acc, roots_jerk, roots_snap, t_1, t_2);
  if (result != InputFeasibilityResult::kInputFeasible) {
    return result;
  }

  if (segment.D() == 4) {
    // Yaw feasibility (assumed independent of translation in the rigid body
    // model).
    // Check the single axis minimum / maximum yaw rate:
    double yaw_rate_limit;
    if (input_constraints_.getConstraint(InputConstraintType::kOmegaZMax,
                                         &yaw_rate_limit)) {
      std::pair<double, double> yaw_rate_min, yaw_rate_max;
      if (!segment[3].computeMinMax(t_1, t_2,
                                    derivative_order::ANGULAR_VELOCITY,
                                    &yaw_rate_min, &yaw_rate_max)) {
        return InputFeasibilityResult::kInputIndeterminable;
      }
      if (std::max(std::abs(yaw_rate_min.second),
                   std::abs(yaw_rate_max.second)) > yaw_rate_limit) {
        return InputFeasibilityResult::kInputInfeasibleYawRates;
      }
    }

    // Check the single axis minimum / maximum yaw acceleration:
    double yaw_acc_limit;
    if (input_constraints_.getConstraint(InputConstraintType::kOmegaZDotMax,
                                         &yaw_acc_limit)) {
      std::pair<double, double> yaw_acc_min, yaw_acc_max;
      if (!segment[3].computeMinMax(t_1, t_2,
                                    derivative_order::ANGULAR_ACCELERATION,
                                    &yaw_acc_min, &yaw_acc_max)) {
        return InputFeasibilityResult::kInputIndeterminable;
      }
      if (std::max(std::abs(yaw_acc_min.second), std::abs(yaw_acc_max.second)) >
          yaw_acc_limit) {
        return InputFeasibilityResult::kInputInfeasibleYawAcc;
      }
    }
  }

  // Segment definitely feasible.
  return InputFeasibilityResult::kInputFeasible;
}

InputFeasibilityResult FeasibilityRecursive::recursiveFeasibility(
    const Segment& segment, const Roots& roots_acc, const Roots& roots_jerk,
    const Roots& roots_snap, double t_1, double t_2) const {
  if (t_2 - t_1 < settings_.getMinSectionTimeS()) {
    return InputFeasibilityResult::kInputIndeterminable;
  }
  // Evaluate the thrust at the boundaries of the section:
  if (input_constraints_.hasConstraint(InputConstraintType::kFMin) ||
      input_constraints_.hasConstraint(InputConstraintType::kFMax)) {
    const double f_t_1 = evaluateThrust(segment, t_1);
    const double f_t_2 = evaluateThrust(segment, t_2);
    double f_min_limit, f_max_limit;
    if (input_constraints_.getConstraint(InputConstraintType::kFMin,
                                         &f_min_limit) &&
        std::min(f_t_1, f_t_2) < f_min_limit) {
      return InputFeasibilityResult::kInputInfeasibleThrustLow;
    }
    if (input_constraints_.getConstraint(InputConstraintType::kFMax,
                                         &f_max_limit) &&
        std::max(f_t_1, f_t_2) > f_max_limit) {
      return InputFeasibilityResult::kInputInfeasibleThrustHigh;
    }
  }

  // Evaluate the velocity at the boundaries of the section:
  double v_max_limit;
  if (input_constraints_.getConstraint(InputConstraintType::kVMax,
                                       &v_max_limit)) {
    const double v_t_1 =
        segment.evaluate(t_1, derivative_order::VELOCITY).head<3>().norm();
    const double v_t_2 =
        segment.evaluate(t_2, derivative_order::VELOCITY).head<3>().norm();
    if (std::max(v_t_1, v_t_2) > v_max_limit) {
      return InputFeasibilityResult::kInputInfeasibleVelocity;
    }
  }

  // Upper and lower bounds on thrust, velocity and jerk for this section.
  double f_min_sqr = 0.0;
  double f_max_sqr = 0.0;
  double v_max_sqr = 0.0;
  double j_max_sqr = 0.0;

  // Compute upper bound for velocity.
  if (input_constraints_.hasConstraint(InputConstraintType::kVMax)) {
    for (size_t i = 0; i < 3; i++) {
      std::pair<double, double> v_min, v_max;
      segment[i].selectMinMaxFromRoots(t_1, t_2, derivative_order::VELOCITY,
                                       roots_acc[i], &v_min, &v_max);
      // Definitly infeasible:
      // The velocity on a single axis is higher than the allowed total
      // velocity.
      if (std::max(std::pow(v_min.second, 2), std::pow(v_max.second, 2)) >
          std::pow(v_max_limit, 2)) {
        return InputFeasibilityResult::kInputInfeasibleVelocity;
      }
      // Add single axis extrema to upper bound on squared velocity.
      v_max_sqr +=
          std::pow(std::max(std::abs(v_min.second), std::abs(v_max.second)), 2);
    }
  }

  // Compute upper and lower bound on thrust.
  if (input_constraints_.hasConstraint(InputConstraintType::kFMin) ||
      input_constraints_.hasConstraint(InputConstraintType::kFMax) ||
      input_constraints_.hasConstraint(InputConstraintType::kOmegaXYMax)) {
    for (size_t i = 0; i < 3; i++) {
      std::pair<double, double> a_min, a_max;
      segment[i].selectMinMaxFromRoots(t_1, t_2, derivative_order::ACCELERATION,
                                       roots_jerk[i], &a_min, &a_max);
      // Distance from zero thrust point in this axis.
      const double f_i_min = a_min.second + gravity_[i];
      const double f_i_max = a_max.second + gravity_[i];

      // Definitly infeasible:
      // The thrust on a single axis is higher than the allowed total thrust.
      double f_max_limit;
      if (input_constraints_.getConstraint(InputConstraintType::kFMax,
                                           &f_max_limit) &&
          std::max(std::abs(f_i_min), std::abs(f_i_max)) > f_max_limit) {
        return InputFeasibilityResult::kInputInfeasibleThrustHigh;
      }

      // Add single axis extrema to lower bound on squared acceleration.
      if (f_i_min * f_i_max < 0.0) {
        // Sign of acceleration changes. The minimum squared value is zero.
        f_min_sqr += 0.0;
      } else {
        f_min_sqr +=
            std::pow(std::min(std::abs(f_i_min), std::abs(f_i_max)), 2);
      }

      // Add single axis extrema to upper bound on squared acceleration.
      f_max_sqr += std::pow(std::max(std::abs(f_i_min), std::abs(f_i_max)), 2);
    }
  }

  // Compute upper limit on jerk.
  if (input_constraints_.hasConstraint(InputConstraintType::kOmegaXYMax)) {
    for (size_t i = 0; i < 3; i++) {
      // Find the minimum / maximum of each axis.
      std::pair<double, double> j_min, j_max;
      segment[i].selectMinMaxFromRoots(t_1, t_2, derivative_order::JERK,
                                       roots_snap[i], &j_min, &j_max);
      // Add single axis extrema to upper bound on squared velocity,
      // acceleration, and jerk.
      j_max_sqr +=
          std::pow(std::max(std::abs(j_min.second), std::abs(j_max.second)), 2);
    }
  }

  double f_lower_bound = std::sqrt(f_min_sqr);
  double f_upper_bound = std::sqrt(f_max_sqr);
  double v_upper_bound = std::sqrt(v_max_sqr);
  double omega_xy_upper_bound;
  // Divide-by-zero protection.
  if (f_min_sqr > 1.0e-6) {
    omega_xy_upper_bound = std::sqrt(j_max_sqr / f_min_sqr);
  } else {
    omega_xy_upper_bound = std::numeric_limits<double>::max();
  }

  // Definitely infeasible:
  double f_min_limit;
  if (input_constraints_.getConstraint(InputConstraintType::kFMin,
                                       &f_min_limit) &&
      f_upper_bound < f_min_limit) {
    return InputFeasibilityResult::kInputInfeasibleThrustLow;
  }

  double f_max_limit;
  if (input_constraints_.getConstraint(InputConstraintType::kFMax,
                                       &f_max_limit) &&
      f_lower_bound > f_max_limit) {
    return InputFeasibilityResult::kInputInfeasibleThrustHigh;
  }

  // Possible infeasible (one of the bounds is exceeding the limits):
  double omega_xy_limit;
  if ((input_constraints_.hasConstraint(InputConstraintType::kFMin) &&
       f_lower_bound < f_min_limit) ||
      (input_constraints_.hasConstraint(InputConstraintType::kFMax) &&
       f_upper_bound > f_max_limit) ||
      (input_constraints_.hasConstraint(InputConstraintType::kVMax) &&
       v_upper_bound > v_max_limit) ||
      (input_constraints_.getConstraint(InputConstraintType::kOmegaXYMax,
                                        &omega_xy_limit) &&
       omega_xy_upper_bound > omega_xy_limit)) {
    // Indeterminate. Must check more closely:
    double t_half = (t_1 + t_2) / 2;
    InputFeasibilityResult result_1 = recursiveFeasibility(
        segment, roots_acc, roots_jerk, roots_snap, t_1, t_half);

    if (result_1 == InputFeasibilityResult::kInputFeasible) {
      // Continue with second half.
      return recursiveFeasibility(segment, roots_acc, roots_jerk, roots_snap,
                                  t_half, t_2);
    } else {
      // First half is already infeasible or inderterminate:
      return result_1;
    }
  }
  // Definitely feasible:
  return InputFeasibilityResult::kInputFeasible;
}

double FeasibilityRecursive::evaluateThrust(const Segment& segment,
                                            double time) const {
  return (segment.evaluate(time, derivative_order::ACCELERATION).head<3>() +
          gravity_)
      .norm();
}
}  // namespace mav_trajectory_generation
