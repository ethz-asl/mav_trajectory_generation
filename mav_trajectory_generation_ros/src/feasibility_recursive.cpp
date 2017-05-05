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

#include <Eigen/Core>

#include <mav_msgs/default_values.h>
#include <mav_trajectory_generation/motion_defines.h>

namespace mav_trajectory_generation {
FeasibilityRecursive::Settings::Settings() : min_section_time_s(0.05) {}

FeasibilityRecursive::FeasibilityRecursive(const Settings& settings)
    : FeasibilityBase(), settings_(settings) {}
FeasibilityRecursive::FeasibilityRecursive(
    const Settings& settings, const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints), settings_(settings) {}

InputFeasibilityResult FeasibilityRecursive::checkInputFeasibility(
    const Segment& segment) {
  // Check user input.
  if (segment.D() != 3 || segment.D() != 4) {
    return InputFeasibilityResult::kInputIndeterminable;
  }

  // Find roots to determine single axis minima / maxima:
  const Roots roots_acc, roots_jerk, roots_snap;
  roots_acc.resize(3);
  roots_jerk.resize(3);
  roots_snap.resize(3);
  for (size_t i = 0; i < 3; i++) {
    roots_acc[i] = findRootsJenkinsTraub(segment[i].getCoefficients(),
                                         derivative_order::ACCELERATION);
    roots_jerk[i] = findRootsJenkinsTraub(segment[i].getCoefficients(),
                                          derivative_order::JERK);
    roots_snap[i] = findRootsJenkinsTraub(segment[i].getCoefficients(),
                                          derivative_order::SNAP);
    if (roots_acc[i].size() == 0 || roots_jerk[i].size() == 0 ||
        roots_snap[i].size() == 0) {
      // Failed root computation.
      return InputFeasibilityResult::kInputIndeterminable;
    }
  }

  // Recursive test for velocity, acceleration, and roll and pitch rate
  // feasibility.
  double t_1 = 0.0;
  double t_2 = segment.getTime();
  InputFeasibilityResult result = recursiveFeasibility(
      segment, roots_acc, roots_jerk, roots_snap, t_1, t_2);

  if (segment.D() == 4) {
    // Yaw feasibility (assumed independent in the rigid body model).
  }
  return result;
}

InputFeasibilityResult FeasibilityRecursive::recursiveFeasibility(
    const Segment& segment, const Roots& roots_acc, const Roots& roots_jerk,
    const Roots& roots_snap, double t_1, double t_2) const {
  // Evaluate the thrust at the boundaries of the section:
  const double f_t_1 = evaluateThrust(segment, t_1);
  const double f_t_2 = evaluateThrust(segment, t_2);
  if (std::max(f_t_1, f_t_2) > f_max) {
    return InputFeasibilityResult::kInputInfeasibleThrustHigh;
  } else if (std::min(f_t_1, f_t_2) < f_min) {
    return InputFeasibilityResult::kInputInfeasibleThrustLow;
  }

  // Evaluate the velocity at the boundaries of the section:
  const double v_t_1 = segment.evaluate(t_1, derivative_order::VELOCITY).norm();
  const double v_t_2 = segment.evaluate(t_2, derivative_order::VELOCITY).norm();
  if (std::max(v_t_1, v_t_2) > v_max) {
    return InputFeasibilityResult::kInputInfeasibleVelocity;
  }

  // Upper and lower boundaries of thrust, velocity and jerk.
  double f_min_sqr = 0.0;
  double f_max_sqr = 0.0;
  double v_max_sqr = 0.0;
  double j_max_sqr = 0.0;

  for (size_t i = 0; i < 3; i++) {
    // Find the minimum / maximum of each axis.
    double v_min, v_max, a_min, a_max, j_min, j_max;

  }
}

double FeasibilityRecursive::evaluateThrust(const Segment& segment,
                                            double time) const {
  return (segment.evaluate(time, derivative_order::ACCELERATION) +
          Eigen::Vector3d(0.0, 0.0, mav_msgs::kGravity))
      .norm();
}
}  // namespace mav_trajectory_generation
