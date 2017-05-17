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

#include "mav_trajectory_generation_ros/feasibility_analytic.h"

#include <limits>

#include <mav_trajectory_generation/extremum.h>
std::vector<int> kPosDim = {0, 1, 2};

namespace mav_trajectory_generation {
FeasibilityAnalytic::Settings::Settings() : min_section_time_s_(0.05) {}

FeasibilityAnalytic::FeasibilityAnalytic(const Settings& settings)
    : FeasibilityBase(), settings_(settings) {}
FeasibilityAnalytic::FeasibilityAnalytic(
    const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints) {}
FeasibilityAnalytic::FeasibilityAnalytic(
    const Settings& settings, const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints), settings_(settings) {}

InputFeasibilityResult FeasibilityAnalytic::checkInputFeasibility(
    const Segment& segment) {
  // Check user input.
  if (!(segment.D() == 3 || segment.D() == 4)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }
  // Check constraints.
  // Thrust:
  std::vector<Extremum> thrust_candidates;
  Segment thrust_segment(segment.N() - 2, kPosDim.size());
  InputFeasibilityResult thrust_result =
      analyticThrustFeasibility(&segment, &thrust_candidates, &thrust_segment);
  if (thrust_result != InputFeasibilityResult::kInputFeasible) {
    return thrust_result;
  }

  // Velocity:
  std::vector<Extremum> velocity_candidates;
  if (!segment.findMinMaxMagnitudeCandidates(derivative_order::VELOCITY,
                                             kPosDim, &velocity_candidates)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }
  const double v_max =
      std::max_element(velocity_candidates.begin(), velocity_candidates.end())
          ->value;
  if (v_max > input_constraints_.getVMax()) {
    return InputFeasibilityResult::kInputInfeasibleVelocity;
  }

  // Roll / Pitch rates using recursive test:
  std::vector<Extremum> jerk_candidates;
  if (!segment.findMinMaxMagnitudeCandidates(derivative_order::JERK, kPosDim,
                                             &jerk_candidates)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }
  InputFeasibilityResult omega_xy_result =
      recursiveRollPitchFeasibility(segment, thrust_segment, thrust_candidates,
                                    jerk_candidates, 0.0, segment.getTime());
  if (omega_xy_result != InputFeasibilityResult::kInputFeasible) {
    return omega_xy_result;
  }

  // Yaw feasibility (assumed independent of translation in the rigid body
  // model)
  if (segment.D() == 4) {
    // Check the single axis minimum / maximum yaw rate:
    std::pair<double, double> yaw_rate_min, yaw_rate_max;
    if (!segment[3].findMinMax(t_1, t_2, derivative_order::ANGULAR_VELOCITY,
                               &yaw_rate_min, &yaw_rate_max)) {
      return InputFeasibilityResult::kInputIndeterminable;
    }
    if (std::max(std::abs(yaw_rate_min.second), std::abs(yaw_rate_max.second)) >
        input_constraints_.getOmegaZMax()) {
      return InputFeasibilityResult::kInputInfeasibleYawRates;
    }
    // Check the single axis minimum / maximum yaw acceleration:
    std::pair<double, double> yaw_acc_min, yaw_acc_max;
    if (!segment[3].findMinMax(t_1, t_2, derivative_order::ANGULAR_ACCELERATION,
                               &yaw_acc_min, &yaw_acc_max)) {
      return InputFeasibilityResult::kInputIndeterminable;
    }
    if (std::max(std::abs(yaw_acc_min.second), std::abs(yaw_acc_max.second)) >
        input_constraints_.getOmegaZDotMax()) {
      return InputFeasibilityResult::kInputInfeasibleYawAcc;
    }
  }

  return InputFeasibilityResult::kInputFeasible;
}

InputFeasibilityResult analyticThrustFeasibility(
    const Segment& segment, std::vector<Extremum>* thrust_candidates,
    Segment* thrust_segment) const {
  // Create thrust segment, f = ddx + g.
  thrust_segment->setTime(segment.getTime());
  for (size_t i = 0; i < thrust_segment->size(); i++) {
    Eigen::VectorXd thrust_coeffs =
        segment[i]
            .getCoefficients(derivative_order::ACCELERATION)
            .head(thrust_segment->N());
    thrust_coeffs(0) += gravity_[i];
    Polynomial p_thrust(thrust_coeffs);
    (*thrust_segment)[i] = p_thrust;
  }

  // Compute the thrust magnitude extrema candidates.
  if (!thrust_segment.findMinMaxMagnitudeCandidates(
          0, 0.0, thrust_segment.getTime(), kPosDim, thrust_candidates)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }
  // Evaluate the candidates.
  const double f_min =
      std::min_element(thrust_candidates->begin(), thrust_candidates->end())
          ->value;
  const double f_max =
      std::max_element(thrust_candidates->begin(), thrust_candidates->end())
          ->value;
  if (f_max > input_constraints_.getFMax()) {
    return InputFeasibilityResult::kInputInfeasibleThrustHigh;
  }
  if (f_min < input_constraints_.getFMin()) {
    return InputFeasibilityResult::kInputInfeasibleThrustLow;
  }

  return InputFeasibilityResult::kInputFeasible;
}

InputFeasibilityResult FeasibilityAnalytic::recursiveRollPitchFeasibility(
    const Segment& pos_segment, const Segment& thrust_segment,
    const std::vector<Extremum>& thrust_candidates,
    const std::vector<Extremum>& jerk_candidates, double t_1,
    double t_2) const {
  if (t_2 - t_1 < settings_.getMinSectionTimeS()) {
    return InputFeasibilityResult::kInputIndeterminable;
  }

  // Evaluate minimum thrust and maximum jerk of this section.
  Extremum f_min, f_max, j_min, j_max;
  thrust_segment.findMinMaxMagnitude(t_1, t_2, 0, kPosDim, jerk_candidates,
                                     &f_min, &f_max);
  pos_segment.findMinMaxMagnitude(t_1, t_2, derivative_order::JERK, kPosDim,
                                  jerk_candidates, &j_min, &j_max);

  // Upper bound on angular rates according to Müller [1].
  double omega_xy_upper_bound;
  // Divide-by-zero protection.
  if (f_min > 1.0e-6) {
    omega_xy_upper_bound = std::sqrt(j_max / f_min);
  } else {
    omega_xy_upper_bound = std::numeric_limits<double>::max();
  }

  // Possible infeasible.
  if (omega_xy_upper_bound > input_constraints_.getOmegaXYMax()) {
    // Indeterminate. Must check more closely:
    double t_half = (t_1 + t_2) / 2;
    InputFeasibilityResult result_1 = recursiveRollPitchFeasibility(
        pos_segment, thrust_segment, thrust_candidates, jerk_candidates, t_1,
        t_half);

    if (result_1 == InputFeasibilityResult::kInputFeasible) {
      // Continue with second half.
      return recursiveFeasibility(pos_segment, thrust_segment,
                                  thrust_candidates, jerk_candidates, t_half,
                                  t_2);
    } else {
      // First half is already infeasible or inderterminate:
      return result_1;
    }
  } else {
    // Definitely feasible:
    return InputFeasibilityResult::kInputFeasible;
  }
}
}  // namespace mav_trajectory_generation
