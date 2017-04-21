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

#include "mav_trajectory_generation_ros/feasibility_sampling.h"

#include <mav_msgs/conversions.h>

#include "mav_trajectory_generation_ros/trajectory_sampling.h"

namespace mav_trajectory_generation {
FeasibilitySampling::Settings::Settings() : sampling_interval(0.01) {}

FeasibilitySampling::FeasibilitySampling(const Settings& settings)
    : FeasibilityBase(), settings_(settings) {}
FeasibilitySampling::FeasibilitySampling(
    const Settings& settings, const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints), settings_(settings) {}

InputFeasibilityResult FeasibilitySampling::checkInputFeasibility(
    const Trajectory& trajectory) {
  if (!sampleTrajectory(trajectory)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }

  for (size_t i = 0; i < states_.size(); i++) {
    double thrust = states_[i].acceleration_B.norm();
    double v = states_[i].velocity_W.norm();
    double omega_xy = states_[i].angular_velocity_B.block<2, 1>(0, 0).norm();
    double omega_z = std::fabs(states_[i].angular_velocity_B(2));
    if (thrust < input_constraints_.f_min) {
      return InputFeasibilityResult::kInputInfeasibleThrustLow;
    } else if (thrust > input_constraints_.f_max) {
      return InputFeasibilityResult::kInputInfeasibleThrustHigh;
    } else if (v > input_constraints_.v_max) {
      return InputFeasibilityResult::kInputInfeasibleVelocity;
    } else if (omega_xy > input_constraints_.omega_xy_max) {
      return InputFeasibilityResult::kInputInfeasibleRollPitchRates;
    } else if (omega_z > input_constraints_.omega_z_max) {
      return InputFeasibilityResult::kInputInfeasibleYawRates;
    }
  }
return InputFeasibilityResult::kInputFeasible;
}

bool FeasibilitySampling::checkHalfPlaneFeasibility(
    const Trajectory& trajectory) {
  if (!sampleTrajectory(trajectory)) {
    return false;
  }
  return true;
}

bool FeasibilitySampling::sampleTrajectory(const Trajectory& trajectory) {
  // Check if already sampled.
  if (trajectory == trajectory_) {
    return true;
  }
  // Try to sample trajectory.
  mav_msgs::EigenTrajectoryPointVector flat_states;
  bool success = sampleWholeTrajectory(trajectory, settings_.sampling_interval,
                                       &flat_states);
  // Recover full state according to Mellinger.
  states_.resize(flat_states.size());
  for (size_t i = 0; i < flat_states.size(); i++) {
    EigenMavStateFromEigenTrajectoryPoint(flat_states[i], &states_[i]);
  }
  // Save trajectory if successful.
  if (!success) {
    return false;
  } else {
    trajectory_ = trajectory;
    return true;
  }
}
}  // namespace mav_trajectory_generation
