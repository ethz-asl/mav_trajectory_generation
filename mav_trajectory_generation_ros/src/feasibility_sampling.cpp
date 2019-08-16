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
#include <mav_msgs/eigen_mav_msgs.h>

namespace mav_trajectory_generation {
const double kNumNanosecondsPerSecond = 1.0e9;

FeasibilitySampling::Settings::Settings() : sampling_interval_s_(0.01) {}

FeasibilitySampling::FeasibilitySampling(const Settings& settings)
    : FeasibilityBase(), settings_(settings) {}
FeasibilitySampling::FeasibilitySampling(
    const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints) {}
FeasibilitySampling::FeasibilitySampling(
    const Settings& settings, const InputConstraints& input_constraints)
    : FeasibilityBase(input_constraints), settings_(settings) {}

InputFeasibilityResult FeasibilitySampling::checkInputFeasibility(
    const Segment& segment) const {
  // Check user input. Feasilbility sampling only valid for 4DOF flat state trajectories.
  if (!(segment.D() == 3 || segment.D() == 4)) {
    return InputFeasibilityResult::kInputIndeterminable;
  }

  double t = 0.0;
  while (t <= segment.getTime()) {
    // Flat state:
    Eigen::VectorXd position = segment.evaluate(t, derivative_order::POSITION);
    Eigen::VectorXd velocity = segment.evaluate(t, derivative_order::VELOCITY);
    Eigen::VectorXd acc = segment.evaluate(t, derivative_order::ACCELERATION);
    Eigen::VectorXd jerk = segment.evaluate(t, derivative_order::JERK);
    Eigen::VectorXd snap = segment.evaluate(t, derivative_order::SNAP);
    int64_t time_from_start_ns =
        static_cast<int64_t>(t * kNumNanosecondsPerSecond);
    mav_msgs::EigenTrajectoryPoint flat_state;
    flat_state.position_W = position.head<3>();
    flat_state.velocity_W = velocity.head<3>();
    flat_state.acceleration_W = acc.head<3>();
    flat_state.jerk_W = jerk.head<3>();
    flat_state.snap_W = snap.head<3>();
    flat_state.time_from_start_ns = time_from_start_ns;
    if (segment.D() == 4) {
      flat_state.setFromYaw(position(3));
      flat_state.setFromYawRate(velocity(3));
      flat_state.setFromYawAcc(acc(3));
    }

    // Full state:
    mav_msgs::EigenMavState state;
    EigenMavStateFromEigenTrajectoryPoint(flat_state, &state);

    // Feasibility check:
    // Thrust.
    if (input_constraints_.hasConstraint(InputConstraintType::kFMin) ||
        input_constraints_.hasConstraint(InputConstraintType::kFMax)) {
      const double thrust = state.acceleration_B.norm();

      double f_min;
      if (input_constraints_.getConstraint(InputConstraintType::kFMin,
                                           &f_min) &&
          thrust < f_min) {
        return InputFeasibilityResult::kInputInfeasibleThrustLow;
      }

      double f_max;
      if (input_constraints_.getConstraint(InputConstraintType::kFMax,
                                           &f_max) &&
          thrust > f_max) {
        return InputFeasibilityResult::kInputInfeasibleThrustHigh;
      }
    }

    // Velocity.
    double v_max;
    if (input_constraints_.getConstraint(InputConstraintType::kVMax, &v_max) &&
        state.velocity_W.norm() > v_max) {
      return InputFeasibilityResult::kInputInfeasibleVelocity;
    }

    // Evaluate roll/pitch rate and yaw rate assuming independency (rigid body
    // model).
    // Roll/Pitch rates.
    double omega_xy_max;
    if (input_constraints_.getConstraint(InputConstraintType::kOmegaXYMax,
                                         &omega_xy_max) &&
        state.angular_velocity_B.head<2>().norm() > omega_xy_max) {
      return InputFeasibilityResult::kInputInfeasibleRollPitchRates;
    }

    // Yaw rates.
    double omega_z_max;
    if (input_constraints_.getConstraint(InputConstraintType::kOmegaZMax,
                                         &omega_z_max) &&
        std::fabs(state.angular_velocity_B(2)) > omega_z_max) {
      return InputFeasibilityResult::kInputInfeasibleYawRates;
    }

    // Yaw acceleration.
    double omega_z_dot_max;
    if (input_constraints_.getConstraint(InputConstraintType::kOmegaZDotMax,
                                         &omega_z_dot_max) &&
        std::fabs(state.angular_acceleration_B(2)) > omega_z_dot_max) {
      return InputFeasibilityResult::kInputInfeasibleYawAcc;
    }

    // Increment time.
    t += settings_.getSamplingIntervalS();
  }
  return InputFeasibilityResult::kInputFeasible;
}
}  // namespace mav_trajectory_generation
