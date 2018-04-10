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

#include "mav_trajectory_generation/trajectory_sampling.h"

namespace mav_trajectory_generation {

const double kNumNanosecondsPerSecond = 1.e9;

bool sampleTrajectoryAtTime(const Trajectory& trajectory, double sample_time,
                            mav_msgs::EigenTrajectoryPoint* state) {
  CHECK_NOTNULL(state);
  if (sample_time < trajectory.getMinTime() ||
      sample_time > trajectory.getMaxTime()) {
    LOG(ERROR) << "Sample time should be within [" << trajectory.getMinTime()
               << " " << trajectory.getMaxTime() << "] but is " << sample_time;
    return false;
  }

  if (trajectory.D() < 3) {
    LOG(ERROR) << "Dimension has to be 3 or 4, but is " << trajectory.D();
    return false;
  }

  return sampleFlatStateAtTime<Trajectory>(trajectory, sample_time, state);
}

bool sampleTrajectoryInRange(const Trajectory& trajectory, double min_time,
                             double max_time, double sampling_interval,
                             mav_msgs::EigenTrajectoryPointVector* states) {
  CHECK_NOTNULL(states);
  if (min_time < trajectory.getMinTime() ||
      max_time > trajectory.getMaxTime()) {
    LOG(ERROR) << "Sample time should be within [" << trajectory.getMinTime()
               << " " << trajectory.getMaxTime() << "] but is [" << min_time
               << " " << max_time << "]";
    return false;
  }

  if (trajectory.D() < 3) {
    LOG(ERROR) << "Dimension has to be 3 or 4, but is " << trajectory.D();
    return false;
  }

  std::vector<Eigen::VectorXd> position, velocity, acceleration, jerk, snap,
      yaw, yaw_rate;

  trajectory.evaluateRange(min_time, max_time, sampling_interval,
                           derivative_order::POSITION, &position);
  trajectory.evaluateRange(min_time, max_time, sampling_interval,
                           derivative_order::VELOCITY, &velocity);
  trajectory.evaluateRange(min_time, max_time, sampling_interval,
                           derivative_order::ACCELERATION, &acceleration);
  trajectory.evaluateRange(min_time, max_time, sampling_interval,
                           derivative_order::JERK, &jerk);
  trajectory.evaluateRange(min_time, max_time, sampling_interval,
                           derivative_order::SNAP, &snap);

  size_t n_samples = position.size();

  states->resize(n_samples);
  for (size_t i = 0; i < n_samples; ++i) {
    mav_msgs::EigenTrajectoryPoint& state = (*states)[i];

    state.position_W = position[i].head<3>();
    state.velocity_W = velocity[i].head<3>();
    state.acceleration_W = acceleration[i].head<3>();
    state.jerk_W = jerk[i].head<3>();
    state.snap_W = snap[i].head<3>();
    state.time_from_start_ns = static_cast<int64_t>(
        (min_time + sampling_interval * i) * kNumNanosecondsPerSecond);
    if (trajectory.D() > 3) {
      state.setFromYaw(position[i](3));
      state.setFromYawRate(velocity[i](3));
      state.setFromYawAcc(acceleration[i](3));
    }
  }
  return true;
}

bool sampleTrajectoryStartDuration(
    const Trajectory& trajectory, double start_time, double duration,
    double sampling_interval, mav_msgs::EigenTrajectoryPointVector* states) {
  return sampleTrajectoryInRange(trajectory, start_time, start_time + duration,
                                 sampling_interval, states);
}

bool sampleWholeTrajectory(const Trajectory& trajectory,
                           double sampling_interval,
                           mav_msgs::EigenTrajectoryPoint::Vector* states) {
  const double min_time = trajectory.getMinTime();
  const double max_time = trajectory.getMaxTime();

  return sampleTrajectoryInRange(trajectory, min_time, max_time,
                                 sampling_interval, states);
}

bool sampleSegmentAtTime(const Segment& segment, double sample_time,
                         mav_msgs::EigenTrajectoryPoint* state) {
  CHECK_NOTNULL(state);
  if (sample_time < 0.0 || sample_time > segment.getTime()) {
    LOG(ERROR) << "Sample time should be within [" << 0.0 << " "
               << segment.getTime() << "] but is " << sample_time;
    return false;
  }

  return sampleFlatStateAtTime<Segment>(segment, sample_time, state);
}

template <class T>
bool sampleFlatStateAtTime(const T& type, double sample_time,
                           mav_msgs::EigenTrajectoryPoint* state) {
  if (type.D() < 3) {
    LOG(ERROR) << "Dimension has to be 3 or 4, but is " << type.D();
    return false;
  }
  type.evaluate(sample_time, derivative_order::POSITION);

  state->position_W =
      type.evaluate(sample_time, derivative_order::POSITION).head(3);
  state->velocity_W =
      type.evaluate(sample_time, derivative_order::VELOCITY).head(3);
  state->acceleration_W =
      type.evaluate(sample_time, derivative_order::ACCELERATION).head(3);
  state->jerk_W = type.evaluate(sample_time, derivative_order::JERK).head(3);
  state->snap_W = type.evaluate(sample_time, derivative_order::SNAP).head(3);

  if (type.D() > 3) {
    state->setFromYaw(
        (type.evaluate(sample_time, derivative_order::POSITION))(3));
    state->setFromYawRate(
        (type.evaluate(sample_time, derivative_order::VELOCITY))(3));
    state->setFromYawAcc(
        (type.evaluate(sample_time, derivative_order::ACCELERATION))(3));
  }
  state->time_from_start_ns =
      static_cast<int64_t>(sample_time * kNumNanosecondsPerSecond);
  return true;
}

}  // namespace mav_trajectory_generation
