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

#ifndef MAV_TRAJECTORY_GENERATION_TRAJECTORY_SAMPLING_H_
#define MAV_TRAJECTORY_GENERATION_TRAJECTORY_SAMPLING_H_

#include <mav_msgs/eigen_mav_msgs.h>
#include "mav_trajectory_generation/trajectory.h"

namespace mav_trajectory_generation {

// All of these functions sample a trajectory at a time, or a range of times,
// into an EigenTrajectoryPoint (or vector). These support 3D or 4D
// trajectories. If the trajectories are 4D, the 4th dimension is assumed to
// be yaw.
// If no yaw is set, then it is simply left at its current value (0 by default).

bool sampleTrajectoryAtTime(const Trajectory& trajectory, double sample_time,
                            mav_msgs::EigenTrajectoryPoint* state);

bool sampleTrajectoryInRange(const Trajectory& trajectory, double min_time,
                             double max_time, double sampling_interval,
                             mav_msgs::EigenTrajectoryPointVector* states);

bool sampleTrajectoryStartDuration(
    const Trajectory& trajectory, double start_time, double duration,
    double sampling_interval, mav_msgs::EigenTrajectoryPointVector* states);

bool sampleWholeTrajectory(const Trajectory& trajectory,
                           double sampling_interval,
                           mav_msgs::EigenTrajectoryPoint::Vector* states);

bool sampleSegmentAtTime(const Segment& segment, double sample_time,
                         mav_msgs::EigenTrajectoryPoint* state);

template<class T>
bool sampleFlatStateAtTime(const T& type, double sample_time,
                           mav_msgs::EigenTrajectoryPoint* state);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_TRAJECTORY_SAMPLING_H_
