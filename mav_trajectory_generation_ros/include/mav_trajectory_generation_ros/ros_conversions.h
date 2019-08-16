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

#ifndef MAV_TRAJECTORY_GENERATION_ROS_ROS_CONVERSIONS_H_
#define MAV_TRAJECTORY_GENERATION_ROS_ROS_CONVERSIONS_H_

#include <mav_planning_msgs/PolynomialTrajectory.h>
#include <mav_planning_msgs/PolynomialTrajectory4D.h>
#include <mav_planning_msgs/conversions.h>

#include <mav_trajectory_generation/trajectory.h>

namespace mav_trajectory_generation {

// Converts a trajectory into a ROS polynomial trajectory msg.
bool trajectoryToPolynomialTrajectoryMsg(
    const Trajectory& trajectory,
    mav_planning_msgs::PolynomialTrajectory* msg);
    
bool trajectoryToPolynomialTrajectoryMsg(
    const Trajectory& trajectory,
    mav_planning_msgs::PolynomialTrajectory4D* msg);

// Converts a ROS polynomial trajectory msg into a Trajectory.
bool polynomialTrajectoryMsgToTrajectory(
    const mav_planning_msgs::PolynomialTrajectory& msg,
    Trajectory* trajectory);
    
bool polynomialTrajectoryMsgToTrajectory(
    const mav_planning_msgs::PolynomialTrajectory4D& msg,
    Trajectory* trajectory);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_ROS_ROS_CONVERSIONS_H_
