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

#ifndef MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_BASE_H_
#define MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_BASE_H_

#include <glog/logging.h>
#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/StdVector>

#include <mav_trajectory_generation/trajectory.h>

#include "input_constraints.h"

namespace mav_trajectory_generation {
enum InputFeasibilityResult {
  kInputFeasible = 0,    // The trajectory is input feasible.
  kInputIndeterminable,  // Cannot determine whether the trajectory is feasible
                         // with respect to the inputs.
  kInputInfeasibleThrustHigh,      // The trajectory is infeasible, failed max.
                                   // thrust test first.
  kInputInfeasibleThrustLow,       // The trajectory is infeasible, failed min.
                                   // thrust test first.
  kInputInfeasibleVelocity,        // The Trajectory is infeasible, failed max.
                                   // velocity test first.
  kInputInfeasibleRollPitchRates,  // The trajectory is infeasible, failed max.
                                   // roll/pitch rates test first.
  kInputInfeasibleYawRates,  // The trajectory is infeasible, faild max. yaw
                             // rates test first.
  kInputInfeasibleYawAcc,    // The trajectory is infeasible, failed max. yaw
                             // acceleration test first.
};

// Human readable InputFeasibilityResult.
std::string getInputFeasibilityResultName(InputFeasibilityResult fr);

// A half plane is defined through a point and a normal.
class HalfPlane {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef std::vector<HalfPlane, Eigen::aligned_allocator<HalfPlane>> Vector;
  HalfPlane(const Eigen::Vector3d& point, const Eigen::Vector3d& normal);
  // Define the half plane from 3 counter-clockwise points (seen from above).
  HalfPlane(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
            const Eigen::Vector3d& c);

  // Create 6 half planes that form a box with inward-facing normals.
  // point = center of the box, bounding_box_size = x, y, z edge length.
  static HalfPlane::Vector createBoundingBox(
      const Eigen::Vector3d& point, const Eigen::Vector3d& bounding_box_size);
  Eigen::Vector3d point;
  Eigen::Vector3d normal;
};

// A base class for different implementations for dynamic and position
// feasibility checks.
class FeasibilityBase {
 public:
  // Default input constraints, no half plane constraints.
  FeasibilityBase();
  // User input constraints, no half plane constraints.
  FeasibilityBase(const InputConstraints& input_constraints);

  // Checks a trajectory for input feasibility.
  InputFeasibilityResult checkInputFeasibilityTrajectory(
      const Trajectory& trajectory) const;
  // Checks a segment for input feasibility.
  inline virtual InputFeasibilityResult checkInputFeasibility(
      const Segment& segment) const {
    ROS_ERROR_STREAM("Input feasibility check not implemented.");
    return InputFeasibilityResult::kInputIndeterminable;
  }
  inline InputConstraints getInputConstraints() const {
    return input_constraints_;
  }

  // Checks if a trajectory stays within a set of half planes.
  bool checkHalfPlaneFeasibility(const Trajectory& trajectory) const;
  // Checks if a segment stays within a set of half planes.
  // This check computes the extrema for each axis and checks whether these lie
  // in the positive half space as in
  // https://github.com/markwmuller/RapidQuadrocopterTrajectories/blob/master/C%2B%2B/RapidTrajectoryGenerator.cpp#L149
  bool checkHalfPlaneFeasibility(const Segment& segment) const;

  // Input constraints.
  InputConstraints input_constraints_;
  // Half plane constraints, e.g., the ground plane or a box.
  HalfPlane::Vector half_plane_constraints_;
  Eigen::Vector3d gravity_;
};
}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_BASE_H_
