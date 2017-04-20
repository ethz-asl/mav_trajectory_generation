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

#ifndef MAV_TRAJECTORY_GENERATION_FEASIBILITY_BASE_H_
#define MAV_TRAJECTORY_GENERATION_FEASIBILITY_BASE_H_

#include <glog/logging.h>
#include <ros/ros.h>
#include <Eigen/Core>
#include <Eigen/StdVector>

namespace mav_trajectory_generation {
const double kGravity = 9.81;  // [m/s/s]

// Dynamic constraints of the MAV.
struct InputConstraints {
  // Reasonable default constraints.
  InputConstraints();
  InputConstraints(double f_min, double f_max, double v_max,
                   double omega_xy_max, double omega_z_max,
                   double omega_dot_z_max);

  // TODO(rikba): Make class and use setter and getter to check user input,
  // e.g., 0 < f_min < f_max?

  // Minimum input acceleration (mass normalized thrust) [m/s/s].
  double f_min;
  // Maximum input acceleration (mass normalized thrust) [m/s/s].
  double f_max;
  // Maximum velocity [m/s].
  double v_max;
  // Maximum roll/pitch input rate [rad/s].
  double omega_xy_max;
  // Maximum yaw rate [rad/s].
  double omega_z_max;
  // Maximum yaw acc [rad/s^2].
  double omega_dot_z_max;
};

// A half plane is defined through a point and a normal.
class HalfPlane {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef std::vector<HalfPlane, Eigen::aligned_allocator<HalfPlane>> Vector;
  HalfPlane(const Eigen::Vector3d& point, const Eigen::Vector3d& normal);
  // Define the half plane from 3 counter-clockwise points (seen from above).
  HalfPlane(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
            const Eigen::Vector3d& c);

 private:
  Eigen::Vector3d point_;
  Eigen::Vector3d normal_;
};

// A base class for different implementations for dynamic and position
// feasibility checks.
class FeasibilityBase {
 public:
  // Default input constraints, no half plane constraints.
  FeasibilityBase() {};
  // User input constraints, no half plane constraints.
  FeasibilityBase(const InputConstraints& input_constraints);
  // Checks a segment for input feasibility.
  inline virtual bool checkInputFeasibility() {
    ROS_ERROR_STREAM("Input feasibility check not implemented.");
    return false;
  }

  // Checks if a segment stays within a set of half planes.
  inline virtual bool checkHalfPlaneFeasibility() {
    ROS_ERROR_STREAM("Half plane feasibility check not implemented.");
    return false;
  }

  // Input constraints.
  InputConstraints input_constraints_;
  // Half plane constraints, e.g., the ground plane or a box.
  HalfPlane::Vector half_plane_constraints_;
};
}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_FEASIBILITY_BASE_H_
