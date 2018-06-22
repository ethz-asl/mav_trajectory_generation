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

#include "mav_trajectory_generation_ros/feasibility_base.h"

#include <cmath>
#include <limits>

#include <Eigen/Geometry>

#include <mav_msgs/default_values.h>

namespace mav_trajectory_generation {

std::string getInputFeasibilityResultName(InputFeasibilityResult fr) {
  switch (fr) {
    case InputFeasibilityResult::kInputFeasible:
      return "Feasible";
    case InputFeasibilityResult::kInputIndeterminable:
      return "Indeterminable";
    case InputFeasibilityResult::kInputInfeasibleThrustHigh:
      return "InfeasibleThrustHigh";
    case InputFeasibilityResult::kInputInfeasibleThrustLow:
      return "InfeasibleThrustLow";
    case InputFeasibilityResult::kInputInfeasibleVelocity:
      return "InfeasibleVelocity";
    case InputFeasibilityResult::kInputInfeasibleRollPitchRates:
      return "InfeasibleRollPitchRates";
    case InputFeasibilityResult::kInputInfeasibleYawRates:
      return "InfeasibleYawRates";
    case InputFeasibilityResult::kInputInfeasibleYawAcc:
      return "InfeasibleYawAcc";
  }
  return "Unknown!";
}

HalfPlane::HalfPlane(const Eigen::Vector3d& point,
                     const Eigen::Vector3d& normal)
    : point(point), normal(normal) {
  CHECK_GT(normal.norm(), 0.0) << "Invalid normal.";
  this->normal.normalize();
}

HalfPlane::HalfPlane(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                     const Eigen::Vector3d& c) {
  CHECK_NE(a, b);
  CHECK_NE(a, c);
  point = a;
  normal = (b - a).cross(c - a).normalized();
}

HalfPlane::Vector HalfPlane::createBoundingBox(
    const Eigen::Vector3d& point, const Eigen::Vector3d& bounding_box_size) {
  HalfPlane::Vector bounding_box;
  bounding_box.reserve(6);

  Eigen::Vector3d bbx_min = point - bounding_box_size / 2.0;
  Eigen::Vector3d bbx_max = point + bounding_box_size / 2.0;

  for (size_t axis = 0; axis < 3; axis++) {
    Eigen::Vector3d normal_min(Eigen::Vector3d::Zero());
    Eigen::Vector3d normal_max(Eigen::Vector3d::Zero());
    normal_min(axis) = 1.0;
    normal_max(axis) = -1.0;
    bounding_box.emplace_back(bbx_min, normal_min);
    bounding_box.emplace_back(bbx_max, normal_max);
  }
  return bounding_box;
}

FeasibilityBase::FeasibilityBase()
    : gravity_((Eigen::Vector3d() << 0.0, 0.0, mav_msgs::kGravity).finished()) {
}

FeasibilityBase::FeasibilityBase(const InputConstraints& input_constraints)
    : input_constraints_(input_constraints),
      gravity_((Eigen::Vector3d() << 0.0, 0.0, mav_msgs::kGravity).finished()) {
}

InputFeasibilityResult FeasibilityBase::checkInputFeasibilityTrajectory(
    const Trajectory& trajectory) const {
  InputFeasibilityResult result = InputFeasibilityResult::kInputIndeterminable;
  for (const Segment segment : trajectory.segments()) {
    result = checkInputFeasibility(segment);
    if (result != InputFeasibilityResult::kInputFeasible) {
      return result;
    }
  }
  return result;
}

bool FeasibilityBase::checkHalfPlaneFeasibility(
    const Trajectory& trajectory) const {
  for (const Segment segment : trajectory.segments()) {
    if (!checkHalfPlaneFeasibility(segment)) {
      return false;
    }
  }
  return true;
}

bool FeasibilityBase::checkHalfPlaneFeasibility(const Segment& segment) const {
  // Check user input.
  if (!(segment.D() == 3 || segment.D() == 4)) {
    LOG(WARNING) << "Feasibility check only implemented for segment dimensions "
                    "3 and 4. Got dimension "
                 << segment.D() << ".";
    return false;
  }

  for (const HalfPlane& half_plane : half_plane_constraints_) {
    // Project the segment position polynomial onto the half plane direction.
    // This polynomial represents the distance of the segment from the origin in
    // half plane direction.
    Polynomial projection(segment.N());
    for (size_t dim = 0; dim < 3; dim++) {
      projection += segment[dim] * half_plane.normal(dim);
    }
    // Find critical times.
    // These are the points in the projected polynomial with zero velocity. The
    // MAV changes directions towards the boundary normal.
    std::vector<double> extrema_candidates;
    projection.computeMinMaxCandidates(0.0, segment.getTime(),
                                       derivative_order::POSITION,
                                       &extrema_candidates);
    // Evaluate position at these critical points.
    for (double t : extrema_candidates) {
      if ((segment.evaluate(t).head(3) - half_plane.point).dot(half_plane.normal) <=
          0.0) {
        // The vector connecting the critical position and the point on the half
        // plane face different / perpendicular direction.
        return false;
      }
    }
  }
  return true;
}

}  // namespace mav_trajectory_generation
