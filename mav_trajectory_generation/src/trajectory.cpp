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

#include "mav_trajectory_generation/trajectory.h"
#include <limits>

// fixes error due to std::iota (has been introduced in c++ standard lately
// and may cause compilation errors depending on compiler)
#if __cplusplus <= 199711L
#include <algorithm>
#else
#include <numeric>
#endif

namespace mav_trajectory_generation {

bool Trajectory::operator==(const Trajectory& rhs) const {
  if (segments_.size() != rhs.segments_.size()) {
    // Different number of segments.
    return false;
  } else {
    for (int i = 0; i < K(); i++) {
      if (segments_ != rhs.segments_) {
        return false;
      }
    }
    return true;
  }
}

Eigen::VectorXd Trajectory::evaluate(double t, int derivative_order) const {
  double accumulated_time = 0.0;

  // Look for the correct segment.
  size_t i = 0;
  for (i = 0; i < segments_.size(); ++i) {
    accumulated_time += segments_[i].getTime();
    // |<--t_accumulated -->|
    // x----------x---------x
    //               ^t_start
    //            |chosen it|
    // in case t_start falls on a vertex, the iterator right of the vertex is
    // chosen, hence accumulated_segment_time_ns > t_start
    if (accumulated_time > t) {
      break;
    }
  }
  if (t > accumulated_time) {
    LOG(ERROR) << "Time out of range of the trajectory!";
    return Eigen::VectorXd::Zero(D(), 1);
  }

  // Make sure we don't go off the end of the segments (can happen if t is
  // equal to trajectory max time).
  if (i >= segments_.size()) {
    i = segments_.size() - 1;
  }
  // Go back to the start of this segment.
  accumulated_time -= segments_[i].getTime();

  return segments_[i].evaluate(t - accumulated_time, derivative_order);
}

void Trajectory::evaluateRange(double t_start, double t_end, double dt,
                               int derivative_order,
                               std::vector<Eigen::VectorXd>* result,
                               std::vector<double>* sampling_times) const {
  const size_t expected_number_of_samples = (t_end - t_start) / dt + 1;

  result->clear();
  result->reserve(expected_number_of_samples);

  if (sampling_times != nullptr) {
    sampling_times->clear();
    sampling_times->reserve(expected_number_of_samples);
  }

  double accumulated_time = 0.0;

  // Look for the correct segment to start.
  size_t i = 0;
  for (i = 0; i < segments_.size(); ++i) {
    accumulated_time += segments_[i].getTime();
    // |<--t_accumulated -->|
    // x----------x---------x
    //               ^t_start
    //            |chosen it|
    // in case t_start falls on a vertex, the iterator right of the vertex is
    // chosen, hence accumulated_segment_time_ns > t_start
    if (accumulated_time > t_start) {
      break;
    }
  }
  if (t_start > accumulated_time) {
    LOG(ERROR) << "Start time out of range of the trajectory!";
    return;
  }

  // Go back to the start of this segment.
  accumulated_time -= segments_[i].getTime();
  double time_in_segment = t_start - accumulated_time;

  // Get all the samples, incrementing the segments as we go.
  while (accumulated_time < t_end) {
    if (time_in_segment > segments_[i].getTime()) {
      time_in_segment = time_in_segment - segments_[i].getTime();
      i++;
      // Make sure we don't access segments that don't exist!
      if (i >= segments_.size()) {
        break;
      }
      continue;
    }

    result->push_back(segments_[i].evaluate(time_in_segment, derivative_order));

    if (sampling_times != nullptr) {
      sampling_times->push_back(accumulated_time);
    }

    time_in_segment += dt;
    accumulated_time += dt;
  }
}

Trajectory Trajectory::getTrajectoryWithSingleDimension(int dimension) const {
  CHECK_LT(dimension, D_);

  // Create a new set of segments with just 1 dimension.
  Segment::Vector segments;
  segments.reserve(segments_.size());

  for (size_t k = 0; k < segments_.size(); ++k) {
    Segment segment(N_, 1);
    segment[0] = (segments_[k])[dimension];
    segments.push_back(segment);
  }

  Trajectory traj;
  traj.setSegments(segments);
  return traj;
}

bool Trajectory::getTrajectoryWithAppendedDimension(
    const Trajectory& trajectory_to_append, Trajectory* new_trajectory) const {
  // Handle the case of one of the trajectories being empty.
  if (N_ == 0 || D_ == 0) {
    *new_trajectory = trajectory_to_append;
    return true;
  }
  if (trajectory_to_append.N() == 0 || trajectory_to_append.D() == 0) {
    *new_trajectory = *this;
    return true;
  }
  CHECK_EQ(static_cast<int>(segments_.size()), trajectory_to_append.K());

  // Create a new set of segments with all of the dimensions.
  Segment::Vector segments;
  segments.reserve(segments_.size());

  for (size_t k = 0; k < segments_.size(); ++k) {
    Segment new_segment(0, 0);
    if (!segments_[k].getSegmentWithAppendedDimension(
            trajectory_to_append.segments()[k], &new_segment)) {
      return false;
    }
    segments.push_back(new_segment);
  }

  new_trajectory->setSegments(segments);
  return true;
}

bool Trajectory::computeMinMaxMagnitude(int derivative,
                                        const std::vector<int>& dimensions,
                                        Extremum* minimum,
                                        Extremum* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  minimum->value = std::numeric_limits<double>::max();
  maximum->value = std::numeric_limits<double>::lowest();

  // For all segments in the trajectory:
  for (size_t segment_idx = 0; segment_idx < segments_.size(); segment_idx++) {
    // Compute candidates.
    std::vector<Extremum> candidates;
    if (!segments_[segment_idx].computeMinMaxMagnitudeCandidates(
            derivative, 0.0, segments_[segment_idx].getTime(), dimensions,
            &candidates)) {
      return false;
    }
    // Evaluate candidates.
    Extremum minimum_candidate, maximum_candidate;
    if (!segments_[segment_idx].selectMinMaxMagnitudeFromCandidates(
            derivative, 0.0, segments_[segment_idx].getTime(), dimensions,
            candidates, &minimum_candidate, &maximum_candidate)) {
      return false;
    }
    // Select minimum / maximum.
    if (minimum_candidate < *minimum) {
      *minimum = minimum_candidate;
      minimum->segment_idx = static_cast<int>(segment_idx);
    }
    if (maximum_candidate > *maximum) {
      *maximum = maximum_candidate;
      maximum->segment_idx = static_cast<int>(segment_idx);
    }
  }
  return true;
}

std::vector<double> Trajectory::getSegmentTimes() const {
  std::vector<double> segment_times(segments_.size());
  for (size_t i = 0; i < segment_times.size(); ++i) {
    segment_times[i] = segments_[i].getTime();
  }
  return segment_times;
}

bool Trajectory::addTrajectories(const std::vector<Trajectory>& trajectories,
                                 Trajectory* merged) const {
  CHECK_NOTNULL(merged);
  merged->clear();
  *merged = *this;

  for (const Trajectory& t : trajectories) {
    // Check dimensions and coefficients.
    // TODO(rikba): Allow different number of coefficients.
    if (t.D() != D_ || t.N() != N_) {
      LOG(WARNING) << "Dimension to append: " << t.D()
                   << " this dimension: " << D_;
      LOG(WARNING) << "Number of coefficients to append: " << t.N()
                   << " this number of coefficients: " << N_;
      return false;
    }
    // Add segments.
    Segment::Vector segments;
    t.getSegments(&segments);
    merged->addSegments(segments);
  }

  return true;
}

bool Trajectory::offsetTrajectory(const Eigen::VectorXd& A_r_B) {
  if (A_r_B.size() < std::min(D_, 3)) {
    LOG(WARNING) << "Offset vector size smaller than trajectory dimension.";
    return false;
  }
  
  for (Segment& s : segments_) {
    // Returns false if dimension check fails at segment level.
    if (!s.offsetSegment(A_r_B)) return false;
  }

  return true;
}

Vertex Trajectory::getVertexAtTime(double t, int max_derivative_order) const {
  Vertex v(D_);
  for (size_t i = 0; i <= max_derivative_order; i++) {
    v.addConstraint(i, evaluate(t, i));
  }
  return v;
}

Vertex Trajectory::getStartVertex(int max_derivative_order) const {
  return getVertexAtTime(0.0, max_derivative_order);
}

Vertex Trajectory::getGoalVertex(int max_derivative_order) const {
  return getVertexAtTime(max_time_, max_derivative_order);
}

bool Trajectory::getVertices(int max_derivative_order_pos,
                             int max_derivative_order_yaw,
                             Vertex::Vector* pos_vertices,
                             Vertex::Vector* yaw_vertices) const {
  CHECK_NOTNULL(pos_vertices);
  CHECK_NOTNULL(yaw_vertices);
  const std::vector<size_t> kPosDimensions = {0, 1, 2};
  const std::vector<size_t> kYawDimensions = {3};
  const int kMaxDerivativeOrder =
      std::max(max_derivative_order_pos, max_derivative_order_yaw);
  pos_vertices->resize(segments_.size() + 1, Vertex(3));
  yaw_vertices->resize(segments_.size() + 1, Vertex(1));

  // Start vertex.
  Vertex temp_vertex(4);
  temp_vertex = getStartVertex(kMaxDerivativeOrder);
  if (!temp_vertex.getSubdimension(kPosDimensions, max_derivative_order_pos,
                                   &pos_vertices->front()))
    return false;
  if (!temp_vertex.getSubdimension(kYawDimensions, max_derivative_order_yaw,
                                   &yaw_vertices->front()))
    return false;

  double t = 0.0;
  for (size_t i = 0; i < segments_.size(); ++i) {
    t += segments_[i].getTime();
    temp_vertex = getVertexAtTime(t, kMaxDerivativeOrder);
    if (!temp_vertex.getSubdimension(kPosDimensions, max_derivative_order_pos,
                                     &(*pos_vertices)[i + 1]))
      return false;
    if (!temp_vertex.getSubdimension(kYawDimensions, max_derivative_order_yaw,
                                     &(*yaw_vertices)[i + 1]))
      return false;
  }
  return true;
}

bool Trajectory::getVertices(int max_derivative_order,
                             Vertex::Vector* vertices) const {
  CHECK_NOTNULL(vertices);
  vertices->resize(segments_.size() + 1, D_);
  vertices->front() = getStartVertex(max_derivative_order);
  
  double t = 0.0;
  for (size_t i = 0; i < segments_.size(); ++i) {
    t += segments_[i].getTime();
    (*vertices)[i + 1] = getVertexAtTime(t, max_derivative_order);
  }
  return true;
}

// Compute max velocity and max acceleration.
bool Trajectory::computeMaxVelocityAndAcceleration(double* v_max,
                                                   double* a_max) const {
  std::vector<int> dimensions(D_);  // Evaluate in whatever dimensions we have.
  std::iota(dimensions.begin(), dimensions.end(), 0);

  Extremum v_min_traj, v_max_traj, a_min_traj, a_max_traj;

  bool success = computeMinMaxMagnitude(
      mav_trajectory_generation::derivative_order::VELOCITY, dimensions,
      &v_min_traj, &v_max_traj);
  success &= computeMinMaxMagnitude(
      mav_trajectory_generation::derivative_order::ACCELERATION, dimensions,
      &a_min_traj, &a_max_traj);

  *v_max = v_max_traj.value;
  *a_max = a_max_traj.value;
  return success;
}

bool Trajectory::scaleSegmentTimes(double scaling) {
  if (scaling < 1.0e-6) return false;

  // Scale the segment times of each segment.
  double new_max_time = 0.0;
  double scaling_inverse = 1.0 / scaling;
  for (size_t i = 0; i < segments_.size(); i++) {
    double new_time = segments_[i].getTime() * scaling;
    for (int d = 0; d < segments_[i].D(); d++) {
      (segments_[i])[d].scalePolynomialInTime(scaling_inverse);
    }
    segments_[i].setTime(new_time);
    new_max_time += new_time;
  }
  max_time_ = new_max_time;

  return true;
}

// This method SCALES the segment times evenly to ensure that the trajectory
// is feasible given the provided v_max and a_max. Does not change the shape
// of the trajectory, and only *increases* segment times.
bool Trajectory::scaleSegmentTimesToMeetConstraints(double v_max,
                                                    double a_max) {
  // In vast majority of cases, this will converge within 1 iteration.
  constexpr size_t kMaxCounter = 20;
  constexpr double kTolerance = 1e-3;

  bool within_range = false;

  for (size_t i = 0; i < kMaxCounter; i++) {
    // From Liu, Sikang, et al. "Planning Dynamically Feasible Trajectories for
    // Quadrotors Using Safe Flight Corridors in 3-D Complex Environments." IEEE
    // Robotics and Automation Letters 2.3 (2017).
    double v_max_actual, a_max_actual;
    computeMaxVelocityAndAcceleration(&v_max_actual, &a_max_actual);

    // Reevaluate constraint/bound violation
    double velocity_violation = v_max_actual / v_max;
    double acceleration_violation = a_max_actual / a_max;

    within_range = velocity_violation <= 1.0 + kTolerance &&
                   acceleration_violation <= 1.0 + kTolerance;
    if (within_range) {
      break;
    }

    double violation_scaling = std::max(
        1.0, std::max(velocity_violation, sqrt(acceleration_violation)));

    // First figure out how to stretch the trajectory in time.
    double violation_scaling_inverse = 1.0 / violation_scaling;

    // Scale the segment times of each segment.
    double new_max_time = 0.0;
    for (size_t i = 0; i < segments_.size(); i++) {
      double new_time = segments_[i].getTime() * violation_scaling;
      for (int d = 0; d < segments_[i].D(); d++) {
        (segments_[i])[d].scalePolynomialInTime(violation_scaling_inverse);
      }
      segments_[i].setTime(new_time);
      new_max_time += new_time;
    }
    max_time_ = new_max_time;
  }
  return within_range;
}

}  // namespace mav_trajectory_generation
