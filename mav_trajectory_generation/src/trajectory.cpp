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
            0.0, segments_[segment_idx].getTime(), derivative, dimensions,
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

}  // namespace mav_trajectory_generation
