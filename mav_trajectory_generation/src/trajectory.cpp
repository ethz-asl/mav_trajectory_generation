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

namespace mav_trajectory_generation {

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
  for (; accumulated_time < t_end; accumulated_time += dt) {
    if (time_in_segment > segments_[i].getTime()) {
      time_in_segment = time_in_segment - segments_[i].getTime();
      i++;
      // Make sure we don't access segments that don't exist!
      if (i >= segments_.size()) {
        break;
      }
    }

    result->push_back(segments_[i].evaluate(time_in_segment, derivative_order));

    if (sampling_times != nullptr) {
      sampling_times->push_back(accumulated_time);
    }

    time_in_segment += dt;
  }
}

}  // namespace mav_trajectory_generation
