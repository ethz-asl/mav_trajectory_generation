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

#ifndef MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_
#define MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_

#include "mav_trajectory_generation/segment.h"

namespace mav_trajectory_generation {

class Trajectory {
 public:
  Trajectory();

  Eigen::VectorXd evaluate(double t, int derivative_order) const;

  int D() const { return D_; }
  int N() const { return N_; }
  int K() const { return segments_.size(); }

  void setSegments(const Segment::Vector& segments) {
    CHECK(!segments_.empty());
    segments_ = segments;
    D_ = segments_.front().D();
    N_ = segments_.front().N();
  }

 private:
  int D_;  ///< Number of dimensions.
  int N_;  ///< Number of coefficients.

  // K is number of segments...
  Segment::Vector segments_;
};

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_
