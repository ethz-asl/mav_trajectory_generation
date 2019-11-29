/*
 * Copyright (c) 2016, Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Helen Oleynikova, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Rik Bähnemann, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Marija Popovic, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2019, Maximilian brunner, ASL, ETH Zurich, Switzerland
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

#ifndef MAV_TRAJECTORY_GENERATION_SQUAD_OPTIMIZATION_H_
#define MAV_TRAJECTORY_GENERATION_SQUAD_OPTIMIZATION_H_

#include <glog/logging.h>
#include <Eigen/Sparse>
#include <tuple>

#include "mav_trajectory_generation/extremum.h"
#include "mav_trajectory_generation/motion_defines.h"
#include "mav_trajectory_generation/polynomial.h"
#include "mav_trajectory_generation/segment.h"
#include "mav_trajectory_generation/trajectory.h"
#include "mav_trajectory_generation/vertex.h"

namespace mav_trajectory_generation {

class SquadOptimization {

 public:

  // Constructor
  SquadOptimization();

  // Sets up the optimization problem from a vector of Vertex objects and
  // a vector of times between the vertices.
  // Input: vertices = Vector containing the vertices defining the support
  // points and constraints of the path.
  // Input: segment_times = Vector containing the time between two vertices.
  // Thus, its size is size(vertices) - 1.
  // Input: derivative_to_optimize = Specifies the derivative of which the
  // cost is optimized.
  bool setupFromVertices(
      const Vertex::Vector& vertices, const std::vector<double>& segment_times);

  // TBD: Sets up the optimization problem from a vector of quaternions.
  bool setupFromRotations(const std::vector<Eigen::Quaterniond>& quaternions,
                         const std::vector<double>& times);
  void addToStates(mav_msgs::EigenTrajectoryPoint::Vector* states);

 private:

  bool getInterpolation(const double &t, Eigen::Quaterniond *result);
  Eigen::Quaterniond getQuaternionControlPoint(const Eigen::Quaterniond &q0,
                                            const Eigen::Quaterniond &q1,
                                            const Eigen::Quaterniond &q2);
  Eigen::Quaterniond quaternionExponential(const Eigen::Quaterniond &q);
  Eigen::Quaterniond quaternionLogarithm(const Eigen::Quaterniond &q);
  Eigen::Quaterniond quaternionSum(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2);
  Eigen::Quaterniond quaternionPow(const Eigen::Quaterniond &q, const double &t);
  Eigen::Quaterniond slerp(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2, const double &t);

  // Original vertices containing the constraints.
  Vertex::Vector vertices_;
  std::vector<Eigen::Quaterniond> quaternions_;

  // The actual segments containing the solution.
  Segment::Vector segments_;

  std::vector<double> segment_times_;
  std::vector<double> waypoint_times_;

  size_t n_vertices_;
  size_t n_segments_;

};

}  // namespace mav_trajectory_generation

#include "mav_trajectory_generation/impl/squad_optimization_impl.h"

#endif  // MAV_TRAJECTORY_GENERATION_SQUAD_OPTIMIZATION_H_
