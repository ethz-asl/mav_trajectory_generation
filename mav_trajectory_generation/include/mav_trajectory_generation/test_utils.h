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

#ifndef MAV_TRAJECTORY_GENERATION_TEST_UTILS_H_
#define MAV_TRAJECTORY_GENERATION_TEST_UTILS_H_

#include <random>
#include <Eigen/Eigen>

#include "mav_trajectory_generation/trajectory.h"

namespace mav_trajectory_generation {

inline double createRandomDouble(double min, double max) {
  return (max - min) * (static_cast<double>(std::rand()) /
                        static_cast<double>(RAND_MAX)) +
         min;
}

template <class T1, class T2>
bool checkMatrices(const Eigen::MatrixBase<T1>& m1,
                   const Eigen::MatrixBase<T2>& m2, double tol) {
  return (m1 - m2).cwiseAbs().maxCoeff() < tol;
}

double getMaximumMagnitude(const Trajectory& trajectory, size_t derivative,
                           double dt = 0.01) {
  double maximum = -1e9;

  for (double ts = 0; ts < trajectory.getMaxTime(); ts += dt) {
    double current_value = trajectory.evaluate(ts, derivative).norm();
    if (current_value > maximum) {
      maximum = current_value;
    }
  }
  return maximum;
}

double computeCostNumeric(const Trajectory& trajectory, size_t derivative,
                          double dt = 0.001) {
  double cost = 0;

  for (double ts = 0; ts < trajectory.getMaxTime(); ts += dt) {
    cost += trajectory.evaluate(ts, derivative).squaredNorm() * dt;
  }
  return cost;
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_TEST_UTILS_H_
