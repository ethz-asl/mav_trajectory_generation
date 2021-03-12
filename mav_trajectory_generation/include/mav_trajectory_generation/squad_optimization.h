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

#include <mav_msgs/common.h>
#include <mav_msgs/conversions.h>

#include <glog/logging.h>
#include <Eigen/Core>

namespace mav_trajectory_generation {

class SquadOptimization {
 public:
  SquadOptimization();
  SquadOptimization(const bool &use_slerp);

  bool setupFromRotations(const std::vector<Eigen::Quaterniond> &quaternions,
                          const std::vector<double> &times);
  bool addToStates(mav_msgs::EigenTrajectoryPoint::Vector &states) const;
  bool addToStatesSmoothed(mav_msgs::EigenTrajectoryPoint::Vector &states,
                           const Eigen::Vector3d &inertia, const Eigen::Vector3d &kp,
                           const Eigen::Vector3d &kd) const;
  bool smoothAttitude(mav_msgs::EigenTrajectoryPoint::Vector &states) const;

 private:
  bool getInterpolation(const double &t, Eigen::Quaterniond *result) const;
  bool addVelAcc(mav_msgs::EigenTrajectoryPoint::Vector &states) const;
  Eigen::Quaterniond getQuaternionControlPoint(const Eigen::Quaterniond &q0,
                                               const Eigen::Quaterniond &q1,
                                               const Eigen::Quaterniond &q2) const;
  Eigen::Quaterniond quaternionExponential(const Eigen::Quaterniond &q) const;
  Eigen::Quaterniond quaternionLogarithm(const Eigen::Quaterniond &q) const;
  Eigen::Quaterniond quaternionSum(const Eigen::Quaterniond &q1,
                                   const Eigen::Quaterniond &q2) const;
  Eigen::Quaterniond quaternionPow(const Eigen::Quaterniond &q, const double &t) const;
  Eigen::Quaterniond slerp(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2,
                           const double &t) const;
  inline Eigen::Matrix<double, 4, 1> quaternionDiff(const Eigen::Quaterniond &q1,
                                                    const Eigen::Quaterniond &q2) const;
  Eigen::Vector3d getOmega(const Eigen::Quaterniond &q1,
                                            const Eigen::Quaterniond &q2,
                                            const Eigen::Quaterniond &qm, const double &dt) const;
  bool addAcc(mav_msgs::EigenTrajectoryPoint::Vector &states) const;

  // Original vertices containing the constraints.
  std::vector<Eigen::Quaterniond> quaternions_;

  std::vector<double> segment_times_;
  std::vector<double> waypoint_times_;

  size_t n_vertices_;
  size_t n_segments_;
  bool verbose_;

  bool use_slerp_;
};

}  // namespace mav_trajectory_generation

#include "mav_trajectory_generation/impl/squad_optimization_impl.h"

#endif  // MAV_TRAJECTORY_GENERATION_SQUAD_OPTIMIZATION_H_
