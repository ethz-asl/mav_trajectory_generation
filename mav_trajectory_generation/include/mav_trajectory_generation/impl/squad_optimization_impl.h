/*
 * Copyright (c) 2016, Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Helen Oleynikova, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Rik Bähnemann, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2016, Marija Popovic, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2019, Maximilian Brunner, ASL, ETH Zurich, Switzerland
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

#ifndef MAV_TRAJECTORY_GENERATION_IMPL_SQUAD_OPTIMIZATION_IMPL_H_
#define MAV_TRAJECTORY_GENERATION_IMPL_SQUAD_OPTIMIZATION_IMPL_H_

// fixes error due to std::iota (has been introduced in c++ standard lately
// and may cause compilation errors depending on compiler)
#if __cplusplus <= 199711L
#include <algorithm>
#else
#include <numeric>
#endif

namespace mav_trajectory_generation {

SquadOptimization::SquadOptimization()
    : n_vertices_(0), n_segments_(0), verbose_(false), use_slerp_(false) {}
SquadOptimization::SquadOptimization(const bool &use_slerp)
    : n_vertices_(0), n_segments_(0), verbose_(false), use_slerp_(use_slerp) {}

bool SquadOptimization::setupFromRotations(const std::vector<Eigen::Quaterniond> &quaternions,
                                           const std::vector<double> &segment_times) {
  if (quaternions.size() != segment_times.size() + 1) {
    std::cout << "Size of times must be one less than quaternions." << std::endl;
    return false;
  }
  n_vertices_ = quaternions.size();
  n_segments_ = n_vertices_ - 1;
  quaternions_ = quaternions;
  segment_times_ = segment_times;
  waypoint_times_.clear();
  waypoint_times_.push_back(0.0);
  for (size_t i = 1; i <= segment_times.size(); i++) {
    waypoint_times_.push_back(waypoint_times_[i - 1] + segment_times[i - 1]);
  }
  if (verbose_) {
    std::cout << "List of reference quaternions and segment times: " << std::endl;
    for (size_t i = 0; i < quaternions.size(); i++) {
      std::cout << i << " = " << quaternions[i].coeffs().transpose();
      if (i < segment_times.size()) {
        std::cout << ", t = " << segment_times[i];
      }
      std::cout << std::endl;
    }
  }
  return true;
}

void SquadOptimization::addToStates(mav_msgs::EigenTrajectoryPoint::Vector *states) const {
  for (size_t i = 0; i < states->size(); i++) {
    Eigen::Quaterniond q;
    double t = static_cast<double>((*states)[i].time_from_start_ns) / 1.e9;
    getInterpolation(t, &q);
    (*states)[i].orientation_W_B = q;
    if (i > 0) {
      Eigen::Vector3d omega, omega_dot;
      Eigen::Quaterniond q_dot;
      double dt = static_cast<double>((*states)[i].time_from_start_ns -
                                      (*states)[i - 1].time_from_start_ns) *
                  1.e-9;
      if (dt > 0) {
        q_dot.coeffs() =
            ((*states)[i].orientation_W_B.coeffs() - (*states)[i - 1].orientation_W_B.coeffs()) /
            dt;
        omega = 2.0 * ((*states)[i].orientation_W_B.conjugate() * q_dot).vec();
        omega_dot = (omega - (*states)[i - 1].angular_velocity_W) / dt;
        (*states)[i].angular_velocity_W = omega;
        (*states)[i].angular_acceleration_W = omega_dot;
      }
    }
  }
  // Set last waypoint's angular velocity and acceleration to zero
  states->back().angular_velocity_W = Eigen::Vector3d::Zero();
  states->back().angular_acceleration_W = Eigen::Vector3d::Zero();
}

bool SquadOptimization::getInterpolation(const double &t, Eigen::Quaterniond *result) const {
  // t is the index describing where we'd like to get the interpolation at.

  if (n_vertices_ == 2) {
    *result = quaternions_[0].slerp(t / segment_times_[0], quaternions_[1]);
    return true;
  }
  if (t > waypoint_times_.back()) {
    ROS_WARN("Requested time after last quaternion. t requested is %f", t);
    return false;
  }
  // Get index of the segment that is requested
  size_t idx = 0;
  while (idx < n_vertices_ - 1 && t >= waypoint_times_[idx]) {
    idx++;
  }
  if (idx > 0) {
    idx--;
  }

  Eigen::Quaterniond q0, q1, q2, q3, s1, s2;
  q0 = quaternions_[std::max(0, static_cast<int>(idx) - 1)];
  q1 = quaternions_[idx];
  q2 = quaternions_[std::min(static_cast<int>(n_vertices_) - 1, static_cast<int>(idx) + 1)];
  q3 = quaternions_[std::min(static_cast<int>(n_vertices_) - 1, static_cast<int>(idx) + 2)];

  // If there's barely any change, don't interpolate
  if (q1.angularDistance(q2) < 0.001) {
    *result = q1;
    return true;
  }

  double tau = (t - waypoint_times_[idx]) / segment_times_[idx];
  if (tau > 1 || tau < 0) {
    ROS_WARN("[SQUAD OPTIMIZATION] tau out of range, this should not happen.");
  }

  if (use_slerp_) {
    *result = q1.slerp(tau, q2);
    return true;
  }

  s1 = getQuaternionControlPoint(q0, q1, q2);
  s2 = getQuaternionControlPoint(q1, q2, q3);

  Eigen::Quaterniond qa, qb;
  qa = slerp(q1, q2, tau);
  qb = slerp(s1, s2, tau);
  *result = slerp(qa, qb, 2.0 * tau * (1.0 - tau));
  return true;
}

Eigen::Quaterniond SquadOptimization::getQuaternionControlPoint(
    const Eigen::Quaterniond &q0, const Eigen::Quaterniond &q1,
    const Eigen::Quaterniond &q2) const {
  Eigen::Quaterniond qa =
      quaternionSum(quaternionLogarithm(q1.inverse() * q2), quaternionLogarithm(q1.inverse() * q0));
  qa.coeffs() = -0.25 * qa.coeffs();
  return q1 * quaternionExponential(qa);
}

Eigen::Quaterniond SquadOptimization::quaternionSum(const Eigen::Quaterniond &q1,
                                                    const Eigen::Quaterniond &q2) const {
  Eigen::Quaterniond q;
  q.coeffs() = q1.coeffs() + q2.coeffs();
  return q;
}

Eigen::Quaterniond SquadOptimization::quaternionPow(const Eigen::Quaterniond &q,
                                                    const double &t) const {
  if (q.vec().norm() < 0.00001) {
    return Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);
  }
  Eigen::Quaterniond q2;
  q2.coeffs() = t * quaternionLogarithm(q).coeffs();
  return quaternionExponential(q2);
}

Eigen::Quaterniond SquadOptimization::slerp(const Eigen::Quaterniond &q0,
                                            const Eigen::Quaterniond &q1, const double &t) const {
  Eigen::Quaterniond res;
  if (t > 1 || t < 0) {
    ROS_WARN("[SQUAD OPTIMIZATION] t out of range, this should not happen.");
  }
  if (t <= 0) {
    return q0;
  } else if (t >= 1.0) {
    return q1;
  }
  res = quaternionPow(q1 * q0.conjugate(), t) * q0;
  return res;
}

Eigen::Quaterniond SquadOptimization::quaternionExponential(const Eigen::Quaterniond &q) const {
  // Norm of the vector part:
  double v_abs = q.vec().norm();
  Eigen::Quaterniond result;
  if (v_abs < 0.00001) {
    result = Eigen::Quaterniond(1.0, 0.0, 0.0, 0.0);
  } else {
    result.w() = exp(q.w()) * cos(v_abs);
    result.vec() = q.vec().normalized() * sin(v_abs);
  }
  return result;
}

Eigen::Quaterniond SquadOptimization::quaternionLogarithm(const Eigen::Quaterniond &q) const {
  Eigen::Quaterniond result;
  double theta = atan2(q.vec().norm(), q.w());
  result.w() = log(q.norm());
  if (q.vec().norm() < 1e-5) {
    result.vec() = Eigen::Vector3d(0, 0, 0);
  } else {
    result.vec() = q.vec().normalized() * theta;
  }
  return result;
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_IMPL_SQUAD_OPTIMIZATION_IMPL_H_
