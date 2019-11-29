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

#include <glog/logging.h>
#include <Eigen/Sparse>
#include <set>
#include <tuple>

// fixes error due to std::iota (has been introduced in c++ standard lately
// and may cause compilation errors depending on compiler)
#if __cplusplus <= 199711L
  #include <algorithm>
#else
  #include <numeric>
#endif

#include "mav_trajectory_generation/convolution.h"

namespace mav_trajectory_generation {

SquadOptimization::SquadOptimization()
    : n_vertices_(0),
      n_segments_(0) {
}

bool SquadOptimization::setupFromVertices(
    const Vertex::Vector& vertices, const std::vector<double>& times) {
  vertices_ = vertices;
  segment_times_ = times;
  return true;
}

bool SquadOptimization::setupFromRotations(
  const std::vector<Eigen::Quaterniond>& quaternions,
                         const std::vector<double>& segment_times) {
  n_vertices_ = quaternions.size();
  n_segments_ = n_vertices_ - 1;
  quaternions_ = quaternions;
  segment_times_ = segment_times;
  waypoint_times_.clear();
  waypoint_times_.push_back(0.0);
  for (int i=1;i<=segment_times.size();i++) {
    waypoint_times_.push_back(waypoint_times_[i-1]+segment_times[i-1]);
  }

  std::cout << "List of reference quaternions: " << std::endl;
  for (int i=0;i<quaternions.size();i++) {
    std::cout << i << " = " << quaternions[i].coeffs().transpose() << std::endl;
  }
}

void SquadOptimization::addToStates(mav_msgs::EigenTrajectoryPoint::Vector* states) {
  for (int i=0;i<states->size();i++) {
    Eigen::Quaterniond q;
    double t = double(states->at(i).time_from_start_ns)/1.e9;
    getInterpolation(t,&q);
    states->at(i).orientation_W_B = q;
    if (i>0) {
      Eigen::Vector3d omega, omega_dot;
      Eigen::Quaterniond q_dot;
      double dt = double(states->at(i).time_from_start_ns - states->at(i-1).time_from_start_ns)*1.e-9;
      q_dot.coeffs() = (states->at(i).orientation_W_B.coeffs() - states->at(i-1).orientation_W_B.coeffs()) / 
                        dt;
      omega = 2.0*( q_dot * states->at(i).orientation_W_B ).vec();
      omega_dot = (omega - states->at(i-1).angular_velocity_W) / dt;
      states->at(i).angular_velocity_W = omega;
      states->at(i).angular_acceleration_W = omega_dot;
    }
  }
}

bool SquadOptimization::getInterpolation(const double &t, Eigen::Quaterniond *result) {
  // t is the index describing where we'd like to get the interpolation at.
  // if(t > n_vertices_-1 || t < 0.0) {
  if (t > waypoint_times_.back()) {
    ROS_WARN("Requested time after last quaternion. t requested is %f", t);
    return false;
  }
  // Get index of the segment that is requested
  int idx = 0;
  while (idx < n_vertices_-1 && t >= waypoint_times_[idx]) {
    idx++;
  }
  idx--;
  // int idx = std::min(int(t),int(n_vertices_)-1);
  Eigen::Quaterniond q0, q1, q2, q3, s1, s2;
  q0 = quaternions_[std::max(0,idx-1)];
  q1 = quaternions_[idx];
  q2 = quaternions_[std::min(int(n_vertices_)-1, idx+1)];
  q3 = quaternions_[std::min(int(n_vertices_)-1, idx+2)];

  double tau = (t - waypoint_times_[idx])/segment_times_[idx];

  s1 = getQuaternionControlPoint(q0,q1,q2);
  s2 = getQuaternionControlPoint(q1,q2,q3);

  Eigen::Quaterniond qa, qb;
  qa = slerp(q1, q2, tau);
  qb = slerp(s1, s2, tau);
  *result = slerp(qa, qb, 2.0*tau*(1.0-tau));
  // std::cout << "Quaternion indices: " << std::max(0,idx-1) << ", " << idx << ", " << std::min(int(n_vertices_)-1, idx+1) << ", " << std::min(int(n_vertices_)-1, idx+2) << std::endl;
  // std::cout << "Result at t = " << t << " is " << (*result).coeffs().transpose() << std::endl;
  // std::cout << "t = " << t << ", idx = " << idx << ", tau = " << tau << "s1 = " << s1.coeffs().transpose() << ", s2 = " << s2.coeffs().transpose() << std::endl;
}

Eigen::Quaterniond SquadOptimization::getQuaternionControlPoint(const Eigen::Quaterniond &q0,
                                            const Eigen::Quaterniond &q1,
                                            const Eigen::Quaterniond &q2) {
  // qa = quaternionLogarithm( q1.inverse() * q2 );
  // qb = quaternionLogarithm( q1.inverse() * q0 );
  Eigen::Quaterniond qa = quaternionSum(quaternionLogarithm( q1.inverse() * q2 ),
                                        quaternionLogarithm( q1.inverse() * q0 ));
  qa.coeffs() = -0.25 * qa.coeffs();
  // q = q1 * quaternionExponential(qa);
  return q1 * quaternionExponential(qa);
}

Eigen::Quaterniond SquadOptimization::quaternionSum(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2) {
  Eigen::Quaterniond q;
  q.coeffs() = q1.coeffs() + q2.coeffs();
  return q;
}

Eigen::Quaterniond SquadOptimization::quaternionPow(const Eigen::Quaterniond &q, const double &t) {
  Eigen::Quaterniond q2;
  q2.coeffs() = t*quaternionLogarithm(q).coeffs();
  return quaternionExponential(q2);
}

Eigen::Quaterniond SquadOptimization::slerp(const Eigen::Quaterniond &q0, const Eigen::Quaterniond &q1, const double &t) {
  // q_interp = tf.quaternion_multiply(quatpow(tf.quaternion_multiply(q1,tf.quaternion_conjugate(q0)),h),q0)
  Eigen::Quaterniond res;
  res = quaternionPow( q1 * q0.conjugate(), t ) * q0;
  return res;
}

Eigen::Quaterniond SquadOptimization::quaternionExponential(const Eigen::Quaterniond &q) {
  // Norm of the vector part:
  double v_abs = q.vec().norm();
  Eigen::Quaterniond result;
  if(v_abs < 0.00001) {
    result = Eigen::Quaterniond(1.0,0.0,0.0,0.0);
  }
  else {
    result.w() = exp(q.w()) * cos(v_abs);
    result.vec() = q.vec().normalized() * sin(v_abs);
  }
  return result;
}

Eigen::Quaterniond SquadOptimization::quaternionLogarithm(const Eigen::Quaterniond &q) {
  Eigen::Quaterniond result;
  double theta = atan2(q.vec().norm(),q.w());
  result.w() = log(q.norm());
  result.vec() = q.vec().normalized() * theta;
  return result;  
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_IMPL_SQUAD_OPTIMIZATION_IMPL_H_
