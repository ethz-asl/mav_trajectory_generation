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

#include "mav_trajectory_generation/extremum.h"
#include "mav_trajectory_generation/segment.h"
#include "mav_trajectory_generation/vertex.h"

namespace mav_trajectory_generation {

// Holder class for trajectories of D dimensions, of K segments, and
// polynomial order N-1. (N=12 -> 11th order polynomial, with 12 coefficients).
class Trajectory {
 public:
  Trajectory() : D_(0), N_(0), max_time_(0.0) {}
  ~Trajectory() {}

  bool operator==(const Trajectory& rhs) const;
  inline bool operator!=(const Trajectory& rhs) const {
    return !operator==(rhs);
  }

  int D() const { return D_; }
  int N() const { return N_; }
  int K() const { return segments_.size(); }

  bool empty() const { return segments_.empty(); }
  void clear() {
    segments_.clear();
    D_ = 0;
    N_ = 0;
    max_time_ = 0.0;
  }

  void setSegments(const Segment::Vector& segments) {
    CHECK(!segments.empty());
    // Reset states.
    D_ = segments.front().D();
    N_ = segments.front().N();
    max_time_ = 0.0;
    segments_.clear();

    addSegments(segments);
  }

  void addSegments(const Segment::Vector& segments) {
    for (const Segment& segment : segments) {
      CHECK_EQ(segment.D(), D_);
      CHECK_EQ(segment.N(), N_);
      max_time_ += segment.getTime();
    }
    segments_.insert(segments_.end(), segments.begin(), segments.end());
  }

  void getSegments(Segment::Vector* segments) const {
    CHECK_NOTNULL(segments);
    *segments = segments_;
  }

  const Segment::Vector& segments() const { return segments_; }

  double getMinTime() const { return 0.0; }
  double getMaxTime() const { return max_time_; }
  std::vector<double> getSegmentTimes() const;

  // Functions to create new trajectories by splitting (getting a NEW trajectory
  // with a single dimension) or compositing (create a new trajectory with
  // another trajectory appended).
  Trajectory getTrajectoryWithSingleDimension(int dimension) const;
  bool getTrajectoryWithAppendedDimension(
      const Trajectory& trajectory_to_append, Trajectory* new_trajectory) const;

  // Add trajectories with same dimensions and coefficients to this trajectory.
  bool addTrajectories(const std::vector<Trajectory>& trajectories,
                       Trajectory* merged) const;

  // Offset this trajectory by vector A_r_B.
  bool offsetTrajectory(const Eigen::VectorXd& A_r_B);

  // Evaluate the vertex constraint at time t.
  Vertex getVertexAtTime(double t, int max_derivative_order) const;
  // Evaluate the vertex constraint at start time.
  Vertex getStartVertex(int max_derivative_order) const;
  // Evaluate the vertex constraint at goal time.
  Vertex getGoalVertex(int max_derivative_order) const;
  // Evaluate all underlying vertices.
  bool getVertices(int max_derivative_order_pos, int max_derivative_order_yaw,
                   Vertex::Vector* pos_vertices,
                   Vertex::Vector* yaw_vertices) const;
  bool getVertices(int max_derivative_order,
                   Vertex::Vector* vertices) const;

  // Evaluation functions.
  // Evaluate at a single time, and a single derivative. Return type of
  // dimension D.
  Eigen::VectorXd evaluate(
      double t, int derivative_order = derivative_order::POSITION) const;

  // Evaluates the trajectory in a specified range and derivative.
  // Outputs are a vector of the sampled values (size of VectorXd is D) by
  // time and optionally the actual sampling times.
  void evaluateRange(double t_start, double t_end, double dt,
                     int derivative_order, std::vector<Eigen::VectorXd>* result,
                     std::vector<double>* sampling_times = nullptr) const;

  // Compute the analytic minimum and maximum of magnitude for a given
  // derivative and dimensions, e.g., [0, 1, 2] for position or [3] for yaw.
  // Returns false in case of extremum calculation failure.
  bool computeMinMaxMagnitude(int derivative,
                              const std::vector<int>& dimensions,
                              Extremum* minimum, Extremum* maximum) const;

  // Compute max velocity and max acceleration. Shorthand for the method above.
  bool computeMaxVelocityAndAcceleration(double* v_max, double* a_max) const;

  // This method SCALES the segment times evenly.
  bool scaleSegmentTimes(double scaling);

  // This method SCALES the segment times evenly to ensure that the trajectory
  // is feasible given the provided v_max and a_max. Does not change the shape
  // of the trajectory, and only *increases* segment times.
  bool scaleSegmentTimesToMeetConstraints(double v_max, double a_max);

 private:
  int D_;            // Number of dimensions.
  int N_;            // Number of coefficients.
  double max_time_;  // Time at the end of the trajectory.

  // K is number of segments...
  Segment::Vector segments_;
};

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_TRAJECTORY_H_
