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

#ifndef MAV_TRAJECTORY_GENERATION_SEGMENT_H_
#define MAV_TRAJECTORY_GENERATION_SEGMENT_H_

#include <glog/logging.h>
#include <Eigen/Core>
#include <chrono>
#include <map>
#include <vector>

#include "mav_trajectory_generation/extremum.h"
#include "mav_trajectory_generation/motion_defines.h"
#include "mav_trajectory_generation/polynomial.h"

namespace mav_trajectory_generation {

constexpr double kNumNSecPerSec = 1.0e9;
constexpr double kNumSecPerNsec = 1.0e-9;

// Class holding the properties of parametric segment of a path.
// Time of the segment and one polynomial for each dimension.
//    X------------X---------------X
//  vertex             segment
class Segment {
 public:
  typedef std::vector<Segment> Vector;

  Segment(int N, int D) : time_(0.0), N_(N), D_(D) {
    polynomials_.resize(D_, Polynomial(N_));
  }
  Segment(const Segment& segment) = default;

  bool operator==(const Segment& rhs) const;
  inline bool operator!=(const Segment& rhs) const { return !operator==(rhs); }

  int D() const { return D_; }
  int N() const { return N_; }
  double getTime() const { return time_; }
  uint64_t getTimeNSec() const {
    return static_cast<uint64_t>(kNumNSecPerSec * time_);
  }

  void setTime(double time_sec) { time_ = time_sec; }
  void setTimeNSec(uint64_t time_ns) { time_ = time_ns * kNumSecPerNsec; }

  Polynomial& operator[](size_t idx);

  const Polynomial& operator[](size_t idx) const;

  const Polynomial::Vector& getPolynomialsRef() const { return polynomials_; }

  Eigen::VectorXd evaluate(
      double t, int derivative_order = derivative_order::POSITION) const;

  // Computes the candidates for the minimum and maximum magnitude of a single
  // segment in the specified derivative. In the 1D case, it simply returns the
  // roots of the derivative of the segment-polynomial. For higher dimensions,
  // e.g. 3D, we need to find the extrema of \sqrt{x(t)^2 + y(t)^2 + z(t)^2}
  // where x, y, z are polynomials describing the position or the derivative,
  // specified by Derivative. Taking the derivative yields
  // 2 x \dot{x} + 2 y \dot{y} + 2 z \dot{z}, which needs to be zero at the
  // extrema. The multiplication of two polynomials is a convolution of their
  // coefficients. Re-ordering by their powers and addition yields a polynomial,
  // for which we can compute the roots.
  // Input: derivative = Derivative of position, in which to find the maxima.
  // Input: t_start = Only maxima >= t_start are returned. Usually set to 0.
  // Input: t_stop = Only maxima <= t_stop are returned. Usually set to segment
  // time.
  // Input: dimensions = Vector containing the dimensions that are evaluated.
  // Usually [0, 1, 2] for position, [3] for yaw.
  // Output: candidates = Vector containing the candidate extrema times.
  // Returns whether the computation succeeded -- false means no candidates
  // were found by Jenkins-Traub.
  bool computeMinMaxMagnitudeCandidateTimes(
      int derivative, double t_start, double t_end,
      const std::vector<int>& dimensions,
      std::vector<double>* candidate_times) const;

  // Convenience function. Additionally evaluates the candidate times.
  bool computeMinMaxMagnitudeCandidates(
      int derivative, double t_start, double t_end,
      const std::vector<int>& dimensions,
      std::vector<Extremum>* candidates) const;

  // Convenience function. Evaluates the magnitudes between t_start and t_end
  // for a set of candidates for given dimensions.
  bool selectMinMaxMagnitudeFromCandidates(
      int derivative, double t_start, double t_end,
      const std::vector<int>& dimensions,
      const std::vector<Extremum>& candidates, Extremum* minimum,
      Extremum* maximum) const;

  // Split a segment to get a segment with the specified dimension.
  bool getSegmentWithSingleDimension(int dimension, Segment* new_segment) const;
  // Compose this segment and another segment to a new segment.
  bool getSegmentWithAppendedDimension(const Segment& segment_to_append,
                                       Segment* new_segment) const;

 // Offset this segment by vector A_r_B.
 bool offsetSegment(const Eigen::VectorXd& A_r_B);

 protected:
  Polynomial::Vector polynomials_;
  double time_;

 private:
  int N_;  // Number of coefficients.
  int D_;  // Number of dimensions.
};

// Prints the properties of the segment.
// Polynomial coefficients are printed with increasing powers,
// i.e. c_0 + c_1*t ... c_{N-1} * t^{N-1}
void printSegment(std::ostream& stream, const Segment& s, int derivative);

std::ostream& operator<<(std::ostream& stream, const Segment& s);

std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Segment>& segments);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_SEGMENT_H_
