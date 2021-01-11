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

#include "mav_trajectory_generation/segment.h"

#include <cmath>
#include <limits>

namespace mav_trajectory_generation {

bool Segment::operator==(const Segment& rhs) const {
  if (D_ != rhs.D_ || time_ != rhs.time_) {
    return false;
  } else {
    for (int i = 0; i < D(); i++) {
      if (polynomials_[i] != rhs[i]) {
        return false;
      }
    }
  }
  return true;
}

Polynomial& Segment::operator[](size_t idx) {
  CHECK_LT(idx, static_cast<size_t>(D_));
  return polynomials_[idx];
}

const Polynomial& Segment::operator[](size_t idx) const {
  CHECK_LT(idx, static_cast<size_t>(D_));
  return polynomials_[idx];
}

Eigen::VectorXd Segment::evaluate(double t, int derivative) const {
  Eigen::VectorXd result(D_);
  result.setZero();
  for (int d = 0; d < D_; ++d) {
    result[d] = polynomials_[d].evaluate(t, derivative);
  }
  return result;
}

void printSegment(std::ostream& stream, const Segment& s, int derivative) {
  CHECK(derivative >= 0 && derivative < s.N());
  stream << "t: " << s.getTime() << std::endl;
  stream << " coefficients for " << positionDerivativeToString(derivative)
         << ": " << std::endl;
  for (int i = 0; i < s.D(); ++i) {
    stream << "dim " << i << ": " << std::endl;
    stream << s[i].getCoefficients(derivative) << std::endl;
  }
}

std::ostream& operator<<(std::ostream& stream, const Segment& s) {
  printSegment(stream, s, derivative_order::POSITION);
  return stream;
}

std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Segment>& segments) {
  for (const Segment& s : segments) stream << s << std::endl;

  return stream;
}

bool Segment::computeMinMaxMagnitudeCandidateTimes(
    int derivative, double t_start, double t_end,
    const std::vector<int>& dimensions,
    std::vector<double>* candidate_times) const {
  CHECK_NOTNULL(candidate_times);
  candidate_times->clear();
  // Compute magnitude derivative roots.
  if (dimensions.empty()) {
    LOG(WARNING) << "No dimensions specified." << std::endl;
    return false;
  } else if (dimensions.size() > 1) {
    const int n_d = N_ - derivative;
    const int n_dd = n_d - 1;
    const int convolved_coefficients_length =
        Polynomial::getConvolutionLength(n_d, n_dd);
    Eigen::VectorXd convolved_coefficients(convolved_coefficients_length);
    convolved_coefficients.setZero();
    for (int dim : dimensions) {
      if (dim < 0 || dim >= D_) {
        LOG(WARNING) << "Specified dimensions " << dim
                     << " are out of bounds [0.." << D_ - 1 << "]."
                     << std::endl;
        return false;
      }
      // Our coefficients are INCREASING, so when you take the derivative,
      // only the lower powers of t have non-zero coefficients.
      // So we take the head.
      Eigen::VectorXd d =
          polynomials_[dim].getCoefficients(derivative).head(n_d);
      Eigen::VectorXd dd =
          polynomials_[dim].getCoefficients(derivative + 1).head(n_dd);
      convolved_coefficients += Polynomial::convolve(d, dd);
    }
    Polynomial polynomial_convolved(convolved_coefficients);

    // derivative = -1 because the convolved polynomial is the derivative
    // already. We wish to find the minimum and maximum candidates for the
    // integral.
    if (!polynomial_convolved.computeMinMaxCandidates(t_start, t_end, -1,
                                                      candidate_times)) {
      return false;
    }
  } else {
    // For dimension.size() == 1  we can simply evaluate the roots of the
    // derivative.
    if (!polynomials_[dimensions[0]].computeMinMaxCandidates(
            t_start, t_end, derivative, candidate_times)) {
      return false;
    }
  }
  return true;
}

bool Segment::computeMinMaxMagnitudeCandidates(
    int derivative, double t_start, double t_end,
    const std::vector<int>& dimensions,
    std::vector<Extremum>* candidates) const {
  CHECK_NOTNULL(candidates);
  // Find candidate times (roots + start + end).
  std::vector<double> candidate_times;
  computeMinMaxMagnitudeCandidateTimes(derivative, t_start, t_end, dimensions,
                                       &candidate_times);

  // Evaluate candidate times.
  candidates->resize(candidate_times.size());
  for (size_t i = 0; i < candidate_times.size(); i++) {
    double magnitude = 0.0;
    for (int dim : dimensions) {
      magnitude += std::pow(
          polynomials_[dim].evaluate(candidate_times[i], derivative), 2);
    }
    magnitude = std::sqrt(magnitude);
    (*candidates)[i] = Extremum(candidate_times[i], magnitude, 0);
  }

  return true;
}

bool Segment::selectMinMaxMagnitudeFromCandidates(
    int derivative, double t_start, double t_end,
    const std::vector<int>& dimensions, const std::vector<Extremum>& candidates,
    Extremum* minimum, Extremum* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  if (t_start > t_end) {
    LOG(WARNING) << "t_start is greater than t_end.";
    return false;
  }

  minimum->value = std::numeric_limits<double>::max();
  maximum->value = std::numeric_limits<double>::lowest();

  // Evaluate passed candidates.
  for (const Extremum& candidate : candidates) {
    if (candidate.time < t_start || candidate.time > t_end) {
      continue;
    }
    *maximum = std::max(*maximum, candidate);
    *minimum = std::min(*minimum, candidate);
  }

  return true;
}

bool Segment::getSegmentWithSingleDimension(int dimension,
                                            Segment* new_segment) const {
  if (dimension < 0 || dimension >= D_) {
    LOG(WARNING)
        << "You shan't ask for a dimension that does not exist in the segment.";
    return false;
  }

  *new_segment = Segment(N_, 1);
  (*new_segment)[0] = polynomials_[dimension];
  new_segment->setTime(time_);
  return true;
}

bool Segment::getSegmentWithAppendedDimension(const Segment& segment_to_append,
                                              Segment* new_segment) const {
  if (N_ == 0 || D_ == 0) {
    *new_segment = segment_to_append;
    return true;
  }
  if (segment_to_append.N() == 0 || segment_to_append.D() == 0) {
    *new_segment = *this;
    return true;
  }

  // Get common polynomial order.
  const int new_N = std::max(segment_to_append.N(), N_);
  const int new_D = D_ + segment_to_append.D();

  // Create temporary segments to scale polynomials if necessary.
  Segment current_segment = *this;
  Segment segment_to_append_temp = segment_to_append;

  // Scale segment polynomials to the longer segment time.
  const double new_time = std::max(time_, segment_to_append.getTime());
  if (time_ < new_time && new_time > 0.0){
    for (int d = 0; d < D_; d++) {
      current_segment[d].scalePolynomialInTime(time_ / new_time);
    }
  } else if (segment_to_append.getTime() < new_time && new_time > 0.0) {
    for (int d = 0; d < segment_to_append.D(); d++) {
      segment_to_append_temp[d].scalePolynomialInTime(segment_to_append.getTime() / new_time);
    }
  }

  *new_segment = Segment(new_N, new_D);

  if (N_ == segment_to_append.N()) {
    for (int i = 0; i < new_D; i++) {
      if (i < D_) {
        (*new_segment)[i] = current_segment[i];
      } else {
        (*new_segment)[i] = segment_to_append_temp[i - D_];
      }
    }
  } else {
    for (int i = 0; i < new_D; i++) {
      Polynomial polynomial_to_append(new_N);
      if (i < D_) {
        if (!polynomials_[i].getPolynomialWithAppendedCoefficients(
                new_N, &polynomial_to_append)) {
          return false;
        }
      } else {
        if (!segment_to_append[i - D_].getPolynomialWithAppendedCoefficients(
                new_N, &polynomial_to_append)) {
          return false;
        }
      }
      (*new_segment)[i] = polynomial_to_append;
    }
  }

  new_segment->setTime(new_time);
  return true;
}

bool Segment::offsetSegment(const Eigen::VectorXd& A_r_B) {
  if (A_r_B.size() < std::min(D_, 3)) {
    LOG(WARNING) << "Offset vector size smaller than segment dimension.";
    return false;
  }

  // Only translate the first three dimensions.
  for (size_t i = 0; i < std::min(D_, 3); ++i) {
    polynomials_[i].offsetPolynomial(A_r_B(i));
  }

  return true;
}

}  // namespace mav_trajectory_generation
