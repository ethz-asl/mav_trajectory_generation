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
    for (size_t i = 0; i < D(); i++) {
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

bool Segment::findMinMaxMagnitudeCandidates(
    int derivative, double t_start, double t_end,
    const std::vector<int>& dimensions, std::vector<double>* candidates) const {
  CHECK_NOTNULL(candidates);

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
      if (dim < 0 || dim > D_ - 1) {
        LOG(WARNING) << "Specified dimension " << dim
                     << " is out of bounds [0.." << D_ - 1 << "]." << std::endl;
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
    if (!polynomial_convolved.findMinMaxCandidates(t_start, t_end, -1,
                                                   candidates)) {
      return false;
    }
  } else {
    // For dimension.size() == 1  we can simply evaluate the roots of the
    // derivative.
    if (!polynomials_[dimensions[0]].findMinMaxCandidates(
            t_start, t_end, derivative, candidates)) {
      return false;
    }
  }
  return true;
}

bool Segment::findMinMaxMagnitude(int derivative, double t_start, double t_end,
                                  const std::vector<int>& dimensions,
                                  Extremum* minimum, Extremum* maximum,
                                  std::vector<Extremum>* candidates) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  CHECK_NOTNULL(candidates);
  candidates->clear();

  // Compute extrema candidates.
  std::vector<double> extrema_candidates;
  extrema_candidates.reserve(N_ - 1);
  if (!findMinMaxMagnitudeCandidates(derivative, t_start, t_end, dimensions,
                                     &extrema_candidates)) {
    return false;
  }

  // Evaluate candidate times.
  candidates->resize(extrema_candidates.size());
  for (size_t i = 0; i < extrema_candidates.size(); i++) {
    double value = 0.0;
    for (int dim : dimensions) {
      value += std::pow(
          polynomials_[dim].evaluate(extrema_candidates[i], derivative), 2);
    }
    value = std::sqrt(value);
    (*candidates)[i] = Extremum(extrema_candidates[i], value, 0);
  }

  // Evaluate candidates.
  findMinMaxMagnitude(t_start, t_end, *candidates, minimum, maximum);

  return true;
}

void Segment::findMinMaxMagnitude(double t_start, double t_end,
                                  const std::vector<Extremum>& candidates,
                                  Extremum* minimum, Extremum* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  minimum->value = std::numeric_limits<double>::max();
  maximum->value = std::numeric_limits<double>::lowest();

  for (const Extremum& candidate : candidates) {
    if (candidate.time < t_start || candidate.time > t_end) {
      continue;
    }
    if (candidate > *maximum) {
      *maximum = candidate;
    }
    if (candidate < *minimum) {
      *minimum = candidate;
    }
  }
}

}  // namespace mav_trajectory_generation
