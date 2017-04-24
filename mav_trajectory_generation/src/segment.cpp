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

namespace mav_trajectory_generation {

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

}  // namespace mav_trajectory_generation
