#include "mav_trajectory_generation/segment.h"

namespace mav_trajectory_generation {

Polynomial& Segment::operator[](size_t idx) {
  CHECK_LT(idx, dimension_);
  return polynomials_[idx];
}

const Polynomial& Segment::operator[](size_t idx) const {
  CHECK_LT(idx, dimension_);
  return polynomials_[idx];
}

Eigen::VectorXd Segment::evaluate(double t, int derivative) const {
  Eigen::VectorXd result(dimension_);
  for (size_t d = 0; d < dimension_; ++d)
    result[d] = polynomials_[d].evaluate(t, derivative);

  return result;
}

void printSegment(std::ostream& stream, const Segment& s, int derivative) {
  CHECK(derivative >= 0 && derivative < s.N);
  stream << "t: " << s.getTime() << std::endl;
  stream << " coefficients for " << positionDerivativeToString(derivative)
         << ": " << std::endl;
  for (size_t i = 0; i < s.getDimension(); ++i)
    stream << s[i].getCoefficients(derivative) << std::endl;
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
