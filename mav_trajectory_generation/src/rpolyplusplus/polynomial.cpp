// Copyright (C) 2015 Chris Sweeney. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of Chris Sweeney nor the names of its contributors may
//       be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Please contact the author of this library if you have any questions.
// Author: Chris Sweeney (cmsweeney@cs.ucsb.edu)

#include "mav_trajectory_generation/rpolyplusplus/polynomial.h"

#include <Eigen/Core>
#include <cmath>

namespace rpoly_plus_plus {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXcd;

// Remove leading terms with zero coefficients.
VectorXd RemoveLeadingZeros(const VectorXd& polynomial_in) {
  int i = 0;
  while (i < (polynomial_in.size() - 1) && polynomial_in(i) == 0) {
    ++i;
  }
  return polynomial_in.tail(polynomial_in.size() - i);
}

VectorXd DifferentiatePolynomial(const VectorXd& polynomial) {
  const int degree = polynomial.rows() - 1;

  // Degree zero polynomials are constants, and their derivative does
  // not result in a smaller degree polynomial, just a degree zero
  // polynomial with value zero.
  if (degree == 0) {
    return VectorXd::Zero(1);
  }

  VectorXd derivative(degree);
  for (int i = 0; i < degree; ++i) {
    derivative(i) = (degree - i) * polynomial(i);
  }

  return derivative;
}

VectorXd MultiplyPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  VectorXd multiplied_poly = VectorXd::Zero(poly1.size() + poly2.size() - 1);;
  for (int i = 0; i < poly1.size(); i++) {
    for (int j = 0; j < poly2.size(); j++) {
      multiplied_poly.reverse()(i + j) +=
          poly1.reverse()(i) * poly2.reverse()(j);
    }
  }
  return multiplied_poly;
}

VectorXd AddPolynomials(const VectorXd& poly1, const VectorXd& poly2) {
  if (poly1.size() > poly2.size()) {
    VectorXd sum = poly1;
    sum.tail(poly2.size()) += poly2;
    return sum;
  } else {
    VectorXd sum = poly2;
    sum.tail(poly1.size()) += poly1;
    return sum;
  }
}

double FindRootIterativeNewton(const Eigen::VectorXd& polynomial,
                               const double x0,
                               const double epsilon,
                               const int max_iterations) {
  double root = x0;
  const Eigen::VectorXd derivative = DifferentiatePolynomial(polynomial);
  double prev;
  for (int i = 0; i < max_iterations; i++) {
    prev = root;
    root -= EvaluatePolynomial(polynomial, root) /
            EvaluatePolynomial(derivative, root);
    if (std::abs(prev - root) < epsilon) {
      break;
    }
  }
  return root;
}

}  // namespace rpoly_plus_plus
