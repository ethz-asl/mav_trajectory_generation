// Copyright (C) 2015 Chris Sweeney
// All rights reserved.
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

#ifndef MAV_TRAJECTORY_GENERATION_RPOLY_PLUS_PLUS_POLYNOMIAL_H_
#define MAV_TRAJECTORY_GENERATION_RPOLY_PLUS_PLUS_POLYNOMIAL_H_

#include <Eigen/Core>

namespace rpoly_plus_plus {

// All polynomials are assumed to be the form
//
//   sum_{i=0}^N polynomial(i) x^{N-i}.
//
// and are given by a vector of coefficients of size N + 1.

// Remove leading terms with zero coefficients.
Eigen::VectorXd RemoveLeadingZeros(const Eigen::VectorXd& polynomial_in);

// Evaluate the polynomial at x using the Horner scheme.
template <typename T>
inline T EvaluatePolynomial(const Eigen::VectorXd& polynomial, const T& x) {
  T v = 0.0;
  for (int i = 0; i < polynomial.size(); ++i) {
    v = v * x + polynomial(i);
  }
  return v;
}

// Return the derivative of the given polynomial. It is assumed that
// the input polynomial is at least of degree zero.
Eigen::VectorXd DifferentiatePolynomial(const Eigen::VectorXd& polynomial);

// Multiplies the two polynoimals together.
Eigen::VectorXd MultiplyPolynomials(const Eigen::VectorXd& poly1,
                                    const Eigen::VectorXd& poly2);

// Adds two polynomials together.
Eigen::VectorXd AddPolynomials(const Eigen::VectorXd& poly1,
                               const Eigen::VectorXd& poly2);

// Find a root from the starting guess using Newton's method.
double FindRootIterativeNewton(const Eigen::VectorXd& polynomial,
                               const double x0,
                               const double epsilon,
                               const int max_iterations);

}  // namespace rpoly_plus_plus

#endif  // MAV_TRAJECTORY_GENERATION_RPOLY_PLUS_PLUS_POLYNOMIAL_H_
