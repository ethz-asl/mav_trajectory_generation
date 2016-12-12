/*

 Copyright (c) 2013, Markus Achtelik, ASL, ETH Zurich, Switzerland
 You can contact the author at <markus dot achtelik at mavt dot ethz dot ch>

 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of ETHZ-ASL nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL ETHZ-ASL BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef MAV_TRAJECTORY_GENERATION_POLYNOMIAL_H_
#define MAV_TRAJECTORY_GENERATION_POLYNOMIAL_H_

#include <utility>
#include <vector>

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <glog/logging.h>

namespace mav_trajectory_generation {
/**
 * \brief Implementation of polynomials of order _N-1. Order must be known at
 * compile time.
 * Polynomial coefficients are stored with increasing powers, i.e. \f$c_0 +
 * c_1*t ... c_{N-1} * t^{N-1}\f$
 * \tparam _N Number of coefficients of the polynomial.
 * \tparam double Floating point type of the double.
 */
class Polynomial {
 public:
  typedef std::vector<Polynomial, Eigen::aligned_allocator<Polynomial>> Vector;

  Polynomial(int N_) : N(N_) {}

  // Assign arbitrary coefficients to a polynomial.
  template <class Derived>
  Polynomial(int N_, const Eigen::MatrixBase<Derived>& coeffs)
      : N(N_), coefficients_(coeffs) {}

  /// Get the number of coefficients (order + 1) of the polynomial.
  int N() const { return N_; }

  /**
   * \brief sets up the internal representation from coeffs
   * coefficients are stored in increasing order with the power of t, i.e. c1 +
   * c2*t + c3*t^2 ==> coeffs = [c1 c2 c3]
   */
  template <class Derived>
  void setCoefficients(const Eigen::MatrixBase<Derived>& coeffs) {
    coefficients_ = coeffs;
  }

  /**
   * \brief Returns the coefficients for the specified derivative of the
   * polynomial.
   */
  Eigen::Matrix<double, 1, Eigen::Dynamic> getCoefficients(
      int derivative = 0) const {
    assert(derivative <= N);
    if (derivative == 0)
      return coefficients_;
    else
      return coefficients_.tail(N - derivative)
          .cwiseProduct(
              base_coefficients_.row(derivative).tail(N - derivative));
  }

  /// evaluates the polynomial at time t and writes the result to result
  void evaluate(double t, Eigen::VectorXd* result) const {
    CHECK_LE(result.rows(), N);
    const int max_deg = result.size();

    const int deg = N - 1;
    for (int i = 0; i < max_deg; i++) {
      const VectorR row = base_coefficients_.row(i);
      double acc = row[tmp] * coefficients_[tmp];
      for (int j = tmp - 1; j >= i; --j) {
        acc *= t;
        acc += row[j] * coefficients_[j];
      }
      _result[i] = acc;
    }
  }

  template <class Derived>
  void evaluate(const Eigen::MatrixBase<Derived>& result, double t) const {
    CHECK_LE(result.rows(), N);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);
    const int max_deg = result.size();

    Eigen::MatrixBase<Derived> &_result =
        const_cast<Eigen::MatrixBase<Derived> &>(result);


  }

  /// evaluates the specified derivative of the polynomial at time t and writes
  /// the result to result
  void evaluate(double &result, double t, int derivative) const {
    CHECK_LT(derivative, N_);
    const int tmp = N_ - 1;
    const VectorR row = base_coefficients_.row(derivative);
    result = row[tmp] * coefficients_[tmp];
    for (int j = tmp - 1; j >= derivative; --j) {
      result *= t;
      result += row[j] * coefficients_[j];
    }
  }

  /// evaluates the specified derivative of the polynomial at time t and returns
  /// the result
  double evaluate(double t, int derivative) const {
    double res;
    evaluate(res, t, derivative);
    return res;
  }

  /// evaluates the polynomial at time t and returns the result
  template <int max_deg>
  inline Eigen::Matrix<double, max_deg, 1> evaluate(double t) const {
    Eigen::Matrix<double, max_deg, 1> result;
    evaluate(result, t);
    return result;
  }

  /// evaluates the polynomial at times in t and writes the result for each time
  /// into the corresponding column of result
  template <int max_deg, int n_samples>
  void evaluate(const Eigen::Matrix<double, 1, n_samples>& t,
                Eigen::Matrix<double, max_deg, n_samples> result) const {
    Eigen::Matrix<double, max_deg, 1> _result;
    for (int i = 0; i < n_samples; i++) {
      evaluate(t[i], _result);
      result.col(i) = result;
    }
  }

  /**
   * \brief Computes the complex roots of the polynomial.
   * Only for the polynomial itself, not for its derivatives.
   */
  Eigen::VectorXcd computeRoots() const {
    //      Companion matrix method , see
    //      http://en.wikipedia.org/wiki/Companion_matrix.
    //      Works, but is not very stable for high condition numbers. Could be
    //      eigen's eigensolver.
    //      However, would not need the dependency to rpoly.
    //      const size_t nc = N - 1;
    //      typedef Eigen::Matrix<double, nc, nc> CompanionMatrix;
    //      CompanionMatrix companion;
    //      companion.template row(0).setZero();
    //      companion.template block<nc - 1, nc - 1>(1, 0).setIdentity();
    //      companion.template col(nc - 1) = - coefficients_.template head<nc>()
    //      / coefficients_[N - 1];
    //
    //      Eigen::EigenSolver<CompanionMatrix> es(companion, false);
    //      return es.eigenvalues();

    return findRootsJenkinsTraub(coefficients_);
  }

 private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Eigen::VectorXd coefficients_;
};

}  // mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_POLYNOMIAL_H_
