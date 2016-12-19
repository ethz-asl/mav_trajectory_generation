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

  // Maximum degree of a polynomial for which the static derivative (basis
  // coefficient) matrix should be evaluated for.
  constexpr int kMaxN = 12;
  // One static shared across all members of the class, computed up to order
  // kMaxN.
  static MatrixXd base_coefficients_;

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
      Eigen::VectorXd row = base_coefficients_.block(i, 0, 1, N);
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

    Eigen::MatrixBase<Derived>& _result =
        const_cast<Eigen::MatrixBase<Derived>&>(result);
  }

  /// evaluates the specified derivative of the polynomial at time t and writes
  /// the result to result
  void evaluate(double& result, double t, int derivative) const {
    CHECK_LT(derivative, N_);
    const int tmp = N_ - 1;
    Eigen::VectorXd row = base_coefficients_.block(derivative, 0, 1, N);
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

// Static functions to compute base coefficients.

/// Computes the base coefficients of the derivatives of the polynomial,
/// up to order N.
Eigen::MatrixXd Polynomial::computeBaseCoefficients(int N) {
  Eigen::MatrixXd base_coefficients(N, N);

  base_coefficients.setZero();
  base_coefficients.row(0).setOnes();

  const int DEG = N - 1;
  int order = DEG;
  for (int n = 1; n < N; n++) {
    for (int i = DEG - order; i < N; i++) {
      base_coefficients(n, i) = (order - DEG + i) * base_coefficients(n - 1, i);
    }
    order--;
  }
  return base_coefficients;
}

Polynomial::base_coefficients_ = computeBaseCoefficients(Polynomial::kMaxN);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_POLYNOMIAL_H_
