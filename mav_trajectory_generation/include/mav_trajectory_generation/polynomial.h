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

#include <Eigen/Eigen>
#include <Eigen/SVD>
#include <glog/logging.h>
#include <utility>
#include <vector>

#include "mav_trajectory_generation/rpoly.h"

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
  static constexpr int kMaxN = 12;
  // One static shared across all members of the class, computed up to order
  // kMaxN.
  static Eigen::MatrixXd base_coefficients_;

  Polynomial(int N) : N_(N), coefficients_(N) { coefficients_.setZero(); }

  // Assign arbitrary coefficients to a polynomial.
  template <class Derived>
  Polynomial(int N, const Eigen::MatrixBase<Derived>& coeffs)
      : N_(N), coefficients_(coeffs) {}

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
    assert(derivative <= N_);
    if (derivative == 0)
      return coefficients_;
    else
      return coefficients_.tail(N_ - derivative)
          .cwiseProduct(base_coefficients_.block(derivative, 0, 1, N_));
  }

  /// evaluates the polynomial at time t and writes the result to result
  void evaluate(double t, Eigen::VectorXd* result) const {
    CHECK_LE(result->size(), N_);
    const int max_deg = result->size();

    const int tmp = N_ - 1;
    for (int i = 0; i < max_deg; i++) {
      Eigen::Matrix<double, 1, Eigen::Dynamic> row =
          base_coefficients_.block(i, 0, 1, N_);
      double acc = row[tmp] * coefficients_[tmp];
      for (int j = tmp - 1; j >= i; --j) {
        acc *= t;
        acc += row[j] * coefficients_[j];
      }
      (*result)[i] = acc;
    }
  }

  /// evaluates the specified derivative of the polynomial at time t and returns
  /// the result
  double evaluate(double t, int derivative) const {
    CHECK_LT(derivative, N_);
    double result;
    const int tmp = N_ - 1;
    Eigen::Matrix<double, 1, Eigen::Dynamic> row =
        base_coefficients_.block(derivative, 0, 1, N_);
    result = row[tmp] * coefficients_[tmp];
    for (int j = tmp - 1; j >= derivative; --j) {
      result *= t;
      result += row[j] * coefficients_[j];
    }
    return result;
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

  /**
    * \brief Computes the base coefficients with the according powers of t, as
    * e.g. needed for computation of (in)equality constraints
    * \param[out] coeffs vector to write the coefficients to
    * \param[in] derivative of the polynomial for which the coefficients have to
    * be computed
    * \param[in] t time of evaluation
    */
  template <class Derived>
  static void baseCoeffsWithTime(int N, int derivative, double t,
                                 const Eigen::MatrixBase<Derived>& coeffs) {
    CHECK_LT(derivative, N);
    CHECK_GE(derivative, 0);
    Eigen::MatrixBase<Derived>& c =
        const_cast<Eigen::MatrixBase<Derived>&>(coeffs);

    c.setZero();
    // first coefficient doesn't get multiplied
    c[derivative] = base_coefficients_(derivative, derivative);

    if (t == 0) return;

    double _t = t;
    // now multiply increasing power of t towards the right
    for (int j = derivative + 1; j < N; j++) {
      c[j] = base_coefficients_(derivative, j) * _t;
      _t = _t * t;
    }
  }

  /**
   * \brief Convenience method to compute the base coefficents with time
   * \sa static void baseCoeffsWithTime(const Eigen::MatrixBase<Derived> &
   * coeffs, int derivative, double t)
   */
  static Eigen::Matrix<double, 1, Eigen::Dynamic> baseCoeffsWithTime(
      int N, int derivative, double t) {
    Eigen::Matrix<double, Eigen::Dynamic, 1> c(N, 1);
    baseCoeffsWithTime(N, derivative, t, c);
    return c;
  }

  /**
   * \brief Computes the Jacobian of the integral over the squared derivative
   * \param[out] C Jacobian matrix to write into. If C is dynamic, the correct
   * size has to be set.
   * \param[in] t time of evaluation
   * \param[in] derivative used to compute the cost
   */
  template <class Derived>
  static void quadraticCostJacobian(int N, int derivative, double t,
                                    const Eigen::MatrixBase<Derived>& C) {
    CHECK_LT(derivative, N);
    if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
        Derived::ColsAtCompileTime == Eigen::Dynamic) {
      CHECK_EQ(C.rows(), N);
      CHECK_EQ(C.cols(), N);
    } else {
      // EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Derived, N, N);
    }

    Eigen::MatrixBase<Derived>& _C = const_cast<Eigen::MatrixBase<Derived>&>(C);
    _C.setZero();

    for (int col = 0; col < N - derivative; col++) {
      for (int row = 0; row < N - derivative; row++) {
        double exp = (N - 1 - derivative) * 2 + 1 - row - col;

        _C(N - 1 - row, N - 1 - col) =
            base_coefficients_(derivative, N - 1 - row) *
            base_coefficients_(derivative, N - 1 - col) * pow(t, exp) * 2 / exp;
      }
    }
  }

  /**
   * \brief convenience method to compute the Jacobian of the quadratic cost
   * \sa static void quadraticCostJacobian(const Eigen::MatrixBase<Derived> & C,
   * Scalar t, int derivative)
   */
  static Eigen::MatrixXd quadraticCostJacobian(int N, int derivative,
                                               double t) {
    Eigen::MatrixXd C(N, N);
    quadraticCostJacobian(N, derivative, t, C);
    return C;
  }

 private:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  int N_;
  Eigen::VectorXd coefficients_;
};

// Static functions to compute base coefficients.

/// Computes the base coefficients of the derivatives of the polynomial,
/// up to order N.
Eigen::MatrixXd computeBaseCoefficients(int N);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_POLYNOMIAL_H_
