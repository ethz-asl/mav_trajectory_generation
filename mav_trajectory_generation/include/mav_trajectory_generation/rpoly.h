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

#ifndef MAV_TRAJECTORY_GENERATION_RPOLY_H_
#define MAV_TRAJECTORY_GENERATION_RPOLY_H_

#include <Eigen/Core>

namespace mav_trajectory_generation {

/**
 * \brief Finds the roots of a polynomial with real coefficients, using the
 * Jenkins-Traub method.
 *        Interface to the original implementation.
 *
 * http://en.wikipedia.org/wiki/Jenkins–Traub_algorithm
 * \param[in] coefficients_decreasing Coefficients of the polynomial in
 * DECREASING! order.
 * \param[in] degree Degree of the polynomial. Usually n_coefficients - 1.
 * \param[out] roots_real Real part of the roots.
 * \param[out] roots_imag Imaginary part of the roots.
 * \param[out] info Contains info about the roots found.
 * \return Number of roots found, or -1 in case of a failure.
 *
 * \note Size of info is not documented. It seems like it just stores the number
 * of iterations per root.
 * --> Recommend to pass NULL and not use it.
 */
int findRootsJenkinsTraub(const double* coefficients_decreasing, int degree,
                          double* roots_real, double* roots_imag, int info[]);

/**
 * \brief Finds the roots of a polynomial with real coefficients, using the
 * Jenkins-Traub method.
 *        Eigen wrapper for convenience.
 *
 * http://en.wikipedia.org/wiki/Jenkins–Traub_algorithm
 * \param[in] coefficients_increasing Coefficients of the polynomial in
 * INDECREASING! order.
 * \return roots Complex roots of the polynomial.
 */
template <typename Derived>
Eigen::VectorXcd findRootsJenkinsTraub(
    const Eigen::MatrixBase<Derived>& coefficients_increasing);

template <typename Derived>
Eigen::VectorXcd findRootsJenkinsTraub(
    const Eigen::MatrixBase<Derived>& coefficients_increasing) {
  EIGEN_STATIC_ASSERT_VECTOR_ONLY(Derived);

  const Derived coefficients_decreasing = coefficients_increasing.reverse();

  const int n_coefficients = coefficients_increasing.size();
  double* roots_real = new double[n_coefficients];
  double* roots_imag = new double[n_coefficients];

  int ret =
      findRootsJenkinsTraub(coefficients_decreasing.data(), n_coefficients - 1,
                            roots_real, roots_imag, NULL);

  Eigen::VectorXcd roots;

  if (ret > 0) {
    roots.resize(ret);
    for (int i = 0; i < ret; ++i) {
      roots[i] = std::complex<double>(roots_real[i], roots_imag[i]);
    }
  }

  delete[] roots_real;
  delete[] roots_imag;
  return roots;
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_RPOLY_H_
