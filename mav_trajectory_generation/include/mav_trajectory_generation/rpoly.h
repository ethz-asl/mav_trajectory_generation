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
#include <iostream>

namespace mav_trajectory_generation {

// Finds the roots of a polynomial with real coefficients, using the
// Jenkins-Traub method. Interface to the original implementation.
// http://en.wikipedia.org/wiki/Jenkins–Traub_algorithm
// Input: coefficients_decreasing = Coefficients of the polynomial in
// DECREASING! order.
// Input: degree = Degree of the polynomial. Usually n_coefficients - 1.
// Output: roots_real = Real part of the roots.
// Output: roots_imag = Imaginary part of the roots.
// Output: info = Contains info about the roots found.
// Output: return = Number of roots found, or -1 in case of a failure.
// Note: Size of info is not documented. It seems like it just stores the number
// of iterations per root. --> Recommend to pass NULL and not use it.
int findRootsJenkinsTraub(const double* coefficients_decreasing, int degree,
                          double* roots_real, double* roots_imag, int info[]);

// Finds the roots of a polynomial with real coefficients, using the
// Jenkins-Traub method. Eigen wrapper for convenience.
// http://en.wikipedia.org/wiki/Jenkins–Traub_algorithm
// Input: coefficients_increasing = Coefficients of the polynomial in
// INDECREASING! order.
// Output: roots = Complex roots of the polynomial.
// Output: return = Root calculation success.
bool findRootsJenkinsTraub(const Eigen::VectorXd& coefficients_increasing,
                           Eigen::VectorXcd* roots);

// Finds the roots of a polynomial with real coefficients, using the
// Jenkins-Traub method. Convenience method.
Eigen::VectorXcd findRootsJenkinsTraub(
    const Eigen::VectorXd& coefficients_increasing);

// Find the last non-zero entry in a vector and returns its index. Returns -1 in
// case of all zeros.
int findLastNonZeroCoeff(const Eigen::VectorXd& coefficients);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_RPOLY_H_
