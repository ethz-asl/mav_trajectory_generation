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

#include "mav_trajectory_generation/polynomial.h"

#include <algorithm>

namespace mav_trajectory_generation {

bool Polynomial::findMinMax(double t_1, double t_2, int order_to_evaluate,
                            const Eigen::VectorXcd& roots_of_derivative,
                            double* t_min, double* t_max, double* min,
                            double* max) {
  // Make sure user input is correct.
  if (t_1 > t_2) {
    const double temp = t_1;
    t_1 = t_2;
    t_2 = temp;
  }

  // Evaluate polynomial at critical points.
  *min = std::numeric_limits<double>::max();
  *max = std::numeric_limits<double>::lowest();
  // Evaluate roots:
  for (size_t i = 0; i < roots_of_derivative.size(); i++) {
    // Only real roots are considered as critical points.
    if (roots_of_derivative(i).imag() != 0.0) {
      continue;
    }
    // Do not evaluate points outside the domain.
    if (roots_of_derivative(i).real() < t_1 ||
        roots_of_derivative(i).real() > t_2) {
      continue;
    }
    const double candidate =
        evaluate(roots_of_derivative(i).real(), order_to_evaluate);
    if (candidate < *min) {
      *min = candidate;
      *t_min = roots_of_derivative(i).real();
    }
    if (candidate > *max) {
      *max = candidate;
      *t_max = roots_of_derivative(i).real();
    }
  }
  // Evaluate interval end points:
  const double candidate_t_1 = evaluate(t_1, order_to_evaluate);
  const double candidate_t_2 = evaluate(t_2, order_to_evaluate);
  if (candidate_t_1 < *min) {
    *min = candidate_t_1;
    *t_min = t_1;
  }
  if (candidate_t_1 > *max) {
    *max = candidate_t_1;
    *t_max = t_1;
  }
  if (candidate_t_2 < *min) {
    *min = candidate_t_2;
    *t_min = t_2;
  }
  if (candidate_t_2 > *max) {
    *max = candidate_t_2;
    *t_max = t_2;
  }

  return true;
}

bool Polynomial::findMinMax(double t_1, double t_2, int order_to_evaluate,
                            double* t_min, double* t_max, double* min,
                            double* max) {
  Eigen::VectorXcd roots_of_derivative;
  if (findRootsJenkinsTraub(getCoefficients(order_to_evaluate + 1),
                            &roots_of_derivative) ||
      N_ - order_to_evaluate < 3) {
    return findMinMax(t_1, t_2, order_to_evaluate, roots_of_derivative, t_min, t_max,
                      min, max);
  } else {
    return false;
  }
}

Eigen::MatrixXd computeBaseCoefficients(int N) {
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

Eigen::MatrixXd Polynomial::base_coefficients_ =
    computeBaseCoefficients(Polynomial::kMaxN);

}  // namespace mav_trajectory_generation
