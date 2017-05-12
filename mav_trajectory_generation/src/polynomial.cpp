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

Eigen::VectorXd Polynomial::convolve(const Eigen::VectorXd& data,
                                            const Eigen::VectorXd& kernel) {
  const int convolution_dimension = data.size() + kernel.size() - 1;
  Eigen::VectorXd convolved = Eigen::VectorXd::Zero(convolution_dimension);
  Eigen::VectorXd kernel_reverse = kernel.reverse();

  for (int i = 0; i < convolution_dimension; i++) {
    const int data_idx = i - kernel.size() + 1;

    int lower_bound = std::max(0, -data_idx);
    int upper_bound = std::min(kernel.size(), data.size() - data_idx);

    for (int kernel_idx = lower_bound; kernel_idx < upper_bound; ++kernel_idx) {
      convolved[i] += kernel_reverse[kernel_idx] * data[data_idx + kernel_idx];
    }
  }
  return convolved;
}

Eigen::MatrixXd Polynomial::base_coefficients_ =
    computeBaseCoefficients(Polynomial::kMaxN);

}  // namespace mav_trajectory_generation
