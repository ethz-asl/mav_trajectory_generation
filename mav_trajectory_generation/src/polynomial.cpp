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
#include <limits>

namespace mav_trajectory_generation {

void Polynomial::findMinMaxCandidates(
    double t_start, double t_end,
    const Eigen::VectorXcd& roots_derivative_of_derivative,
    std::vector<double>* candidates) const {
  CHECK_NOTNULL(candidates);
  candidates->clear();
  candidates->reserve(roots_derivative_of_derivative.size() + 2);
  candidates->push_back(t_start);
  candidates->push_back(t_end);
  for (size_t i = 0; i < roots_derivative_of_derivative.size(); i++) {
    // Only real roots are considered as critical points.
    if (std::abs(roots_derivative_of_derivative[i].imag()) >
        std::numeric_limits<double>::epsilon()) {
      continue;
    }
    const double candidate = roots_derivative_of_derivative[i].real();
    // Do not evaluate points outside the domain.
    if (candidate < std::min(t_start, t_end) ||
        candidate > std::max(t_start, t_end)) {
      continue;
    } else {
      candidates->push_back(candidate);
    }
  }
}

bool Polynomial::findMinMaxCandidates(double t_start, double t_end,
                                      int derivative,
                                      std::vector<double>* candidates) const {
  CHECK_NOTNULL(candidates);
  candidates->clear();
  if (N_ - derivative - 1 < 0) {
    LOG(WARNING) << "N - derivative - 1 has to be at least 0.";
    return false;
  }
  Eigen::VectorXcd roots_derivative_of_derivative;
  if (!findRootsJenkinsTraub(getCoefficients(derivative + 1),
                             &roots_derivative_of_derivative)) {
    return false;
  } else {
    findMinMaxCandidates(t_start, t_end, roots_derivative_of_derivative,
                         candidates);
    return true;
  }
}

bool Polynomial::findMinMax(
    double t_start, double t_end, int derivative,
    const Eigen::VectorXcd& roots_derivative_of_derivative,
    std::pair<double, double>* min, std::pair<double, double>* max) const {
  CHECK_NOTNULL(min);
  CHECK_NOTNULL(max);
  // Find candidates in interval t_start to t_end computing the roots.
  std::vector<double> candidates;
  findMinMaxCandidates(t_start, t_end, roots_derivative_of_derivative,
                       &candidates);
  // Evaluate minimum and maximum.
  return findMinMax(candidates, derivative, min, max);
}

bool Polynomial::findMinMax(double t_start, double t_end, int derivative,
                            std::pair<double, double>* min,
                            std::pair<double, double>* max) const {
  CHECK_NOTNULL(min);
  CHECK_NOTNULL(max);
  // Find candidates in interval t_start to t_end by computing the roots.
  std::vector<double> candidates;
  if (!findMinMaxCandidates(t_start, t_end, derivative, &candidates)) {
    return false;
  }
  // Evaluate minimum and maximum.
  return findMinMax(candidates, derivative, min, max);
}

bool Polynomial::findMinMax(const std::vector<double>& candidates,
                            int derivative, std::pair<double, double>* min,
                            std::pair<double, double>* max) const {
  CHECK_NOTNULL(min);
  CHECK_NOTNULL(max);
  if (candidates.empty()) {
    LOG(WARNING) << "Cannot find extrema from an empty candidates vector.";
    return false;
  }
  min->first = candidates[0];
  min->second = std::numeric_limits<double>::max();
  max->first = candidates[0];
  max->second = std::numeric_limits<double>::lowest();

  for (const double& t : candidates) {
    const double value = evaluate(t, derivative);
    if (value < min->second) {
      min->first = t;
      min->second = value;
    }
    if (value > max->second) {
      max->first = t;
      max->second = value;
    }
  }
  return true;
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

Eigen::VectorXd Polynomial::convolve(const Eigen::VectorXd& data,
                                     const Eigen::VectorXd& kernel) {
  const int convolution_dimension =
      getConvolutionLength(data.size(), kernel.size());
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
