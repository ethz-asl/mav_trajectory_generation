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
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */
#include "mav_trajectory_generation/polynomial.h"
#include "mav_trajectory_generation/rpoly/rpoly_ak1.h"

#include <algorithm>
#include <limits>

namespace mav_trajectory_generation {

bool Polynomial::getRoots(int derivative, Eigen::VectorXcd* roots) const {
  return findRootsJenkinsTraub(getCoefficients(derivative), roots);
}

bool Polynomial::selectMinMaxCandidatesFromRoots(
    double t_start, double t_end,
    const Eigen::VectorXcd& roots_derivative_of_derivative,
    std::vector<double>* candidates) {
  CHECK_NOTNULL(candidates);
  if (t_start > t_end) {
    LOG(WARNING) << "t_start is greater than t_end.";
    return false;
  }
  candidates->clear();
  candidates->reserve(roots_derivative_of_derivative.size() + 2);
  // Put start and end in, as they are valid candidates.
  candidates->push_back(t_start);
  candidates->push_back(t_end);
  for (size_t i = 0;
       i < static_cast<size_t>(roots_derivative_of_derivative.size()); i++) {
    // Only real roots are considered as critical points.
    if (std::abs(roots_derivative_of_derivative[i].imag()) >
        std::numeric_limits<double>::epsilon()) {
      continue;
    }
    const double candidate = roots_derivative_of_derivative[i].real();

    // Do not evaluate points outside the domain.
    if (candidate < t_start || candidate > t_end) {
      continue;
    } else {
      candidates->push_back(candidate);
    }
  }
  return true;
}

bool Polynomial::computeMinMaxCandidates(
    double t_start, double t_end, int derivative,
    std::vector<double>* candidates) const {
  CHECK_NOTNULL(candidates);
  candidates->clear();
  if (N_ - derivative - 1 < 0) {
    LOG(WARNING) << "N - derivative - 1 has to be at least 0.";
    return false;
  }
  Eigen::VectorXcd roots;
  bool success = getRoots(derivative + 1, &roots);
  if (!success) {
    VLOG(1) << "Couldn't find roots, polynomial may be constant.";
  }
  if (!selectMinMaxCandidatesFromRoots(t_start, t_end, roots, candidates)) {
    return false;
  }
  return true;
}

bool Polynomial::selectMinMaxFromRoots(
    double t_start, double t_end, int derivative,
    const Eigen::VectorXcd& roots_derivative_of_derivative,
    std::pair<double, double>* minimum,
    std::pair<double, double>* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  // Find candidates in interval t_start to t_end computing the roots.
  std::vector<double> candidates;
  if (!selectMinMaxCandidatesFromRoots(
          t_start, t_end, roots_derivative_of_derivative, &candidates)) {
    return false;
  }
  // Evaluate minimum and maximum.
  return selectMinMaxFromCandidates(candidates, derivative, minimum, maximum);
}

bool Polynomial::computeMinMax(double t_start, double t_end, int derivative,
                               std::pair<double, double>* minimum,
                               std::pair<double, double>* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  // Find candidates in interval t_start to t_end by computing the roots.
  std::vector<double> candidates;
  if (!computeMinMaxCandidates(t_start, t_end, derivative, &candidates)) {
    return false;
  }
  // Evaluate minimum and maximum.
  return selectMinMaxFromCandidates(candidates, derivative, minimum, maximum);
}

bool Polynomial::selectMinMaxFromCandidates(
    const std::vector<double>& candidates, int derivative,
    std::pair<double, double>* minimum,
    std::pair<double, double>* maximum) const {
  CHECK_NOTNULL(minimum);
  CHECK_NOTNULL(maximum);
  if (candidates.empty()) {
    LOG(WARNING) << "Cannot find extrema from an empty candidates vector.";
    return false;
  }
  minimum->first = candidates[0];
  minimum->second = std::numeric_limits<double>::max();
  maximum->first = candidates[0];
  maximum->second = std::numeric_limits<double>::lowest();

  for (const double& t : candidates) {
    const double value = evaluate(t, derivative);
    if (value < minimum->second) {
      minimum->first = t;
      minimum->second = value;
    }
    if (value > maximum->second) {
      maximum->first = t;
      maximum->second = value;
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

bool Polynomial::getPolynomialWithAppendedCoefficients(
    int new_N, Polynomial* new_polynomial) const {
  if (new_N == N_) {
    *new_polynomial = *this;
    return true;
  } else if (new_N < N_) {
    LOG(WARNING) << "You shan't decrease the number of coefficients.";
    *new_polynomial = *this;
    return false;
  } else {
    Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(new_N);
    coeffs.head(N_) = coefficients_;
    *new_polynomial = Polynomial(coeffs);
    return true;
  }
}

void Polynomial::scalePolynomialInTime(double scaling_factor) {
  double scale = 1.0;
  for (int n = 0; n < N_; n++) {
    coefficients_[n] *= scale;
    scale *= scaling_factor;
  }
}

void Polynomial::offsetPolynomial(const double offset) {
  if (coefficients_.size() == 0) return;

  coefficients_[0] += offset;
}

Eigen::MatrixXd Polynomial::base_coefficients_ =
    computeBaseCoefficients(Polynomial::kMaxConvolutionSize);

}  // namespace mav_trajectory_generation
