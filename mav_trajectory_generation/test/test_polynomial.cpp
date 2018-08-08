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

#include <cmath>
#include <iostream>
#include <random>

#include <eigen-checks/entrypoint.h>
#include <eigen-checks/glog.h>
#include <eigen-checks/gtest.h>

#include "mav_trajectory_generation/motion_defines.h"
#include "mav_trajectory_generation/polynomial.h"
#include "mav_trajectory_generation/test_utils.h"
#include "mav_trajectory_generation/timing.h"

using namespace mav_trajectory_generation;

const double kSamplingInterval = 1.0e-3;
const double kEqualityResolution = 1.0e-2;
const int kDerivative = derivative_order::POSITION;

void findMinMaxBySampling(const Polynomial& polynomial, int derivative,
                          double t_start, double t_end,
                          std::pair<double, double>* min,
                          std::pair<double, double>* max) {
  min->first = t_start;
  min->second = std::numeric_limits<double>::max();
  max->first = t_start;
  max->second = std::numeric_limits<double>::lowest();
  double t = t_start;
  while (t <= t_end) {
    const double value = polynomial.evaluate(t, derivative);
    if (value < min->second) {
      min->first = t;
      min->second = value;
    }
    if (value > max->second) {
      max->first = t;
      max->second = value;
    }
    t += kSamplingInterval;
  }
}

bool approxEqual(double x_1, double x_2) {
  double dist = std::abs(x_1 - x_2);
  return dist < kEqualityResolution;
}

TEST(PolynomialTest, Convolution) {
  Eigen::VectorXd coeffs_1(2), coeffs_2(2);
  coeffs_1 << 1.0, 2.0;
  coeffs_2 << -1.0, 3.0;
  Polynomial p(coeffs_1), q(coeffs_2);
  Polynomial convolution = p * q;
  Eigen::VectorXd expected_convolution(3);
  expected_convolution << coeffs_1(0) * coeffs_2(0),
      coeffs_1(0) * coeffs_2(1) + coeffs_1(1) * coeffs_2(0),
      coeffs_1(1) * coeffs_2(1);
  CHECK_EIGEN_MATRIX_EQUAL(expected_convolution, convolution.getCoefficients());
}

TEST(PolynomialTest, FindMinMax) {
  const double kTMin = -100.0;
  const double kTMax = 100.0;
  const double kCoeffMin = -100.0;
  const double kCoeffMax = 100.0;

  std::srand(1234567);
  static int num_failures = 0;
  const int kNumPolynomials = 1e2;
  std::vector<int> derivatives_to_test = {derivative_order::POSITION,
                                          derivative_order::VELOCITY,
                                          derivative_order::ACCELERATION};
  for (int derivative : derivatives_to_test) {
    for (size_t i = 0; i < kNumPolynomials; i++) {
      // Create random polynomial.
      int num_coeffs = std::rand() % (Polynomial::kMaxN - 1) + 1 + derivative;
      Eigen::VectorXd coeffs(num_coeffs);
      for (size_t i = 0; i < num_coeffs; i++) {
        coeffs[i] = createRandomDouble(kCoeffMin, kCoeffMax);
      }
      Polynomial p(coeffs);

      // Calculate minimum and maximum.
      std::pair<double, double> min_sampling, max_sampling, min_computing,
          max_computing;
      const double t_start = createRandomDouble(kTMin, kTMax);
      const double t_end = createRandomDouble(t_start, kTMax);
      timing::Timer timer_sampling("find_min_max_sampling");
      findMinMaxBySampling(p, derivative, t_start, t_end, &min_sampling,
                           &max_sampling);
      timer_sampling.Stop();
      timing::Timer timer_analytic("find_min_max_analytic");
      bool success = p.computeMinMax(t_start, t_end, derivative, &min_computing,
                                     &max_computing);
      timer_analytic.Stop();
      if (!success) {
        std::cout << "Failed to compute roots of derivative of polynomial: "
                  << coeffs.transpose() << std::endl;
        num_failures++;
        continue;
      }
      EXPECT_TRUE(approxEqual(max_sampling.first, max_computing.first))
          << "t_max_sampling: " << max_sampling.first << std::endl
          << "t_max_computing: " << max_computing.first << std::endl
          << "max_sampling: " << max_sampling.second << std::endl
          << "max_computing: " << max_computing.second << std::endl;
      EXPECT_TRUE(approxEqual(min_sampling.first, min_computing.first))
          << "t_min_sampling: " << min_sampling.first << std::endl
          << "t_min_computing: " << min_computing.first << std::endl
          << "min_sampling: " << min_sampling.second << std::endl
          << "min_computing: " << min_computing.second << std::endl;
    }
  }
  EXPECT_EQ(num_failures, 0);
  std::cout << "Failed to compute minimum for " << num_failures << " / "
            << kNumPolynomials << " polynomials." << std::endl;
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  timing::Timing::Print(std::cout);

  return result;
}
