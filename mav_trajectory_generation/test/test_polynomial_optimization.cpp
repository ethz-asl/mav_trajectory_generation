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

#include <iostream>
#include <limits>
#include <random>

#include <eigen-checks/entrypoint.h>
#include <eigen-checks/glog.h>
#include <eigen-checks/gtest.h>

#include "mav_trajectory_generation/polynomial_optimization_linear.h"
#include "mav_trajectory_generation/polynomial_optimization_nonlinear.h"
#include "mav_trajectory_generation/timing.h"

using namespace mav_trajectory_generation;

const int N = 10;
const int max_derivative = derivative_order::SNAP;
const size_t derivative_to_optimize = derivative_order::SNAP;

Eigen::IOFormat matlab_format(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "[",
                              "]");

template <class T1, class T2>
bool checkMatrices(const Eigen::MatrixBase<T1>& m1,
                   const Eigen::MatrixBase<T2>& m2, double tol) {
  return (m1 - m2).cwiseAbs().maxCoeff() < tol;
}

double getMaximumMagnitude(const std::vector<Segment>& segments,
                           size_t derivative, double dt = 0.01) {
  double maximum = -1e9;

  for (const Segment& s : segments) {
    for (double ts = 0; ts < s.getTime(); ts += dt) {
      double current_value = s.evaluate(ts, derivative).norm();
      if (current_value > maximum) maximum = current_value;
    }
  }
  return maximum;
}

double computeCostNumeric(const std::vector<Segment>& segments,
                          size_t derivative, double dt = 0.001) {
  double cost = 0;

  for (const Segment& s : segments) {
    for (double ts = 0; ts < s.getTime(); ts += dt) {
      cost += s.evaluate(ts, derivative).squaredNorm() * dt;
    }
  }
  return cost;
}

void checkPath(const Vertex::Vector& vertices,
               const std::vector<Segment>& segments) {
  const double tol = 1e-6;
  size_t n_vertices = vertices.size();
  size_t n_segments = segments.size();
  EXPECT_EQ(n_segments, n_vertices - 1);

  for (size_t i = 0; i < n_segments; ++i) {
    const Vertex& v_begin = vertices[i];
    const Vertex& v_end = vertices[i + 1];
    const Segment& segment = segments[i];

    // check if fixed constraints are met.
    for (Vertex::Constraints::const_iterator it = v_begin.cBegin();
         it != v_begin.cEnd(); ++it) {
      const int derivative = it->first;
      const Vertex::ConstraintValue desired = it->second;
      const Vertex::ConstraintValue actual = segment.evaluate(0, derivative);
      Eigen::Matrix3d m, n;
      std::stringstream segment_derivative;
      printSegment(segment_derivative, segment, derivative);
      EXPECT_TRUE(EIGEN_MATRIX_NEAR(desired, actual, tol))
          << "[fixed constraint check t=0] at vertex " << i
          << " and constraint " << positionDerivativeToString(derivative)
          << "\nsegment:\n" << segment << segment_derivative.str();
    }
    for (Vertex::Constraints::const_iterator it = v_end.cBegin();
         it != v_end.cEnd(); ++it) {
      const int derivative = it->first;
      const Vertex::ConstraintValue desired = it->second;
      const Vertex::ConstraintValue actual =
          segment.evaluate(segment.getTime(), derivative);
      std::stringstream segment_derivative;
      printSegment(segment_derivative, segment, derivative);
      EXPECT_TRUE(EIGEN_MATRIX_NEAR(desired, actual, tol))
          << "[fixed constraint check] at vertex " << i + 1
          << " and constraint " << positionDerivativeToString(derivative)
          << "\nsegment:\n" << segment << segment_derivative.str();
    }

    // Check if values at vertices are continuous.
    if (i > 0) {
      const Segment& last_segment = segments[i - 1];
      for (size_t derivative = 0; derivative < N / 2; ++derivative) {
        const Vertex::ConstraintValue last_segment_value =
            last_segment.evaluate(last_segment.getTime(), derivative);
        const Vertex::ConstraintValue current_segment_value =
            segment.evaluate(0, derivative);
        std::stringstream segment_derivative;
        printSegment(segment_derivative, segment, derivative);
        EXPECT_TRUE(
            EIGEN_MATRIX_NEAR(last_segment_value, current_segment_value, tol))
            << "[vertex continuity] at vertex " << i << " and constraint "
            << positionDerivativeToString(derivative) << "\nsegment:\n"
            << segment << segment_derivative.str();
      }
    }
  }
}

bool checkCost(double cost_to_check, const std::vector<Segment>& segments,
               size_t derivative, double relative_tolerance) {
  CHECK_GE(derivative, size_t(0));
  CHECK(relative_tolerance >= 0.0 && relative_tolerance <= 1.0);
  const double sampling_interval = 0.001;
  double cost_numeric =
      computeCostNumeric(segments, derivative, sampling_interval);

  if (std::abs(cost_numeric - cost_to_check) >
      cost_numeric * relative_tolerance) {
    std::cout << "[FAIL]: tested cost is" << cost_to_check
              << " but real cost is " << cost_numeric
              << ". Difference vs. allowed difference: "
              << std::abs(cost_numeric - cost_to_check) << "/"
              << cost_numeric * relative_tolerance << std::endl;
    return false;
  }

  return true;
}

TEST(MavTrajectoryGeneration, PathPlanning_TestVertexGeneration1D) {
  Vertex::Vector vertices;
  const double p_min = -50;
  const double p_max = 50;
  vertices = createRandomVertices1D(max_derivative, 100, p_min, p_max, 0);

  EXPECT_EQ(vertices.front().getNumberOfConstraints(), N / 2);
  EXPECT_EQ(vertices.back().getNumberOfConstraints(), N / 2);

  for (const Vertex& v : vertices) {
    EXPECT_TRUE(v.hasConstraint(derivative_order::POSITION));
    Eigen::VectorXd c;
    v.getConstraint(derivative_order::POSITION, &c);
    EXPECT_LE(c[0], p_max);
    EXPECT_GE(c[0], p_min);
  }
}

TEST(MavTrajectoryGeneration, PathPlanning_TestVertexGeneration3D) {
  Eigen::VectorXd pos_min(3), pos_max(3);
  pos_min << -10.0, -20.0, -10.0;
  pos_max << 10.0, 20.0, 10.0;
  Vertex::Vector vertices;

  vertices = createRandomVertices(max_derivative, 100, pos_min, pos_max, 12345);

  EXPECT_EQ(vertices.front().getNumberOfConstraints(), N / 2);
  EXPECT_EQ(vertices.back().getNumberOfConstraints(), N / 2);

  for (const Vertex& v : vertices) {
    EXPECT_TRUE(v.hasConstraint(derivative_order::POSITION));
    Eigen::VectorXd c;
    v.getConstraint(derivative_order::POSITION, &c);
    for (int i = 0; i < 3; ++i) {
      EXPECT_LE(c[i], pos_max[i]);
      EXPECT_GE(c[i], pos_min[i]);
    }
  }
}

TEST(MavTrajectoryGeneration, PathPlanning_A_matrix_inversion) {
  const double max_time = 60;
  for (double t = 1; t <= max_time; t += 1) {
    Eigen::Matrix<double, N, N> A, Ai, Ai_eigen;
    PolynomialOptimization<N>::setupMappingMatrix(t, &A);
    PolynomialOptimization<N>::invertMappingMatrix(A, &Ai);
    Ai_eigen = A.inverse();
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(Ai, Ai_eigen, 1.0e-10)) << "time was " << t
                                                          << std::endl;
  }
}

TEST(MavTrajectoryGeneration, PathPlanningUnconstrained_1D_10_segments) {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 10, -10, 10, 12);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  timing::Timer timer_setup("setup_1D_10s");
  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear_1D_10s");
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);

  std::cout << "Base coefficients: "
            << Polynomial::base_coefficients_.block(3, 0, 1, N) << std::endl;

  checkPath(vertices, segments);
  double v_max = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max << " a_max: " << a_max << std::endl;

  timing::Timer timer_cost("cost_1D_10s");
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  EXPECT_LT(v_max, approximate_v_max * 2.0);
  EXPECT_LT(a_max, approximate_a_max * 2.0);
  EXPECT_TRUE(
      checkCost(opt.computeCost(), segments, derivative_to_optimize, 0.1));
}

TEST(MavTrajectoryGeneration, PathPlanningUnconstrained_1D_50_segments) {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 50, -10, 10, 123);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  timing::Timer timer_setup("setup_1D_50s");
  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear_1D_50s");
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);

  checkPath(vertices, segments);
  double v_max = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max << " a_max: " << a_max << std::endl;

  timing::Timer timer_cost("cost_1D_50s");
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  EXPECT_LT(v_max, approximate_v_max * 2.0);
  EXPECT_LT(a_max, approximate_a_max);
  EXPECT_TRUE(
      checkCost(opt.computeCost(), segments, derivative_to_optimize, 0.1));
}

TEST(MavTrajectoryGeneration, PathPlanningUnconstrained_1D_100_segments) {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 100, -10, 10, 1234);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  timing::Timer timer_setup("setup_1D_100s");
  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear_1D_100s");
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);

  checkPath(vertices, segments);
  double v_max = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max << " a_max: " << a_max << std::endl;

  timing::Timer timer_cost("cost_1D_100s");
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  EXPECT_LT(v_max, approximate_v_max * 5.0);
  EXPECT_LT(a_max, approximate_a_max * 2.0);
  EXPECT_TRUE(
      checkCost(opt.computeCost(), segments, derivative_to_optimize, 0.1));
}

TEST(MavTrajectoryGeneration,
     PathPlanningUnconstrained_1D_100_segments_high_segment_times) {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 100, -50, 50, 12345);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  timing::Timer timer_setup("setup_1D_100s_long");
  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear_1D_100s_long");
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);

  checkPath(vertices, segments);
  double v_max = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max << " a_max: " << a_max << std::endl;

  timing::Timer timer_cost("cost_1D_100s_long");
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  // max. velocity estimation  beforehand is not accurate here anymore
  EXPECT_LT(v_max, approximate_v_max * 5.0);
  EXPECT_LT(a_max, approximate_a_max * 2.0);
  EXPECT_TRUE(
      checkCost(opt.computeCost(), segments, derivative_to_optimize, 0.1));
}

TEST(MavTrajectoryGeneration,
     PathPlanningUnconstrained_3D_100_segments_high_segment_times) {
  Eigen::VectorXd pos_min(3), pos_max(3);
  pos_min << -10.0, -20.0, -10.0;
  pos_max << 10.0, 20.0, 10.0;
  Vertex::Vector vertices;
  vertices = createRandomVertices(max_derivative, 100, pos_min, pos_max, 12345);

  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  timing::Timer timer_setup("setup_3D_100s_long");
  PolynomialOptimization<N> opt(3);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear_3D_100s_long");
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);

  checkPath(vertices, segments);
  double v_max = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max << " a_max: " << a_max << std::endl;

  timing::Timer timer_cost("cost_3D_100s_long");
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  /// max. velocity estimation beforehand is not accurate here anymore
  EXPECT_LT(v_max, approximate_v_max * 5.0);
  EXPECT_LT(a_max, approximate_a_max * 2.0);
  EXPECT_TRUE(
      checkCost(opt.computeCost(), segments, derivative_to_optimize, 0.1));
}

bool checkExtrema(const std::vector<double>& testee,
                  const std::vector<double>& reference, double tol = 0.01) {
  for (double t : testee) {
    bool found_match = false;
    for (double r : reference) {
      if (std::abs(t - r) < tol) {
        found_match = true;
        break;
      }
    }
    if (!found_match) {
      std::cout << "[ERROR]: did not find matching extremum for " << t
                << " in ";
      for (double r : reference) {
        std::cout << r << " ";
      }
      std::cout << std::endl;
      return false;
    }
  }
  return true;
}

TEST(MavTrajectoryGeneration,
     PathOptimization_1D_segment_extrema_of_magnitude) {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 100, -10, 10, 1234);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt.solveLinear();

  Segment::Vector segments;
  opt.getSegments(&segments);

  timing::Timer time_analytic("time_extrema_analytic_1", false);
  timing::Timer time_analytic_template_free(
      "time_extrema_analytic_1_template_free", false);
  timing::Timer time_sampling("time_extrema_sampling_1", false);
  int segment_idx = 0;
  for (const Segment& s : segments) {
    std::vector<double> res;
    time_analytic.Start();
    opt.computeSegmentMaximumMagnitudeCandidates<1>(s, 0, s.getTime(), &res);
    time_analytic.Stop();

    std::vector<double> res_template_free;
    std::vector<int> dimensions = {0};
    time_analytic_template_free.Start();
    s.computeMinMaxMagnitudeCandidateTimes(1, 0.0, s.getTime(), dimensions,
                                           &res_template_free);
    time_analytic_template_free.Stop();

    std::vector<double> res_sampling;
    time_sampling.Start();
    opt.computeSegmentMaximumMagnitudeCandidatesBySampling<1>(
        s, 0, s.getTime(), 0.001, &res_sampling);
    time_sampling.Stop();

    constexpr double check_tolerance = 0.01;
    bool success = checkExtrema(res, res_sampling, check_tolerance);
    if (!success) {
      std::cout << "############CHECK XTREMA FAILED: \n";
      std::cout << "segment idx: " << segment_idx << "/" << segments.size()
                << std::endl;

      std::cout << "real eigs for segment with time " << s.getTime() << "s : ";

      for (const double& t : res) std::cout << t << " ";
      std::cout << std::endl;
      std::cout << "sampling: ";
      for (const double& t : res_sampling) std::cout << t << " ";
      std::cout << std::endl;

      std::cout << "vx = [ "
                << s[0].getCoefficients(derivative_order::VELOCITY)
                       .reverse()
                       .format(matlab_format)
                << "];\n";
      std::cout << "t = 0:0.001:" << s.getTime() << "; \n";
    }
    EXPECT_TRUE(success);

    EXPECT_EQ(res.size(), res_template_free.size() - 2);
    for (size_t i = 0; i < res.size(); i++) {
      EXPECT_EQ(res[i], res_template_free[i + 2]);
    }

    ++segment_idx;
  }

  double v_max_ref = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max_ref =
      getMaximumMagnitude(segments, derivative_order::ACCELERATION);

  timing::Timer time_analytic_v("time_extrema_analytic_v");
  Extremum v_max =
      opt.computeMaximumOfMagnitude<derivative_order::VELOCITY>(nullptr);
  time_analytic_v.Stop();
  timing::Timer time_analytic_a("time_extrema_analytic_a");
  Extremum a_max =
      opt.computeMaximumOfMagnitude<derivative_order::ACCELERATION>(nullptr);
  time_analytic_a.Stop();

  std::cout << "v_max " << v_max << std::endl << "a_max " << a_max << std::endl;

  EXPECT_NEAR(v_max_ref, v_max.value, 0.01);
  EXPECT_NEAR(a_max_ref, a_max.value, 0.01);
}

TEST(MavTrajectoryGeneration, PathOptimization3D_segment_extrema_of_magnitude) {
  Eigen::VectorXd pos_min(3), pos_max(3);
  pos_min << -10.0, -9.0, -8.0;
  pos_max << 8.0, 9.0, 10.0;
  Vertex::Vector vertices;
  vertices = createRandomVertices(max_derivative, 100, pos_min, pos_max, 978);

  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  PolynomialOptimization<N> opt(3);
  opt.setupFromVertices(vertices, segment_times);
  opt.solveLinear();

  Segment::Vector segments;
  opt.getSegments(&segments);

  timing::Timer time_analytic("time_extrema_analytic_3", false);
  timing::Timer time_analytic_template_free(
      "time_extrema_analytic_3_template_free", false);
  timing::Timer time_sampling("time_extrema_sampling_3", false);
  int segment_idx = 0;
  for (const Segment& s : segments) {
    std::vector<double> res;
    time_analytic.Start();
    opt.computeSegmentMaximumMagnitudeCandidates<1>(s, 0, s.getTime(), &res);
    time_analytic.Stop();

    std::vector<double> res_template_free;
    std::vector<int> dimensions = {0, 1, 2};
    time_analytic_template_free.Start();
    s.computeMinMaxMagnitudeCandidateTimes(1, 0.0, s.getTime(), dimensions,
                                           &res_template_free);
    time_analytic_template_free.Stop();

    std::vector<double> res_sampling;
    time_sampling.Start();
    opt.computeSegmentMaximumMagnitudeCandidatesBySampling<1>(
        s, 0, s.getTime(), 0.001, &res_sampling);
    time_sampling.Stop();

    constexpr double check_tolerance = 0.01;
    bool success = checkExtrema(res, res_sampling, check_tolerance);
    if (!success) {
      std::cout << "############CHECK XTREMA FAILED: \n";
      std::cout << "segment idx: " << segment_idx << "/" << segments.size()
                << std::endl;

      std::cout << "real eigs for segment with time " << s.getTime() << "s : ";

      for (const double& t : res) std::cout << t << " ";
      std::cout << std::endl;
      std::cout << "sampling: ";
      for (const double& t : res_sampling) std::cout << t << " ";
      std::cout << std::endl;

      std::cout << "vx = [ "
                << s[0].getCoefficients(derivative_order::VELOCITY)
                       .reverse()
                       .format(matlab_format)
                << "];\n";
      std::cout << "vy = [ "
                << s[1].getCoefficients(derivative_order::VELOCITY)
                       .reverse()
                       .format(matlab_format)
                << "];\n";
      std::cout << "vz = [ "
                << s[2].getCoefficients(derivative_order::VELOCITY)
                       .reverse()
                       .format(matlab_format)
                << "];\n";
      std::cout << "t = 0:0.001:" << s.getTime() << "; \n";
    }
    EXPECT_TRUE(success);

    EXPECT_EQ(res.size(), res_template_free.size() - 2);
    for (size_t i = 0; i < res.size(); i++) {
      EXPECT_EQ(res[i], res_template_free[i + 2]);
    }
    ++segment_idx;
  }

  double v_max_ref = getMaximumMagnitude(segments, derivative_order::VELOCITY);
  double a_max_ref =
      getMaximumMagnitude(segments, derivative_order::ACCELERATION);

  timing::Timer time_analytic_v("time_extrema_analytic_v");
  Extremum v_max =
      opt.computeMaximumOfMagnitude<derivative_order::VELOCITY>(nullptr);
  time_analytic_v.Stop();
  timing::Timer time_analytic_a("time_extrema_analytic_a");
  Extremum a_max =
      opt.computeMaximumOfMagnitude<derivative_order::ACCELERATION>(nullptr);
  time_analytic_a.Stop();

  std::cout << "v_max " << v_max << std::endl << "a_max " << a_max << std::endl;

  EXPECT_NEAR(v_max_ref, v_max.value, 0.01);
  EXPECT_NEAR(a_max_ref, a_max.value, 0.01);
}

TEST(MavTrajectoryGeneration,
     PathPlanningUnconstrained_3D_10_segments_nonlinear) {
  Eigen::VectorXd pos_min(3), pos_max(3);
  pos_min << -10.0, -20.0, -10.0;
  pos_max << 10.0, 20.0, 10.0;
  Vertex::Vector vertices;
  vertices = createRandomVertices(max_derivative, 10, pos_min * 0.2,
                                  pos_max * 0.2, 12345);

  const double approximate_v_max = 4.0;
  const double approximate_a_max = 5.0;
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max / 4, approximate_a_max);

  NonlinearOptimizationParameters parameters;
  parameters.max_iterations = 1000;
  parameters.f_rel = 0.05;
  parameters.x_rel = 0.1;
  parameters.time_penalty = 500.0;
  parameters.initial_stepsize_rel = 0.1;
  parameters.inequality_constraint_tolerance = 0.1;
  //  parameters.algorithm = nlopt::GN_ORIG_DIRECT;
  //  parameters.algorithm = nlopt::GN_ORIG_DIRECT_L;
  parameters.algorithm = nlopt::GN_ISRES;
  //  parameters.algorithm = nlopt::LN_COBYLA;

  parameters.random_seed = 12345678;

  int ret;

  timing::Timer timer_setup("setup_3D_10s_nonlinear_time_only");
  PolynomialOptimizationNonLinear<N> opt(3, parameters, true);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt.addMaximumMagnitudeConstraint(derivative_order::VELOCITY,
                                    approximate_v_max);
  opt.addMaximumMagnitudeConstraint(derivative_order::ACCELERATION,
                                    approximate_a_max);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_3D_10s_nonlinear_time_only");
  ret = opt.optimize();
  timer_solve.Stop();

  std::cout << "nlopt1 stopped for reason: " << nlopt::returnValueToString(ret)
            << std::endl;

  timing::Timer timer_setup2("setup_3D_10s_nonlinear_time_and_derivatives");
  PolynomialOptimizationNonLinear<N> opt2(3, parameters, false);
  opt2.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt2.addMaximumMagnitudeConstraint(derivative_order::VELOCITY,
                                     approximate_v_max);
  opt2.addMaximumMagnitudeConstraint(derivative_order::ACCELERATION,
                                     approximate_a_max);
  timer_setup2.Stop();
  timing::Timer timer_solve2("solve_3D_10s_nonlinear_time_and_derivatives");
  ret = opt2.optimize();
  timer_solve2.Stop();

  std::cout << "nlopt2 stopped for reason: " << nlopt::returnValueToString(ret)
            << std::endl;

  Segment::Vector segments1, segments2;
  opt.getPolynomialOptimizationRef().getSegments(&segments1);
  opt2.getPolynomialOptimizationRef().getSegments(&segments2);

  checkPath(vertices, segments1);
  checkPath(vertices, segments2);
  double v_max = getMaximumMagnitude(segments1, derivative_order::VELOCITY);
  double a_max = getMaximumMagnitude(segments1, derivative_order::ACCELERATION);
  std::cout << "v_max 1: " << v_max << " a_max 1: " << a_max << std::endl;
  EXPECT_LT(v_max, approximate_v_max * 1.5);
  EXPECT_LT(a_max, approximate_a_max * 1.5);

  EXPECT_TRUE(checkCost(opt.getPolynomialOptimizationRef().computeCost(),
                        segments1, derivative_to_optimize, 0.1));

  v_max = getMaximumMagnitude(segments2, derivative_order::VELOCITY);
  a_max = getMaximumMagnitude(segments2, derivative_order::ACCELERATION);
  std::cout << "v_max 2: " << v_max << " a_max 2: " << a_max << std::endl;

  EXPECT_LT(v_max, approximate_v_max * 1.5);
  EXPECT_LT(a_max, approximate_a_max * 1.5);

  EXPECT_TRUE(checkCost(opt2.getPolynomialOptimizationRef().computeCost(),
                        segments2, derivative_to_optimize, 0.1));
}

TEST(MavTrajectoryGeneration, 2_vertices_setup) {
  const int kDim = 1;
  // Create a known 1D spline.
  Vertex start_vertex(kDim);
  const double kStartX = 0.0;
  start_vertex.addConstraint(derivative_order::POSITION, kStartX);
  start_vertex.addConstraint(derivative_order::VELOCITY, 0.0);
  start_vertex.addConstraint(derivative_order::ACCELERATION, 0.0);
  start_vertex.addConstraint(derivative_order::JERK, 0.0);
  start_vertex.addConstraint(derivative_order::SNAP, 0.0);

  Vertex goal_vertex = start_vertex;
  const double kGoalX = 5.0;
  goal_vertex.addConstraint(derivative_order::POSITION, kGoalX);

  const double kVMax = 2.0;
  const double kSegmentTime = std::fabs(kGoalX - kStartX) * 2.0 / kVMax;

  const int kDerivativeToOptimize = derivative_order::SNAP;
  const int kNumCoefficients = 10;

  // Setup optimization with two vertices.
  PolynomialOptimization<kNumCoefficients> opt(kDim);
  Vertex::Vector vertices{start_vertex, goal_vertex};
  std::vector<double> segment_times{kSegmentTime};
  opt.setupFromVertices(vertices, segment_times, kDerivativeToOptimize);
  opt.solveLinear();

  Segment::Vector segments;
  opt.getSegments(&segments);
  checkPath(vertices, segments);

  // Matlab solution:
  Eigen::VectorXd matlab_coeffs(kNumCoefficients);
  matlab_coeffs << -0.000000000000004, 0.000000000000004, -0.000000000000006,
      0.000000000000003, -0.000000000000001, 0.201600000000015,
      -0.134400000000012, 0.034560000000004, -0.004032000000000,
      0.000179200000000;
  // Solution with two vertices:
  Eigen::VectorXd coeffs = segments[0].getPolynomialsRef()[0].getCoefficients();

  // The matlab solution is only approximately valid, e.g., the first
  // coefficient is not equal 0.0.
  CHECK_EIGEN_MATRIX_EQUAL_DOUBLE(matlab_coeffs, coeffs);
}

// Test 2 vertices setup
TEST(MavTrajectoryGeneration, 2_vertices_rand) {
  const int kMaxDerivative = derivative_order::ACCELERATION;
  const size_t kNumSegments = 1;
  Eigen::VectorXd min_pos, max_pos;
  min_pos = Eigen::Vector3d::Constant(-50.0);
  max_pos = -min_pos;
  const int kSeed = 12345;
  const size_t kNumSetups = 1e2;
  for (size_t i = 0; i < kNumSetups; i++) {
    Vertex::Vector vertices;
    vertices = createRandomVertices(kMaxDerivative, kNumSegments, min_pos,
                                    max_pos, kSeed + i);

    const double kApproxVMax = 3.0;
    const double kApproxAMax = 5.0;
    std::vector<double> segment_times =
        estimateSegmentTimes(vertices, kApproxVMax, kApproxAMax);

    PolynomialOptimization<N> opt(3);
    opt.setupFromVertices(vertices, segment_times);
    opt.solveLinear();

    Segment::Vector segments;
    opt.getSegments(&segments);

    checkPath(vertices, segments);
  }
}

// Test unpacking and repacking constraints between [d_f; d_p] and p.
TEST(MavTrajectoryGeneration, ConstraintPacking) {
  const int kMaxDerivative = derivative_order::JERK;
  const size_t kNumSegments = 5;
  Eigen::VectorXd min_pos, max_pos;
  min_pos = Eigen::Vector3d::Constant(-50.0);
  max_pos = -min_pos;
  const int kSeed = 12345;
  const size_t kNumSetups = 100;

  for (size_t i = 0; i < kNumSetups; i++) {
    Vertex::Vector vertices;
    vertices = createRandomVertices(kMaxDerivative, kNumSegments, min_pos,
                                    max_pos, kSeed + i);

    const double kApproxVMax = 3.0;
    const double kApproxAMax = 5.0;
    std::vector<double> segment_times =
        estimateSegmentTimes(vertices, kApproxVMax, kApproxAMax);

    PolynomialOptimization<N> opt(3);
    opt.setupFromVertices(vertices, segment_times);
    opt.solveLinear();

    Segment::Vector segments;
    opt.getSegments(&segments);

    std::vector<Eigen::VectorXd> fixed_constraints;
    std::vector<Eigen::VectorXd> free_constraints;

    opt.getFixedConstraints(&fixed_constraints);
    opt.getFreeConstraints(&free_constraints);

    // Get the mapping matrices.
    Eigen::MatrixXd M, A_inv, A, M_pinv;
    opt.getM(&M);
    opt.getAInverse(&A_inv);
    opt.getA(&A);
    opt.getMpinv(&M_pinv);

    ASSERT_EQ(fixed_constraints.size(), 3);
    ASSERT_EQ(free_constraints.size(), 3);

    // For each dimension... We can test that [df;dp] -> p -> [df;dp].
    for (int i = 0; i < 3; ++i) {
      Eigen::VectorXd d_all_ordered(fixed_constraints[i].size() +
                                    free_constraints[i].size());
      d_all_ordered << fixed_constraints[i], free_constraints[i];
      Eigen::VectorXd p = A_inv * M * d_all_ordered;
      Eigen::VectorXd d_unordered = A * p;
      Eigen::VectorXd d_reordered = M_pinv * d_unordered;
      EXPECT_TRUE(EIGEN_MATRIX_NEAR(d_all_ordered, d_reordered, 1e-6));

      // Now also check for each segment.
      for (size_t j = 0; j < segments.size(); ++j) {
        Eigen::VectorXd p_seg = segments[j][i].getCoefficients(0);
        EXPECT_TRUE(EIGEN_MATRIX_NEAR(p_seg, p.segment<N>(j * N), 1e-6));
      }
    }
  }
}

void createTestPolynomials() {
  Vertex::Vector vertices;
  vertices = createRandomVertices1D(max_derivative, 100, -50, 50, 12345);
  const double approximate_v_max = 3.0;
  const double approximate_a_max = 5.0;

  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  PolynomialOptimization<N> opt(1);
  opt.setupFromVertices(vertices, segment_times);
  opt.solveLinear();

  Segment::Vector segments;
  opt.getSegments(&segments);

  int s_idx = 1;
  for (const Segment& s : segments) {
    std::cout << "s(" << s_idx << ", :) = "
              << s[0].getCoefficients(derivative_order::POSITION)
                     .reverse()
                     .format(matlab_format)
              << ";\n";
    ++s_idx;
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  timing::Timing::Print(std::cout);
  //  createTestPolynomials();

  return result;
}
