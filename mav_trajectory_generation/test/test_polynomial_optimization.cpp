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
#include "mav_trajectory_generation/test_utils.h"
#include "mav_trajectory_generation/timing.h"

using namespace mav_trajectory_generation;

Eigen::IOFormat matlab_format(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "[",
                              "]");
const int N = 10;

struct OptimizationParams {
  int D;
  int max_derivative;
  int num_segments;
  int seed;
  double pos_bounds;
  double v_max;
  double a_max;
};

std::ostream& operator<<(std::ostream& stream, const OptimizationParams& val) {
  stream << "D: " << val.D << " max_deriv: " << val.max_derivative
         << " num_seg: " << val.num_segments << " seed: " << val.seed
         << " pos_bounds: " << val.pos_bounds << " v_max: " << val.v_max
         << " a_max: " << val.a_max;
  return stream;
}

// Parameters to test: N, max_derivative, num_segments, seed
class PolynomialOptimizationTests
    : public ::testing::TestWithParam<OptimizationParams> {
 public:
  void SetUp() override {
    params_ = GetParam();

    // Unpack params.
    D = params_.D;
    max_derivative = params_.max_derivative;
    v_max = params_.v_max;
    a_max = params_.a_max;

    Eigen::VectorXd pos_min(D), pos_max(D);
    pos_min.setConstant(-params_.pos_bounds);
    pos_max.setConstant(params_.pos_bounds);

    vertices_ =
        createRandomVertices(getHighestDerivativeFromN(N), params_.num_segments,
                             pos_min, pos_max, params_.seed);
  }

  // Helper checking functions.
  void checkPath(const Vertex::Vector& vertices,
                 const std::vector<Segment>& segments) const;
  bool checkCost(double cost_to_check, const Trajectory& trajectory,
                 size_t derivative, double relative_tolerance) const;
  void getMaxVelocityAndAccelerationAnalytical(const Trajectory& trajectory,
                                               double* v_max,
                                               double* a_max) const;
  void getMaxVelocityAndAccelerationNumerical(const Trajectory& trajectory,
                                              double* v_max,
                                              double* a_max) const;
  bool checkExtrema(const std::vector<double>& testee,
                    const std::vector<double>& reference,
                    double tol = 0.01) const;

  std::string getSuffix() const {
    std::ostringstream sstream;
    sstream << "_" << D << "D_" << params_.num_segments << "s";
    return sstream.str();
  }

 protected:
  OptimizationParams params_;

  // Unfold params to be a bit simpler.
  int D;
  int max_derivative;
  double v_max;
  double a_max;

  Vertex::Vector vertices_;
};

void PolynomialOptimizationTests::checkPath(
    const Vertex::Vector& vertices,
    const std::vector<Segment>& segments) const {
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
          << "\nsegment:\n"
          << segment << segment_derivative.str();
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
          << "\nsegment:\n"
          << segment << segment_derivative.str();
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

bool PolynomialOptimizationTests::checkCost(double cost_to_check,
                                            const Trajectory& trajectory,
                                            size_t derivative,
                                            double relative_tolerance) const {
  CHECK_GE(derivative, size_t(0));
  CHECK(relative_tolerance >= 0.0 && relative_tolerance <= 1.0);
  const double sampling_interval = 0.001;
  double cost_numeric =
      computeCostNumeric(trajectory, derivative, sampling_interval);

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

void PolynomialOptimizationTests::getMaxVelocityAndAccelerationAnalytical(
    const Trajectory& trajectory, double* v_max, double* a_max) const {
  std::vector<int> dimensions(D);  // Evaluate in whatever dimensions we have.
  std::iota(dimensions.begin(), dimensions.end(), 0);

  mav_trajectory_generation::Extremum v_min_traj, v_max_traj, a_min_traj,
      a_max_traj;

  trajectory.computeMinMaxMagnitude(
      mav_trajectory_generation::derivative_order::VELOCITY, dimensions,
      &v_min_traj, &v_max_traj);
  trajectory.computeMinMaxMagnitude(
      mav_trajectory_generation::derivative_order::ACCELERATION, dimensions,
      &a_min_traj, &a_max_traj);

  *v_max = v_max_traj.value;
  *a_max = a_max_traj.value;
}

void PolynomialOptimizationTests::getMaxVelocityAndAccelerationNumerical(
    const Trajectory& trajectory, double* v_max, double* a_max) const {
  *v_max = getMaximumMagnitude(
      trajectory, mav_trajectory_generation::derivative_order::VELOCITY);
  *a_max = getMaximumMagnitude(
      trajectory, mav_trajectory_generation::derivative_order::ACCELERATION);
}

bool PolynomialOptimizationTests::checkExtrema(
    const std::vector<double>& testee, const std::vector<double>& reference,
    double tol) const {
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

TEST_P(PolynomialOptimizationTests, VertexGeneration) {
  Eigen::VectorXd pos_min(D), pos_max(D);
  pos_min.setConstant(-params_.pos_bounds);
  pos_max.setConstant(params_.pos_bounds);

  EXPECT_EQ(vertices_.front().getNumberOfConstraints(),
            getHighestDerivativeFromN(N) + 1);
  EXPECT_EQ(vertices_.back().getNumberOfConstraints(),
            getHighestDerivativeFromN(N) + 1);

  for (const Vertex& v : vertices_) {
    EXPECT_TRUE(v.hasConstraint(derivative_order::POSITION));
    Eigen::VectorXd c;
    v.getConstraint(derivative_order::POSITION, &c);
    for (int i = 0; i < D; ++i) {
      EXPECT_LE(c[i], pos_max[i]);
      EXPECT_GE(c[i], pos_min[i]);
    }
  }
}

TEST_P(PolynomialOptimizationTests, UnconstrainedLinearEstimateSegmentTimes) {
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices_, v_max, a_max);

  timing::Timer timer_setup("setup" + getSuffix());
  PolynomialOptimization<N> opt(D);
  opt.setupFromVertices(vertices_, segment_times, max_derivative);
  timer_setup.Stop();
  timing::Timer timer_solve("solve_linear" + getSuffix());
  opt.solveLinear();
  timer_solve.Stop();

  Segment::Vector segments;
  opt.getSegments(&segments);
  Trajectory trajectory;
  opt.getTrajectory(&trajectory);

  std::cout << "Base coefficients: "
            << Polynomial::base_coefficients_.block(3, 0, 1, N) << std::endl;

  checkPath(vertices_, segments);
  double v_max_traj =
      getMaximumMagnitude(trajectory, derivative_order::VELOCITY);
  double a_max_traj =
      getMaximumMagnitude(trajectory, derivative_order::ACCELERATION);
  std::cout << "v_max: " << v_max_traj << " a_max: " << a_max_traj << std::endl;

  timing::Timer timer_cost("cost" + getSuffix());
  double cost = opt.computeCost();
  timer_cost.Stop();
  std::cout << "cost: " << cost << std::endl;

  EXPECT_LT(v_max_traj, v_max * 2.5);
  EXPECT_LT(a_max_traj, a_max * 2.5);
  EXPECT_TRUE(checkCost(opt.computeCost(), trajectory, max_derivative, 0.1));
}

TEST_P(PolynomialOptimizationTests, ExtremaOfMagnitude) {
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices_, v_max, a_max);
  constexpr int kDerivative = derivative_order::VELOCITY;

  PolynomialOptimization<N> opt(D);
  opt.setupFromVertices(vertices_, segment_times, max_derivative);
  opt.solveLinear();

  Segment::Vector segments;
  opt.getSegments(&segments);

  Trajectory trajectory;
  opt.getTrajectory(&trajectory);

  std::vector<int> dimensions(D);
  std::iota(dimensions.begin(), dimensions.end(), 0);

  timing::Timer time_analytic("time_extrema_analytic" + getSuffix(), false);
  timing::Timer time_sampling("time_extrema_sampling" + getSuffix(), false);
  int segment_idx = 0;
  for (const Segment& s : segments) {
    std::vector<double> res;
    time_analytic.Start();
    EXPECT_TRUE(opt.computeSegmentMaximumMagnitudeCandidates(
        kDerivative, s, 0, s.getTime(), &res));
    time_analytic.Stop();

    std::vector<double> res_sampling;
    time_sampling.Start();
    const double dt = 0.01;
    opt.computeSegmentMaximumMagnitudeCandidatesBySampling<kDerivative>(
        s, 0, s.getTime(), dt, &res_sampling);
    time_sampling.Stop();

    const double check_tolerance = 2 * dt; // Nyquist resolution
    bool success = checkExtrema(res_sampling, res, check_tolerance);
    if (!success) {
      std::cout << "############CHECK XTREMA FAILED: \n";
      std::cout << "segment idx: " << segment_idx << "/" << segments.size()
                << " time: " << s.getTime() << std::endl
                << "segment: " << s << std::endl;

      std::cout << "analytically found: ";
      for (const double& t : res) {
        std::cout << t << " ";
      }
      std::cout << std::endl;
      std::cout << "sampling: ";
      for (const double& t : res_sampling) {
        std::cout << t << " ";
      }
      std::cout << std::endl;
      for (int i = 0; i < D; i++) {
        std::cout
            << "vx" << i << " = "
            << s[i].getCoefficients(kDerivative).reverse().format(matlab_format)
            << ";\n";
      }
      std::cout << "t = 0:0.001:" << s.getTime() << "; \n";
    }
    EXPECT_TRUE(success);

    ++segment_idx;
  }

  double v_max_ref =
      getMaximumMagnitude(trajectory, derivative_order::VELOCITY);
  double a_max_ref =
      getMaximumMagnitude(trajectory, derivative_order::ACCELERATION);

  timing::Timer time_analytic_v("time_extrema_analytic_v" + getSuffix());
  Extremum v_max_opt =
      opt.computeMaximumOfMagnitude<derivative_order::VELOCITY>(nullptr);
  time_analytic_v.Stop();
  timing::Timer time_analytic_a("time_extrema_analytic_a" + getSuffix());
  Extremum a_max_opt =
      opt.computeMaximumOfMagnitude<derivative_order::ACCELERATION>(nullptr);
  time_analytic_a.Stop();

  Extremum v_min_traj, v_max_traj, a_min_traj, a_max_traj;

  trajectory.computeMinMaxMagnitude(derivative_order::VELOCITY, dimensions,
                                    &v_min_traj, &v_max_traj);
  trajectory.computeMinMaxMagnitude(derivative_order::ACCELERATION, dimensions,
                                    &a_min_traj, &a_max_traj);

  EXPECT_NEAR(v_max_ref, v_max_opt.value, 0.01);
  EXPECT_NEAR(a_max_ref, a_max_opt.value, 0.01);

  EXPECT_NEAR(v_max_ref, v_max_traj.value, 0.01);
  EXPECT_NEAR(a_max_ref, a_max_traj.value, 0.01);
}

TEST_P(PolynomialOptimizationTests, UnconstrainedNonlinear) {
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices_, v_max, a_max);

  NonlinearOptimizationParameters parameters;
  parameters.max_iterations = 3000;
  parameters.f_rel = 0.05;
  parameters.x_rel = 0.1;
  parameters.time_penalty = 500.0;
  parameters.initial_stepsize_rel = 0.1;
  parameters.inequality_constraint_tolerance = 0.1;
  // For global methods (GN_): Non-inf boundary conditions for all
  // optimization parameters needed. Change bounds to non-inf values.
  // Otherwise infinite compile time or "error: nlopt invalid argument".
  //  parameters.algorithm = nlopt::GN_ORIG_DIRECT;
  //  parameters.algorithm = nlopt::GN_ORIG_DIRECT_L;
  //  parameters.algorithm = nlopt::GN_ISRES;
  //  parameters.algorithm = nlopt::LN_COBYLA;
  parameters.algorithm = nlopt::LN_BOBYQA;
  // parameters.algorithm = nlopt::LN_SBPLX;

  parameters.random_seed = 12345678;

  int ret;

  timing::Timer timer_setup("setup_nonlinear_time_only" + getSuffix());
  parameters.time_alloc_method = NonlinearOptimizationParameters::kSquaredTime;
  PolynomialOptimizationNonLinear<N> opt(D, parameters);
  opt.setupFromVertices(vertices_, segment_times, max_derivative);
  opt.addMaximumMagnitudeConstraint(derivative_order::VELOCITY, v_max);
  opt.addMaximumMagnitudeConstraint(derivative_order::ACCELERATION, a_max);
  opt.solveLinear();
  double initial_cost = opt.getTotalCostWithSoftConstraints();
  timer_setup.Stop();
  timing::Timer timer_solve("solve_nonlinear_time_only" + getSuffix());
  ret = opt.optimize();
  timer_solve.Stop();
  double final_cost = opt.getTotalCostWithSoftConstraints();

  EXPECT_NE(ret, nlopt::FAILURE);
  EXPECT_NE(ret, nlopt::INVALID_ARGS);

  EXPECT_LE(final_cost, initial_cost * 1.1);

  std::cout << "nlopt1 stopped for reason: " << nlopt::returnValueToString(ret)
            << std::endl;

  timing::Timer timer_setup2("setup_nonlinear_time_and_derivatives" +
                             getSuffix());
  parameters.time_alloc_method =
      NonlinearOptimizationParameters::kSquaredTimeAndConstraints;
  PolynomialOptimizationNonLinear<N> opt2(D, parameters);
  opt2.setupFromVertices(vertices_, segment_times, max_derivative);
  opt2.addMaximumMagnitudeConstraint(derivative_order::VELOCITY, v_max);
  opt2.addMaximumMagnitudeConstraint(derivative_order::ACCELERATION, a_max);
  timer_setup2.Stop();
  timing::Timer timer_solve2("solve_nonlinear_time_and_derivatives" +
                             getSuffix());
  ret = opt2.optimize();
  timer_solve2.Stop();
  double final_cost2 = opt2.getTotalCostWithSoftConstraints();

  EXPECT_NE(ret, nlopt::FAILURE);
  EXPECT_NE(ret, nlopt::INVALID_ARGS);

  EXPECT_LE(final_cost2, initial_cost * 1.1);
  EXPECT_LE(final_cost2, final_cost * 1.5);

  std::cout << "nlopt2 stopped for reason: " << nlopt::returnValueToString(ret)
            << std::endl;

  Segment::Vector segments1, segments2;
  Trajectory trajectory1, trajectory2;
  opt.getPolynomialOptimizationRef().getSegments(&segments1);
  opt2.getPolynomialOptimizationRef().getSegments(&segments2);
  opt.getPolynomialOptimizationRef().getTrajectory(&trajectory1);
  opt2.getPolynomialOptimizationRef().getTrajectory(&trajectory2);

  checkPath(vertices_, segments1);
  checkPath(vertices_, segments2);
  double v_max_1 = getMaximumMagnitude(trajectory1, derivative_order::VELOCITY);
  double a_max_1 =
      getMaximumMagnitude(trajectory1, derivative_order::ACCELERATION);
  std::cout << "v_max_1: " << v_max_1 << " a_max_1: " << a_max_1 << std::endl;
  EXPECT_LE(v_max_1, v_max * 1.5);
  EXPECT_LE(a_max_1, a_max * 1.5);

  EXPECT_TRUE(checkCost(opt.getPolynomialOptimizationRef().computeCost(),
                        trajectory1, max_derivative, 0.1));

  double v_max_2 = getMaximumMagnitude(trajectory2, derivative_order::VELOCITY);
  double a_max_2 =
      getMaximumMagnitude(trajectory2, derivative_order::ACCELERATION);
  std::cout << "v_max_2: " << v_max_2 << " a_max_2: " << a_max_2 << std::endl;

  EXPECT_LE(v_max_2, v_max * 1.5);
  EXPECT_LE(a_max_2, a_max * 1.5);

  EXPECT_TRUE(checkCost(opt2.getPolynomialOptimizationRef().computeCost(),
                        trajectory2, max_derivative, 0.1));
}

// Test unpacking and repacking constraints between [d_f; d_p] and p.
TEST_P(PolynomialOptimizationTests, ConstraintPacking) {
  const int kMaxDerivative = max_derivative;
  const size_t kNumSegments = params_.num_segments;
  Eigen::VectorXd min_pos, max_pos;
  min_pos = Eigen::VectorXd::Constant(D, -50.0);
  max_pos = -min_pos;
  const int kSeed = 12345;
  const size_t kNumSetups = 5;

  for (size_t i = 0; i < kNumSetups; i++) {
    Vertex::Vector vertices;
    vertices = createRandomVertices(kMaxDerivative, kNumSegments, min_pos,
                                    max_pos, kSeed + i);

    const double kApproxVMax = 3.0;
    const double kApproxAMax = 5.0;
    std::vector<double> segment_times =
        estimateSegmentTimes(vertices, kApproxVMax, kApproxAMax);

    PolynomialOptimization<N> opt(D);
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

    ASSERT_EQ(fixed_constraints.size(), D);
    ASSERT_EQ(free_constraints.size(), D);

    // For each dimension... We can test that [df;dp] -> p -> [df;dp].
    for (int i = 0; i < D; ++i) {
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

TEST_P(PolynomialOptimizationTests, TimeAllocation) {
  std::vector<double> segment_times_ramp, segment_times_nfabian;
  double time_factor = 1.0;
  segment_times_ramp =
      estimateSegmentTimesVelocityRamp(vertices_, v_max, a_max);
  segment_times_nfabian = estimateSegmentTimesNfabian(vertices_, v_max, a_max);

  EXPECT_EQ(segment_times_ramp.size(), segment_times_nfabian.size());
  EXPECT_EQ(params_.num_segments, segment_times_nfabian.size());

  for (size_t i = 0; i < params_.num_segments; i++) {
    EXPECT_LT(0.0, segment_times_ramp[i]);
    EXPECT_LT(0.0, segment_times_nfabian[i]);
    EXPECT_GT(1e5, segment_times_ramp[i]);
    EXPECT_GT(1e5, segment_times_nfabian[i]);
  }

  PolynomialOptimization<N> opt_ramp(D);
  opt_ramp.setupFromVertices(vertices_, segment_times_ramp, max_derivative);
  opt_ramp.solveLinear();
  Trajectory trajectory_ramp;
  opt_ramp.getTrajectory(&trajectory_ramp);

  PolynomialOptimization<N> opt_nfabian(D);
  opt_nfabian.setupFromVertices(vertices_, segment_times_nfabian,
                                max_derivative);
  opt_nfabian.solveLinear();
  Trajectory trajectory_nfabian;
  opt_nfabian.getTrajectory(&trajectory_nfabian);

  double v_max_ramp, a_max_ramp, v_max_nfabian, a_max_nfabian;
  getMaxVelocityAndAccelerationAnalytical(trajectory_ramp, &v_max_ramp,
                                          &a_max_ramp);
  getMaxVelocityAndAccelerationAnalytical(trajectory_nfabian, &v_max_nfabian,
                                          &a_max_nfabian);

  EXPECT_LT(v_max_ramp, v_max * 2.5);
  EXPECT_LT(a_max_ramp, a_max * 2.5);
  EXPECT_LT(v_max_nfabian, v_max * 2.5);
  EXPECT_LT(a_max_nfabian, a_max * 2.5);
}

TEST_P(PolynomialOptimizationTests, TimeScaling) {
  std::vector<double> segment_times;
  double time_factor = 1.0;
  constexpr double kTolerance = 1e-3;

  // Allocate using the ramp.
  segment_times =
      estimateSegmentTimesVelocityRamp(vertices_, v_max, a_max, time_factor);
  EXPECT_EQ(segment_times.size(), params_.num_segments);

  // Now re-scale the times, using various non-linear optimization techniques.
  mav_trajectory_generation::NonlinearOptimizationParameters nlopt_parameters;
  nlopt_parameters.algorithm = nlopt::LD_LBFGS;
  nlopt_parameters.time_alloc_method = mav_trajectory_generation::
      NonlinearOptimizationParameters::kMellingerOuterLoop;
  nlopt_parameters.print_debug_info_time_allocation = false;
  nlopt_parameters.random_seed = 12345678;
  mav_trajectory_generation::PolynomialOptimizationNonLinear<N> nlopt(
      D, nlopt_parameters);

  nlopt.setupFromVertices(vertices_, segment_times, max_derivative);
  nlopt.addMaximumMagnitudeConstraint(
      mav_trajectory_generation::derivative_order::VELOCITY, v_max);
  nlopt.addMaximumMagnitudeConstraint(
      mav_trajectory_generation::derivative_order::ACCELERATION, a_max);
  nlopt.solveLinear();

  Trajectory trajectory;
  double v_max_traj, a_max_traj, v_max_num, a_max_num;
  nlopt.getTrajectory(&trajectory);

  double initial_cost = nlopt.getCost();
  getMaxVelocityAndAccelerationAnalytical(trajectory, &v_max_traj, &a_max_traj);
  getMaxVelocityAndAccelerationNumerical(trajectory, &v_max_num, &a_max_num);

  EXPECT_LE(v_max_traj, v_max * 2.5);
  EXPECT_LE(a_max_traj, a_max * 2.5);

  std::cout << "Starting v max: " << v_max_traj << " (" << v_max_num
            << ") a max: " << a_max_traj << " (" << a_max_num << ")"
            << std::endl;

  EXPECT_NEAR(v_max_traj, v_max_num, kTolerance);
  EXPECT_NEAR(a_max_traj, a_max_num, kTolerance);

  // Scaling of segment times
  nlopt.scaleSegmentTimesWithViolation();
  nlopt.getTrajectory(&trajectory);
  double scaled_cost = nlopt.getCost();
  getMaxVelocityAndAccelerationAnalytical(trajectory, &v_max_traj, &a_max_traj);
  getMaxVelocityAndAccelerationNumerical(trajectory, &v_max_num, &a_max_num);

  std::cout << "Scaled v max: " << v_max_traj << " a max: " << a_max_traj
            << std::endl;

  EXPECT_NEAR(v_max_traj, v_max_num, kTolerance);
  EXPECT_NEAR(a_max_traj, a_max_num, kTolerance);

  EXPECT_LE(v_max_traj, v_max + kTolerance);
  EXPECT_LE(a_max_traj, a_max + kTolerance);

  nlopt.optimize();
  double mellinger_cost = nlopt.getCost();
  nlopt.getTrajectory(&trajectory);
  getMaxVelocityAndAccelerationAnalytical(trajectory, &v_max_traj, &a_max_traj);
  getMaxVelocityAndAccelerationNumerical(trajectory, &v_max_num, &a_max_num);

  std::cout << "Mellinger v max: " << v_max_traj << " a max: " << a_max_traj
            << std::endl;

  EXPECT_NEAR(v_max_traj, v_max_num, kTolerance);
  EXPECT_NEAR(a_max_traj, a_max_num, kTolerance);

  EXPECT_LE(v_max_traj, v_max + kTolerance);
  EXPECT_LE(a_max_traj, a_max + kTolerance);

  EXPECT_LE(scaled_cost, initial_cost);
  EXPECT_LE(mellinger_cost, scaled_cost);
}

TEST_P(PolynomialOptimizationTests, TimeScalingInTrajectory) {
  constexpr double kTolerance = 1e-3;

  std::vector<double> segment_times =
      estimateSegmentTimesVelocityRamp(vertices_, v_max * 2.0, a_max * 2.0);
  EXPECT_EQ(segment_times.size(), params_.num_segments);

  PolynomialOptimization<N> opt(D);
  opt.setupFromVertices(vertices_, segment_times, max_derivative);
  opt.solveLinear();

  Trajectory trajectory;
  opt.getTrajectory(&trajectory);

  double v_max_traj, a_max_traj, v_max_num, a_max_num;
  getMaxVelocityAndAccelerationAnalytical(trajectory, &v_max_traj, &a_max_traj);
  getMaxVelocityAndAccelerationNumerical(trajectory, &v_max_num, &a_max_num);

  std::cout << "Starting v max: " << v_max_traj << " (" << v_max_num
            << ") a max: " << a_max_traj << " (" << a_max_num << ")"
            << std::endl;

  EXPECT_LE(v_max_traj, v_max * 5.0);
  EXPECT_LE(a_max_traj, a_max * 5.0);

  EXPECT_NEAR(v_max_traj, v_max_num, kTolerance);
  EXPECT_NEAR(a_max_traj, a_max_num, kTolerance);

  // Scaling of segment times
  trajectory.scaleSegmentTimesToMeetConstraints(v_max, a_max);
  getMaxVelocityAndAccelerationAnalytical(trajectory, &v_max_traj, &a_max_traj);
  getMaxVelocityAndAccelerationNumerical(trajectory, &v_max_num, &a_max_num);
  std::cout << "Scaled v max: " << v_max_traj << " (" << v_max_num
            << ") a max: " << a_max_traj << " (" << a_max_num << ")"
            << std::endl;

  EXPECT_NEAR(v_max_traj, v_max_num, kTolerance);
  EXPECT_NEAR(a_max_traj, a_max_num, kTolerance);

  EXPECT_LE(v_max_traj, v_max + kTolerance);
  EXPECT_LE(a_max_traj, a_max + kTolerance);
}

TEST_P(PolynomialOptimizationTests, AMatrixInversion) {
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

TEST_P(PolynomialOptimizationTests, TwoVerticesSetup) {
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

// Set up some common cases.
OptimizationParams segment_1_dim_1 = {1 /* D */,
                                      derivative_order::SNAP,
                                      1 /* num_segments*/,
                                      100 /* seed */,
                                      10.0 /* pos_max */,
                                      3.0 /* v_max */,
                                      5.0 /* a_mav */};

OptimizationParams segment_10_dim_1 = {1 /* D */,
                                       derivative_order::SNAP,
                                       10 /* num_segments*/,
                                       102 /* seed */,
                                       10.0 /* pos_max */,
                                       3.0 /* v_max */,
                                       5.0 /* a_mav */};

OptimizationParams segment_50_dim_1 = {1 /* D */,
                                       derivative_order::SNAP,
                                       50 /* num_segments*/,
                                       103 /* seed */,
                                       10.0 /* pos_max */,
                                       3.0 /* v_max */,
                                       5.0 /* a_mav */};

OptimizationParams segment_1_dim_3 = {3 /* D */,
                                      derivative_order::SNAP,
                                      1 /* num_segments*/,
                                      104 /* seed */,
                                      10.0 /* pos_max */,
                                      3.0 /* v_max */,
                                      5.0 /* a_mav */};

OptimizationParams segment_10_dim_3 = {3 /* D */,
                                       derivative_order::SNAP,
                                       10 /* num_segments*/,
                                       105 /* seed */,
                                       10.0 /* pos_max */,
                                       3.0 /* v_max */,
                                       5.0 /* a_mav */};

OptimizationParams segment_50_dim_3 = {3 /* D */,
                                       derivative_order::SNAP,
                                       50 /* num_segments*/,
                                       106 /* seed */,
                                       10.0 /* pos_max */,
                                       3.0 /* v_max */,
                                       5.0 /* a_mav */};

OptimizationParams deriv_accel_1 = {1 /* D */,
                                    derivative_order::ACCELERATION,
                                    5 /* num_segments*/,
                                    107 /* seed */,
                                    10.0 /* pos_max */,
                                    1.0 /* v_max */,
                                    2.0 /* a_mav */};

OptimizationParams deriv_accel_3_1 = {3 /* D */,
                                      derivative_order::ACCELERATION,
                                      1 /* num_segments*/,
                                      108 /* seed */,
                                      10.0 /* pos_max */,
                                      1.0 /* v_max */,
                                      2.0 /* a_mav */};

OptimizationParams deriv_accel = {3 /* D */,
                                  derivative_order::ACCELERATION,
                                  5 /* num_segments*/,
                                  109 /* seed */,
                                  10.0 /* pos_max */,
                                  1.0 /* v_max */,
                                  2.0 /* a_mav */};
OptimizationParams deriv_jerk = {3 /* D */,
                                 derivative_order::JERK,
                                 5 /* num_segments*/,
                                 110 /* seed */,
                                 10.0 /* pos_max */,
                                 1.0 /* v_max */,
                                 2.0 /* a_mav */};

INSTANTIATE_TEST_CASE_P(OneDimension, PolynomialOptimizationTests,
                        ::testing::Values(segment_1_dim_1, segment_10_dim_1,
                                          segment_50_dim_1));

INSTANTIATE_TEST_CASE_P(ThreeDimensions, PolynomialOptimizationTests,
                        ::testing::Values(segment_1_dim_3, segment_10_dim_3,
                                          segment_50_dim_3));

INSTANTIATE_TEST_CASE_P(Derivatives, PolynomialOptimizationTests,
                        ::testing::Values(deriv_accel_1, deriv_accel_3_1,
                                          deriv_accel, deriv_jerk));

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  timing::Timing::Print(std::cout);

  return result;
}
