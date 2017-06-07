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
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include <eigen-checks/gtest.h>

#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation/segment.h>
#include <mav_trajectory_generation/test_utils.h>
#include <mav_trajectory_generation/timing.h>
#include <mav_trajectory_generation/vertex.h>

#include "mav_trajectory_generation_ros/feasibility_analytic.h"
#include "mav_trajectory_generation_ros/feasibility_recursive.h"
#include "mav_trajectory_generation_ros/feasibility_sampling.h"
#include "mav_trajectory_generation_ros/feasibility_base.h"

using namespace mav_trajectory_generation;

void writeResultsToFile(const std::string& file_name,
                        const std::vector<InputFeasibilityResult>& results) {
  std::ofstream file;
  file.open(file_name, std::ofstream::out | std::ofstream::trunc);
  if (file.is_open()) {
    for (const InputFeasibilityResult& result : results) {
      file << std::to_string(result) + "\n";
    }
    file.close();
  } else {
    std::cout << "Unable to open file " << file_name << "." << std::endl;
  }
}

TEST(FeasibilityTest, CompareFeasibilityTests) {
  std::srand(1234567);
  // Create random segments.
  const int kNumSegments = 1e3;
  const int kN = 12;
  const int kD = 4;
  Segment::Vector segments(kNumSegments, Segment(kN, kD));
  for (size_t i = 0; i < kNumSegments; i++) {
    // Position segment.
    const Eigen::VectorXd kMinPos = Eigen::VectorXd::Constant(3, -5.0);
    const Eigen::VectorXd kMaxPos = -kMinPos;
    Vertex::Vector vertices_pos =
        createRandomVertices(derivative_order::SNAP, 1, kMinPos, kMaxPos);

    // Add some random start and end velocity to make segment more dynamic.
    const double kAvgVMin = 0.5;
    const double kAvgVMax = 2.0;
    Eigen::Vector3d direction_start, direction_end;
    direction_start.x() = createRandomDouble(0.01, 1.0);
    direction_start.y() = createRandomDouble(0.01, 1.0);
    direction_start.z() = createRandomDouble(0.01, 1.0);
    direction_start.normalize();
    vertices_pos[0].addConstraint(
        derivative_order::VELOCITY,
        direction_start * createRandomDouble(0.0, kAvgVMax));

    direction_end.x() = createRandomDouble(0.01, 1.0);
    direction_end.y() = createRandomDouble(0.01, 1.0);
    direction_end.z() = createRandomDouble(0.01, 1.0);
    direction_end.normalize();
    vertices_pos[1].addConstraint(
        derivative_order::VELOCITY,
        direction_end * createRandomDouble(0.0, kAvgVMax));

    // Segment time.
    double v_avg_segment = createRandomDouble(kAvgVMin, kAvgVMax);
    Eigen::VectorXd pos_start, pos_end;
    vertices_pos.front().getConstraint(derivative_order::POSITION, &pos_start);
    vertices_pos.back().getConstraint(derivative_order::POSITION, &pos_end);
    double distance = (pos_end - pos_start).norm();
    std::vector<double> segment_times = {distance / v_avg_segment};

    // Optimization
    PolynomialOptimization<kN> opt_pos(3);
    opt_pos.setupFromVertices(vertices_pos, segment_times,
                              derivative_order::SNAP);
    opt_pos.solveLinear();
    Segment::Vector pos_segments;
    opt_pos.getSegments(&pos_segments);
    Trajectory pos_trajectory;
    pos_trajectory.setSegments(pos_segments);

    // Yaw segment.
    const double kMinYaw = -3 * M_PI;
    const double kMaxYaw = -kMinYaw;
    Vertex::Vector vertices_yaw = createRandomVertices1D(
        derivative_order::ANGULAR_VELOCITY, 1, kMinYaw, kMaxYaw);
    PolynomialOptimization<kN> opt_yaw(1);
    opt_yaw.setupFromVertices(vertices_yaw, segment_times,
                              derivative_order::ANGULAR_ACCELERATION);
    opt_yaw.solveLinear();
    Segment::Vector yaw_segments;
    opt_yaw.getSegments(&yaw_segments);
    Trajectory yaw_trajectory;
    yaw_trajectory.setSegments(yaw_segments);

    Trajectory trajectory;
    EXPECT_TRUE(pos_trajectory.getTrajectoryWithAppendedDimension(
        yaw_trajectory, &trajectory));

    EXPECT_EQ(trajectory.segments().size(), 1);
    segments[i] = trajectory.segments()[0];
  }
  std::cout << "Created " << segments.size() << " random segments."
            << std::endl;

  // Set regular input constraints.
  InputConstraints input_constraints;
  input_constraints.setDefaultValues();

  // Create feasibility checks with different user settings.
  FeasibilitySampling feasibility_sampling_01(input_constraints);
  feasibility_sampling_01.settings_.setSamplingIntervalS(0.01);
  FeasibilitySampling feasibility_sampling_05(input_constraints);
  feasibility_sampling_05.settings_.setSamplingIntervalS(0.05);
  FeasibilitySampling feasibility_sampling_10(input_constraints);
  feasibility_sampling_10.settings_.setSamplingIntervalS(0.10);

  FeasibilityRecursive feasibility_recursive_01(input_constraints);
  feasibility_recursive_01.settings_.setMinSectionTimeS(0.01);
  FeasibilityRecursive feasibility_recursive_05(input_constraints);
  feasibility_recursive_05.settings_.setMinSectionTimeS(0.05);
  FeasibilityRecursive feasibility_recursive_10(input_constraints);
  feasibility_recursive_10.settings_.setMinSectionTimeS(0.10);

  FeasibilityAnalytic feasibility_analytic_01(input_constraints);
  feasibility_analytic_01.settings_.setMinSectionTimeS(0.01);
  FeasibilityAnalytic feasibility_analytic_05(input_constraints);
  feasibility_analytic_05.settings_.setMinSectionTimeS(0.05);
  FeasibilityAnalytic feasibility_analytic_10(input_constraints);
  feasibility_analytic_10.settings_.setMinSectionTimeS(0.10);

  // Some statistics.
  std::vector<InputFeasibilityResult> result_sampling_01(kNumSegments);
  std::vector<InputFeasibilityResult> result_sampling_05(kNumSegments);
  std::vector<InputFeasibilityResult> result_sampling_10(kNumSegments);

  std::vector<InputFeasibilityResult> result_recursive_01(kNumSegments);
  std::vector<InputFeasibilityResult> result_recursive_05(kNumSegments);
  std::vector<InputFeasibilityResult> result_recursive_10(kNumSegments);

  std::vector<InputFeasibilityResult> result_analytic_01(kNumSegments);
  std::vector<InputFeasibilityResult> result_analytic_05(kNumSegments);
  std::vector<InputFeasibilityResult> result_analytic_10(kNumSegments);

  timing::Timer time_sampling_01("time_sampling_01", false);
  timing::Timer time_sampling_05("time_sampling_05", false);
  timing::Timer time_sampling_10("time_sampling_10", false);

  timing::Timer time_recursive_01("time_recursive_01", false);
  timing::Timer time_recursive_05("time_recursive_05", false);
  timing::Timer time_recursive_10("time_recursive_10", false);

  timing::Timer time_analytic_01("time_analytic_01", false);
  timing::Timer time_analytic_05("time_analytic_05", false);
  timing::Timer time_analytic_10("time_analytic_10", false);
  for (size_t i = 0; i < segments.size(); i++) {
    // Sampling.
    time_sampling_01.Start();
    result_sampling_01[i] =
        feasibility_sampling_01.checkInputFeasibility(segments[i]);
    time_sampling_01.Stop();

    time_sampling_05.Start();
    result_sampling_05[i] =
        feasibility_sampling_05.checkInputFeasibility(segments[i]);
    time_sampling_05.Stop();

    time_sampling_10.Start();
    result_sampling_10[i] =
        feasibility_sampling_10.checkInputFeasibility(segments[i]);
    time_sampling_10.Stop();

    // Recursive.
    time_recursive_01.Start();
    result_recursive_01[i] =
        feasibility_recursive_01.checkInputFeasibility(segments[i]);
    time_recursive_01.Stop();

    time_recursive_05.Start();
    result_recursive_05[i] =
        feasibility_recursive_05.checkInputFeasibility(segments[i]);
    time_recursive_05.Stop();

    time_recursive_10.Start();
    result_recursive_10[i] =
        feasibility_recursive_10.checkInputFeasibility(segments[i]);
    time_recursive_10.Stop();

    // Analtic.
    time_analytic_01.Start();
    result_analytic_01[i] =
        feasibility_analytic_01.checkInputFeasibility(segments[i]);
    time_analytic_01.Stop();

    time_analytic_05.Start();
    result_analytic_05[i] =
        feasibility_analytic_05.checkInputFeasibility(segments[i]);
    time_analytic_05.Stop();

    time_analytic_10.Start();
    result_analytic_10[i] =
        feasibility_analytic_10.checkInputFeasibility(segments[i]);
    time_analytic_10.Stop();

    // If the recursive test shows feasibility, the sampling test also needs to
    // show feasibility.
    if (result_recursive_01[i] == InputFeasibilityResult::kInputFeasible ||
        result_recursive_05[i] == InputFeasibilityResult::kInputFeasible ||
        result_recursive_10[i] == InputFeasibilityResult::kInputFeasible) {
      EXPECT_EQ(result_sampling_01[i], InputFeasibilityResult::kInputFeasible);
    }

    // If the analytic test shows feasibility all other tests should show
    // feasibility or indeterminability.
    if (result_analytic_01[i] == InputFeasibilityResult::kInputFeasible) {
      EXPECT_TRUE(result_analytic_05[i] ==
                      InputFeasibilityResult::kInputFeasible ||
                  result_analytic_05[i] ==
                      InputFeasibilityResult::kInputIndeterminable);
      EXPECT_TRUE(result_analytic_10[i] ==
                      InputFeasibilityResult::kInputFeasible ||
                  result_analytic_10[i] ==
                      InputFeasibilityResult::kInputIndeterminable);

      EXPECT_TRUE(result_recursive_01[i] ==
                      InputFeasibilityResult::kInputFeasible ||
                  result_recursive_01[i] ==
                      InputFeasibilityResult::kInputIndeterminable);
      EXPECT_TRUE(result_recursive_05[i] ==
                      InputFeasibilityResult::kInputFeasible ||
                  result_recursive_05[i] ==
                      InputFeasibilityResult::kInputIndeterminable);
      EXPECT_TRUE(result_recursive_10[i] ==
                      InputFeasibilityResult::kInputFeasible ||
                  result_recursive_10[i] ==
                      InputFeasibilityResult::kInputIndeterminable);

      // The sampling test needs to show feasibility even.
      EXPECT_TRUE(result_sampling_01[i] ==
                  InputFeasibilityResult::kInputFeasible);
      EXPECT_TRUE(result_sampling_05[i] ==
                  InputFeasibilityResult::kInputFeasible);
      EXPECT_TRUE(result_sampling_10[i] ==
                  InputFeasibilityResult::kInputFeasible);
    }
  }

  // Write results to txt.
  writeResultsToFile("result_sampling_01.txt", result_sampling_01);
  writeResultsToFile("result_sampling_05.txt", result_sampling_05);
  writeResultsToFile("result_sampling_10.txt", result_sampling_10);
  writeResultsToFile("result_recursive_01.txt", result_recursive_01);
  writeResultsToFile("result_recursive_05.txt", result_recursive_05);
  writeResultsToFile("result_recursive_10.txt", result_recursive_10);
  writeResultsToFile("result_analytic_01.txt", result_recursive_01);
  writeResultsToFile("result_analytic_05.txt", result_recursive_05);
  writeResultsToFile("result_analytic_10.txt", result_recursive_10);

  // Write timing to txt.
  std::ofstream time_file;
  std::string time_file_name("feasibility_times");
  time_file.open("feasibility_times.txt",
                 std::ofstream::out | std::ofstream::trunc);
  if (time_file.is_open()) {
    time_file << timing::Timing::Print() + "\n";
    time_file.close();
  } else {
    std::cout << "Unable to open file " << time_file_name << "." << std::endl;
  }
}

TEST(FeasibilityTest, HalfPlaneFeasibility) {
  FeasibilityBase half_space_check;
  Eigen::VectorXd coeffs_x(Eigen::Vector3d::Zero());
  Eigen::VectorXd coeffs_y(Eigen::Vector3d::Zero());
  Eigen::VectorXd coeffs_z(Eigen::Vector3d::Zero());
  // Parabola.
  coeffs_x(1) = 1;
  coeffs_z(2) = 1;

  Segment segment(3, 3);
  segment[0] = Polynomial(coeffs_x);
  segment[1] = Polynomial(coeffs_y);
  segment[2] = Polynomial(coeffs_z);
  segment.setTime(1.0);

  Eigen::Vector3d point(0.0, 0.0, 0.0);
  Eigen::Vector3d normal(-1.0, 0.0, 1.0);
  // Shift boundary down.
  while (point.z() > -1.0) {
    half_space_check.half_plane_constraints_.emplace_back(point, normal);
    bool feasible = half_space_check.checkHalfPlaneFeasibility(segment);
    if (point.z() >= -0.25) {
      EXPECT_FALSE(feasible) << point.transpose();
    } else {
      EXPECT_TRUE(feasible) << point.transpose();;
    }
    half_space_check.half_plane_constraints_.clear();
    point.z() -= 0.05;
  }

  FeasibilityBase box_check;
  Eigen::Vector3d box_center(0.0, 0.0, 0.0);

  // Grow box.
  double l = 0.0;
  while(l < 4.0) {
    Eigen::Vector3d box_size(Eigen::Vector3d::Constant(l));
    box_check.half_plane_constraints_ = HalfPlane::createBoundingBox(box_center, box_size);
    bool feasible = box_check.checkHalfPlaneFeasibility(segment);
    if(l <= 2.0) {
      EXPECT_FALSE(feasible);
    }
    else {
      EXPECT_TRUE(feasible);
    }
    l += 0.05;
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  timing::Timing::Print(std::cout);

  return result;
}
