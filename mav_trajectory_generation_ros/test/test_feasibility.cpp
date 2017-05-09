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
#include <eigen-checks/gtest.h>

#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation/segment.h>
#include <mav_trajectory_generation/timing.h>
#include <mav_trajectory_generation/vertex.h>

#include "mav_trajectory_generation_ros/feasibility_recursive.h"
#include "mav_trajectory_generation_ros/feasibility_sampling.h"

using namespace mav_trajectory_generation;

double createRandomDouble(double min, double max) {
  // No seed for repeatability.
  return (max - min) * (static_cast<double>(std::rand()) /
                        static_cast<double>(RAND_MAX)) +
         min;
}

TEST(PolynomialTest, CompareFeasibilityTests) {
  // Create random segments.
  const int kNumSegments = 1e3;
  const int kN = 12;
  const int kD = 4;
  Segment::Vector segments(kNumSegments, Segment(kN, kD));
  for (size_t i = 0; i < kNumSegments; i++) {
    // Position segment.
    const Eigen::VectorXd kMinPos = Eigen::VectorXd::Constant(3, -100.0);
    const Eigen::VectorXd kMaxPos = -kMinPos;
    Vertex::Vector vertices_pos =
        createRandomVertices(derivative_order::SNAP, 1, kMinPos, kMaxPos);
    const double kApproxVMin = 2.0;
    const double kApproxVMax = 5.0;
    double kApproxV = createRandomDouble(kApproxVMin, kApproxVMax);
    const double kApproxAMax = 5.0;
    std::vector<double> segment_times =
        estimateSegmentTimes(vertices_pos, kApproxV, kApproxAMax);
    PolynomialOptimization<kN> opt_pos(3);
    opt_pos.setupFromVertices(vertices_pos, segment_times,
                              derivative_order::SNAP);
    opt_pos.solveLinear();
    Segment::Vector pos_segments;
    opt_pos.getSegments(&pos_segments);
    Trajectory pos_trajectory;
    pos_trajectory.setSegments(pos_segments);

    // Yaw segment.
    const double kMinYaw = -M_PI;
    const double kMaxYaw = -kMinYaw;
    Vertex::Vector vertices_yaw = createRandomVertices1D(
        derivative_order::ANGULAR_ACCELERATION, 1, kMinYaw, kMaxYaw);
    PolynomialOptimization<kN> opt_yaw(1);
    opt_yaw.setupFromVertices(vertices_yaw, segment_times,
                              derivative_order::ANGULAR_ACCELERATION);
    opt_yaw.solveLinear();
    Segment::Vector yaw_segments;
    opt_yaw.getSegments(&yaw_segments);
    Trajectory yaw_trajectory;
    yaw_trajectory.setSegments(yaw_segments);

    Trajectory trajectory =
        pos_trajectory.getTrajectoryWithAppendedDimension(yaw_trajectory);

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

  // Some statistics.
  std::vector<InputFeasibilityResult> result_sampling_01(kNumSegments);
  std::vector<InputFeasibilityResult> result_sampling_05(kNumSegments);
  std::vector<InputFeasibilityResult> result_sampling_10(kNumSegments);

  std::vector<InputFeasibilityResult> result_recursive_01(kNumSegments);
  std::vector<InputFeasibilityResult> result_recursive_05(kNumSegments);
  std::vector<InputFeasibilityResult> result_recursive_10(kNumSegments);

  timing::Timer time_sampling_01("time_sampling_01", false);
  timing::Timer time_sampling_05("time_sampling_05", false);
  timing::Timer time_sampling_10("time_sampling_10", false);

  timing::Timer time_recursive_01("time_recursive_01", false);
  timing::Timer time_recursive_05("time_recursive_05", false);
  timing::Timer time_recursive_10("time_recursive_10", false);
  for (size_t i = 0; i < segments.size(); i++) {
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

    // If the recursive test shows feasibility, the sampling test also needs to
    // show feasibility.
    if (result_recursive_01[i] == InputFeasibilityResult::kInputFeasible ||
        result_recursive_05[i] == InputFeasibilityResult::kInputFeasible ||
        result_recursive_10[i] == InputFeasibilityResult::kInputFeasible) {
      EXPECT_EQ(result_sampling_01[i], InputFeasibilityResult::kInputFeasible);
    }
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  timing::Timing::Print(std::cout);

  return result;
}
