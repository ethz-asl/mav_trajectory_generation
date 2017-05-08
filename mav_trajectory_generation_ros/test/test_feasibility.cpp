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
#include <mav_trajectory_generation/vertex.h>

using namespace mav_trajectory_generation;

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
    const double kApproxVMax = 3.0;
    const double kApproxAMax = 5.0;
    std::vector<double> segment_times =
        estimateSegmentTimes(vertices_pos, kApproxVMax, kApproxAMax);
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
                              derivative_order::ANGULAR_VELOCITY);
    opt_yaw.solveLinear();
    Segment::Vector yaw_segments;
    opt_yaw.getSegments(&yaw_segments);
    Trajectory yaw_trajectory;
    yaw_trajectory.setSegments(yaw_segments);

    Trajectory trajectory =
        pos_trajectory.getTrajectoryWithAppendedDimension(yaw_trajectory);

    segments[i] = trajectory.segments()[0];
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
