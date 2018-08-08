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

#include <chrono>
#include <iostream>
#include <limits>
#include <random>

#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation/timing.h>

const int N = 10;
const int max_derivative = mav_trajectory_generation::derivative_order::SNAP;
const size_t derivative_to_optimize =
    mav_trajectory_generation::derivative_order::SNAP;

mav_trajectory_generation::Vertex::Vector createRandomVerticesPath(
    int dimension, size_t n_segments, double average_distance,
    int maximum_derivative, size_t seed) {
  CHECK_GE(static_cast<int>(n_segments), 1);

  CHECK_GT(maximum_derivative, 0);

  mav_trajectory_generation::Vertex::Vector vertices;
  std::mt19937 generator(seed);
  std::vector<std::uniform_real_distribution<double> > distribution;
  std::uniform_real_distribution<double> random_distance(0,
                                                         2 * average_distance);

  distribution.resize(dimension);

  for (int i = 0; i < dimension; ++i) {
    distribution[i] = std::uniform_real_distribution<double>(-1, 1);
  }

  const double min_distance = 0.2;
  const int n_vertices = n_segments + 1;

  Eigen::VectorXd last_position(dimension);
  for (int i = 0; i < dimension; ++i) {
    last_position[i] = distribution[i](generator);
  }

  vertices.reserve(n_segments + 1);
  vertices.push_back(mav_trajectory_generation::Vertex(dimension));

  vertices.front().makeStartOrEnd(last_position, maximum_derivative);

  double distance_accumulated = 0;

  for (int i = 1; i < n_vertices; ++i) {
    Eigen::VectorXd position_sample(dimension);

    while (true) {
      for (int d = 0; d < dimension; ++d) {
        position_sample[d] = distribution[d](generator);
      }
      if (position_sample.norm() > min_distance) break;
    }

    position_sample = position_sample.normalized() * random_distance(generator);

    distance_accumulated += position_sample.norm();

    mav_trajectory_generation::Vertex v(dimension);
    v.addConstraint(mav_trajectory_generation::derivative_order::POSITION,
                    position_sample + last_position);
    vertices.push_back(v);
    last_position = position_sample;
  }
  vertices.back().makeStartOrEnd(last_position, maximum_derivative);

  return vertices;
}

bool timeEval(int n_segments, double average_distance, size_t seed) {
  mav_trajectory_generation::Vertex::Vector vertices;
  vertices = createRandomVerticesPath(3, n_segments, average_distance,
                                      max_derivative, seed);

  const double approximate_v_max = 2.0;
  const double approximate_a_max = 2.0;
  const double magic_fabian_constant = 6.5;
  std::vector<double> segment_times =
      estimateSegmentTimes(vertices, approximate_v_max, approximate_a_max);

  mav_trajectory_generation::timing::Timer timer_solve(
      "polynomial_optimization_template_" + std::to_string(n_segments));
  mav_trajectory_generation::PolynomialOptimization<N> opt(3);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);

  opt.solveLinear();
  timer_solve.Stop();
  return true;
}

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  int n_segments_to_test[4] = {2, 10, 50, 100};
  double average_distance = 5;
  unsigned long seed = 1;

  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 1000; ++i) {
      int n_segments = n_segments_to_test[j];
      timeEval(n_segments, average_distance, seed);
    }
  }
  mav_trajectory_generation::timing::Timing::Print(std::cout);
}
