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

#ifndef MAV_TRAJECTORY_GENERATION_IMPL_VERTEX_IMPL_H_
#define MAV_TRAJECTORY_GENERATION_IMPL_VERTEX_IMPL_H_

#include <random>

namespace mav_trajectory_generation {

Vertex::Vector createRandomVertices(int maximum_derivative, size_t n_segments,
                                    const Eigen::VectorXd& pos_min,
                                    const Eigen::VectorXd& pos_max,
                                    size_t seed) {
  CHECK_GE(static_cast<int>(n_segments), 1);
  CHECK_EQ(pos_min.size(), pos_min.size());
  CHECK_GT(maximum_derivative, 0);

  Vertex::Vector vertices;
  std::mt19937 generator(seed);
  std::vector<std::uniform_real_distribution<double> > distribution;

  const size_t dimension = pos_min.size();

  distribution.resize(dimension);

  for (size_t i = 0; i < dimension; ++i) {
    distribution[i] =
        std::uniform_real_distribution<double>(pos_min[i], pos_max[i]);
  }

  const double min_distance = 0.2;
  const size_t n_vertices = n_segments + 1;

  Eigen::VectorXd last_pos(dimension);
  for (size_t i = 0; i < dimension; ++i) {
    last_pos[i] = distribution[i](generator);
  }

  vertices.reserve(n_segments + 1);
  vertices.push_back(Vertex(dimension));

  vertices.front().makeStartOrEnd(last_pos, maximum_derivative);

  for (size_t i = 1; i < n_vertices; ++i) {
    Eigen::VectorXd pos(dimension);

    while (true) {
      for (size_t d = 0; d < dimension; ++d) {
        pos[d] = distribution[d](generator);
      }
      if ((pos - last_pos).norm() > min_distance) break;
    }

    Vertex v(dimension);
    v.addConstraint(derivative_order::POSITION, pos);
    vertices.push_back(v);
    last_pos = pos;
  }
  vertices.back().makeStartOrEnd(last_pos, maximum_derivative);

  return vertices;
}

inline Vertex::Vector createRandomVertices1D(int maximum_derivative,
                                             size_t n_segments, double pos_min,
                                             double pos_max, size_t seed) {
  return createRandomVertices(maximum_derivative, n_segments,
                              Eigen::VectorXd::Constant(1, pos_min),
                              Eigen::VectorXd::Constant(1, pos_max), seed);
}

void Vertex::addConstraint(int derivative_order,
                           const Eigen::VectorXd& constraint) {
  CHECK_EQ(constraint.rows(), static_cast<long>(D_));
  constraints_[derivative_order] = constraint;
}

void Vertex::makeStartOrEnd(const Eigen::VectorXd& constraint,
                            int up_to_derivative) {
  addConstraint(derivative_order::POSITION, constraint);
  for (int i = 1; i <= up_to_derivative; ++i) {
    constraints_[i] = ConstraintValue::Zero(static_cast<int>(D_));
  }
}

bool Vertex::getConstraint(int derivative_order, Eigen::VectorXd* value) const {
  CHECK_NOTNULL(value);
  typename Constraints::const_iterator it = constraints_.find(derivative_order);
  if (it != constraints_.end()) {
    *value = it->second;
    return true;
  } else
    return false;
}
}

#endif  // MAV_TRAJECTORY_GENERATION_IMPL_VERTEX_IMPL_H_