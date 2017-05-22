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

#include <random>

#include "mav_trajectory_generation/vertex.h"

namespace mav_trajectory_generation {

Vertex::Vector createRandomVertices(int maximum_derivative, size_t n_segments,
                                    const Eigen::VectorXd& pos_min,
                                    const Eigen::VectorXd& pos_max,
                                    size_t seed) {
  CHECK_GE(static_cast<int>(n_segments), 1);
  CHECK_EQ(pos_min.size(), pos_max.size());
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

Vertex::Vector createRandomVertices1D(int maximum_derivative, size_t n_segments,
                                      double pos_min, double pos_max,
                                      size_t seed) {
  return createRandomVertices(maximum_derivative, n_segments,
                              Eigen::VectorXd::Constant(1, pos_min),
                              Eigen::VectorXd::Constant(1, pos_max), seed);
}

void Vertex::addConstraint(int derivative_order,
                           const Eigen::VectorXd& constraint) {
  CHECK_EQ(constraint.rows(), static_cast<long>(D_));
  constraints_[derivative_order] = constraint;
}

bool Vertex::removeConstraint(int type) {
  Constraints::const_iterator it = constraints_.find(type);
  if (it != constraints_.end()) {
    constraints_.erase(it);
    return true;
  } else {
    // Constraint not found.
    return false;
  }
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

bool Vertex::hasConstraint(int derivative_order) const {
  typename Constraints::const_iterator it = constraints_.find(derivative_order);
  return it != constraints_.end();
}

bool Vertex::isEqualTol(const Vertex& rhs, double tol) const {
  if (constraints_.size() != rhs.constraints_.size()) return false;
  // loop through lhs constraint map
  for (typename Constraints::const_iterator it = cBegin(); it != cEnd(); ++it) {
    // look for matching key
    typename Constraints::const_iterator rhs_it =
        rhs.constraints_.find(it->first);
    if (rhs_it == rhs.constraints_.end()) return false;
    // check value
    if (!((it->second - rhs_it->second).isZero(tol))) return false;
  }
  return true;
}

std::ostream& operator<<(std::ostream& stream, const Vertex& v) {
  stream << "constraints: " << std::endl;
  Eigen::IOFormat format(4, 0, ", ", "\n", "[", "]");
  for (typename Vertex::Constraints::const_iterator it = v.cBegin();
       it != v.cEnd(); ++it) {
    stream << "  type: " << positionDerivativeToString(it->first);
    stream << "  value: " << it->second.transpose().format(format) << std::endl;
  }
  return stream;
}

std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Vertex>& vertices) {
  for (const Vertex& v : vertices) {
    stream << v << std::endl;
  }
  return stream;
}

std::vector<double> estimateSegmentTimes(const Vertex::Vector& vertices,
                                         double v_max, double a_max,
                                         double magic_fabian_constant) {
  std::vector<double> segment_times;
  segment_times.reserve(vertices.size() - 1);
  for (size_t i = 0; i < vertices.size() - 1; ++i) {
    Eigen::VectorXd start, end;
    vertices[i].getConstraint(derivative_order::POSITION, &start);
    vertices[i + 1].getConstraint(derivative_order::POSITION, &end);
    double distance = (end - start).norm();
    double t = distance / v_max * 2 *
               (1.0 + magic_fabian_constant * v_max / a_max *
                          exp(-distance / v_max * 2));
    segment_times.push_back(t);
  }
  return segment_times;
}

}  // namespace mav_trajectory_generation
