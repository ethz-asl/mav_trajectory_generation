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

#ifndef MAV_TRAJECTORY_GENERATION_VERTEX_H_
#define MAV_TRAJECTORY_GENERATION_VERTEX_H_

#include <chrono>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <glog/logging.h>
#include <map>
#include <vector>

#include "mav_trajectory_generation/motion_defines.h"
#include "mav_trajectory_generation/polynomial.h"

namespace mav_trajectory_generation {

/**
 * \brief Vertex describes the properties of a support point of a path.
 *
 * A vertex has a set of constraints, the derivative of position a value, that
 * have to be matched
 * during optimization procedures. In case of a multi-dimensional vertex (mostly
 * 3D), the constraint
 * for a derivative exists in all dimensions, but can have different values in
 * each dimension.
 *
 * \code
 *     X------------X---------------X
 *   vertex             segment
 * \endcode
 */
class Vertex {
 public:
  typedef std::vector<Vertex> Vector;
  typedef Eigen::VectorXd ConstraintValue;
  typedef std::pair<int, ConstraintValue> Constraint;
  typedef std::map<int, ConstraintValue> Constraints;

  /// constructs an empty vertex and sets time_to_next and
  /// derivative_to_optimize to zero
  Vertex(size_t dimension) : D_(dimension) {}

  int D() const { return D_; }

  /**
   * \brief Adds a constraint for the specified derivative order with the given
   * value.
   * If this is a multi-dimensional vertex, all dimensions are set to the same
   * value.
   */
  void addConstraint(int derivative_order, double value) {
    constraints_[derivative_order] = ConstraintValue::Constant(D_, value);
  }

  /// adds a constraint for the derivative specified in type with the given
  /// values in c
  /**
   * the constraints for the specified derivative get set to the values in the
   * vector c. the dimension has to match the
   * dimension of the vertex
   */
  void addConstraint(int type, const Eigen::VectorXd& constraint);

  /**
   * \brief Sets a constraint for position and sets all derivatives up to
   * (including) up_to_derivative to zero.
   * Convenience method for beginning / end vertices.
   */
  void makeStartOrEnd(const Eigen::VectorXd& constraint,
                      int up_to_derivative);

  void makeStartOrEnd(double value, int up_to_derivative) {
    makeStartOrEnd(Eigen::VectorXd::Constant(D_, value), up_to_derivative);
  }

  /// Returns whether the vertex has a constraint for the specified derivative
  /// order.
  bool hasConstraint(int derivative_order) const;

  /// Passes the value of the constraint for derivative order to *value,
  /// and returns whether the
  ///        constraint is set.
  bool getConstraint(int derivative_order,
                     Eigen::VectorXd* constraint) const;

  /// Returns a const iterator to the first constraint.
  typename Constraints::const_iterator cBegin() const {
    return constraints_.begin();
  }

  /// Returns a const iterator to the end of the constraints, i.e. one after the
  /// last element.
  typename Constraints::const_iterator cEnd() const {
    return constraints_.end();
  }

  /// Returns the number of constraints.
  size_t getNumberOfConstraints() const { return constraints_.size(); }

  /// Checks if both lhs and rhs are equal up to tol in case of double values.
  bool isEqualTol(const Vertex& rhs, double tol) const;

 private:
  int D_;
  Constraints constraints_;
};

std::ostream& operator<<(std::ostream& stream, const Vertex& v);

std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Vertex>& vertices);

/**
 * \brief Makes a rough estimate based on v_max and a_max about the time
 * required to get from one vertex to the next.
 *
 * t_est = 2 * distance/v_max * (1 + magic_fabian_constant * v_max/a_max * exp(-
 * distance/v_max * 2);
 * magic_fabian_constant was determined to 6.5 in a student project ...
 */
std::vector<double> estimateSegmentTimes(const Vertex::Vector& vertices,
                                         double v_max, double a_max,
                                         double magic_fabian_constant = 6.5);

/**
 * \brief Creates random vertices for position within minimum_position and
 * maximum_position.
 *
 * Vertices at the beginning and end have only fixed constraints with their
 * derivative set to zero, while
 * all vertices in between have position as fixed constraint and the derivatives
 * are left free.
 *
 * \param[in] maximum_derivative The maximum derivative to be set to zero for
 * beginning and end.
 * \param[in] s_segments Number of segments of the resulting trajectory. Number
 * of vertices is n_segments + 1.
 * \param[in] minimum_position Minimum position of the space to sample.
 * \param[in] maximum_position Maximum position of the space to sample.
 * \param[in] seed Initial seed for random number generation.
 * \return Vector containing n_segment + 1 vertices.
 */
template <class Derived1, class Derived2>
Vertex::Vector createRandomVertices(
    int maximum_derivative, size_t n_segments,
    const Eigen::VectorXd& minimum_position,
    const Eigen::VectorXd& maximum_position, size_t seed = 0);

/**
 * \brief Conveninence function to create 1D vertices.
 *
 * \sa createRandomVertices
 */
inline Vertex::Vector createRandomVertices1D(int maximum_derivative,
                                             size_t n_segments,
                                             double minimum_position,
                                             double maximum_position,
                                             size_t seed = 0);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_VERTEX_H_

#include "mav_trajectory_generation/impl/vertex_impl.h"
