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
  Vertex(size_t dimension) : dimension_(dimension) {}

  size_t getDimension() const { return dimension_; }

  /**
   * \brief Adds a constraint for the specified derivative order with the given
   * value.
   * If this is a multi-dimensional vertex, all dimensions are set to the same
   * value.
   */
  void addConstraint(int derivative_order, double value) {
    constraints_[derivative_order] =
        ConstraintValue::Constant(dimension_, value);
  }

  /// adds a constraint for the derivative specified in type with the given
  /// values in c
  /**
   * the constraints for the specified derivative get set to the values in the
   * vector c. the dimension has to match the
   * dimension of the vertex
   */
  template <class Derived>
  void addConstraint(int type, const Eigen::MatrixBase<Derived>& c);

  /**
   * \brief Sets a constraint for position and sets all derivatives up to
   * (including) up_to_derivative to zero.
   * Convenience method for beginning / end vertices.
   */
  template <class Derived>
  void makeStartOrEnd(const Eigen::MatrixBase<Derived>& c,
                      int up_to_derivative);

  void makeStartOrEnd(double value, int up_to_derivative) {
    makeStartOrEnd(Eigen::VectorXd::Constant(dimension_, value),
                   up_to_derivative);
  }

  /// Returns whether the vertex has a constraint for the specified derivative
  /// order.
  bool hasConstraint(int derivative_order) const;

  /// Passes the value of the constraint for derivative order to *value,
  /// and returns whether the
  ///        constraint is set.
  template <class Derived>
  bool getConstraint(int derivative_order,
                     Eigen::MatrixBase<Derived>* value) const;

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
  size_t dimension_;
  Constraints constraints_;
};

std::ostream& operator<<(std::ostream& stream, const Vertex& v);

std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Vertex>& vertices);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_VERTEX_H_

#include "mav_trajectory_generation/impl/vertex_impl.h"
