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
#include <map>
#include <vector>

#include <glog/logging.h>
#include <Eigen/Core>
#include <Eigen/StdVector>

#include <mav_trajectory_generation/motion_defines.h>
#include <mav_trajectory_generation/polynomial.h>

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
 * \endcodbop<int, ConstraintValue> Constraints;

 private:
  int D_;  ///< Number of dimensions.
  Constraints constraints_;

 public:
  /**
   * \brief Constructs an empty vertex and sets time_to_next and
   * derivative_to_optimize to zero.
   */
  Vertex(int D) : D_(D) {}

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

  /**
   * \brief Adds a constraint for the derivative specified in type with the
   * given values in c. The dimension has to match the dimension of the vertex.
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
    makeStartOrEnd(Eigen::VectorXd::Constant(D_, value), up_to_derivative);
  }

  /**
   * \brief Returns whether the vertex has a constraint for the specified
   * derivative order.
   */
  bool hasConstraint(int derivative_order) const;

  /**
   * \brief Passes the value of the constraint for derivative order to *value,
   * and returns whether the
   *        constraint is set.
   */
  template <class Derived>
  bool getConstraint(int derivative_order,
                     Eigen::MatrixBase<Derived>* value) const;

  /**
   * \brief Returns a const iterator to the first constraint.
   */
  typename Constraints::const_iterator cBegin() const {
    return constraints_.begin();
  }

  /**
   * \brief Returns a const iterator to the end of the constraints,
      i.e. one after the last element.
   */
  typename Constraints::const_iterator cEnd() const {
    return constraints_.end();
  }

  /**
   * \brief Returns the number of constraints.
   */
  int getNumberOfConstraints() const { return constraints_.size(); }

  /**
   * \brief Checks if both lhs and rhs are equal up to tol in case of double
   * values.
   */
  bool isEqualTol(const Vertex& rhs, double tol) const;
};

std::ostream& operator<<(std::ostream& stream, const Vertex& v);
std::ostream& operator<<(std::ostream& stream,
                         const std::vector<Vertex>& vertices);

}  // end namespace mav_trajectory_generation
#endif // MAV_TRAJECTORY_GENERATION_VERTEX_H_
