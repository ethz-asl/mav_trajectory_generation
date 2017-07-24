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

#include <cmath>

#include <glog/logging.h>

#include <mav_msgs/default_values.h>

#include "mav_trajectory_generation_ros/input_constraints.h"

namespace mav_trajectory_generation {
typedef InputConstraintType ICT;

void InputConstraints::addConstraint(int constraint_type, double value) {
  // Correct user input.
  value = std::abs(value);
  if (constraint_type == ICT::kFMin && hasConstraint(ICT::kFMax)) {
    constraints_[ICT::kFMax] =
        value > constraints_[ICT::kFMax] ? value : constraints_[ICT::kFMax];
  } else if (constraint_type == ICT::kFMax && hasConstraint(ICT::kFMin)) {
    constraints_[ICT::kFMin] =
        value < constraints_[ICT::kFMin] ? value : constraints_[ICT::kFMin];
  }
  constraints_[constraint_type] = value;
}

void InputConstraints::setDefaultValues() {
  constraints_[ICT::kFMin] = 0.5 * mav_msgs::kGravity;
  constraints_[ICT::kFMax] = 1.5 * mav_msgs::kGravity;
  constraints_[ICT::kVMax] = 3.0;
  constraints_[ICT::kOmegaXYMax] = M_PI / 2.0;
  constraints_[ICT::kOmegaZMax] = M_PI / 2.0;
  constraints_[ICT::kOmegaZDotMax] = 2.0 * M_PI;
}

bool InputConstraints::getConstraint(int constraint_type, double* value) const {
  CHECK_NOTNULL(value);
  std::map<int, double>::const_iterator it = constraints_.find(constraint_type);
  if (it != constraints_.end()) {
    *value = it->second;
    return true;
  } else {
    return false;
  }
}

bool InputConstraints::hasConstraint(int constraint_type) const {
  std::map<int, double>::const_iterator it = constraints_.find(constraint_type);
  return it != constraints_.end();
}

bool InputConstraints::removeConstraint(int constraint_type) {
  std::map<int, double>::const_iterator it = constraints_.find(constraint_type);
  if (it != constraints_.end()) {
    constraints_.erase(it);
    return true;
  } else {
    return false;
  }
}

}  // namespace mav_trajectory_generation
