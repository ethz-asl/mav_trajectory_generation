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

#ifndef MAV_TRAJECTORY_GENERATION_ROS_INPUT_CONSTRAINTS_H_
#define MAV_TRAJECTORY_GENERATION_ROS_INPUT_CONSTRAINTS_H_

#include <yaml-cpp/yaml.h>
#include <map>

namespace mav_trajectory_generation {

enum InputConstraintType {
  kFMin = 0,     // Minimum acceleration (normalized thrust) in [m/s/s].
  kFMax,         // Maximum acceleration (normalized thrust) in [m/s/s].
  kVMax,         // Maximum velocity in [m/s].
  kOmegaXYMax,   // Maximum roll/pitch rate in [rad/s].
  kOmegaZMax,    // Maximum yaw rate in [rad/s].
  kOmegaZDotMax  // Maximum yaw acceleration in [rad/s/s].
};

std::string getInputConstraintName(InputConstraintType type);

// Dynamic constraints of the MAV.
class InputConstraints {
 public:
  // Empty constraints object.
  InputConstraints() {}

  // Set a constraint given type and value.
  void addConstraint(int constraint_type, double value);

  // Sets all constraints to reasonable default values.
  void setDefaultValues();

  // Return a constraint. Returns false if constraint is not set.
  bool getConstraint(int constraint_type, double* value) const;

  // Check if a specific constraint type is set.
  bool hasConstraint(int constraint_type) const;

  // Remove a specific constraint type. Returns false if constraint was not set.
  bool removeConstraint(int constraint_type);

  // Save this to a YAML node.
  YAML::Node toYaml() const;

  // Load this from a YAML node.
  void fromYaml(const YAML::Node& node);

 private:
  std::map<int, double> constraints_;
};
}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_ROS_INPUT_CONSTRAINTS_H_
