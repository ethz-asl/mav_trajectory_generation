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

#include "mav_trajectory_generation/motion_defines.h"

namespace mav_trajectory_generation {

std::string positionDerivativeToString(int derivative) {
  if (derivative >= derivative_order::POSITION &&
      derivative <= derivative_order::SNAP) {
    static constexpr const char* text[] = {"position", "velocity",
                                           "acceleration", "jerk", "snap"};
    return std::string(text[derivative]);
  } else {
    return std::string("invalid");
  }
}

int positionDerivativeToInt(const std::string& string) {
  using namespace derivative_order;
  if (string == "position") {
    return POSITION;
  } else if (string == "velocity") {
    return VELOCITY;
  } else if (string == "acceleration") {
    return ACCELERATION;
  } else if (string == "jerk") {
    return JERK;
  } else if (string == "snap") {
    return SNAP;
  } else {
    return INVALID;
  }
}

std::string orintationDerivativeToString(int derivative) {
  if (derivative >= derivative_order::ORIENTATION &&
      derivative <= derivative_order::ANGULAR_ACCELERATION) {
    static constexpr const char* text[] = {"orientation", "angular_velocity",
                                           "angular_acceleration"};
    return std::string(text[derivative]);
  } else {
    return std::string("invalid");
  }
}

int orientationDerivativeToInt(const std::string& string) {
  using namespace derivative_order;
  if (string == "orientation") {
    return ORIENTATION;
  } else if (string == "angular_velocity") {
    return ANGULAR_VELOCITY;
  } else if (string == "angular_acceleration") {
    return ANGULAR_ACCELERATION;
  } else {
    return INVALID;
  }
}

}  // namespace mav_trajectory_generation
