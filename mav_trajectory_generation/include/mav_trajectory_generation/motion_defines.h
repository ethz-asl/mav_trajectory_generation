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

#ifndef MAV_TRAJECTORY_MOTION_DEFINES_H_
#define MAV_TRAJECTORY_MOTION_DEFINES_H_

#include <string>

namespace mav_trajectory_generation {

namespace derivative_order {
static constexpr int POSITION = 0;
static constexpr int VELOCITY = 1;
static constexpr int ACCELERATION = 2;
static constexpr int JERK = 3;
static constexpr int SNAP = 4;

static constexpr int ORIENTATION = 0;
static constexpr int ANGULAR_VELOCITY = 1;
static constexpr int ANGULAR_ACCELERATION = 2;

static constexpr int INVALID = -1;
}

std::string positionDerivativeToString(int derivative);
int positionDerivativeToInt(const std::string& string);

std::string orintationDerivativeToString(int derivative);
int orientationDerivativeToInt(const std::string& string);

}  // namespace mav_trajectory_generation

#endif
