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

#ifndef MAV_TRAJECTORY_GENERATION_YAML_IO_H_
#define MAV_TRAJECTORY_GENERATION_YAML_IO_H_

#include "mav_trajectory_generation/segment.h"
#include "mav_trajectory_generation/trajectory.h"

namespace mav_trajectory_generation {

bool segmentsToFile(const std::string& filename,
                    const mav_trajectory_generation::Segment::Vector& segments);

inline bool trajectoryToFile(
    const std::string& filename,
    const mav_trajectory_generation::Trajectory& trajectory) {
  mav_trajectory_generation::Segment::Vector segments;
  trajectory.getSegments(&segments);
  return segmentsToFile(filename, segments);
}

bool segmentsFromFile(const std::string& filename,
                      mav_trajectory_generation::Segment::Vector* segments);

inline bool trajectoryFromFile(
    const std::string& filename,
    mav_trajectory_generation::Trajectory* trajectory) {
  mav_trajectory_generation::Segment::Vector segments;
  bool success = segmentsFromFile(filename, &segments);
  trajectory->setSegments(segments);
  return success;
}

bool sampledTrajectoryStatesToFile(const std::string& filename,
                                   const Trajectory& trajectory);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_YAML_IO_H_
