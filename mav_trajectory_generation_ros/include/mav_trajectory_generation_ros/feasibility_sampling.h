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

#ifndef MAV_TRAJECTORY_GENERATION_FEASIBILITY_SAMPLING_H_
#define MAV_TRAJECTORY_GENERATION_FEASIBILITY_SAMPLING_H_

#include <mav_msgs/eigen_mav_msgs.h>
#include <mav_trajectory_generation/trajectory.h>

#include "mav_trajectory_generation_ros/feasibility_base.h"

namespace mav_trajectory_generation {

// Sampling based input and position feasibility checks.
class FeasibilitySampling : public FeasibilityBase {
 public:
  struct Settings {
    Settings();
    double sampling_interval;
  };

  FeasibilitySampling();
  FeasibilitySampling(const InputConstraints& input_constraints);
  // Checks a trajectory for input feasibility.
  virtual bool checkInputFeasibility(const Trajectory& trajectory);

  // Checks if a trajectory stays within a set of half planes.
  virtual bool checkHalfPlaneFeasibility(const Trajectory& trajectory);

  Settings settings_;

 private:
  // Sets and samples the whole trajectory if not done already.
  bool sampleTrajectory(const Trajectory& trajectory);
  // The current trajectory.
  Trajectory trajectory_;
  // The sampled states along the trajectory (differentially flat assumption).
  mav_msgs::EigenMavState::Vector states_;
};
}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_FEASIBILITY_SAMPLING_H_
