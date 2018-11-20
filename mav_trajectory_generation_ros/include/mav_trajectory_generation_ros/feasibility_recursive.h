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

#ifndef MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_RECURSIVE_H_
#define MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_RECURSIVE_H_

#include <cmath>

#include <mav_trajectory_generation/segment.h>

#include "mav_trajectory_generation_ros/feasibility_base.h"

namespace mav_trajectory_generation {

/* Recursive input feasibility checks.
 * This implementation is based on [1] and extended to test yaw rates and
 * higher order polynomials.
 * The general idea is to check a segment for lower and upper bounds that can be
 * found easily by evaluating the single axis minima and maxima.
 * We extend checking for maximum velocity constraints and yaw rate and
 * acceleration constraints.
 *
 * [1] Mueller, Mark W., Markus Hehn, and Raffaello D'Andrea. "A
 * Computationally Efficient Motion Primitive for Quadrocopter
 * Trajectory Generation." Robotics, IEEE Transactions on 31.6
 * (2015): 1294-1310.
 *
 * See also https://github.com/markwmuller/RapidQuadrocopterTrajectories
 */
class FeasibilityRecursive : public FeasibilityBase {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef std::vector<Eigen::VectorXcd,
                      Eigen::aligned_allocator<Eigen::VectorXcd>>
      Roots;

  class Settings {
   public:
    Settings();

    inline void setMinSectionTimeS(double min_section_time_s) {
      min_section_time_s_ = std::abs(min_section_time_s);
    }
    inline double getMinSectionTimeS() const { return min_section_time_s_; }

   private:
    // The minimum section length the binary search is going to check. If the
    // trajectory is feasible with respect to an upper bound for this section
    // length it is considered overall feasible.
    double min_section_time_s_;
  };

  FeasibilityRecursive() {}
  FeasibilityRecursive(const Settings& settings);
  FeasibilityRecursive(const InputConstraints& input_constraints);
  FeasibilityRecursive(const Settings& settings,
                       const InputConstraints& input_constraints);

  // Checks a segment for input feasibility.
  virtual InputFeasibilityResult checkInputFeasibility(
      const Segment& segment) const;

  // The user settings.
  Settings settings_;

 private:
  // Recursive test to determine velocity, acceleration, and angular rate (roll,
  // pitch) feasibility between t_1 and t_2.
  InputFeasibilityResult recursiveFeasibility(const Segment& segment,
                                              const Roots& roots_acc,
                                              const Roots& roots_jerk,
                                              const Roots& roots_snap,
                                              double t_1, double t_2) const;

  double evaluateThrust(const Segment& segment, double time) const;
};
}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_ROS_FEASIBILITY_RECURSIVE_H_
