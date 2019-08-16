/*
 * Copyright (c) 2017, Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2017, Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2017, Helen Oleynikova, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2017, Rik Bähnemann, ASL, ETH Zurich, Switzerland
 * Copyright (c) 2017, Marija Popovic, ASL, ETH Zurich, Switzerland
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

#ifndef TRAJECTORY_SAMPLER_NODE_H
#define TRAJECTORY_SAMPLER_NODE_H

#include <mav_msgs/conversions.h>
#include <mav_msgs/default_topics.h>
#include <mav_msgs/eigen_mav_msgs.h>
#include <mav_planning_msgs/PolynomialSegment.h>
#include <mav_planning_msgs/PolynomialTrajectory.h>
#include <ros/ros.h>
#include <std_srvs/Empty.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

#include <mav_trajectory_generation/polynomial.h>
#include <mav_trajectory_generation/trajectory_sampling.h>
#include <mav_trajectory_generation_ros/ros_conversions.h>

// deprecated
#include <mav_planning_msgs/PolynomialSegment4D.h>
#include <mav_planning_msgs/PolynomialTrajectory4D.h>

class TrajectorySamplerNode {
 public:
  TrajectorySamplerNode(const ros::NodeHandle& nh,
                        const ros::NodeHandle& nh_private);
  ~TrajectorySamplerNode();

 private:
  void pathSegmentsCallback(
      const mav_planning_msgs::PolynomialTrajectory& segments_message);
  void pathSegments4DCallback(
      const mav_planning_msgs::PolynomialTrajectory4D& segments_message);
  bool stopSamplingCallback(std_srvs::Empty::Request& request,
                            std_srvs::Empty::Response& response);
  void commandTimerCallback(const ros::TimerEvent&);
  void processTrajectory();

  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;

  ros::Timer publish_timer_;
  ros::Subscriber trajectory_sub_;
  ros::Subscriber trajectory4D_sub_;
  ros::Publisher command_pub_;
  ros::ServiceServer stop_srv_;
  ros::Time start_time_;

  // Service client for getting the MAV interface to listen to our sent
  // commands.
  ros::ServiceClient position_hold_client_;

  // Flag whether to publish entire trajectory at once or not.
  bool publish_whole_trajectory_;
  // Trajectory sampling interval.
  double dt_;
  // Time at currently published trajectory sample.
  double current_sample_time_;

  // The trajectory to sub-sample.
  mav_trajectory_generation::Trajectory trajectory_;
};

#endif  // TRAJECTORY_SAMPLER_NODE_H
