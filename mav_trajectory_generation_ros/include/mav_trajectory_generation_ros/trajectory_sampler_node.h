/*
 * Copyright 2017 Marija Popovic, ASL, ETH Zurich, Switzerland
 *
 * This code can be used only for academic, non-commercial use.
 * This code cannot be redistributed under any license, open source or
 * otherwise.
 *
 */

#ifndef TRAJECTORY_SAMPLER_NODE_H
#define TRAJECTORY_SAMPLER_NODE_H

#include <mav_msgs/conversions.h>
#include <mav_msgs/default_topics.h>
#include <mav_msgs/eigen_mav_msgs.h>
#include <planning_msgs/PolynomialSegment4D.h>
#include <planning_msgs/PolynomialTrajectory4D.h>
#include <ros/ros.h>
#include <std_srvs/Empty.h>
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

#include <mav_trajectory_generation/polynomial.h>
#include <mav_trajectory_generation_ros/ros_conversions.h>
#include <mav_trajectory_generation_ros/trajectory_sampling.h>

class TrajectorySamplerNode {
 public:
  TrajectorySamplerNode(const ros::NodeHandle& nh,
                        const ros::NodeHandle& nh_private);
  ~TrajectorySamplerNode();

 private:
  void pathSegmentsCallback(
      const planning_msgs::PolynomialTrajectory4D& segments_message);
  bool stopSamplingCallback(std_srvs::Empty::Request& request,
                            std_srvs::Empty::Response& response);
  void commandTimerCallback(const ros::TimerEvent&);

  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;

  ros::Timer publish_timer_;
  ros::Subscriber trajectory_sub_;
  ros::Publisher command_pub_;
  ros::ServiceServer stop_srv_;
  ros::Time start_time_;

  // Flag whether to publish entire trajectory at once or not.
  bool publish_whole_trajectory_ = false;
  // Trajectory sampling interval.
  double dt_ = 0;
  // Time at currently published trajectory sample.
  double current_sample_time_ = 0;

  // The trajectory to sub-sample.
  mav_trajectory_generation::Trajectory trajectory_;
};

#endif  // TRAJECTORY_SAMPLER_NODE_H