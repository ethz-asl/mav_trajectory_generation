/*
 * Simple example that shows a trajectory planner using
 *  mav_trajectory_generation.
 *
 *
 * Launch via
 *   roslaunch mav_trajectory_generation_example example.launch
 *
 * Wait for console to run through all gazebo/rviz messages and then
 * you should see the example below
 *  - After Enter, it receives the current uav position
 *  - After second enter, publishes trajectory information
 *  - After third enter, executes trajectory (sends it to the sampler)
 */

#include  "ros/ros.h"
#include <mav_trajectory_generation_example/example_planner.h>

#include <iostream>

int main(int argc, char** argv) {

  ros::init(argc, argv, "simple_planner");

  ros::NodeHandle n;
  ExamplePlanner planner(n);
  ROS_WARN_STREAM("SLEEPING FOR 5s TO WAIT FOR CLEAR CONSOLE");
  ros::Duration(5.0).sleep();
  ROS_WARN_STREAM("WARNING: CONSOLE INPUT/OUTPUT ONLY FOR DEMONSTRATION!");

  // THIS SHOULD NORMALLY RUN INSIDE ROS::SPIN!!! JUST FOR DEMO PURPOSES LIKE THIS.
  ROS_WARN_STREAM("PRESS ENTER TO UPDATE CURRENT POSITION");
  std::cin.get();
  for (int i = 0; i < 10; i++) {
    ros::spinOnce();  // process a few messages in the background - causes the uavPoseCallback to happen
  }
  Eigen::Affine3d pose;
  planner.getCurrentPose(&pose);
  ROS_WARN_STREAM("CURRENT POSITION: " << pose.translation());

  ROS_WARN_STREAM("PRESS ENTER TO VISUALIZE TRAJECTORY....");
  std::cin.get();
  planner.planTakeOffTrajectory(5.0);
  planner.visualizeCachedTrajectory();

  ROS_WARN_STREAM("SENT RVIZ MARKERS. PRESS ENTER TO EXECUTE.");
  std::cin.get();
  planner.executeCachedTrajectory();
  ROS_WARN_STREAM("DONE. GOODBYE.");

  return 0;
}