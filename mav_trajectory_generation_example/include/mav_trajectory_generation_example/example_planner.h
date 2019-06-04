#include <iostream>
#include <ros/ros.h>
#include <Eigen/Dense>
#include <eigen_conversions/eigen_msg.h>
#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation_ros/ros_visualization.h>
#include <mav_trajectory_generation_ros/ros_conversions.h>

class ExamplePlanner {
 public:
  ExamplePlanner(ros::NodeHandle& nh);

  void uavPoseCallback(const geometry_msgs::Pose::ConstPtr& pose);

  void getCurrentPose(Eigen::Affine3d* current);

  void setMaxSpeed(double max_v);

  // Plans a trajectory to take off from the current position and
  // fly to the given altitude (while maintaining x,y, and yaw).
  bool planTakeOffTrajectory(double altitude);

  // Visualize the last planned trajectory
  void visualizeCachedTrajectory();

  // Output last planned trajectory
  void executeCachedTrajectory();

 private:
  ros::Publisher pub_markers_;
  ros::Publisher pub_trajectory_;
  ros::Subscriber sub_pose_;

  ros::NodeHandle& nh_;
  Eigen::Affine3d current_position_;
  double max_v_; // m/s
  double max_a_; // m/s^2
  double safety_altitude_; // m above take-off height.
  mav_trajectory_generation::Trajectory trajectory_;

};