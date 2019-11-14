#include <mav_trajectory_generation_example/example_planner.h>

ExamplePlanner::ExamplePlanner(ros::NodeHandle& nh)
    : nh_(nh),
      max_v_(2.0),
      max_a_(2.0),
      max_ang_v_(2.0),
      max_ang_a_(2.0),
      current_velocity_(Eigen::Vector3d::Zero()),
      current_angular_velocity_(Eigen::Vector3d::Zero()),
      current_pose_(Eigen::Affine3d::Identity()) {
        
  // Load params
  if (!nh_.getParam(ros::this_node::getName() + "/max_v", max_v_)){
    ROS_WARN("[example_planner] param max_v not found");
  }
  if (!nh_.getParam(ros::this_node::getName() + "/max_a", max_a_)){
    ROS_WARN("[example_planner] param max_a not found");
  }
  if (!nh_.getParam(ros::this_node::getName() + "/max_ang_v", max_ang_v_)){
    ROS_WARN("[example_planner] param max_ang_v not found");
  }
  if (!nh_.getParam(ros::this_node::getName() + "/max_ang_a", max_ang_a_)){
    ROS_WARN("[example_planner] param max_ang_a not found");
  }
        
  // create publisher for RVIZ markers
  pub_markers_ =
      nh.advertise<visualization_msgs::MarkerArray>("trajectory_markers", 0);

  pub_trajectory_ =
      nh.advertise<mav_planning_msgs::PolynomialTrajectory>("trajectory",
                                                              0);

  // subscriber for Odometry
  sub_odom_ =
      nh.subscribe("uav_pose", 1, &ExamplePlanner::uavOdomCallback, this);
}

// Callback to get current Pose of UAV
void ExamplePlanner::uavOdomCallback(const nav_msgs::Odometry::ConstPtr& odom) {

  // store current position in our planner
  tf::poseMsgToEigen(odom->pose.pose, current_pose_);

  // store current velocity
  tf::vectorMsgToEigen(odom->twist.twist.linear, current_velocity_);
  tf::vectorMsgToEigen(odom->twist.twist.angular, current_angular_velocity_);
}

// Method to set maximum speed.
void ExamplePlanner::setMaxSpeed(const double max_v) {
  max_v_ = max_v;
}

// Plans a trajectory from the current position to the a goal position and velocity
bool ExamplePlanner::planTrajectory(
    const Eigen::VectorXd& goal_pos, const Eigen::VectorXd& goal_vel,
    mav_trajectory_generation::Trajectory* trajectory) {
  assert(trajectory);
  trajectory->clear();

  // 3 Dimensional trajectory => 3D position
  // 4 Dimensional trajectory => 3D position + yaw
  // 6 Dimensional trajectory => through SE(3) space, position and orientation
  const int dimension = goal_pos.size();
  bool success = false;

  if (dimension == 6) 
  {
    mav_trajectory_generation::Trajectory trajectory_trans, trajectory_rot;

    // Translation trajectory.
    Eigen::Vector3d goal_position = goal_pos.head(3);
    Eigen::Vector3d goal_lin_vel = goal_vel.head(3);
    success = planTrajectory(
        goal_position, goal_lin_vel, current_pose_.translation(),
        current_velocity_, max_v_, max_a_, &trajectory_trans);

    // Rotation trajectory.
    Eigen::Vector3d goal_rotation = goal_pos.tail(3);
    Eigen::Vector3d goal_ang_vel = goal_vel.tail(3);
    Eigen::Vector3d current_rot_vec;
    mav_msgs::vectorFromRotationMatrix(
        current_pose_.rotation(), &current_rot_vec);
    success &= planTrajectory(
        goal_rotation, goal_ang_vel, current_rot_vec, current_angular_velocity_,
        max_ang_v_, max_ang_a_, &trajectory_rot);

    // Combine trajectories.
    success &= trajectory_trans.getTrajectoryWithAppendedDimension(
            trajectory_rot, &(*trajectory));
    return success;
  } 
  else if (dimension == 3) 
  {
    success = planTrajectory(
        goal_pos, goal_vel, current_pose_.translation(), current_velocity_,
        max_v_, max_a_, &(*trajectory));
    return success;
  } 
  else if (dimension == 4) 
  {
    Eigen::Vector4d start_pos_4d, start_vel_4d;
    start_pos_4d << current_pose_.translation(),
        mav_msgs::yawFromQuaternion(
            (Eigen::Quaterniond)current_pose_.rotation());
    start_vel_4d << current_velocity_, 0.0;
    success = planTrajectory(
        goal_pos, goal_vel, start_pos_4d, start_vel_4d, max_v_, max_a_,
        &(*trajectory));
    return success;
  } 
  else 
  {
    LOG(WARNING) << "Dimension must be 3, 4 or 6 to be valid.";
    return false;
  }
}

// Plans a trajectory from a start position and velocity to a goal position and velocity
bool ExamplePlanner::planTrajectory(const Eigen::VectorXd& goal_pos,
                                    const Eigen::VectorXd& goal_vel,
                                    const Eigen::VectorXd& start_pos,
                                    const Eigen::VectorXd& start_vel,
                                    double v_max, double a_max,
                                    mav_trajectory_generation::Trajectory* trajectory) {
  assert(trajectory);
  const int dimension = goal_pos.size();
  // Array for all waypoints and their constraints
  mav_trajectory_generation::Vertex::Vector vertices;

  // Optimze up to 4th order derivative (SNAP)
  const int derivative_to_optimize =
      mav_trajectory_generation::derivative_order::SNAP;

  // we have 2 vertices:
  // start = desired start vector
  // end = desired end vector
  mav_trajectory_generation::Vertex start(dimension), end(dimension);

  /******* Configure start point *******/
  start.makeStartOrEnd(start_pos, derivative_to_optimize);
  start.addConstraint(mav_trajectory_generation::derivative_order::VELOCITY,
                      start_vel);
  vertices.push_back(start);

  /******* Configure end point *******/
  // set end point constraints to desired position and set all derivatives to zero
  end.makeStartOrEnd(goal_pos, derivative_to_optimize);
  end.addConstraint(mav_trajectory_generation::derivative_order::VELOCITY,
                    goal_vel);
vertices.push_back(end);

  // setimate initial segment times
  std::vector<double> segment_times;
  segment_times = estimateSegmentTimes(vertices, v_max, a_max);

  // Set up polynomial solver with default params
  mav_trajectory_generation::NonlinearOptimizationParameters parameters;

  // set up optimization problem
  const int N = 10;
  mav_trajectory_generation::PolynomialOptimizationNonLinear<N> opt(dimension, parameters);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);

  // constrain velocity and acceleration
  opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::VELOCITY, v_max);
  opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::ACCELERATION, a_max);

  // solve trajectory
  opt.optimize();

  // get trajectory as polynomial parameters
  opt.getTrajectory(&(*trajectory));
  trajectory->scaleSegmentTimesToMeetConstraints(v_max, a_max);
  
  return true;
}
                                    

bool ExamplePlanner::publishTrajectory(const mav_trajectory_generation::Trajectory& trajectory){
  // send trajectory as markers to display them in RVIZ
  visualization_msgs::MarkerArray markers;
  double distance =
      0.2; // Distance by which to seperate additional markers. Set 0.0 to disable.
  std::string frame_id = "world";

  mav_trajectory_generation::drawMavTrajectory(trajectory,
                                               distance,
                                               frame_id,
                                               &markers);
  pub_markers_.publish(markers);

  // send trajectory to be executed on UAV
  mav_planning_msgs::PolynomialTrajectory msg;
  mav_trajectory_generation::trajectoryToPolynomialTrajectoryMsg(trajectory,
                                                                 &msg);
  msg.header.frame_id = "world";
  pub_trajectory_.publish(msg);

  return true;
}

