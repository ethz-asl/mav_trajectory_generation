#include <mav_trajectory_generation_example/example_planner.h>

ExamplePlanner::ExamplePlanner(ros::NodeHandle& nh) :
    nh_(nh),
    max_v_(2.0),
    max_a_(2.0),
    safety_altitude_(0.5),
    current_position_(Eigen::Affine3d::Identity()) {

  // create publisher for RVIZ markers
  pub_markers_ =
      nh.advertise<visualization_msgs::MarkerArray>("trajectory_markers", 0);

  pub_trajectory_ =
      nh.advertise<mav_planning_msgs::PolynomialTrajectory4D>("trajectory",
                                                              0);

  // subscriber for pose
  sub_pose_ =
      nh.subscribe("uav_pose", 1, &ExamplePlanner::uavPoseCallback, this);
}

void ExamplePlanner::uavPoseCallback(const geometry_msgs::Pose::ConstPtr& pose) {
  tf::poseMsgToEigen(*pose, current_position_);
}

void ExamplePlanner::getCurrentPose(Eigen::Affine3d* current) {
  *current = current_position_;
}

void ExamplePlanner::setMaxSpeed(const double max_v) {
  max_v_ = max_v;
}

// Plans a trajectory to take off from the current position and
// fly to the given altitude (while maintaining x,y, and yaw).
bool ExamplePlanner::planTakeOffTrajectory(const double altitude) {

  // check if actually a take off
  if (altitude < 0) {
    ROS_WARN("Not a takeoff - not executing");
    return false;
  }

  const double
      target_altitude = current_position_.translation().z() + altitude;

  bool above_safety = altitude > safety_altitude_;

  mav_trajectory_generation::Vertex::Vector vertices;
  const int dimension = 3;
  const int derivative_to_optimize =
      mav_trajectory_generation::derivative_order::SNAP;

  // we have 3 vertices:
  // Start = current position
  // middle = a point 50 cm above ground that we approach with slower velocity.
  // end = Final point that is approach with max velocity.
  mav_trajectory_generation::Vertex start(dimension), middle(dimension),
      end(dimension);

  // set start point constraints
  // (current position, and everything else zero)
  start.makeStartOrEnd(current_position_.translation(),
                       derivative_to_optimize);
  vertices.push_back(start);

  // set middle point constraints
  Eigen::Vector3d middle_point_position = current_position_.translation();

  // if we are above safety threshold, this is only an intermediate position
  if (above_safety) {
    middle_point_position.z() += safety_altitude_;
    middle.addConstraint(mav_trajectory_generation::derivative_order::POSITION,
                         middle_point_position);

    Eigen::Vector3d middle_point_velocity;
    middle_point_velocity << 0.0, 0.0, max_v_ * 0.1;
    middle.addConstraint(mav_trajectory_generation::derivative_order::VELOCITY,
                         middle_point_velocity);
  } else {
    // if not, this is the final waypoint.
    middle_point_position.z() = target_altitude;
    middle.makeStartOrEnd(middle_point_position, derivative_to_optimize);
  }
  vertices.push_back(middle);

  // plan final point if needed (to end position at rest).
  if (above_safety) {
    Eigen::Vector3d end_point_position = current_position_.translation();
    end_point_position.z() = target_altitude;
    end.makeStartOrEnd(end_point_position, derivative_to_optimize);
    vertices.push_back(end);
  }


  // compute segment times
  std::vector<double> segment_times;
  segment_times = estimateSegmentTimes(vertices, max_v_, max_a_);


  // compute segment times at mach slower speed / acceleration and use for
  // first segment.
  // a bit hacky....
  std::vector<double> segment_times_slow;
  segment_times_slow =
      estimateSegmentTimes(vertices, max_v_ * 0.2, max_a_ * 0.2);
  segment_times[0] = segment_times_slow[0];


  // solve trajectory
  const int N = 10;
  mav_trajectory_generation::PolynomialOptimization<N> opt(dimension);
  opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
  opt.solveLinear();

  // get trajectory
  trajectory_.clear();
  opt.getTrajectory(&trajectory_);

  return true;
}

void ExamplePlanner::visualizeCachedTrajectory() {
  visualization_msgs::MarkerArray markers;
  double distance =
      0.2; // Distance by which to seperate additional markers. Set 0.0 to disable.
  std::string frame_id = "world";

  mav_trajectory_generation::drawMavTrajectory(trajectory_,
                                               distance,
                                               frame_id,
                                               &markers);
  pub_markers_.publish(markers);
}

void ExamplePlanner::executeCachedTrajectory() {
  mav_planning_msgs::PolynomialTrajectory4D msg;
  mav_trajectory_generation::trajectoryToPolynomialTrajectoryMsg(trajectory_,
                                                                 &msg);
  msg.header.frame_id = "world";
  pub_trajectory_.publish(msg);
}
