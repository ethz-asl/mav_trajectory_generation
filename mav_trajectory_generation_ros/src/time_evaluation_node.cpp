#include <ros/ros.h>

#include <mav_visualization/helpers.h>
#include <mav_trajectory_generation/polynomial_optimization_linear.h>
#include <mav_trajectory_generation/polynomial_optimization_nonlinear.h>
#include <mav_trajectory_generation/timing.h>

#include "mav_trajectory_generation_ros/ros_conversions.h"
#include "mav_trajectory_generation_ros/ros_visualization.h"

namespace mav_trajectory_generation {

// Benchmarking utilities to evaluate different methods of time allocation for
// polynomial trajectories.

struct TimeAllocationBenchmarkResult {
  TimeAllocationBenchmarkResult()
      : trial_number(-1),
        method_name("none"),
        num_segments(0),
        nominal_length(0.0),
        optimization_success(false),
        bounds_violated(false),
        trajectory_time(0.0),
        trajectory_length(0.0),
        computation_time(0.0) {}

  // Evaluation settings
  int trial_number;
  std::string method_name;

  // Trajectory settings
  int num_segments;
  double nominal_length;

  // Evaluation results
  bool optimization_success;
  bool bounds_violated;
  double trajectory_time;
  double trajectory_length;
  double computation_time;
  double a_max_actual;
  double v_max_actual;

  // More to come: convex hull/bounding box, etc.
};

class TimeEvaluationNode {
 public:
  TimeEvaluationNode(const ros::NodeHandle& nh,
                     const ros::NodeHandle& nh_private);

  // Number of Coefficients
  const static int kN = 10;  // has to be even !!
  // Dimension
  const static int kDim = 3;

  // Running the actual benchmark, one trial at a time (so that it can be
  // paused between for visualization).
  void runBenchmark(int trial_number, int num_segments);

  // Generate trajectories with different methods.
  void runNfabian(const Vertex::Vector& vertices, Trajectory* trajectory) const;
  void runNonlinear(const Vertex::Vector& vertices,
                    Trajectory* trajectory) const;

  void evaluateTrajectory(const std::string& method_name,
                          const Trajectory& traj,
                          TimeAllocationBenchmarkResult* result) const;

  // Accessors.
  bool visualize() const { return visualize_; }

 private:
  ros::NodeHandle nh_;
  ros::NodeHandle nh_private_;

  // General settings.
  std::string frame_id_;
  bool visualize_;

  // Dynamic constraints.
  double v_max_;
  double a_max_;

  // General trajectory settings.
  int max_derivative_order_;

  // Store all the results.
  std::vector<TimeAllocationBenchmarkResult> results_;

  // ROS stuff.
  ros::Publisher path_marker_pub_;
};

TimeEvaluationNode::TimeEvaluationNode(const ros::NodeHandle& nh,
                                       const ros::NodeHandle& nh_private)
    : nh_(nh),
      nh_private_(nh_private),
      frame_id_("world"),
      visualize_(true),
      v_max_(1.0),
      a_max_(2.0),
      max_derivative_order_(derivative_order::JERK) {
  nh_private_.param("frame_id", frame_id_, frame_id_);
  nh_private_.param("visualize", visualize_, visualize_);
  nh_private_.param("v_max", v_max_, v_max_);
  nh_private_.param("a_max", a_max_, a_max_);

  path_marker_pub_ =
      nh_private_.advertise<visualization_msgs::MarkerArray>("path", 1, true);
}

void TimeEvaluationNode::runBenchmark(int trial_number, int num_segments) {
  srand(trial_number);

  const Eigen::VectorXd min_pos = Eigen::VectorXd::Constant(kDim, -5.0);
  const Eigen::VectorXd max_pos = -min_pos;

  // Use trial number as seed to create the trajectory.
  Vertex::Vector vertices;
  vertices = createRandomVertices(max_derivative_order_, num_segments, min_pos,
                                  max_pos, trial_number);

  TimeAllocationBenchmarkResult result;
  // Fill in all the basics in the results that are shared between all the
  // evaluations.
  result.trial_number = trial_number;
  result.num_segments = num_segments;
  double nominal_length = 0.0;
  // TODO(helenol): compute nominal length from the vertices...

  // Run all the evaluations.
  Trajectory trajectory_nfabian;
  runNfabian(vertices, &trajectory_nfabian);
  evaluateTrajectory("nfabian", trajectory_nfabian, &result);
  results_.push_back(result);

  Trajectory trajectory_nonlinear;
  runNonlinear(vertices, &trajectory_nonlinear);
  evaluateTrajectory("nonlinear", trajectory_nonlinear, &result);
  results_.push_back(result);
}

void TimeEvaluationNode::runNfabian(const Vertex::Vector& vertices,
                                    Trajectory* trajectory) const {
  std::vector<double> segment_times;
  segment_times =
      mav_trajectory_generation::estimateSegmentTimes(vertices, v_max_, a_max_);

  mav_trajectory_generation::PolynomialOptimization<kN> linopt(kDim);
  linopt.setupFromVertices(vertices, segment_times, max_derivative_order_);
  linopt.solveLinear();
  linopt.getTrajectory(trajectory);
}

void TimeEvaluationNode::runNonlinear(const Vertex::Vector& vertices,
                                      Trajectory* trajectory) const {
  std::vector<double> segment_times;
  segment_times =
      mav_trajectory_generation::estimateSegmentTimes(vertices, v_max_, a_max_);

  mav_trajectory_generation::NonlinearOptimizationParameters nlopt_parameters;
  mav_trajectory_generation::PolynomialOptimizationNonLinear<kN> nlopt(
      kDim, nlopt_parameters, false);
  nlopt.setupFromVertices(vertices, segment_times, max_derivative_order_);
  nlopt.addMaximumMagnitudeConstraint(derivative_order::VELOCITY, v_max_);
  nlopt.addMaximumMagnitudeConstraint(derivative_order::ACCELERATION, a_max_);
  nlopt.optimize();
  nlopt.getTrajectory(trajectory);
}

void TimeEvaluationNode::evaluateTrajectory(
    const std::string& method_name, const Trajectory& traj,
    TimeAllocationBenchmarkResult* result) const {
  result->method_name = method_name;

  result->trajectory_time = traj.getMaxTime();

  // TODO(helenol): evaluate the path length, min/max extrema, etc.
}

}  // namespace mav_trajectory_generation

int main(int argc, char** argv) {
  ros::init(argc, argv, "time_evaluation_node");
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  ros::NodeHandle nh("");
  ros::NodeHandle nh_private("~");

  mav_trajectory_generation::TimeEvaluationNode time_eval_node(nh, nh_private);

  ROS_INFO("Initialized time evaluation node.");

  ros::spin();
  return 0;
}
