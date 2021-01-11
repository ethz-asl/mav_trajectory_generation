# mav_trajectory_generation
This repository contains tools for polynomial trajectory generation and optimization based on methods described in [1].
These techniques are especially suitable for rotary-wing micro aerial vehicles (MAVs).
This README provides a brief overview of our trajectory generation utilities with some examples.

**Authors**: Markus Achtelik, Michael Burri, Helen Oleynikova, Rik Bähnemann, Marija Popović  
**Maintainer**: Rik Bähnemann, brik@ethz.ch  
**Affiliation**: Autonomous Systems Lab, ETH Zurich  

## Bibliography
This implementation is largely based on the work of C. Richter *et al*, who should be cited if this is used in a scientific publication (or the preceding conference papers):  
[1] C. Richter, A. Bry, and N. Roy, “**Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments,**” in *International Journal of Robotics Research*, Springer, 2016.
```
@incollection{richter2016polynomial,
  title={Polynomial trajectory planning for aggressive quadrotor flight in dense indoor environments},
  author={Richter, Charles and Bry, Adam and Roy, Nicholas},
  booktitle={Robotics Research},
  pages={649--666},
  year={2016},
  publisher={Springer}
}
```

Furthermore, the nonlinear optimization features our own extensions, described in:  

Michael Burri, Helen Oleynikova, Markus Achtelik, and Roland Siegwart, “**Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Previously Unknown Environments**”. In *IEEE Int. Conf. on Intelligent Robots and Systems* (IROS), September 2015.
```
@inproceedings{burri2015real-time,
  author={Burri, Michael and Oleynikova, Helen and  and Achtelik, Markus W. and Siegwart, Roland},
  booktitle={Intelligent Robots and Systems (IROS 2015), 2015 IEEE/RSJ International Conference on},
  title={Real-Time Visual-Inertial Mapping, Re-localization and Planning Onboard MAVs in Unknown Environments},
  year={2015},
  month={Sept}
}
```

## Installation Instructions (Ubuntu)
To install this package with [ROS Indigo](http://wiki.ros.org/indigo/Installation/Ubuntu) or [ROS Kinetic](http://wiki.ros.org/kinetic/Installation/Ubuntu):

1. Install additional system dependencies (swap indigo for kinetic or melodic as necessary):

** Note: ROS melodic requires libyaml-cpp-dev and does not build with yaml_cpp_catkin in your catkin workspace!

```
sudo apt-get install python-wstool python-catkin-tools ros-indigo-cmake-modules libyaml-cpp-dev
```

2. Set up a catkin workspace (if not already done):

```
mkdir -p ~/catkin_ws/src
cd ~/catkin_ws
catkin init
catkin config --extend /opt/ros/indigo
catkin config --cmake-args -DCMAKE_BUILD_TYPE=Release
catkin config --merge-devel
```

3. Install the repository and its dependencies (with rosinstall):

```
cd src
wstool init
wstool set --git mav_trajectory_generation git@github.com:ethz-asl/mav_trajectory_generation.git -y
wstool update
wstool merge mav_trajectory_generation/install/mav_trajectory_generation_https.rosinstall
wstool update -j8
echo "source ~/catkin_ws/devel/setup.bash" >> ~/.bashrc
source /opt/ros/indigo/setup.bash
```

In case you have your SSH keys for github set up, feel free to use the ssh rosinstall instead:
```
wstool merge mav_trajectory_generation/install/mav_trajectory_generation_ssh.rosinstall
```

4. Use [catkin_build](http://catkin-tools.readthedocs.io/en/latest/verbs/catkin_build.html) to build the repository:

```
catkin build mav_trajectory_generation_ros
```


## Basics
A **vertex** describes the properties of a support point of a **polynomial** path. Pairs of vertices are connected together to form **segments**.
Each vertex has a set of constraints: the values of position derivatives that must be matched during optimization procedures.
Typically, only the positions are specified for vertices along a path, while start and end vertices have all derivatives of position set to zero.
```
  x----------x-----------------x
vertex            segment
```
In the case of multi-dimensional vertices, the derivative constraints exist in all dimensions, with possibly different values.

## Linear Optimization
In this section, we consider how to generate polynomial segments passing through a set of arbitrary vertices using the unconstrained **linear** optimization approach described in [1].

Necessary includes:
```c++
#include <mav_trajectory_generation/polynomial_optimization_linear.h>
```

1. Create a list of three (x,y,z) vertices to fly through, e.g. (0,0,1) -> (1,2,3) -> (2,1,5), and define some parameters. The ``dimension`` variable denotes the spatial dimension of the path (here, 3D). The ``derivative_to_optimize`` should usually be set to the last derivative that should be continuous (here, snap).

```c++
mav_trajectory_generation::Vertex::Vector vertices;
const int dimension = 3;
const int derivative_to_optimize = mav_trajectory_generation::derivative_order::SNAP;
mav_trajectory_generation::Vertex start(dimension), middle(dimension), end(dimension);
```

2. Add constraints to the vertices.

```c++
start.makeStartOrEnd(Eigen::Vector3d(0,0,1), derivative_to_optimize);
vertices.push_back(start);

middle.addConstraint(mav_trajectory_generation::derivative_order::POSITION, Eigen::Vector3d(1,2,3));
vertices.push_back(middle);

end.makeStartOrEnd(Eigen::Vector3d(2,1,5), derivative_to_optimize);
vertices.push_back(end);
```

3. Compute the segment times.

```c++
std::vector<double> segment_times;
const double v_max = 2.0;
const double a_max = 2.0;
segment_times = estimateSegmentTimes(vertices, v_max, a_max);
```

4. Create an optimizer object and solve. The template parameter (N) denotes the number of coefficients of the underlying polynomial, which has to be even. If we want the trajectories to be snap-continuous, N needs to be at least 10; for minimizing jerk, 8.

```c++
const int N = 10;
mav_trajectory_generation::PolynomialOptimization<N> opt(dimension);
opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
opt.solveLinear();
```

5. Obtain the polynomial segments.

```c++
mav_trajectory_generation::Segment::Vector segments;
opt.getSegments(&segments);
```

## Nonlinear Optimization
In this section, we consider how to generate polynomial segments passing through a set of arbitrary vertices using the unconstrained **nonlinear** optimization approach described in [1]. The same approach is followed as in the previous section.

Necessary includes:

```c++
#include <mav_trajectory_generation/polynomial_optimization_nonlinear.h>
```

1. Set up the problem by following Steps 1-3 in the **Linear Optimization** section.

2. Set the parameters for nonlinear optimization. Below is an example, but the default parameters should be reasonable enough to use without fine-tuning.

```c++
NonlinearOptimizationParameters parameters;
parameters.max_iterations = 1000;
parameters.f_rel = 0.05;
parameters.x_rel = 0.1;
parameters.time_penalty = 500.0;
parameters.initial_stepsize_rel = 0.1;
parameters.inequality_constraint_tolerance = 0.1;
```

3. Create an optimizer object and solve. The third argument of the optimization object (true/false) specifies whether the optimization is run over the segment times only.

```c++
const int N = 10;
PolynomialOptimizationNonLinear<N> opt(dimension, parameters, false);
opt.setupFromVertices(vertices, segment_times, derivative_to_optimize);
opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::VELOCITY, v_max);                                opt.addMaximumMagnitudeConstraint(mav_trajectory_generation::derivative_order::ACCELERATION, a_max);
opt.optimize();
```

4. Obtain the polynomial segments.

```c++
mav_trajectory_generation::Segment::Vector segments;
opt.getPolynomialOptimizationRef().getSegments(&segments);
```
## Creating Trajectories
In this section, we consider how to use our trajectory optimization results. We first need to convert our optimization object into the Trajectory class:

```c++
#include <mav_trajectory_generation/trajectory.h>

mav_trajectory_generation::Trajectory trajectory;
opt.getTrajectory(&trajectory);
```

We can also create new trajectories by splitting (getting a trajectory with a single dimension) or compositing (getting a trajectory by appending another trajectory:

```c++
// Splitting:
mav_trajectory_generation::Trajectory x_trajectory = trajectory.getTrajectoryWithSingleDimension(1);

// Compositing:
mav_trajectory_generation::Trajectory trajectory_with_yaw; trajectory.getTrajectoryWithAppendedDimension(yaw_trajectory, &trajectory_with_yaw);
```

## Sampling Trajectories
In this section, we consider methods of evaluating the trajectory at particular instances of time. There are two methods of doing this.

1. By evaluating directly from the Trajectory class.

```c++
// Single sample:
double sampling_time = 2.0;
int derivative_order = mav_trajectory_generation::derivative_order::POSITION;
Eigen::VectorXd sample = trajectory.evaluate(sampling_time, derivative_order);

// Sample range:
double t_start = 2.0;
double t_end = 10.0;
double dt = 0.01;
std::vector<Eigen::VectorXd> result;
std::vector<double> sampling_times; // Optional.
trajectory.evaluateRange(t_start, t_end, dt, derivative_order, &result, &sampling_times);
```

2. By conversion to ```mav_msgs::EigenTrajectoryPoint``` state(s). These functions support 3D or 4D trajectories (the 4th dimension is assumed to be yaw if it exists).

```c++
#include <mav_trajectory_generation_ros/trajectory_sampling.h>

mav_msgs::EigenTrajectoryPoint state;
mav_msgs::EigenTrajectoryPoint::Vector states;

// Single sample:
double sampling_time = 2.0;
bool success = mav_trajectory_generation::sampleTrajectoryAtTime(trajectory, sample_time, &state);

// Sample range:
double t_start = 2.0;
double duration = 10.0;
double dt = 0.01;
success = mav_trajectory_generation::sampleTrajectoryInRange(trajectory, t_start, duration, dt, &states);

// Whole trajectory:
double sampling_interval = 0.01;
success = mav_trajectory_generation::sampleWholeTrajectory(trajectory, sampling_interval, &states);
```

## Visualizing Trajectories
In this section, we describe how to visualize trajectories in [rviz](http://wiki.ros.org/rviz).

<p align="center"><img src="http://i68.tinypic.com/2cp5oxs.png" height="350"/></p>

For a simple visualization:

```c++
#include <mav_trajectory_generation_ros/ros_visualization.h>

visualization_msgs::MarkerArray markers;
double distance = 1.0; // Distance by which to seperate additional markers. Set 0.0 to disable.
std::string frame_id = "world";

// From Trajectory class:
mav_trajectory_generation::drawMavTrajectory(trajectory, distance, frame_id, &markers);

// From mav_msgs::EigenTrajectoryPoint::Vector states:
mav_trajectory_generation::drawMavSampledTrajectory(states, distance, frame_id, &markers)
```

For a visualization including an additional marker at a set distance (e.g. hexacopter marker):


```c++
mav_visualization::HexacopterMarker hex(simple);

// From Trajectory class:
mav_trajectory_generation::drawMavTrajectoryWithMavMarker(trajectory, distance, frame_id, hex &markers);

// From mav_msgs::EigenTrajectoryPoint::Vector states:
mav_trajectory_generation::drawMavSampledTrajectoryWithMavMarker(states, distance, frame_id, hex, &markers)
```

## Checking Input Feasibility
The package contains three implementations to check generated trajectories for
input feasibility. The checks are based on the rigid-body model assumption and
flat state characteristics presented in [Mellinger2011](http://www-personal.acfr.usyd.edu.au/spns/cdm/papers/Mellinger.pdf).

```
@inproceedings{mellinger2011minimum,
  title={Minimum snap trajectory generation and control for quadrotors},
  author={Mellinger, Daniel and Kumar, Vijay},
  booktitle={Robotics and Automation (ICRA), 2011 IEEE International Conference on},
  pages={2520--2525},
  year={2011},
  organization={IEEE}
}
```

The trajectories are checked for low and high thrust, high velocities, high roll
and pitch rates, high yaw rates and high yaw angular accelerations.

`FeasibilitySampling` implements a naive sampling-based check.
`FeasibilityRecursive` implements a slightly adapted recursive feasibility test
presented in [Müller2015](http://flyingmachinearena.org/wp-content/publications/2015/mueTRO15.pdf).
```
@article{mueller2015computationally,
  title={A computationally efficient motion primitive for quadrocopter trajectory generation},
  author={Mueller, Mark W and Hehn, Markus and D'Andrea, Raffaello},
  journal={IEEE Transactions on Robotics},
  volume={31},
  number={6},
  pages={1294--1310},
  year={2015},
  publisher={IEEE}
}
```
We adapted [RapidQuadrocopterTrajectories](https://github.com/markwmuller/RapidQuadrocopterTrajectories) to
check arbitrary polynomial order trajectories for yaw and velocity feasibility.

`FeasibilityAnalytic` analytically checks the magnitudes except for the roll and
pitch rates, where it runs the recursive test (recommended: low in computation
time, no false positives).

### Example
```c++
// Create input constraints.
typedef InputConstraintType ICT;
InputConstraints input_constraints;
input_constraints.addConstraint(ICT::kFMin, 0.5 * 9.81); // minimum acceleration in [m/s/s].
input_constraints.addConstraint(ICT::kFMax, 1.5 * 9.81); // maximum acceleration in [m/s/s].
input_constraints.addConstraint(ICT::kVMax, 3.5); // maximum velocity in [m/s].
input_constraints.addConstraint(ICT::kOmegaXYMax, M_PI / 2.0); // maximum roll/pitch rates in [rad/s].
input_constraints.addConstraint(ICT::kOmegaZMax, M_PI / 2.0); // maximum yaw rates in [rad/s].
input_constraints.addConstraint(ICT::kOmegaZDotMax, M_PI); // maximum yaw acceleration in [rad/s/s].

// Create feasibility object of choice (FeasibilityAnalytic,
// FeasibilitySampling, FeasibilityRecursive).
FeasibilityAnalytic feasibility_check(input_constraints);
feasibility_check.settings_.setMinSectionTimeS(0.01);

// Check feasibility.
Segment dummy_segment;
InputFeasibilityResult result =
    feasibility_check.checkInputFeasibility(dummy_segment);
std::cout << "The segment input is " << getInputFeasibilityResultName(result);
<< "." << std::endl;
```

### Benchmarking
Both recursive and analytic checks are comparably fast.
The recursive check may have a couple more false negatives, i.e., segments, that
can not be determined feasible although they are. But this is neglectable.
The sampling based check is both slow and may have false positives, i.e.,
consider segments feasible although they are not. We do not recommend using
this.

Here are the computational results over 1000 random segments with different
parameter settings:
<p align="center"><img src="https://cloud.githubusercontent.com/assets/11293852/26199226/903dea9a-3bc9-11e7-945e-91e5a119e63f.png" /></p>

## Checking Half-Space Feasibility
The package also contains a method to check if a trajectory or segment is inside
an arbitrary set of half spaces based on [RapidQuadrocopterTrajectories](https://github.com/markwmuller/RapidQuadrocopterTrajectories).
This is useful to check if a segment is inside a box or above ground.

Example ground plane feasibility:
```c++
// Create feasibility check.
FeasibilityBase feasibility_check;
// Create ground plane.
Eigen::Vector3d point(0.0, 0.0, 0.0);
Eigen::Vector3d normal(0.0, 0.0, 1.0);
feasibility_check.half_plane_constraints_.emplace_back(point, normal);
// Check feasibility.
Segment dummy_segment;
if(!feasibility_check.checkHalfPlaneFeasibility(segment)) {
  std::cout << "The segment is in collision with the ground plane." << std::endl;
}
```

Example box feasibility:
```c++
// Create feasibility check.
FeasibilityBase feasibility_check;
// Create box constraints.
Eigen::Vector3d box_center(0.0, 0.0, 0.0);
Eigen::Vector3d box_size(1.0, 1.0, 1.0);
feasibility_check.half_plane_constraints_ =
    HalfPlane::createBoundingBox(box_center, box_size);
// Check feasibility.
Segment dummy_segment;
if(!feasibility_check.checkHalfPlaneFeasibility(segment)) {
  std::cout << "The segment is not inside the box." << std::endl;
}
```
