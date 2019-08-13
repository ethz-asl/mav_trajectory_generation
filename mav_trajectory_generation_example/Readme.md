# Overview
This is a minimal example for a trajectory generator that outputs a trajectory to a specified goal state (position and velocity).
The initial start state is based on the current odometry of the UAV in the simulation. In reality, this can be current reference or the output of a state estimater too.

# How to use
After building, launch the included launch file:
roslaunch mav_trajectory_generation_example example.launch

The following things are launched:
* Gazebo for the UAV simulation
* The trajectory generator called ExamplePlanner that, based on console input, publishes visualization markers and a parametrized trajectory (not a sampled one yet). 
* The trajectory sampler that subscribes to the parametrized trajectory and samples the resulting setpoints (position, orientation, velocity etc) at 100 hz and sends them to the controller
* RVIZ, that visualizes the trajectory markers published by the generator.

After starting the launch file, wait for the console to clear and show this message:
![console](https://raw.githubusercontent.com/ethz-asl/mav_trajectory_generation/feature/example_planner/mav_trajectory_generation_example/img/console.png)


Once you press enter in that console, the trajectory is generated and sent to the controller. You should now see the UAV move in Gazebo. Simultaenously, RVIZ should display the planned trajectory:
![rviz](https://raw.githubusercontent.com/ethz-asl/mav_trajectory_generation/feature/example_planner/mav_trajectory_generation_example/img/traject_rviz.png)


# Node Graph
The following image visualizes the node graph of the simulation:
![nodes](https://raw.githubusercontent.com/ethz-asl/mav_trajectory_generation/feature/example_planner/mav_trajectory_generation_example/img/nodes.png)



* Red is the planner that waits on the console and publishes the markers (orange) and the parametrized trajectory (blue)
* Green is the sampler that subscribes to the parametrized trajectory (blue) and publishes sampled trajectories to the controller (cyan / light blue).

# Note
Please note that waiting for console input is not the proper way to trigger trajectory generation. Depending on your use case, the planner should be triggered by services or based on information of other messages that it subscribes to etc.
