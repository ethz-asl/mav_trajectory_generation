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

#ifndef MAV_TRAJECTORY_GENERATION_ROS_ROS_VISUALIZATION_H_
#define MAV_TRAJECTORY_GENERATION_ROS_ROS_VISUALIZATION_H_

#include <mav_msgs/eigen_mav_msgs.h>
#include <mav_visualization/marker_group.h>
#include <visualization_msgs/MarkerArray.h>

#include <mav_trajectory_generation/trajectory.h>
#include <mav_trajectory_generation/vertex.h>

namespace mav_trajectory_generation {

// Draws the trajectory of the MAV, with additional markers spaced by distance.
// If distance = 0.0, then these additional markers are disabled.
void drawMavTrajectory(const Trajectory& trajectory, double distance,
                       const std::string& frame_id,
                       visualization_msgs::MarkerArray* marker_array);

// Draw an eigen trajectory with additional markers spaced by distance (0.0 to
// disable).
void drawMavSampledTrajectory(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double distance,
    const std::string& frame_id, visualization_msgs::MarkerArray* marker_array);

// Draw an eigen trajectory with markers that are a minimum of dt apart in time.
void drawMavSampledTrajectorybyTime(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double dt,
    const std::string& frame_id, visualization_msgs::MarkerArray* marker_array);

// Same as drawMavTrajectory, but also draws an additional marker at a set
// distance.
void drawMavTrajectoryWithMavMarker(
    const Trajectory& trajectory, double distance, const std::string& frame_id,
    const mav_visualization::MarkerGroup& additional_marker,
    visualization_msgs::MarkerArray* marker_array);

// Draw a eigen trajectory with additional marker.
void drawMavSampledTrajectoryWithMavMarker(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double distance,
    const std::string& frame_id,
    const mav_visualization::MarkerGroup& additional_marker,
    visualization_msgs::MarkerArray* marker_array);

// Visualize original vertices.
void drawVertices(const Vertex::Vector& vertices, const std::string& frame_id,
                  visualization_msgs::MarkerArray* marker_array);

// Visualize orginal vertices from trajectory.
void drawVerticesFromTrajectory(const Trajectory& trajectory,
                                const std::string& frame_id,
                                visualization_msgs::MarkerArray* marker_array);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_ROS_ROS_CONVERSIONS_H_
