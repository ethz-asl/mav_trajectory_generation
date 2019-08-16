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

#include <eigen_conversions/eigen_msg.h>
#include <mav_msgs/conversions.h>
#include <mav_visualization/helpers.h>

#include "mav_trajectory_generation/trajectory_sampling.h"
#include "mav_trajectory_generation_ros/ros_visualization.h"

namespace mav_trajectory_generation {

// Declare/define some internal functions that we don't want to expose.
namespace internal {

void appendMarkers(const visualization_msgs::MarkerArray& markers_to_insert,
                   const std::string& marker_namespace,
                   visualization_msgs::MarkerArray* marker_array) {
  marker_array->markers.reserve(marker_array->markers.size() +
                                markers_to_insert.markers.size());
  for (const visualization_msgs::Marker& marker : markers_to_insert.markers) {
    marker_array->markers.push_back(marker);
    if (!marker_namespace.empty()) {
      marker_array->markers.back().ns = marker_namespace;
    }
  }
}

// Overwrites the given properties of the marker array.
void setMarkerProperties(const std_msgs::Header& header, double life_time,
                         const visualization_msgs::Marker::_action_type& action,
                         visualization_msgs::MarkerArray* markers) {
  CHECK_NOTNULL(markers);
  int count = 0;
  for (visualization_msgs::Marker& marker : markers->markers) {
    marker.header = header;
    marker.action = action;
    marker.id = count;
    marker.lifetime = ros::Duration(life_time);
    ++count;
  }
}

}  // end namespace internal

static constexpr double kDefaultSamplingTime = 0.1;

void drawMavTrajectory(const Trajectory& trajectory, double distance,
                       const std::string& frame_id,
                       visualization_msgs::MarkerArray* marker_array) {
  // This is just an empty extra marker that doesn't draw anything.
  mav_visualization::MarkerGroup dummy_marker;
  return drawMavTrajectoryWithMavMarker(trajectory, distance, frame_id,
                                        dummy_marker, marker_array);
}

void drawMavSampledTrajectorybyTime(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double dt,
    const std::string& frame_id, visualization_msgs::MarkerArray* marker_array){

      mav_msgs::EigenTrajectoryPointVector filtered_vector;
      uint32_t dt_ns = ((uint32_t)(dt * 1e9));
      uint32_t last_time_ns = 0;
      filtered_vector.push_back(trajectory_points.front());
      for(size_t i = 1; i < trajectory_points.size() -1; i++){
        uint32_t current_time_ns = trajectory_points[i].time_from_start_ns;

        if(current_time_ns - last_time_ns >= dt_ns){
          filtered_vector.push_back(trajectory_points[i]);
          last_time_ns = current_time_ns;
        }

      }

      filtered_vector.push_back(trajectory_points.back());

      // filter by time.
      mav_visualization::MarkerGroup dummy_marker;
      return drawMavSampledTrajectoryWithMavMarker(filtered_vector, 0.0, frame_id,
                                               dummy_marker, marker_array);
    }


void drawMavSampledTrajectory(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double distance,
    const std::string& frame_id,
    visualization_msgs::MarkerArray* marker_array) {
  // This is just an empty extra marker that doesn't draw anything.
  mav_visualization::MarkerGroup dummy_marker;
  return drawMavSampledTrajectoryWithMavMarker(trajectory_points, distance, frame_id,
                                               dummy_marker, marker_array);
}

void drawMavTrajectoryWithMavMarker(
    const Trajectory& trajectory, double distance, const std::string& frame_id,
    const mav_visualization::MarkerGroup& additional_marker,
    visualization_msgs::MarkerArray* marker_array) {
  // Sample the trajectory.
  mav_msgs::EigenTrajectoryPoint::Vector trajectory_points;

  bool success =
      sampleWholeTrajectory(trajectory, kDefaultSamplingTime, &trajectory_points);
  if (!success) {
    return;
  }
  // Draw the trajectory.
  drawMavSampledTrajectoryWithMavMarker(trajectory_points, distance, frame_id,
                                        additional_marker, marker_array);
}

void drawMavSampledTrajectoryWithMavMarker(
    const mav_msgs::EigenTrajectoryPoint::Vector& trajectory_points, double distance,
    const std::string& frame_id,
    const mav_visualization::MarkerGroup& additional_marker,
    visualization_msgs::MarkerArray* marker_array) {
  CHECK_NOTNULL(marker_array);
  marker_array->markers.clear();

  visualization_msgs::Marker line_strip;
  line_strip.type = visualization_msgs::Marker::LINE_STRIP;
  line_strip.color = mav_visualization::Color::Orange();
  line_strip.scale.x = 0.01;
  line_strip.ns = "path";

  double accumulated_distance = 0.0;
  Eigen::Vector3d last_position = Eigen::Vector3d::Zero();
  for (size_t i = 0; i < trajectory_points.size(); ++i) {
    const mav_msgs::EigenTrajectoryPoint& trajectory_point = trajectory_points[i];

    accumulated_distance += (last_position - trajectory_point.position_W).norm();
    if (accumulated_distance > distance) {
      accumulated_distance = 0.0;
      mav_msgs::EigenMavState mav_state;
      mav_msgs::EigenMavStateFromEigenTrajectoryPoint(trajectory_point, &mav_state);
      mav_state.orientation_W_B = trajectory_point.orientation_W_B;

      visualization_msgs::MarkerArray axes_arrows;
      mav_visualization::drawAxesArrows(mav_state.position_W,
                                        mav_state.orientation_W_B, 0.3, 0.3,
                                        &axes_arrows);
      internal::appendMarkers(axes_arrows, "pose", marker_array);

      visualization_msgs::Marker arrow;
      mav_visualization::drawArrowPoints(
          trajectory_point.position_W,
          trajectory_point.position_W + trajectory_point.acceleration_W,
          mav_visualization::Color((190.0 / 255.0), (81.0 / 255.0),
                                   (80.0 / 255.0)),
          0.3, &arrow);
      arrow.ns = positionDerivativeToString(derivative_order::ACCELERATION);
      marker_array->markers.push_back(arrow);

      mav_visualization::drawArrowPoints(
          trajectory_point.position_W, trajectory_point.position_W + trajectory_point.velocity_W,
          mav_visualization::Color((80.0 / 255.0), (172.0 / 255.0),
                                   (196.0 / 255.0)),
          0.3, &arrow);
      arrow.ns = positionDerivativeToString(derivative_order::VELOCITY);
      marker_array->markers.push_back(arrow);

      mav_visualization::MarkerGroup tmp_marker(additional_marker);
      tmp_marker.transform(mav_state.position_W, mav_state.orientation_W_B);
      tmp_marker.getMarkers(marker_array->markers, 1.0, true);
    }
    last_position = trajectory_point.position_W;
    geometry_msgs::Point last_position_msg;
    tf::pointEigenToMsg(last_position, last_position_msg);
    line_strip.points.push_back(last_position_msg);
  }
  marker_array->markers.push_back(line_strip);

  std_msgs::Header header;
  header.frame_id = frame_id;
  header.stamp = ros::Time::now();
  internal::setMarkerProperties(header, 0.0, visualization_msgs::Marker::ADD,
                                marker_array);
}

void drawVertices(const Vertex::Vector& vertices, const std::string& frame_id,
                  visualization_msgs::MarkerArray* marker_array) {
  CHECK_NOTNULL(marker_array);
  marker_array->markers.resize(1);
  visualization_msgs::Marker& marker = marker_array->markers.front();

  marker.type = visualization_msgs::Marker::LINE_STRIP;
  marker.color = mav_visualization::Color::Green();
  marker.scale.x = 0.01;
  marker.ns = "straight_path";

  for (const Vertex& vertex : vertices) {
    if (vertex.D() != 3) {
      ROS_ERROR("Vertex has dimension %d but should have dimension 3.",
                vertex.D());
      return;
    }

    if (vertex.hasConstraint(derivative_order::POSITION)) {
      Eigen::VectorXd position = Eigen::Vector3d::Zero();
      vertex.getConstraint(derivative_order::POSITION, &position);
      geometry_msgs::Point constraint_msg;
      tf::pointEigenToMsg(position, constraint_msg);
      marker.points.push_back(constraint_msg);
    } else
      ROS_WARN("Vertex does not have a position constraint, skipping.");
  }

  std_msgs::Header header;
  header.frame_id = frame_id;
  header.stamp = ros::Time::now();
  internal::setMarkerProperties(header, 0.0, visualization_msgs::Marker::ADD,
                                marker_array);
}

void drawVerticesFromTrajectory(const Trajectory& trajectory,
                                const std::string& frame_id,
                                visualization_msgs::MarkerArray* marker_array) {
  CHECK_NOTNULL(marker_array);

  std::vector<double> segment_times = trajectory.getSegmentTimes();
  if (segment_times.empty()) {
    ROS_WARN("Empty trajectory.");
    return;
  }
  Vertex::Vector vertices(segment_times.size() + 1, trajectory.D());

  double trajectory_time = 0.0;
  for (size_t i = 0; i < segment_times.size(); ++i) {
    vertices[i] =
        trajectory.getVertexAtTime(trajectory_time, derivative_order::POSITION);
    trajectory_time += segment_times[i];
  }

  // Add final vertex manually to not end up in numeric errors.
  vertices.back() = trajectory.getGoalVertex(derivative_order::POSITION);

  drawVertices(vertices, frame_id, marker_array);
}

}  // namespace mav_trajectory_generation
