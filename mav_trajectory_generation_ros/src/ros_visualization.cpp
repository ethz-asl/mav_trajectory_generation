/*
* Copyright (c) 2015, Markus Achtelik, ASL, ETH Zurich, Switzerland
* You can contact the author at <markus dot achtelik at mavt dot ethz dot ch>
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

#include "mav_trajectory_generation_ros/ros_visualization.h"
#include "mav_trajectory_generation_ros/trajectory_sampling.h"

namespace mav_trajectory_generation {

// Declare/define some internal functions that we don't want to expose.
namespace internal {

void appendMarkers(const visualization_msgs::MarkerArray& markers_to_insert,
                   const std::string& marker_namespace,
                   visualization_msgs::MarkerArray* marker_array) {
  marker_array->markers.reserve(marker_array->markers.size() +
                                markers_to_insert.markers.size());
  for (const auto marker : markers_to_insert.markers) {
    marker_array->markers.push_back(marker);
    if (!marker_namespace.empty()) {
      marker_array->markers.back().ns = marker_namespace;
    }
  }
}

/**
 * \brief Overwrites the given properties of the marker array.
 */
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

void drawMavTrajectoryWithMavMarker(
    const Trajectory& trajectory, double distance, const std::string& frame_id,
    const mav_visualization::MarkerGroup& additional_marker,
    visualization_msgs::MarkerArray* marker_array) {
  CHECK_NOTNULL(marker_array);
  marker_array->markers.clear();

  // Sample the trajectory.
  mav_msgs::EigenTrajectoryPoint::Vector flat_states;
  bool success =
      sampleWholeTrajectory(trajectory, kDefaultSamplingTime, &flat_states);
  if (!success) {
    return;
  }

  visualization_msgs::Marker line_strip;
  line_strip.type = visualization_msgs::Marker::LINE_STRIP;
  line_strip.color = mav_visualization::createColorRGBA(1, 0.5, 0, 1);
  line_strip.scale.x = 0.01;
  line_strip.ns = "path";

  double accumulated_distance = 0;
  Eigen::Vector3d last_position = Eigen::Vector3d::Zero();
  for (size_t i = 0; i < flat_states.size(); ++i) {
    const mav_msgs::EigenTrajectoryPoint& flat_state = flat_states[i];

    accumulated_distance += (last_position - flat_state.position_W).norm();
    if (accumulated_distance > distance) {
      accumulated_distance = 0.0;
      mav_msgs::EigenMavState mav_state;
      mav_msgs::EigenMavStateFromEigenTrajectoryPoint(flat_state, &mav_state);

      visualization_msgs::MarkerArray axes_arrows;
      mav_visualization::drawAxesArrows(mav_state.position_W,
                                        mav_state.orientation_W_B, 0.3, 0.3,
                                        &axes_arrows);
      internal::appendMarkers(axes_arrows, "pose", marker_array);

      visualization_msgs::Marker arrow;
      mav_visualization::drawArrowPoints(
          flat_state.position_W,
          flat_state.position_W + flat_state.acceleration_W,
          mav_visualization::createColorRGBA((190.0 / 255.0), (81.0 / 255.0),
                                             (80.0 / 255.0), 1),
          0.3, &arrow);
      arrow.ns = positionDerivativeToString(derivative_order::ACCELERATION);
      marker_array->markers.push_back(arrow);

      mav_visualization::drawArrowPoints(
          flat_state.position_W, flat_state.position_W + flat_state.velocity_W,
          mav_visualization::createColorRGBA((80.0 / 255.0), (172.0 / 255.0),
                                             (196.0 / 255.0), 1),
          0.3, &arrow);
      arrow.ns = positionDerivativeToString(derivative_order::VELOCITY);
      marker_array->markers.push_back(arrow);

      mav_visualization::MarkerGroup tmp_marker(additional_marker);
      tmp_marker.transform(mav_state.position_W, mav_state.orientation_W_B);
      tmp_marker.getMarkers(marker_array->markers, 1.0, true);
    }
    last_position = flat_state.position_W;
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
  marker.color = mav_visualization::createColorRGBA(0.5, 1.0, 0, 1);
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

}  // namespace mav_trajectory_generation
