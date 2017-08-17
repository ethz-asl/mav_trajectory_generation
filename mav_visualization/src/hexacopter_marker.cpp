/*
 * Copyright (c) 2016, Markus Achtelik, ASL, ETH Zurich, Switzerland
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

#include "mav_visualization/hexacopter_marker.h"

namespace mav_visualization {

HexacopterMarker::HexacopterMarker(bool simple) : MarkerGroup() {
  createHexacopter(simple);
}

void HexacopterMarker::createHexacopter(bool simple) {
  const double sqrt2_2 = sqrt(2) / 2;

  visualization_msgs::Marker marker;
  marker.id = 0;

  // Make rotors.
  marker.type = visualization_msgs::Marker::CYLINDER;
  marker.ns = "hexacopter";
  marker.scale.x = 0.2;
  marker.scale.y = 0.2;
  marker.scale.z = 0.01;
  marker.color.r = 0.8;
  marker.color.g = 0.5;
  marker.color.b = 0.0;
  marker.color.a = 0.5;
  marker.pose.position.z = 0;
  marker.pose.orientation.w = 1;
  marker.pose.orientation.x = 0;
  marker.pose.orientation.y = 0;
  marker.pose.orientation.z = 0;

  // front left/right
  marker.pose.position.x = 0.19;
  marker.pose.position.y = 0.11;
  marker.id++;
  markers_.push_back(marker);

  marker.pose.position.x = 0.19;
  marker.pose.position.y = -0.11;
  marker.id++;
  markers_.push_back(marker);

  marker.color.g = 0.8;
  marker.color.b = 0.8;

  // left
  marker.pose.position.x = 0;
  marker.pose.position.y = 0.22;
  marker.id++;
  markers_.push_back(marker);

  // right
  marker.pose.position.x = 0;
  marker.pose.position.y = -0.22;
  marker.id++;
  markers_.push_back(marker);

  // back left/right
  marker.pose.position.x = -0.19;
  marker.pose.position.y = 0.11;
  marker.id++;
  markers_.push_back(marker);

  marker.pose.position.x = -0.19;
  marker.pose.position.y = -0.11;
  marker.id++;
  markers_.push_back(marker);

  if (simple) {
    // Make arms.
    marker.type = visualization_msgs::Marker::CUBE;
    marker.scale.x = 0.44;
    marker.scale.y = 0.02;
    marker.scale.z = 0.01;
    marker.color.r = 0.3;
    marker.color.g = 0.3;
    marker.color.b = 0.3;
    marker.color.a = 1;

    marker.pose.position.x = 0;
    marker.pose.position.y = 0;
    marker.pose.position.z = -0.015;
    marker.pose.orientation.x = 0;
    marker.pose.orientation.y = 0;

    marker.pose.orientation.w = sqrt2_2;
    marker.pose.orientation.z = sqrt2_2;
    marker.id++;
    markers_.push_back(marker);

    // 30 deg rotation  0.9659  0  0  0.2588

    marker.pose.orientation.w = 0.9659;
    marker.pose.orientation.z = 0.2588;
    marker.id++;
    markers_.push_back(marker);

    marker.pose.orientation.w = 0.9659;
    marker.pose.orientation.z = -0.2588;
    marker.id++;
    markers_.push_back(marker);
  } else {
    int id = marker.id;
    marker = visualization_msgs::Marker();
    marker.id = id + 1;
    marker.ns = "hexacopter";
    marker.mesh_resource =
        "package://mav_visualization/meshes/firefly_carbon.dae";
    marker.mesh_use_embedded_materials = false;
    marker.scale.x = 1;
    marker.scale.y = 1;
    marker.scale.z = 1;
    marker.color.a = 1.0;
    marker.color.r = 0.3;
    marker.color.g = 0.3;
    marker.color.b = 0.3;
    marker.pose.position.z = -0.03;
    marker.type = visualization_msgs::Marker::MESH_RESOURCE;
    markers_.push_back(marker);

    marker.id++;
    marker.mesh_resource =
        "package://mav_visualization/meshes/firefly_cowl.dae";
    marker.mesh_use_embedded_materials = false;
    marker.scale.x = 1;
    marker.scale.y = 1;
    marker.scale.z = 1;
    marker.color.a = 1.0;
    marker.color.r = 1.0;
    marker.color.g = 1.0;
    marker.color.b = 1.0;
    markers_.push_back(marker);
  }

  this->setFrameLocked(true);
}

}  // namespace mav_visualization
