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

#include "mav_visualization/leica_marker.h"

int main(int argc, char** argv) {
  ros::init(argc, argv, "leica_publisher");
  ros::NodeHandle nh;
  ros::NodeHandle nh_private("~");

  ros::Publisher marker_pub =
      nh_private.advertise<visualization_msgs::MarkerArray>("marker_array", 10,
                                                            true);

  std::string frame_id("leica");
  double scale = 1.0;

  nh_private.param("frame_id", frame_id, frame_id);
  nh_private.param("scale", scale, scale);

  mav_visualization::LeicaMarker leica;
  visualization_msgs::MarkerArray markers;

  leica.setLifetime(0.0);
  leica.setAction(visualization_msgs::Marker::ADD);

  std_msgs::Header header;
  header.frame_id = frame_id;

  while (ros::ok()) {
    header.stamp = ros::Time::now();
    leica.setHeader(header);
    leica.getMarkers(markers, scale, false);
    marker_pub.publish(markers);
    ++header.seq;

    ros::Duration(5.0).sleep();
  }
}
