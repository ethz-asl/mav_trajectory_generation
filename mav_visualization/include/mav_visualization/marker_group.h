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

#ifndef MAV_VISUALIZATION_MARKER_GROUP_H_
#define MAV_VISUALIZATION_MARKER_GROUP_H_

#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <string>
#include <vector>

namespace mav_visualization {

typedef std::vector<visualization_msgs::Marker> MarkerVector;

class MarkerGroup {
 public:
  MarkerGroup();
  virtual ~MarkerGroup();
  void getMarkers(visualization_msgs::MarkerArray& marker_array,
                  const double& scale = 1, bool append = false) const;
  void getMarkers(MarkerVector& markers, const double& scale = 1,
                  bool append = false) const;
  void setNamespace(const std::string& ns);
  void setHeader(const std_msgs::Header& header);
  void setHeaderAndNamespace(const std_msgs::Header& header,
                             const std::string& ns);
  void setAction(const int32_t& action);
  void setLifetime(double lifetime);
  void setFrameLocked(bool locked);
  void transform(const Eigen::Vector3d& t, const Eigen::Quaterniond& q);
  void publish(ros::Publisher& pub);

 protected:
  std::string name_;
  std::string description_;
  MarkerVector markers_;
  static void transformMarker(visualization_msgs::Marker& marker,
                              const Eigen::Vector3d& t,
                              const Eigen::Quaterniond& q);
};

}  // namespace mav_visualization

#endif  // MAV_VISUALIZATION_MARKER_GROUP_H_
