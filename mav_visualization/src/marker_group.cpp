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

#include "mav_visualization/marker_group.h"

namespace mav_visualization {

MarkerGroup::MarkerGroup() {}

MarkerGroup::~MarkerGroup() {}

void MarkerGroup::getMarkers(MarkerVector& markers, const double& scale,
                             bool append) const {
  if (!append && scale == 1) {
    markers = markers_;
    return;
  } else if (!append && scale != 1)
    markers.clear();

  markers.reserve(markers_.size());

  for (MarkerVector::const_iterator it = markers_.begin(); it < markers_.end();
       it++) {
    if (scale == 1) {
      markers.push_back(*it);
    } else {
      visualization_msgs::Marker m = *it;
      m.pose.position.x *= scale;
      m.pose.position.y *= scale;
      m.pose.position.z *= scale;
      m.scale.x *= scale;
      m.scale.y *= scale;
      m.scale.z *= scale;
      markers.push_back(m);
    }
  }
}

void MarkerGroup::getMarkers(visualization_msgs::MarkerArray& marker_array,
                             const double& scale, bool append) const {
  getMarkers(marker_array.markers, scale, append);
}

void MarkerGroup::setNamespace(const std::string& ns) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->ns = ns;
  }
}

void MarkerGroup::setHeader(const std_msgs::Header& header) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->header = header;
  }
}

void MarkerGroup::setHeaderAndNamespace(const std_msgs::Header& header,
                                        const std::string& ns) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->ns = ns;
    it->header = header;
  }
}

void MarkerGroup::setAction(const int32_t& action) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->action = action;
  }
}

void MarkerGroup::setLifetime(double lifetime) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->lifetime = ros::Duration(lifetime);
  }
}

void MarkerGroup::setFrameLocked(bool locked) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    it->frame_locked = locked;
  }
}

void MarkerGroup::publish(ros::Publisher& pub) {
  for (MarkerVector::iterator it = markers_.begin(); it < markers_.end();
       it++) {
    pub.publish(*it);
  };
}

void MarkerGroup::transformMarker(visualization_msgs::Marker& marker,
                                  const Eigen::Vector3d& t,
                                  const Eigen::Quaterniond& q) {
  geometry_msgs::Pose::_orientation_type& mq = marker.pose.orientation;
  geometry_msgs::Pose::_position_type& mp = marker.pose.position;

  Eigen::Quaterniond e_mq(mq.w, mq.x, mq.y, mq.z);
  Eigen::Vector3d e_mp(mp.x, mp.y, mp.z);

  e_mp = q * e_mp + t;
  e_mq = q * e_mq;

  mq.w = e_mq.w();
  mq.x = e_mq.x();
  mq.y = e_mq.y();
  mq.z = e_mq.z();

  mp.x = e_mp[0];
  mp.y = e_mp[1];
  mp.z = e_mp[2];
}

void MarkerGroup::transform(const Eigen::Vector3d& t,
                            const Eigen::Quaterniond& q) {
  for (MarkerVector::iterator it = markers_.begin(); it != markers_.end(); ++it)
    transformMarker(*it, t, q);
}

}  // namespace mav_visualization
