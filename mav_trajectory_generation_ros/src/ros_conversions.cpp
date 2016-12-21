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

#include "mav_trajectory_generation_ros/ros_conversions.h"

namespace mav_trajectory_generation {

bool trajectoryToPolynomialTrajectoryMsg(
    const Trajectory& trajectory, planning_msgs::PolynomialTrajectory4D* msg) {
  CHECK_NOTNULL(msg);
  msg->segments.clear();

  bool success = true;

  Segment::Vector segments;
  trajectory.getSegments(&segments);

  msg->segments.reserve(segments.size());
  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment& segment = segments[i];

    if (segment.D() < 3) {
      LOG(ERROR) << "Dimension of position segment has to be 3 or 4, but is "
                 << segment.D();
      success = false;
      break;
    }

    planning_msgs::PolynomialSegment4D segment_msg;
    planning_msgs::EigenPolynomialSegment eigen_segment;
    eigen_segment.x = segment[0].getCoefficients();
    eigen_segment.y = segment[1].getCoefficients();
    eigen_segment.z = segment[2].getCoefficients();
    if (segment.D() > 3) {
      eigen_segment.yaw = segment[3].getCoefficients();
    }
    eigen_segment.num_coeffs = segment.N();
    eigen_segment.segment_time_ns = segment.getTimeNSec();

    planning_msgs::polynomialSegmentMsgFromEigen(eigen_segment, &segment_msg);
    msg->segments.push_back(segment_msg);
  }

  if (!success) msg->segments.clear();
  return success;
}

// Converts a ROS polynomial trajectory msg into a Trajectory.
bool polynomialTrajectoryMsgToTrajectory(
    const planning_msgs::PolynomialTrajectory4D& msg, Trajectory* trajectory) {
  planning_msgs::EigenPolynomialTrajectory eigen_trajectory_msg;
  planning_msgs::eigenPolynomialTrajectoryFromMsg(msg, &eigen_trajectory_msg);
  // TODO(helenol): maybe add more error checking here.
  Segment::Vector segment_vector;
  for (const planning_msgs::EigenPolynomialSegment& msg_segment :
       eigen_trajectory_msg) {
    int D = 3;
    if (msg_segment.yaw.size() > 0) {
      D = 4;
    }

    Segment segment(msg_segment.x.size(), D);
    segment[0].setCoefficients(msg_segment.x);
    segment[1].setCoefficients(msg_segment.y);
    segment[2].setCoefficients(msg_segment.z);
    if (D > 3) {
      segment[3].setCoefficients(msg_segment.yaw);
    }
    segment.setTimeNSec(msg_segment.segment_time_ns);
    segment_vector.push_back(segment);
  }

  trajectory->setSegments(segment_vector);
  return true;
}

}  // namespace mav_trajectory_generation
