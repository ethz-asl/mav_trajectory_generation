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

#include "mav_trajectory_generation/io.h"
#include "mav_trajectory_generation/trajectory_sampling.h"

#include <yaml-cpp/yaml.h>
#include <fstream>

const std::string kSegmentsKey = "segments";
const std::string kNumCoefficientsKey = "N";
const std::string kDimKey = "D";
const std::string kSegmentTimeKey = "time";
const std::string kCoefficientsKey = "coefficients";

namespace mav_trajectory_generation {

YAML::Node coefficientsToYaml(const Eigen::VectorXd& coefficients) {
  YAML::Node node(YAML::NodeType::Sequence);
  for (size_t i = 0; i < coefficients.size(); ++i)
    node.push_back(coefficients(i));
  node.SetStyle(YAML::EmitterStyle::Flow);
  return node;
}

YAML::Node segmentToYaml(const Segment& segment) {
  YAML::Node node;
  node[kNumCoefficientsKey] = segment.N();
  node[kDimKey] = segment.D();
  node[kSegmentTimeKey] = segment.getTimeNSec();

  for (size_t i = 0; i < segment.D(); ++i)
    node[kCoefficientsKey].push_back(
        coefficientsToYaml(segment[i].getCoefficients()));

  return node;
}

YAML::Node segmentsToYaml(const Segment::Vector& segments) {
  YAML::Node node;
  for (const mav_trajectory_generation::Segment& segment : segments)
    node[kSegmentsKey].push_back(segmentToYaml(segment));

  return node;
}

YAML::Node trajectoryToYaml(const Trajectory& trajectory) {
  Segment::Vector segments;
  trajectory.getSegments(&segments);
  return segmentsToYaml(segments);
}

bool coefficientsFromYaml(const YAML::Node& node,
                          Eigen::VectorXd* coefficients) {
  CHECK_NOTNULL(coefficients);
  if (!node.IsSequence()) return false;
  *coefficients = Eigen::VectorXd(node.size());
  for (std::size_t i = 0; i < node.size(); ++i) {
    (*coefficients)(i) = node[i].as<double>();
  }
  return true;
}

bool segmentFromYaml(const YAML::Node& node, Segment* segment) {
  CHECK_NOTNULL(segment);

  if (!node[kNumCoefficientsKey]) return false;
  if (!node[kDimKey]) return false;
  if (!node[kSegmentTimeKey]) return false;
  if (!node[kCoefficientsKey]) return false;
  if (!node[kCoefficientsKey].IsSequence()) return false;

  *segment =
      Segment(node[kNumCoefficientsKey].as<int>(), node[kDimKey].as<int>());

  for (size_t i = 0; i < segment->D(); ++i) {
    Eigen::VectorXd coeffs;
    if (!coefficientsFromYaml(node[kCoefficientsKey][i], &coeffs)) return false;
    (*segment)[i] = coeffs;
  }

  segment->setTimeNSec(node[kSegmentTimeKey].as<uint64_t>());

  return true;
}

bool segmentsFromYaml(const YAML::Node& node, Segment::Vector* segments) {
  CHECK_NOTNULL(segments);
  if (!node.IsSequence()) return false;

  segments->resize(node.size(), Segment(0, 0));
  for (size_t i = 0; i < node.size(); ++i) {
    if (!segmentFromYaml(node[i], &(*segments)[i])) return false;
  }

  return true;
}

bool trajectoryFromYaml(const YAML::Node& node, Trajectory* trajectory) {
  CHECK_NOTNULL(trajectory);

  Segment::Vector segments;
  if (!segmentsFromYaml(node[kSegmentsKey], &segments)) return false;
  trajectory->setSegments(segments);

  return true;
}

bool segmentsToFile(
    const std::string& filename,
    const mav_trajectory_generation::Segment::Vector& segments) {
  YAML::Emitter out;
  out << YAML::BeginMap;
  out << YAML::Key << kSegmentsKey;
  out << YAML::BeginSeq;
  for (const mav_trajectory_generation::Segment& segment : segments) {
    out << YAML::BeginMap;
    // Header.
    out << YAML::Key << kNumCoefficientsKey << YAML::Value << segment.N();
    out << YAML::Key << kDimKey << YAML::Value << segment.D();
    out << YAML::Key << kSegmentTimeKey << YAML::Value << segment.getTimeNSec()
        << YAML::Comment("[ns]");
    // Coefficients.
    out << YAML::Key << kCoefficientsKey;
    out << YAML::BeginSeq;
    for (size_t i = 0; i < segment.D(); i++) {
      out << YAML::Flow;  // List output format.
      out << YAML::BeginSeq;
      Eigen::VectorXd coefficients = segment[i].getCoefficients();
      for (size_t j = 0; j < segment.N(); j++) {
        out << coefficients(j);
      }
      out << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
  }
  out << YAML::EndSeq;
  out << YAML::EndMap;

  // Write to file.
  std::ofstream fout(filename);
  if (!fout) {
    return false;
  }
  fout << out.c_str();
  fout.close();

  return true;
}

bool segmentsFromFile(const std::string& filename,
                      mav_trajectory_generation::Segment::Vector* segments) {
  CHECK_NOTNULL(segments);

  // Check file exists and is readable.
  std::ifstream in(filename);
  if (!in.good()) {
    return false;
  }
  segments->clear();

  // Parse YAML.
  YAML::Node node = YAML::LoadFile(filename);

  if (node[kSegmentsKey]) {
    const YAML::Node& segments_yaml = node[kSegmentsKey];
    for (size_t i = 0; i < segments_yaml.size(); i++) {
      if (segments_yaml[i][kNumCoefficientsKey] && segments_yaml[i][kDimKey] &&
          segments_yaml[i][kSegmentTimeKey] &&
          segments_yaml[i][kCoefficientsKey]) {
        // Header.
        int N = segments_yaml[i][kNumCoefficientsKey].as<int>();
        int D = segments_yaml[i][kDimKey].as<int>();
        mav_trajectory_generation::Segment segment(N, D);
        uint64_t t = segments_yaml[i][kSegmentTimeKey].as<uint64_t>();
        segment.setTimeNSec(t);
        // Coefficients.
        if (segments_yaml[i][kCoefficientsKey].size() != D) {
          return false;  // Coefficients and dimensions do not coincide.
        }
        for (size_t j = 0; j < D; j++) {
          if (segments_yaml[i][kCoefficientsKey][j].size() != N) {
            return false;  // Number of coefficients does no coincide.
          }
          Eigen::VectorXd coeffs(N);
          for (size_t k = 0; k < N; k++) {
            coeffs(k) = segments_yaml[i][kCoefficientsKey][j][k].as<double>();
          }
          segment[j] = coeffs;
        }
        segments->push_back(segment);
      } else {
        return false;  // Wrong format, missing elements.
      }
    }
  } else {
    return false;  // No segments element.
  }

  return true;
}

bool sampledTrajectoryStatesToFile(const std::string& filename,
                                   const Trajectory& trajectory) {
  // Print to file for matlab
  const double sampling_time = 0.01;
  mav_msgs::EigenTrajectoryPoint::Vector trajectory_points;
  bool success =
      sampleWholeTrajectory(trajectory, sampling_time, &trajectory_points);
  if (!success) {
    return false;
  }

  // Layout: [t, x, y, z, vx, vy, vz, jx, jy, jz, sx, sy, sz, qw, qx, qy, qz,
  // wx, wy, wz, ax, ay, az, tm]
  const unsigned int dim = 3;
  Eigen::MatrixXd output(trajectory_points.size(), 8 * dim + 3);
  output.setZero();
  for (int i = 0; i < trajectory_points.size(); ++i) {
    const mav_msgs::EigenTrajectoryPoint state = trajectory_points[i];

    if (trajectory.D() == 4) {
      double yaw = state.getYaw();
      double yaw_rate = state.getYawRate();
      double yaw_acc = state.getYawAcc();
    }

    if (i < output.rows()) {
      output(i, 0) = state.time_from_start_ns;
      output.row(i).segment(1, dim) = state.position_W;
      output.row(i).segment(1 + dim, dim) = state.velocity_W;
      output.row(i).segment(1 + 2 * dim, dim) = state.acceleration_W;
      output.row(i).segment(1 + 3 * dim, dim) = state.jerk_W;
      output.row(i).segment(1 + 4 * dim, dim) = state.snap_W;
      output.row(i).segment(1 + 5 * dim, dim + 1) =
          Eigen::Vector4d(state.orientation_W_B.w(), state.orientation_W_B.x(),
                          state.orientation_W_B.y(), state.orientation_W_B.z());
      output.row(i).segment(2 + 6 * dim, dim) = state.angular_velocity_W;
      output.row(i).segment(2 + 7 * dim, dim) = state.angular_acceleration_W;
    }
  }

  // Set the segment times
  mav_trajectory_generation::Segment::Vector segments;
  trajectory.getSegments(&segments);
  double current_segment_time = 0.0;
  for (int j = 0; j < segments.size(); ++j) {
    double segment_time = segments[j].getTime();
    current_segment_time += segment_time;
    output(j, 2 + 8 * dim) = current_segment_time;
  }

  std::fstream fs;
  fs.open(filename, std::fstream::out);
  if (!fs) {
    return false;
  }
  fs << output;
  fs.close();
  return true;
}

}  // namespace mav_trajectory_generation
