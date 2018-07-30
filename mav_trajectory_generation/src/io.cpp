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

const std::string kSegments = "segments";
const std::string kNumCoefficients = "N";
const std::string kDim = "D";
const std::string kSegmentTime = "time";
const std::string kCoefficients = "coefficients";

namespace mav_trajectory_generation {
bool segmentsToFile(
    const std::string& filename,
    const mav_trajectory_generation::Segment::Vector& segments) {
  YAML::Emitter out;
  out << YAML::BeginMap;
  out << YAML::Key << kSegments;
  out << YAML::BeginSeq;
  for (const mav_trajectory_generation::Segment& segment : segments) {
    out << YAML::BeginMap;
    // Header.
    out << YAML::Key << kNumCoefficients << YAML::Value << segment.N();
    out << YAML::Key << kDim << YAML::Value << segment.D();
    out << YAML::Key << kSegmentTime << YAML::Value << segment.getTimeNSec()
        << YAML::Comment("[ns]");
    // Coefficients.
    out << YAML::Key << kCoefficients;
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

  if (node[kSegments]) {
    const YAML::Node& segments_yaml = node[kSegments];
    for (size_t i = 0; i < segments_yaml.size(); i++) {
      if (segments_yaml[i][kNumCoefficients] && segments_yaml[i][kDim] &&
          segments_yaml[i][kSegmentTime] && segments_yaml[i][kCoefficients]) {
        // Header.
        int N = segments_yaml[i][kNumCoefficients].as<int>();
        int D = segments_yaml[i][kDim].as<int>();
        mav_trajectory_generation::Segment segment(N, D);
        uint64_t t = segments_yaml[i][kSegmentTime].as<uint64_t>();
        segment.setTimeNSec(t);
        // Coefficients.
        if (segments_yaml[i][kCoefficients].size() != D) {
          return false;  // Coefficients and dimensions do not coincide.
        }
        for (size_t j = 0; j < D; j++) {
          if (segments_yaml[i][kCoefficients][j].size() != N) {
            return false;  // Number of coefficients does no coincide.
          }
          Eigen::VectorXd coeffs(N);
          for (size_t k = 0; k < N; k++) {
            coeffs(k) = segments_yaml[i][kCoefficients][j][k].as<double>();
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
  mav_msgs::EigenTrajectoryPoint::Vector flat_states;
  bool success = sampleWholeTrajectory(trajectory, sampling_time, &flat_states);
  if (!success) {
    return false;
  }

  // Layout: [t, x, y, z, vx, vy, vz, jx, jy, jz, sx, sy, sz, tm]
  const unsigned int dim = trajectory.D();
  Eigen::MatrixXd output(flat_states.size(), 5 * dim + 2);
  output.setZero();
  for (int i = 0; i < flat_states.size(); ++i) {
    const mav_msgs::EigenTrajectoryPoint state = flat_states[i];

    if (trajectory.D() > 3) {
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
    }
  }

  // Set the segment times
  mav_trajectory_generation::Segment::Vector segments;
  trajectory.getSegments(&segments);
  double current_segment_time = 0.0;
  for (int j = 0; j < segments.size(); ++j) {
    double segment_time = segments[j].getTime();
    current_segment_time += segment_time;
    output(j, 1 + 5 * dim) = current_segment_time;
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
