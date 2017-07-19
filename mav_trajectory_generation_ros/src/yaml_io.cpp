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

#include "mav_trajectory_generation_ros/yaml_io.h"

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
    for (size_t i = 0; i < segment.D(); i++) {
      out << YAML::Flow;  // List output format.
      out << YAML::BeginSeq;
      Eigen::VectorXd coefficients = segment[i].getCoefficients();
      for (size_t j = 0; j < segment.N(); j++) {
        out << coefficients(j);
      }
      out << YAML::EndSeq;
    }
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

//mav_trajectory_generation::Segment segment()
const YAML::Node& segments_yaml = node[kSegments];


}

}  // namespace mav_trajectory_generation
