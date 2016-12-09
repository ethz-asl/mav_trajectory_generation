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

#ifndef TRAJECTORY_TYPES_TEMPLATELESS_H_
#define TRAJECTORY_TYPES_TEMPLATELESS_H_

#include <chrono>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/StdVector>
#include <glog/logging.h>

#include <mav_planning_utils/motion_defines.h>
#include <mav_planning_utils/polynomial_templateless.h>

namespace mav_planning_utils {
/***
 * \brief Container holding the properties (time, value, segment where it
 * occurred) of an extremum.
 */
class Extremum {
 public:
  Extremum() : time(0), value(0), segment_idx(0) {}

  Extremum(double _time, double _value, int _segment_idx)
      : time(_time), value(_value), segment_idx(_segment_idx) {}

  bool operator<(const Extremum& rhs) const { return value < rhs.value; }

  double
      time;  ///< Time where the extremum occurs, relative to the segment start.
  double value;     ///< Value of the extremum at time.
  int segment_idx;  ///< Index of the segment, where the extremum occurs.
};

std::ostream& operator<<(std::ostream& stream, const Extremum& e);

}  // end namespace mav_planning_utils

#endif

