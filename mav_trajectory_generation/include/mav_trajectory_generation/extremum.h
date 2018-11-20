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

#ifndef MAV_TRAJECTORY_GENERATION_EXTREMUM_H_
#define MAV_TRAJECTORY_GENERATION_EXTREMUM_H_

#include <iostream>

namespace mav_trajectory_generation {

// Container holding the properties of an extremum (time, value,
// segment where it occurred).
struct Extremum {
 public:
  Extremum() : time(0.0), value(0.0), segment_idx(0) {}

  Extremum(double _time, double _value, int _segment_idx)
      : time(_time), value(_value), segment_idx(_segment_idx) {}

  bool operator<(const Extremum& rhs) const { return value < rhs.value; }
  bool operator>(const Extremum& rhs) const { return value > rhs.value; }

  double
      time;  // Time where the extremum occurs, relative to the segment start.
  double value;     // Value of the extremum at time.
  int segment_idx;  // Index of the segment where the extremum occurs.
};

inline std::ostream& operator<<(std::ostream& stream, const Extremum& e) {
  stream << "time: " << e.time << ", value: " << e.value
         << ", segment idx: " << e.segment_idx << std::endl;
  return stream;
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_EXTREMUM_H_
