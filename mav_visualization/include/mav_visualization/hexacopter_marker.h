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

#ifndef MAV_VISUALIZATION_HEXACOPTER_MARKER_H_
#define MAV_VISUALIZATION_HEXACOPTER_MARKER_H_

#include "mav_visualization/marker_group.h"

namespace mav_visualization {

class HexacopterMarker : public MarkerGroup {
 public:
  HexacopterMarker(bool simple = false);

 private:
  void createHexacopter(bool simple = false);
};

}  // namespace mav_visualization

#endif  // MAV_VISUALIZATION_HEXACOPTER_MARKER_H_
