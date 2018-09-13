/*
 * Copyright (c) 2018, Helen Oleynikova, ASL, ETH Zurich, Switzerland
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

#ifndef MAV_TRAJECTORY_GENERATION_RPOLY_RPOLY_AK1_H_
#define MAV_TRAJECTORY_GENERATION_RPOLY_RPOLY_AK1_H_

#include <Eigen/Eigen>

namespace mav_trajectory_generation {

int findLastNonZeroCoeff(const Eigen::VectorXd& coefficients);

bool findRootsJenkinsTraub(const Eigen::VectorXd& coefficients_increasing,
                           Eigen::VectorXcd* roots);

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_RPOLY_RPOLY_AK1_H_
