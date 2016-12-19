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

#ifndef MAV_TRAJECTORY_GENERATION_CONVOLUTION_H_
#define MAV_TRAJECTORY_GENERATION_CONVOLUTION_H_

namespace mav_trajectory_generation {

template <int D, int K>
struct ConvolutionDimension {
  enum { length = D + K - 1 };
};

template <int DataDimension_, int KernelDimension_>
Eigen::Matrix<double,
              ConvolutionDimension<DataDimension_, KernelDimension_>::length, 1>
convolve(const Eigen::Matrix<double, DataDimension_, 1>& data,
         const Eigen::Matrix<double, KernelDimension_, 1>& kernel) {
  const int convolution_dimension =
      ConvolutionDimension<DataDimension_, KernelDimension_>::length;
  Eigen::Matrix<double, convolution_dimension, 1> convolved;
  convolved.setZero();
  Eigen::Matrix<double, KernelDimension_, 1> kernel_reverse(kernel.reverse());

  for (int output_idx = 0; output_idx < convolution_dimension; ++output_idx) {
    const int data_idx = output_idx - KernelDimension_ + 1;

    int lower_bound = std::max(0, -data_idx);
    int upper_bound = std::min(KernelDimension_, DataDimension_ - data_idx);

    for (int kernel_idx = lower_bound; kernel_idx < upper_bound; ++kernel_idx) {
      convolved[output_idx] +=
          kernel_reverse[kernel_idx] * data[data_idx + kernel_idx];
    }
  }

  return convolved;
}

template <typename T>
int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_CONVOLUTION_H_
