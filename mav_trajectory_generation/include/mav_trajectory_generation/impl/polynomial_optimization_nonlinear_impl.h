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

#ifndef MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_NONLINEAR_IMPL_H_
#define MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_NONLINEAR_IMPL_H_

#include <chrono>

#include "mav_trajectory_generation/polynomial_optimization_linear.h"
#include "mav_trajectory_generation/timing.h"

namespace mav_trajectory_generation {

inline void OptimizationInfo::print(std::ostream& stream) const {
  stream << "--- optimization info ---" << std::endl;
  stream << "  optimization time:     " << optimization_time << std::endl;
  stream << "  n_iterations:          " << n_iterations << std::endl;
  stream << "  stopping reason:       "
         << nlopt::returnValueToString(stopping_reason) << std::endl;
  stream << "  cost trajectory:       " << cost_trajectory << std::endl;
  stream << "  cost time:             " << cost_time << std::endl;
  stream << "  cost soft constraints: " << cost_soft_constraints << std::endl;
  stream << "  maxima: " << std::endl;
  for (const std::pair<int, Extremum>& m : maxima) {
    stream << "    " << positionDerivativeToString(m.first) << ": "
           << m.second.value << " in segment " << m.second.segment_idx
           << " and segment time " << m.second.time << std::endl;
  }
}

template <int _N>
PolynomialOptimizationNonLinear<_N>::PolynomialOptimizationNonLinear(
    size_t dimension, const NonlinearOptimizationParameters& parameters,
    bool optimize_time_only)
    : poly_opt_(dimension),
      optimization_parameters_(parameters),
      optimize_time_only_(optimize_time_only) {}

template <int _N>
bool PolynomialOptimizationNonLinear<_N>::setupFromVertices(
    const Vertex::Vector& vertices, const std::vector<double>& segment_times,
    int derivative_to_optimize) {
  bool ret = poly_opt_.setupFromVertices(vertices, segment_times,
                                         derivative_to_optimize);

  size_t n_optimization_parameters;
  if (optimize_time_only_) {
    n_optimization_parameters = segment_times.size();
  } else {
    n_optimization_parameters =
        segment_times.size() +
        poly_opt_.getNumberFreeConstraints() * poly_opt_.getDimension();
  }

  nlopt_.reset(new nlopt::opt(optimization_parameters_.algorithm,
                              n_optimization_parameters));
  nlopt_->set_ftol_rel(optimization_parameters_.f_rel);
  nlopt_->set_ftol_abs(optimization_parameters_.f_abs);
  nlopt_->set_xtol_rel(optimization_parameters_.x_rel);
  nlopt_->set_xtol_abs(optimization_parameters_.x_abs);
  nlopt_->set_maxeval(optimization_parameters_.max_iterations);

  if (optimization_parameters_.random_seed < 0)
    nlopt_srand_time();
  else
    nlopt_srand(optimization_parameters_.random_seed);

  return ret;
}

template <int _N>
bool PolynomialOptimizationNonLinear<_N>::solveLinear() {
  return poly_opt_.solveLinear();
}

template <int _N>
int PolynomialOptimizationNonLinear<_N>::optimize() {
  optimization_info_ = OptimizationInfo();
  int result = nlopt::FAILURE;

  const std::chrono::high_resolution_clock::time_point t_start =
      std::chrono::high_resolution_clock::now();

  if (optimize_time_only_) {
    result = optimizeTime();
  } else {
    result = optimizeTimeAndFreeConstraints();
  }

  const std::chrono::high_resolution_clock::time_point t_stop =
      std::chrono::high_resolution_clock::now();
  optimization_info_.optimization_time =
      std::chrono::duration_cast<std::chrono::duration<double> >(t_stop -
                                                                 t_start)
          .count();

  optimization_info_.stopping_reason = result;

  return result;
}

template <int _N>
int PolynomialOptimizationNonLinear<_N>::optimizeTime() {
  std::vector<double> initial_step, segment_times, upper_bounds;

  poly_opt_.getSegmentTimes(&segment_times);
  const size_t n_segments = segment_times.size();

  initial_step.reserve(n_segments);
  for (double t : segment_times) {
    initial_step.push_back(optimization_parameters_.initial_stepsize_rel * t);
  }

  for (double t : segment_times) {
    upper_bounds.push_back(t * 2.0);
  }

  try {
    // Set a lower bound on the segment time per segment to avoid numerical
    // issues.
    constexpr double kOptimizationTimeLowerBound = 0.1;
    nlopt_->set_initial_step(initial_step);
    nlopt_->set_upper_bounds(upper_bounds);
    nlopt_->set_lower_bounds(kOptimizationTimeLowerBound);
    nlopt_->set_min_objective(
        &PolynomialOptimizationNonLinear<N>::objectiveFunctionTime, this);
  } catch (std::exception& e) {
    LOG(ERROR) << "error while setting up nlopt: " << e.what() << std::endl;
    return nlopt::FAILURE;
  }

  double final_cost = std::numeric_limits<double>::max();
  int result;

  try {
    result = nlopt_->optimize(segment_times, final_cost);
  } catch (std::exception& e) {
    LOG(ERROR) << "error while running nlopt: " << e.what() << std::endl;
    return nlopt::FAILURE;
  }

  return result;
}

template <int _N>
int PolynomialOptimizationNonLinear<_N>::optimizeTimeAndFreeConstraints() {
  std::vector<double> initial_step, initial_solution, segment_times,
      lower_bounds, upper_bounds;

  poly_opt_.getSegmentTimes(&segment_times);
  const size_t n_segments = segment_times.size();

  // compute initial solution
  poly_opt_.solveLinear();
  std::vector<Eigen::VectorXd> free_constraints;
  poly_opt_.getFreeConstraints(&free_constraints);
  CHECK(free_constraints.size() > 0);
  CHECK(free_constraints.front().size() > 0);

  const size_t n_optmization_variables =
      n_segments + free_constraints.size() * free_constraints.front().size();

  initial_solution.reserve(n_optmization_variables);
  initial_step.reserve(n_optmization_variables);
  lower_bounds.reserve(n_optmization_variables);
  upper_bounds.reserve(n_optmization_variables);

  // copy all constraints into one vector:
  for (double t : segment_times) {
    initial_solution.push_back(t);
  }

  for (const Eigen::VectorXd& c : free_constraints) {
    for (int i = 0; i < c.size(); ++i) {
      initial_solution.push_back(c[i]);
    }
  }

  initial_step.reserve(n_optmization_variables);
  for (double x : initial_solution) {
    const double abs_x = std::abs(x);
    initial_step.push_back(optimization_parameters_.initial_stepsize_rel *
                           abs_x);
    lower_bounds.push_back(-abs_x * 2);
    upper_bounds.push_back(abs_x * 2);
  }

  for (size_t i = 0; i < n_segments; ++i) {
    lower_bounds[i] = 0.1;
  }

  try {
    nlopt_->set_initial_step(initial_step);
    nlopt_->set_lower_bounds(lower_bounds);
    nlopt_->set_upper_bounds(upper_bounds);
    nlopt_->set_min_objective(&PolynomialOptimizationNonLinear<
                                  N>::objectiveFunctionTimeAndConstraints,
                              this);
  } catch (std::exception& e) {
    LOG(ERROR) << "error while setting up nlopt: " << e.what() << std::endl;
    return nlopt::FAILURE;
  }

  double final_cost = std::numeric_limits<double>::max();
  int result;

  try {
    timing::Timer timer_solve("optimize_nonlinear_full_total_time");
    result = nlopt_->optimize(initial_solution, final_cost);
    timer_solve.Stop();
  } catch (std::exception& e) {
    LOG(ERROR) << "error while running nlopt: " << e.what() << std::endl;
    return nlopt::FAILURE;
  }

  return result;
}

template <int _N>
bool PolynomialOptimizationNonLinear<_N>::addMaximumMagnitudeConstraint(
    int derivative, double maximum_value) {
  CHECK_GE(derivative, 0);
  CHECK_GE(maximum_value, 0.0);

  std::shared_ptr<ConstraintData> constraint_data(new ConstraintData);
  constraint_data->derivative = derivative;
  constraint_data->value = maximum_value;
  constraint_data->this_object = this;

  // Store the shared_ptrs such that their data will be destroyed later.
  inequality_constraints_.push_back(constraint_data);

  if (!optimization_parameters_.use_soft_constraints) {
    try {
      nlopt_->add_inequality_constraint(
          &PolynomialOptimizationNonLinear<
              N>::evaluateMaximumMagnitudeConstraint,
          constraint_data.get(),
          optimization_parameters_.inequality_constraint_tolerance);
    } catch (std::exception& e) {
      LOG(ERROR) << "ERROR while setting inequality constraint " << e.what()
                 << std::endl;
      return false;
    }
  }
  return true;
}

template <int _N>
double PolynomialOptimizationNonLinear<_N>::objectiveFunctionTime(
    const std::vector<double>& segment_times, std::vector<double>& gradient,
    void* data) {
  CHECK(gradient.empty())
      << "computing gradient not possible, choose a gradient free method";
  CHECK_NOTNULL(data);

  PolynomialOptimizationNonLinear<N>* optimization_data =
      static_cast<PolynomialOptimizationNonLinear<N>*>(data);  // wheee ...

  CHECK_EQ(segment_times.size(),
           optimization_data->poly_opt_.getNumberSegments());

  optimization_data->poly_opt_.updateSegmentTimes(segment_times);
  optimization_data->poly_opt_.solveLinear();
  double cost_trajectory = optimization_data->poly_opt_.computeCost();
  double cost_time = 0;
  double cost_constraints = 0;
  const double total_time = computeTotalTrajectoryTime(segment_times);
  cost_time = total_time * total_time *
              optimization_data->optimization_parameters_.time_penalty;

  if (optimization_data->optimization_parameters_.print_debug_info) {
    std::cout << "---- cost at iteration "
              << optimization_data->optimization_info_.n_iterations << "---- "
              << std::endl;
    std::cout << "  trajectory: " << cost_trajectory << std::endl;
    std::cout << "  time: " << cost_time << std::endl;
  }

  if (optimization_data->optimization_parameters_.use_soft_constraints) {
    cost_constraints =
        optimization_data->evaluateMaximumMagnitudeAsSoftConstraint(
            optimization_data->inequality_constraints_,
            optimization_data->optimization_parameters_.soft_constraint_weight);
  }

  if (optimization_data->optimization_parameters_.print_debug_info) {
    std::cout << "  sum: " << cost_trajectory + cost_time + cost_constraints
              << std::endl;
    std::cout << "  total time: " << total_time << std::endl;
  }

  optimization_data->optimization_info_.n_iterations++;
  optimization_data->optimization_info_.cost_trajectory = cost_trajectory;
  optimization_data->optimization_info_.cost_time = cost_time;
  optimization_data->optimization_info_.cost_soft_constraints =
      cost_constraints;

  return cost_trajectory + cost_time + cost_constraints;
}

template <int _N>
double PolynomialOptimizationNonLinear<_N>::objectiveFunctionTimeAndConstraints(
    const std::vector<double>& x, std::vector<double>& gradient, void* data) {
  CHECK(gradient.empty())
      << "computing gradient not possible, choose a gradient free method";
  CHECK_NOTNULL(data);

  PolynomialOptimizationNonLinear<N>* optimization_data =
      static_cast<PolynomialOptimizationNonLinear<N>*>(data);  // wheee ...

  const size_t n_segments = optimization_data->poly_opt_.getNumberSegments();
  const size_t n_free_constraints =
      optimization_data->poly_opt_.getNumberFreeConstraints();
  const size_t dim = optimization_data->poly_opt_.getDimension();

  CHECK_EQ(x.size(), n_segments + n_free_constraints * dim);

  std::vector<Eigen::VectorXd> free_constraints;
  free_constraints.resize(dim);
  std::vector<double> segment_times;
  segment_times.reserve(n_segments);

  for (size_t i = 0; i < n_segments; ++i) segment_times.push_back(x[i]);

  for (size_t d = 0; d < dim; ++d) {
    const size_t idx_start = n_segments + d * n_free_constraints;

    Eigen::VectorXd& free_constraints_dim = free_constraints[d];
    free_constraints_dim.resize(n_free_constraints, Eigen::NoChange);
    for (size_t i = 0; i < n_free_constraints; ++i) {
      free_constraints_dim[i] = x[idx_start + i];
    }
  }

  optimization_data->poly_opt_.updateSegmentTimes(segment_times);
  optimization_data->poly_opt_.setFreeConstraints(free_constraints);

  double cost_trajectory = optimization_data->poly_opt_.computeCost();
  double cost_time = 0;
  double cost_constraints = 0;

  const double total_time = computeTotalTrajectoryTime(segment_times);
  cost_time = total_time * total_time *
              optimization_data->optimization_parameters_.time_penalty;

  if (optimization_data->optimization_parameters_.print_debug_info) {
    std::cout << "---- cost at iteration "
              << optimization_data->optimization_info_.n_iterations << "---- "
              << std::endl;
    std::cout << "  trajectory: " << cost_trajectory << std::endl;
    std::cout << "  time: " << cost_time << std::endl;
  }

  if (optimization_data->optimization_parameters_.use_soft_constraints) {
    cost_constraints =
        optimization_data->evaluateMaximumMagnitudeAsSoftConstraint(
            optimization_data->inequality_constraints_,
            optimization_data->optimization_parameters_.soft_constraint_weight);
  }

  if (optimization_data->optimization_parameters_.print_debug_info) {
    std::cout << "  sum: " << cost_trajectory + cost_time + cost_constraints
              << std::endl;
    std::cout << "  total time: " << total_time << std::endl;
  }

  optimization_data->optimization_info_.n_iterations++;
  optimization_data->optimization_info_.cost_trajectory = cost_trajectory;
  optimization_data->optimization_info_.cost_time = cost_time;
  optimization_data->optimization_info_.cost_soft_constraints =
      cost_constraints;

  return cost_trajectory + cost_time + cost_constraints;
}

template <int _N>
double PolynomialOptimizationNonLinear<_N>::evaluateMaximumMagnitudeConstraint(
    const std::vector<double>& segment_times, std::vector<double>& gradient,
    void* data) {
  CHECK(gradient.empty())
      << "computing gradient not possible, choose a gradient free method";
  ConstraintData* constraint_data =
      static_cast<ConstraintData*>(data);  // wheee ...
  PolynomialOptimizationNonLinear<N>* optimization_data =
      constraint_data->this_object;

  Extremum max;
  // for now, let's assume that the optimization has been done
  switch (constraint_data->derivative) {
    case derivative_order::POSITION:
      max = optimization_data->poly_opt_
                .template computeMaximumOfMagnitude<derivative_order::POSITION>(
                    nullptr);
      break;
    case derivative_order::VELOCITY:
      max = optimization_data->poly_opt_
                .template computeMaximumOfMagnitude<derivative_order::VELOCITY>(
                    nullptr);
      break;
    case derivative_order::ACCELERATION:
      max = optimization_data->poly_opt_.template computeMaximumOfMagnitude<
          derivative_order::ACCELERATION>(nullptr);
      break;
    case derivative_order::JERK:
      max = optimization_data->poly_opt_
                .template computeMaximumOfMagnitude<derivative_order::JERK>(
                    nullptr);
      break;
    case derivative_order::SNAP:
      max = optimization_data->poly_opt_
                .template computeMaximumOfMagnitude<derivative_order::SNAP>(
                    nullptr);
      break;
    default:
      LOG(WARNING) << "[Nonlinear inequality constraint evaluation]: no "
                      "implementation for derivative: "
                   << constraint_data->derivative;
      return 0;
  }

  optimization_data->optimization_info_.maxima[constraint_data->derivative] =
      max;

  return max.value - constraint_data->value;
}

template <int _N>
double
PolynomialOptimizationNonLinear<_N>::evaluateMaximumMagnitudeAsSoftConstraint(
    const std::vector<std::shared_ptr<ConstraintData> >& inequality_constraints,
    double weight, double maximum_cost) const {
  std::vector<double> dummy;
  double cost = 0;

  if (optimization_parameters_.print_debug_info)
    std::cout << "  soft_constraints: " << std::endl;

  for (std::shared_ptr<const ConstraintData> constraint :
       inequality_constraints_) {
    // need to call the c-style callback function here, thus the ugly cast to
    // void*.
    double abs_violation = evaluateMaximumMagnitudeConstraint(
        dummy, dummy, (void*)constraint.get());

    double relative_violation = abs_violation / constraint->value;
    const double current_cost =
        std::min(maximum_cost, exp(relative_violation * weight));
    cost += current_cost;
    if (optimization_parameters_.print_debug_info) {
      std::cout << "    derivative " << constraint->derivative
                << " abs violation: " << abs_violation
                << " : relative violation: " << relative_violation
                << " cost: " << current_cost << std::endl;
    }
  }
  return cost;
}

template <int _N>
double PolynomialOptimizationNonLinear<_N>::computeTotalTrajectoryTime(
    const std::vector<double>& segment_times) {
  double total_time = 0;
  for (double t : segment_times) total_time += t;
  return total_time;
}

}  // namespace mav_trajectory_generation

namespace nlopt {

inline std::string returnValueToString(int return_value) {
  switch (return_value) {
    case nlopt::SUCCESS:
      return std::string("SUCCESS");
    case nlopt::FAILURE:
      return std::string("FAILURE");
    case nlopt::INVALID_ARGS:
      return std::string("INVALID_ARGS");
    case nlopt::OUT_OF_MEMORY:
      return std::string("OUT_OF_MEMORY");
    case nlopt::ROUNDOFF_LIMITED:
      return std::string("ROUNDOFF_LIMITED");
    case nlopt::FORCED_STOP:
      return std::string("FORCED_STOP");
    case nlopt::STOPVAL_REACHED:
      return std::string("STOPVAL_REACHED");
    case nlopt::FTOL_REACHED:
      return std::string("FTOL_REACHED");
    case nlopt::XTOL_REACHED:
      return std::string("XTOL_REACHED");
    case nlopt::MAXEVAL_REACHED:
      return std::string("MAXEVAL_REACHED");
    case nlopt::MAXTIME_REACHED:
      return std::string("MAXTIME_REACHED");
    default:
      return std::string("ERROR CODE UNKNOWN");
  }
}
}  // namespace nlopt

#endif  // MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_NONLINEAR_IMPL_H_
