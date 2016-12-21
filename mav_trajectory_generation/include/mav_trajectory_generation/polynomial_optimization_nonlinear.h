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

#ifndef MAV_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_NONLINEAR_H_
#define MAV_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_NONLINEAR_H_

#include <memory>
#include <nlopt.hpp>

#include "mav_trajectory_generation/polynomial_optimization_linear.h"

namespace mav_trajectory_generation {

/**
 * \brief Class holding all important parameters for nonlinear optimization.
 */
struct NonlinearOptimizationParameters {
  NonlinearOptimizationParameters()
      // Default parameters should be reasonable enough to use without further
      // fine-tuning.
      : f_abs(-1),
        f_rel(0.05),
        x_rel(-1),
        x_abs(-1),
        initial_stepsize_rel(0.1),
        equality_constraint_tolerance(1.0e-3),
        inequality_constraint_tolerance(0.1),
        max_iterations(3000),
        time_penalty(500.0),
        algorithm(nlopt::LN_SBPLX),
        random_seed(0),
        use_soft_constraints(true),
        soft_constraint_weight(100.0),
        print_debug_info(false) {}

  /// Stopping criteria, if objective function changes less than absolute value.
  /// Disabled if negative.
  double f_abs;

  /// Stopping criteria, if objective function changes less than relative value.
  /// Disabled if negative.
  double f_rel;

  /// Stopping criteria, if state changes less than absolute value. Disabled if
  /// negative.
  double x_rel;

  /// Stopping criteria, if state changes less than relative value. Disabled if
  /// negative.
  double x_abs;

  /// Determines a fraction of the initial guess as initial step size. Heuristic
  /// value if negative.
  double initial_stepsize_rel;

  /// Absolute tolerance, within an equality constraint is considered as met.
  double equality_constraint_tolerance;

  /// Absolute tolerance, within an inequality constraint is considered as met.
  double inequality_constraint_tolerance;

  /// Maximum number of iterations. Disabled if negative.
  int max_iterations;

  /// Penalty for the segment time.
  double time_penalty;

  /// Optimization algorithm used by nlopt, see
  /// http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
  nlopt::algorithm algorithm;

  /**
   * \brief Random seed, if an optimization algorithm involving random numbers
   * is used (e.g. nlopt::GN_ISRES).
   *        If set to a value < 0, a "random" (getTimeofday) value for the seed
   * is chosen.
   */
  int random_seed;

  /// Decide whether to use soft constraints.
  bool use_soft_constraints;

  /// Weights the relative violation of a soft constraint.
  double soft_constraint_weight;

  bool print_debug_info;
};

class OptimizationInfo {
 public:
  OptimizationInfo()
      : n_iterations(0),
        stopping_reason(nlopt::FAILURE),
        cost_trajectory(0),
        cost_time(0),
        cost_soft_constraints(0),
        optimization_time(0) {}
  void print(std::ostream& stream) const;
  int n_iterations;
  int stopping_reason;
  double cost_trajectory;
  double cost_time;
  double cost_soft_constraints;
  double optimization_time;
  std::map<int, Extremum> maxima;
};

/**
 * \brief Implements a nonlinear optimization of the unconstrained optimization
 * of paths
 *        consisting of polynomial segments as described in [1]
 *
 * [1]: Polynomial Trajectory Planning for Aggressive Quadrotor Flight in Dense
 * Indoor Environments.
 *      Charles Richter, Adam Bry, and Nicholas Roy. In ISRR 2013
 *
 * \template _N specifies the number of coefficients for the underlying
 * polynomials. \sa Polynomial
 */
template <int _N = 10>
class PolynomialOptimizationNonLinear {
  static_assert(_N % 2 == 0, "The number of coefficients has to be even.");

 public:
  enum { N = _N };

  /**
   * \brief Sets up the nonlinear optimization problem.
   *
   * \param[in] dimension Spatial dimension of the problem. Usually 1 or 3.
   * \param[in] parameters Parameters for the optimization problem.
   * \param[in] optimize_time_only Specifies whether the optimization is run
   * over the segment times only.
   *                               If true, only the segment times are
   * optimization parameters, and the remaining
   *                               free parameters are found by solving the
   * linear optimization problem with the given
   *                               segment times in every iteration. If false,
   * both segment times and free derivatives
   *                               become optimization variables. The latter
   * case is theoretically correct, but may end
   *                               up in more iterations.
   */
  PolynomialOptimizationNonLinear(
      size_t dimension, const NonlinearOptimizationParameters& parameters,
      bool optimize_time_only);

  /**
   * \brief Sets up the optimization problem from a vector of Vertex objects and
   * a vector of times between the vertices.
   *
   * \param[in] vertices Vector containing the vertices defining the support
   * points and constraints of the path.
   * \param[in] segment_times Vector containing an initial guess of the time
   * between two vertices.
   *                          Thus, its size is size(vertices) - 1.
   * \param[in] derivative_to_optimize Specifies the derivative of which the
   * cost is optimized.
   */
  bool setupFromVertices(
      const Vertex::Vector& vertices, const std::vector<double>& segment_times,
      int derivative_to_optimize =
          PolynomialOptimization<N>::kHighestDerivativeToOptimize);

  /**
   * \brief Adds a constraint for the maximum of magnitude to the optimization
   * problem.
   *
   * \param[in] derivative_order Order of the derivative, for which the
   * constraint should be checked.
   *            Usually velocity (=1) or acceleration (=2).
   * \param[in] maximum_value Maximum magnitude of the specified derivative.
   */
  bool addMaximumMagnitudeConstraint(int derivative_order,
                                     double maximum_value);

  /**
   * \brief Solves the linear optimization problem according to [1].
   *
   * The solver is re-used for every dimension, which means:
   *  - segment times are equal for each dimension.
   *  - each dimension has the same type/set of constraints. Their values can of
   * course differ.
   */
  bool solveLinear();

  /**
   * \brief Runs the optimization until one of the stopping criteria in
   * NonlinearOptimizationParameters
   *        and the constraints are met.
   */
  int optimize();

  // Get the resulting trajectory out -- prefer this as the main method
  // to get the results of the optimization, over getting the reference
  // to the linear optimizer.
  void getTrajectory(Trajectory* trajectory) const {
    poly_opt_.getTrajectory(trajectory);
  }

  /**
   * \brief Returns a const reference to the underlying linear optimization
   * object.
   */
  const PolynomialOptimization<N>& getPolynomialOptimizationRef() const {
    return poly_opt_;
  }

  /**
   * \brief Returns a non-const reference to the underlying linear optimization
   * object.
   */
  PolynomialOptimization<N>& getPolynomialOptimizationRef() {
    return poly_opt_;
  }

  OptimizationInfo getOptimizationInfo() const { return optimization_info_; }

 private:
  /**
   * \brief Holds the data for constraint evaluation, since these methods are
   * static.
   */
  struct ConstraintData {
    PolynomialOptimizationNonLinear<N>* this_object;
    int derivative;
    double value;
  };

  /**
   * \brief Objective function for the time-only version.
   *
   * \param[in] segment_times Segment times in the current iteration.
   * \param[in] gradient Gradient of the objective function w.r.t. changes of
   * parameters. We can't compute the
   *            gradient analytically here. Thus, only gradient free
   * optimization methods are possible.
   * \param[in] data Custum data pointer. In our case, it's an ConstraintData
   * object.
   * \return Cost based on the parameters passed in.
   *
   * \sa PolynomialOptimizationNonLinear
   */
  static double objectiveFunctionTime(const std::vector<double>& segment_times,
                                      std::vector<double>& gradient,
                                      void* data);

  /**
   * \brief Objective function for the version optimizing segment times and free
   * derivatives.
   *
   * \param[in] optimization_variables Optimization variables times in the
   * current iteration. The variables
   *            (time, derivatives) are stacked as follows: [segment_times
   * derivatives_dim_0 ... derivatives_dim_N]
   * \param[in] gradient Gradient of the objective function w.r.t. changes of
   * parameters. We can't compute the
   *            gradient analytically here. Thus, only gradient free
   * optimization methods are possible.
   * \param[in] data Custum data pointer. In our case, it's an ConstraintData
   * object.
   * \return Cost based on the parameters passed in.
   *
   * \sa PolynomialOptimizationNonLinear
   */
  static double objectiveFunctionTimeAndConstraints(
      const std::vector<double>& optimization_variables,
      std::vector<double>& gradient, void* data);

  /**
   * \brief Evaluates the maximum magnitude constraint at the current value of
   * the optimization variables.
   *
   * All input parameters are ignored, all information is contained in data.
   */
  static double evaluateMaximumMagnitudeConstraint(
      const std::vector<double>& optimization_variables,
      std::vector<double>& gradient, void* data);

  /**
   * \brief Does the actual optimization work for the time-only version.
   *
   * \sa PolynomialOptimizationNonLinear
   */
  int optimizeTime();

  /**
   * \brief Does the actual optimization work for the full optimization version.
   *
   * \sa PolynomialOptimizationNonLinear
   */
  int optimizeTimeAndFreeConstraints();

  /**
   * \brief Evaluates the maximum magnitude constraints as soft constraints and
   * returns a cost, depending
   *        on the violation of the constraints.
   *
   *  cost_i = min(maximum_cost, exp(abs_violation_i / max_allowed_i * weight))
   *  cost = sum(cost_i)
   *
   * \param[in] inequality_constraints Vector of ConstraintData shared_ptrs,
   * describing the constraints.
   * \param[in] weight Multiplicative weight of the constraint violation.
   * \param[in] maximum_cost Upper bound of the cost. Necessary, since exp of a
   * high violation can end up in inf.
   * \return Sum of the costs per constraint.
   */
  double evaluateMaximumMagnitudeAsSoftConstraint(
      const std::vector<std::shared_ptr<ConstraintData> >&
          inequality_constraints,
      double weight, double maximum_cost = 1.0e12) const;

  /**
   * \brief Computes the total trajectory time.
   */
  static double computeTotalTrajectoryTime(
      const std::vector<double>& segment_times);

  /// nlopt optimization object.
  std::shared_ptr<nlopt::opt> nlopt_;

  /// Underlying linear optimization object. \sa PolynomialOptimization
  PolynomialOptimization<N> poly_opt_;

  /// Parameters for the nonlinear optimzation.
  NonlinearOptimizationParameters optimization_parameters_;

  /// Holds the data for evaluating inequality constraints.
  std::vector<std::shared_ptr<ConstraintData> > inequality_constraints_;

  /// Specifies whether to run the time only, or the full optimization.
  bool optimize_time_only_;

  OptimizationInfo optimization_info_;
};

}  // namespace mav_trajectory_generation

namespace nlopt {
/// Convenience function that turns nlopt's return values into something
/// readable.
std::string returnValueToString(int return_value);
}  // namespace nlopt

#endif  // MAV_TRAJECTORY_GENERATION_POLYNOMIAL_OPTIMIZATION_NONLINEAR_H_

#include "mav_trajectory_generation/impl/polynomial_optimization_nonlinear_impl.h"
