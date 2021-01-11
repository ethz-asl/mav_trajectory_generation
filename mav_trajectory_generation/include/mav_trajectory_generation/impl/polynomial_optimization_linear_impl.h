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

#ifndef MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_LINEAR_IMPL_H_
#define MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_LINEAR_IMPL_H_

#include <glog/logging.h>
#include <Eigen/Sparse>
#include <set>
#include <tuple>

// fixes error due to std::iota (has been introduced in c++ standard lately
// and may cause compilation errors depending on compiler)
#if __cplusplus <= 199711L
  #include <algorithm>
#else
  #include <numeric>
#endif

#include "mav_trajectory_generation/convolution.h"



namespace mav_trajectory_generation {

template <int _N>
PolynomialOptimization<_N>::PolynomialOptimization(size_t dimension)
    : dimension_(dimension),
      derivative_to_optimize_(derivative_order::INVALID),
      n_vertices_(0),
      n_segments_(0),
      n_all_constraints_(0),
      n_fixed_constraints_(0),
      n_free_constraints_(0) {
  fixed_constraints_compact_.resize(dimension_);
  free_constraints_compact_.resize(dimension_);
}

template <int _N>
bool PolynomialOptimization<_N>::setupFromVertices(
    const Vertex::Vector& vertices, const std::vector<double>& times,
    int derivative_to_optimize) {
  CHECK(derivative_to_optimize >= 0 &&
        derivative_to_optimize <= kHighestDerivativeToOptimize)
      << "You tried to optimize the " << derivative_to_optimize
      << "th derivative of position on a " << N
      << "th order polynomial. This is not possible, you either need a higher "
         "order polynomial or a smaller derivative to optimize.";

  derivative_to_optimize_ = derivative_to_optimize;
  vertices_ = vertices;
  segment_times_ = times;

  n_vertices_ = vertices.size();
  n_segments_ = n_vertices_ - 1;

  segments_.resize(n_segments_, Segment(N, dimension_));

  CHECK(n_vertices_ == times.size() + 1)
      << "Size of times must be one less than positions.";

  inverse_mapping_matrices_.resize(n_segments_);
  cost_matrices_.resize(n_segments_);

  // Iterate through all vertices and remove invalid constraints (order too
  // high).
  for (size_t vertex_idx = 0; vertex_idx < n_vertices_; ++vertex_idx) {
    Vertex& vertex = vertices_[vertex_idx];

    // Check if we have valid constraints.
    bool vertex_valid = true;
    Vertex vertex_tmp(dimension_);
    for (Vertex::Constraints::const_iterator it = vertex.cBegin();
         it != vertex.cEnd(); ++it) {
      if (it->first > kHighestDerivativeToOptimize) {
        vertex_valid = false;
        LOG(WARNING) << "Invalid constraint on vertex " << vertex_idx
                     << ": maximum possible derivative is "
                     << kHighestDerivativeToOptimize << ", but was set to "
                     << it->first << ". Ignoring constraint";
      } else {
        vertex_tmp.addConstraint(it->first, it->second);
      }
    }
    if (!vertex_valid) {
      vertex = vertex_tmp;
    }
  }
  updateSegmentTimes(times);
  setupConstraintReorderingMatrix();
  return true;
}

template <int _N>
void PolynomialOptimization<_N>::setupMappingMatrix(double segment_time,
                                                    SquareMatrix* A) {
  // The sum of fixed/free variables has to be equal on both ends of the
  // segment.
  // Thus, A is created as [A(t=0); A(t=segment_time)].
  for (int i = 0; i < N / 2; ++i) {
    A->row(i) = Polynomial::baseCoeffsWithTime(N, i, 0.0);
    A->row(i + N / 2) = Polynomial::baseCoeffsWithTime(N, i, segment_time);
  }
}

template <int _N>
double PolynomialOptimization<_N>::computeCost() const {
  CHECK(n_segments_ == segments_.size() &&
        n_segments_ == cost_matrices_.size());
  double cost = 0;
  for (size_t segment_idx = 0; segment_idx < n_segments_; ++segment_idx) {
    const SquareMatrix& Q = cost_matrices_[segment_idx];
    const Segment& segment = segments_[segment_idx];
    for (size_t dimension_idx = 0; dimension_idx < dimension_;
         ++dimension_idx) {
      const Eigen::VectorXd c =
          segment[dimension_idx].getCoefficients(derivative_order::POSITION);
      const double partial_cost = c.transpose() * Q * c;
      cost += partial_cost;
    }
  }
  return 0.5 * cost;  // cost = 0.5 * c^T * Q * c
}

template <int _N>
void PolynomialOptimization<_N>::invertMappingMatrix(
    const SquareMatrix& mapping_matrix, SquareMatrix* inverse_mapping_matrix) {
  // The mapping matrix has the following structure:
  // [ x 0 0 0 0 0 ]
  // [ 0 x 0 0 0 0 ]
  // [ 0 0 x 0 0 0 ]
  // [ x x x x x x ]
  // [ 0 x x x x x ]
  // [ 0 0 x x x x ]
  // ==>
  // [ A_diag B=0 ]
  // [ C      D   ]
  // We make use of the Schur-complement, so the inverse is:
  // [ inv(A_diag)               0      ]
  // [ -inv(D) * C * inv(A_diag) inv(D) ]
  const int half_n = N / 2;

  // "template" keyword required below as half_n is dependent on the template
  // parameter.
  const Eigen::Matrix<double, half_n, 1> A_diag =
      mapping_matrix.template block<half_n, half_n>(0, 0).diagonal();
  const Eigen::Matrix<double, half_n, half_n> A_inv =
      A_diag.cwiseInverse().asDiagonal();

  const Eigen::Matrix<double, half_n, half_n> C =
      mapping_matrix.template block<half_n, half_n>(half_n, 0);

  const Eigen::Matrix<double, half_n, half_n> D_inv =
      mapping_matrix.template block<half_n, half_n>(half_n, half_n).inverse();

  inverse_mapping_matrix->template block<half_n, half_n>(0, 0) = A_inv;
  inverse_mapping_matrix->template block<half_n, half_n>(0, half_n).setZero();
  inverse_mapping_matrix->template block<half_n, half_n>(half_n, 0) =
      -D_inv * C * A_inv;
  inverse_mapping_matrix->template block<half_n, half_n>(half_n, half_n) =
      D_inv;
}

template <int _N>
void PolynomialOptimization<_N>::setupConstraintReorderingMatrix() {
  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> reordering_list;

  const size_t n_vertices = vertices_.size();

  std::vector<Constraint> all_constraints;
  std::set<Constraint> fixed_constraints;
  std::set<Constraint> free_constraints;

  all_constraints.reserve(
      n_vertices_ * N /
      2);  // Will have exactly this number of elements in the end.

  for (size_t vertex_idx = 0; vertex_idx < n_vertices; ++vertex_idx) {
    const Vertex& vertex = vertices_[vertex_idx];

    // Extract constraints and sort them to fixed and free. For the start and
    // end Vertex, we need to do this once, while we need to do it twice for the
    // other vertices, since constraints are shared and enforce continuity.
    int n_constraint_occurence = 2;
    if (vertex_idx == 0 || vertex_idx == (n_segments_))
      n_constraint_occurence = 1;
    for (int co = 0; co < n_constraint_occurence; ++co) {
      for (size_t constraint_idx = 0; constraint_idx < N / 2;
           ++constraint_idx) {
        Constraint constraint;
        constraint.vertex_idx = vertex_idx;
        constraint.constraint_idx = constraint_idx;
        bool has_constraint =
            vertex.getConstraint(constraint_idx, &(constraint.value));
        if (has_constraint) {
          all_constraints.push_back(constraint);
          fixed_constraints.insert(constraint);
        } else {
          constraint.value = Vertex::ConstraintValue::Constant(dimension_, 0);
          all_constraints.push_back(constraint);
          free_constraints.insert(constraint);
        }
      }
    }
  }

  n_all_constraints_ = all_constraints.size();
  n_fixed_constraints_ = fixed_constraints.size();
  n_free_constraints_ = free_constraints.size();

  reordering_list.reserve(n_all_constraints_);
  constraint_reordering_ = Eigen::SparseMatrix<double>(
      n_all_constraints_, n_fixed_constraints_ + n_free_constraints_);

  for (Eigen::VectorXd& df : fixed_constraints_compact_)
    df.resize(n_fixed_constraints_, Eigen::NoChange);

  int row = 0;
  int col = 0;
  for (const Constraint& ca : all_constraints) {
    for (const Constraint& cf : fixed_constraints) {
      if (ca == cf) {
        reordering_list.emplace_back(Triplet(row, col, 1.0));
        for (size_t d = 0; d < dimension_; ++d) {
          Eigen::VectorXd& df = fixed_constraints_compact_[d];
          const Eigen::VectorXd constraint_all_dimensions = cf.value;
          df[col] = constraint_all_dimensions[d];
        }
      }
      ++col;
    }
    for (const Constraint& cp : free_constraints) {
      if (ca == cp) reordering_list.emplace_back(Triplet(row, col, 1.0));
      ++col;
    }
    col = 0;
    ++row;
  }

  constraint_reordering_.setFromTriplets(reordering_list.begin(),
                                         reordering_list.end());
}

template <int _N>
void PolynomialOptimization<_N>::updateSegmentsFromCompactConstraints() {
  const size_t n_all_constraints = n_fixed_constraints_ + n_free_constraints_;

  for (size_t dimension_idx = 0; dimension_idx < dimension_; ++dimension_idx) {
    const Eigen::VectorXd& df = fixed_constraints_compact_[dimension_idx];
    const Eigen::VectorXd& dp_opt = free_constraints_compact_[dimension_idx];

    Eigen::VectorXd d_all(n_all_constraints);
    d_all << df, dp_opt;

    for (size_t i = 0; i < n_segments_; ++i) {
      const Eigen::Matrix<double, N, 1> new_d =
          constraint_reordering_.block(i * N, 0, N, n_all_constraints) * d_all;
      const Eigen::Matrix<double, N, 1> coeffs =
          inverse_mapping_matrices_[i] * new_d;
      Segment& segment = segments_[i];
      segment.setTime(segment_times_[i]);
      segment[dimension_idx] = Polynomial(N, coeffs);
    }
  }
}

template <int _N>
void PolynomialOptimization<_N>::updateSegmentTimes(
    const std::vector<double>& segment_times) {
  const size_t n_segment_times = segment_times.size();
  CHECK(n_segment_times == n_segments_)
      << "Number of segment times (" << n_segment_times
      << ") does not match number of segments (" << n_segments_ << ")";

  segment_times_ = segment_times;

  for (size_t i = 0; i < n_segments_; ++i) {
    const double segment_time = segment_times[i];
    CHECK_GT(segment_time, 0) << "Segment times need to be greater than zero";

    computeQuadraticCostJacobian(derivative_to_optimize_, segment_time,
                                 &cost_matrices_[i]);
    SquareMatrix A;
    setupMappingMatrix(segment_time, &A);
    invertMappingMatrix(A, &inverse_mapping_matrices_[i]);
  };
}

template <int _N>
void PolynomialOptimization<_N>::constructR(
    Eigen::SparseMatrix<double>* R) const {
  CHECK_NOTNULL(R);
  typedef Eigen::Triplet<double> Triplet;
  std::vector<Triplet> cost_unconstrained_triplets;
  cost_unconstrained_triplets.reserve(N * N * n_segments_);

  for (size_t i = 0; i < n_segments_; ++i) {
    const SquareMatrix& Ai = inverse_mapping_matrices_[i];
    const SquareMatrix& Q = cost_matrices_[i];
    const SquareMatrix H = Ai.transpose() * Q * Ai;
    const int start_pos = i * N;
    for (int row = 0; row < N; ++row) {
      for (int col = 0; col < N; ++col) {
        cost_unconstrained_triplets.emplace_back(
            Triplet(start_pos + row, start_pos + col, H(row, col)));
      }
    }
  }
  Eigen::SparseMatrix<double> cost_unconstrained(N * n_segments_,
                                                 N * n_segments_);
  cost_unconstrained.setFromTriplets(cost_unconstrained_triplets.begin(),
                                     cost_unconstrained_triplets.end());

  // [1]: R = C^T * H * C. C: constraint_reodering_ ; H: cost_unconstrained,
  // assembled from the block-H above.
  *R = constraint_reordering_.transpose() * cost_unconstrained *
       constraint_reordering_;
}

template <int _N>
bool PolynomialOptimization<_N>::solveLinear() {
  CHECK(derivative_to_optimize_ >= 0 &&
        derivative_to_optimize_ <= kHighestDerivativeToOptimize);
  // Catch the fully constrained case:
  if (n_free_constraints_ == 0) {
    DLOG(WARNING)
        << "No free constraints set in the vertices. Polynomial can "
           "not be optimized. Outputting fully constrained polynomial.";
    updateSegmentsFromCompactConstraints();
    return true;
  }

  // TODO(acmarkus): figure out if sparse becomes less efficient for small
  // problems, and switch back to dense in case.

  // Compute cost matrix for the unconstrained optimization problem.
  // Block-wise H = A^{-T}QA^{-1} according to [1]
  Eigen::SparseMatrix<double> R;
  constructR(&R);

  // Extract block matrices and prepare solver.
  Eigen::SparseMatrix<double> Rpf = R.block(
      n_fixed_constraints_, 0, n_free_constraints_, n_fixed_constraints_);
  Eigen::SparseMatrix<double> Rpp =
      R.block(n_fixed_constraints_, n_fixed_constraints_, n_free_constraints_,
              n_free_constraints_);
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      solver;
  solver.compute(Rpp);

  // Compute dp_opt for every dimension.
  for (size_t dimension_idx = 0; dimension_idx < dimension_; ++dimension_idx) {
    Eigen::VectorXd df =
        -Rpf * fixed_constraints_compact_[dimension_idx];  // Rpf = Rfp^T
    free_constraints_compact_[dimension_idx] =
        solver.solve(df);  // dp = -Rpp^-1 * Rpf * df
  }

  updateSegmentsFromCompactConstraints();
  return true;
}

template <int _N>
void PolynomialOptimization<_N>::printReorderingMatrix(
    std::ostream& stream) const {
  stream << "Mapping matrix:\n" << constraint_reordering_ << std::endl;
}

template <int _N>
template <int Derivative>
bool PolynomialOptimization<_N>::computeSegmentMaximumMagnitudeCandidates(
    const Segment& segment, double t_start, double t_stop,
    std::vector<double>* candidates) {
  return computeSegmentMaximumMagnitudeCandidates(Derivative, segment, t_start,
                                                  t_stop, candidates);
}

template <int _N>
bool PolynomialOptimization<_N>::computeSegmentMaximumMagnitudeCandidates(
    int derivative, const Segment& segment, double t_start, double t_stop,
    std::vector<double>* candidates) {
  CHECK(candidates);
  CHECK(N - derivative - 1 > 0) << "N-Derivative-1 has to be greater 0";

  // Use the implementation of this in the segment (template-free) as it's
  // actually faster.
  std::vector<int> dimensions(segment.D());
  std::iota(dimensions.begin(), dimensions.end(), 0);
  return segment.computeMinMaxMagnitudeCandidateTimes(
      derivative, t_start, t_stop, dimensions, candidates);
}

template <int _N>
template <int Derivative>
void PolynomialOptimization<_N>::
    computeSegmentMaximumMagnitudeCandidatesBySampling(
        const Segment& segment, double t_start, double t_stop, double dt,
        std::vector<double>* candidates) {
  CHECK_NOTNULL(candidates);
  // Start is candidate.
  candidates->push_back(t_start);

  // Determine initial direction from t_start to t_start + dt.
  auto t_old = t_start + dt;
  auto value_new = segment.evaluate(t_old, Derivative);
  auto value_old = segment.evaluate(t_start, Derivative);
  auto direction = value_new.norm() - value_old.norm();

  // Continue with direction from t_start + dt to t_start + 2 dt until t_stop.
  bool last_sample = false;
  for (double t = t_start + dt + dt; t <= t_stop; t += dt) {
    // Update direction.
    value_old = value_new;
    value_new = segment.evaluate(t, Derivative);
    auto direction_new = value_new.norm() - value_old.norm();

    if (std::signbit(direction) != std::signbit(direction_new)) {
      auto value_deriv = segment.evaluate(t_old, Derivative + 1);
      if (value_deriv.norm() < 1e-2) {
        candidates->push_back(t_old);  // extremum was at last dt
      }
    }

    direction = direction_new;
    t_old = t;

    // Check last sample before t_stop.
    if ((t + dt) > t_stop && !last_sample) {
      t = t_stop - dt;
      last_sample = true;
    }
  }

  // End is candidates.
  if (candidates->back() != t_stop) {
    candidates->push_back(t_stop);
  }
}

template <int _N>
template <int Derivative>
Extremum PolynomialOptimization<_N>::computeMaximumOfMagnitude(
    std::vector<Extremum>* candidates) const {
  return computeMaximumOfMagnitude(Derivative, candidates);
}

template <int _N>
Extremum PolynomialOptimization<_N>::computeMaximumOfMagnitude(
    int derivative, std::vector<Extremum>* candidates) const {
  if (candidates != nullptr) candidates->clear();

  int segment_idx = 0;
  Extremum extremum;
  for (const Segment& s : segments_) {
    std::vector<double> extrema_times;
    extrema_times.reserve(N - 1);
    // Add the beginning as well. Call below appends its extrema.
    extrema_times.push_back(0.0);
    computeSegmentMaximumMagnitudeCandidates(derivative, s, 0.0, s.getTime(),
                                             &extrema_times);

    for (double t : extrema_times) {
      const Extremum candidate(t, s.evaluate(t, derivative).norm(),
                               segment_idx);
      if (extremum < candidate) extremum = candidate;
      if (candidates != nullptr) candidates->emplace_back(candidate);
    }
    ++segment_idx;
  }
  // Check last time at last segment.
  const Extremum candidate(
      segments_.back().getTime(),
      segments_.back().evaluate(segments_.back().getTime(), derivative).norm(),
      n_segments_ - 1);
  if (extremum < candidate) extremum = candidate;
  if (candidates != nullptr) candidates->emplace_back(candidate);

  return extremum;
}

template <int _N>
void PolynomialOptimization<_N>::setFreeConstraints(
    const std::vector<Eigen::VectorXd>& free_constraints) {
  CHECK(free_constraints.size() == dimension_);
  for (const Eigen::VectorXd& v : free_constraints)
    CHECK(static_cast<size_t>(v.size()) == n_free_constraints_);

  free_constraints_compact_ = free_constraints;
  updateSegmentsFromCompactConstraints();
}

template <int _N>
void PolynomialOptimization<_N>::getAInverse(Eigen::MatrixXd* A_inv) const {
  CHECK_NOTNULL(A_inv);

  A_inv->resize(N * n_segments_, N * n_segments_);
  A_inv->setZero();

  for (size_t i = 0; i < n_segments_; ++i) {
    (*A_inv).block<N, N>(N * i, N * i) = inverse_mapping_matrices_[i];
  }
}

template <int _N>
void PolynomialOptimization<_N>::getM(Eigen::MatrixXd* M) const {
  CHECK_NOTNULL(M);
  *M = constraint_reordering_;
}

template <int _N>
void PolynomialOptimization<_N>::getR(Eigen::MatrixXd* R) const {
  CHECK_NOTNULL(R);

  Eigen::SparseMatrix<double> R_sparse;
  constructR(&R_sparse);

  *R = R_sparse;
}

template <int _N>
void PolynomialOptimization<_N>::getA(Eigen::MatrixXd* A) const {
  CHECK_NOTNULL(A);
  A->resize(N * n_segments_, N * n_segments_);
  A->setZero();

  // Create a mapping matrix per segment and append them together.
  for (size_t i = 0; i < n_segments_; ++i) {
    const double segment_time = segment_times_[i];
    CHECK_GT(segment_time, 0) << "Segment times need to be greater than zero";

    SquareMatrix A_segment;
    setupMappingMatrix(segment_time, &A_segment);

    (*A).block<N, N>(N * i, N * i) = A_segment;
  }
}

template <int _N>
void PolynomialOptimization<_N>::getMpinv(Eigen::MatrixXd* M_pinv) const {
  CHECK_NOTNULL(M_pinv);

  // Pseudoinverse implementation by @SebastianInd.
  *M_pinv = constraint_reordering_.transpose();
  for (int M_row = 0; M_row < M_pinv->rows(); M_row++) {
    M_pinv->row(M_row) = M_pinv->row(M_row) / M_pinv->row(M_row).sum();
  }
}

template <int _N>
void PolynomialOptimization<_N>::computeQuadraticCostJacobian(
    int derivative, double t, SquareMatrix* cost_jacobian) {
  CHECK_LT(derivative, N);

  cost_jacobian->setZero();
  for (int col = 0; col < N - derivative; col++) {
    for (int row = 0; row < N - derivative; row++) {
      double exponent = (N - 1 - derivative) * 2 + 1 - row - col;

      (*cost_jacobian)(N - 1 - row, N - 1 - col) =
          Polynomial::base_coefficients_(derivative, N - 1 - row) *
          Polynomial::base_coefficients_(derivative, N - 1 - col) *
          pow(t, exponent) * 2.0 / exponent;
    }
  }
}

}  // namespace mav_trajectory_generation

#endif  // MAV_TRAJECTORY_GENERATION_IMPL_POLYNOMIAL_OPTIMIZATION_LINEAR_IMPL_H_
