#include "mav_trajectory_generation/polynomial.h"

namespace mav_trajectory_generation {

Eigen::MatrixXd computeBaseCoefficients(int N) {
  Eigen::MatrixXd base_coefficients(N, N);

  base_coefficients.setZero();
  base_coefficients.row(0).setOnes();

  const int DEG = N - 1;
  int order = DEG;
  for (int n = 1; n < N; n++) {
    for (int i = DEG - order; i < N; i++) {
      base_coefficients(n, i) = (order - DEG + i) * base_coefficients(n - 1, i);
    }
    order--;
  }
  return base_coefficients;
}

Eigen::MatrixXd Polynomial::base_coefficients_ =
    computeBaseCoefficients(Polynomial::kMaxN);

}  // namespace mav_trajectory_generation
