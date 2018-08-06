/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002
 *                       Future Team Aps
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Future Team Software License Agreement which restricts the manner
 *   in which it may be used.
 *
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :   newr.cpp
 * Module ID Nbr   :
 * Description     :   Solve n degree polynominal using Newton's (Madsen)
 *methode
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version  Author/Date   Description of changes
 * -------  -----------   ----------------------
 * 01.01    HVE/021018    Initial release
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <complex>
#include <iostream>

#include "mav_trajectory_generation/cpoly/newr.h"
// For findLastNonZeroCoeff()
#include "mav_trajectory_generation/rpoly/rpoly_ak1.h"

namespace mav_trajectory_generation {

int newtonWrapper(int degree, const double* coefficients,
                  std::complex<double>* roots);

bool findRootsNewtonsMethod(const Eigen::VectorXd& coefficients_increasing,
                            Eigen::VectorXcd* roots) {
  // Remove trailing zeros.
  const int last_non_zero_coefficient =
      findLastNonZeroCoeff(coefficients_increasing);
  if (last_non_zero_coefficient == -1) {
    // The polynomial has all zero coefficients and has no roots.
    roots->resize(0);
    return true;
  }

  // Reverse coefficients in descending order.
  Eigen::VectorXd coefficients_decreasing =
      coefficients_increasing.head(last_non_zero_coefficient + 1).reverse();

  const int n_coefficients = coefficients_decreasing.size();
  if (n_coefficients < 2) {
    // The polynomial is 0th order and has no roots.
    roots->resize(0);
    return true;
  }
  std::complex<double>* roots_arr = new std::complex<double>[n_coefficients];

  int ret = newtonWrapper(n_coefficients - 1, coefficients_decreasing.data(),
                          roots_arr);

  if (ret > -1) {
    roots->resize(n_coefficients);
    for (int i = 0; i < n_coefficients; ++i) {
      (*roots)[i] = roots_arr[i];
    }
  }

  delete[] roots_arr;

  if (ret > -1) {
    return true;
  } else {
    return false;
  }
}

namespace newr_impl {
using namespace std;

#define MAXITER 50

/* Solve linear or quadratic equation
 */
static void quadratic(const int n, const double a[], complex<double> res[]) {
  double r;

  if (n == 2) {
    if (a[1] == 0) {
      r = -a[2] / a[0];
      if (r < 0) {
        res[1] = complex<double>(0, sqrt(-r));
        res[2] = complex<double>(0, -res[1].imag());
      } else {
        res[1] = complex<double>(sqrt(r), 0);
        res[2] = complex<double>(-res[1].real(), 0);
      }
    } else {
      r = 1 - 4 * a[0] * a[2] / (a[1] * a[1]);
      if (r < 0) {
        res[1] =
            complex<double>(-a[1] / (2 * a[0]), a[1] * sqrt(-r) / (2 * a[0]));
        res[2] = complex<double>(res[1].real(), -res[1].imag());
      } else {
        res[1] = complex<double>((-1 - sqrt(r)) * a[1] / (2 * a[0]), 0);
        res[2] = complex<double>(a[2] / (a[0] * res[1].real()), 0);
      }
    }
  } else if (n == 1)
    res[1] = complex<double>(-a[1] / a[0], 0);
}

/* Performed function evaluation. Horners algorithm.
 */
static double feval(const int n, const double a[], const complex<double> z,
                    complex<double>* fz) {
  int i;
  double p, q, r, s, t;

  p = -2.0 * z.real();
  q = z.real() * z.real() + z.imag() * z.imag();
  s = 0;
  r = a[0];
  for (i = 1; i < n; i++) {
    t = a[i] - p * r - q * s;
    s = r;
    r = t;
  }
  *fz = complex<double>(a[n] + z.real() * r - q * s, z.imag() * r);

  return fz->real() * fz->real() + fz->imag() * fz->imag();
}

static double startpoint(const int n, const double a[]) {
  int i;
  double r, u, min;

  /* Determine starting point */
  r = log(fabs(a[n]));
  min = exp((r - log(fabs(a[0]))) / n);
  for (i = 1; i < n; i++)
    if (a[i] != 0) {
      u = exp((r - log(fabs(a[i]))) / (n - i));
      if (u < min) min = u;
    }

  return min;
}

/* Calculate a upper bound for the rounding errors performed in a
   polynomial at a complex point.
   ( Adam's test )
 */
static double upperbound(const int n, const double a[],
                         const complex<double> z) {
  int i;
  double p, q, r, s, t, u, e;

  p = -2.0 * z.real();
  q = z.real() * z.real() + z.imag() * z.imag();
  u = sqrt(q);
  s = 0.0;
  r = a[0];
  e = fabs(r) * (3.5 / 4.5);
  for (i = 1; i < n; i++) {
    t = a[i] - p * r - q * s;
    s = r;
    r = t;
    e = u * e + fabs(t);
  }
  t = a[n] + z.real() * r - q * s;
  e = u * e + fabs(t);
  e = (9.0 * e - 7.0 * (fabs(t) + fabs(r) * u) +
       fabs(z.real()) * fabs(r) * 2.0) *
      std::numeric_limits<double>::epsilon();

  return e * e;
}

static void alterdirection(complex<double>* dz, const double m) {
  double x, y;

  x = (dz->real() * 0.6 - dz->imag() * 0.8) * m;
  y = (dz->real() * 0.8 + dz->imag() * 0.6) * m;
  *dz = complex<double>(x, y);
}

/* Real root forward deflation.
 */
static int realdeflation(const int n, double a[], const double x) {
  int i;
  double r;

  for (r = 0, i = 0; i < n; i++) a[i] = r = r * x + a[i];
  return n - 1;
}

/* Complex root forward deflation.
 */
static int complexdeflation(const int n, double a[], const complex<double> z) {
  int i;
  double r, u;

  r = -2.0 * z.real();
  u = z.real() * z.real() + z.imag() * z.imag();
  a[1] -= r * a[0];
  for (i = 2; i < n - 1; i++) a[i] = a[i] - r * a[i - 1] - u * a[i - 2];

  return n - 2;
}

// Find all root of a polynomial of n degree with real coeeficient using the
// modified Newton by Madsen
//
int newton_real(int n, const double coeff[], complex<double> res[]) {
  int i, itercnt;
  int stage1, div2;
  int err;
  double *a, *a1;
  double u, r, r0, eps;
  double f, f0, ff, f2, fw;
  complex<double> z, z0, dz, fz, fwz, wz, fz0, fz1;

  a = new double[n + 1];
  for (i = 0; i <= n; i++) a[i] = coeff[i];

  err = 0;
  for (; a[n] == 0.0; n--) res[n] = complex<double>(0);

  a1 = new double[n];
  while (n > 2) {
    /* Calculate coefficients of f'(x) */
    for (i = 0; i < n; i++) a1[i] = a[i] * (n - i);

    u = startpoint(n, a);
    z0 = complex<double>(0);
    f0 = ff = 2.0 * a[n] * a[n];
    fz0 = complex<double>(a[n - 1]);
    z = complex<double>(a[n - 1] == 0.0 ? 1 : -a[n] / a[n - 1], 0);
    z = complex<double>(z.real() / (double)fabs(z.real()) * u * 0.5, 0);
    dz = z;
    f = feval(n, a, z, &fz);
    r0 = 2.5 * u;
    r = sqrt(dz.real() * dz.real() + dz.imag() * dz.imag());
    eps = 4 * n * n * f0 * std::numeric_limits<double>::epsilon();

    // Start iteration
    for (stage1 = 1, itercnt = 0; (z.real() + dz.real() != z.real() ||
                                   z.imag() + dz.imag() != z.imag()) &&
                                  f > eps && itercnt < MAXITER;
         itercnt++) { /* Iterativ loop */
      u = feval(n - 1, a1, z, &fz1);
      if (u == 0.0) /* True saddelpoint */
        alterdirection(&dz, 5.0);
      else {
        dz = complex<double>(
            (fz.real() * fz1.real() + fz.imag() * fz1.imag()) / u,
            (fz.imag() * fz1.real() - fz.real() * fz1.imag()) / u);

        /* Which stage are we on */
        fwz = fz0 - fz1;
        wz = z0 - z;
        f2 = (fwz.real() * fwz.real() + fwz.imag() * fwz.imag()) /
             (wz.real() * wz.real() + wz.imag() * wz.imag());
        stage1 = f2 / u > u / f / 4 || f != ff;
        r = sqrt(dz.real() * dz.real() + dz.imag() * dz.imag());
        if (r > r0) alterdirection(&dz, r0 / r);
        r0 = r * 5.0;
      }

      z0 = z;
      f0 = f;
      fz0 = fz;

    iter2:
      z = z0 - dz;
      ff = f = feval(n, a, z, &fz);
      if (stage1) {
        wz = z;
        for (i = 1, div2 = f > f0; i <= n; i++) {
          if (div2) {
            dz *= complex<double>(0.5);
            wz = z0 - dz;
          } else {
            wz -= dz;
          }
          fw = feval(n, a, wz, &fwz);
          if (fw >= f) break;
          f = fw;
          fz = fwz;
          z = wz;
          if (div2 && i == 2) {
            alterdirection(&dz, 0.5);
            z = z0 - dz;
            f = feval(n, a, z, &fz);
            break;
          }
        }
      } else {
        /* calculate the upper bound of erros using Adam's test */
        eps = upperbound(n, a, z);
      }

      if (r < sqrt(z.real() * z.real() + z.imag() * z.imag()) *
                  sqrt(std::numeric_limits<double>::epsilon()) &&
          f >= f0) {
        /* Domain rounding errors */
        z = z0;
        alterdirection(&dz, 0.5);
        if (z + dz != z) goto iter2;
      }
    }

    if (itercnt >= MAXITER) err--;

    z0 = complex<double>(z.real(), 0.0);
    if (feval(n, a, z0, &fz) <= f) {
      /* Real root */
      res[n] = complex<double>(z.real(), 0);
      n = realdeflation(n, a, z.real());
    } else {
      /* Complex root */
      res[n] = z;
      res[n - 1] = complex<double>(z.real(), -z.imag());
      n = complexdeflation(n, a, z);
    }
  }

  quadratic(n, a, res);

  delete[] a1;
  delete[] a;

  return (err);
}

}  // namespace newr_impl;

// Find all root of a polynomial of n degree with real coeeficient using the
// modified Newton by Madsen
//

int newtonWrapper(int degree, const double* coefficients,
                  std::complex<double>* roots) {
  return newr_impl::newton_real(degree, coefficients, roots);
}

}  // namespace mav_trajectory_generation
