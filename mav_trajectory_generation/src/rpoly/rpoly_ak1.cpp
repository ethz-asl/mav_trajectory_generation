// rpoly_ak1.cpp - Program for calculating the roots of a polynomial of real
// coefficients.
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
// 27 May 2014
//
// The sub-routines listed below are translations of the FORTRAN routines
// included in RPOLY.FOR,
// posted off the NETLIB site as TOMS/493:
//
// http://www.netlib.org/toms/493
//
// TOMS/493 is based on the Jenkins-Traub algorithm.
//
// To distinguish the routines posted below from others, an _ak1 suffix has been
// appended to them.
//
// Following is a list of the major changes made in the course of translating
// the TOMS/493 routines
// to the C++ versions posted below:
// 1) All global variables have been eliminated.
// 2) The "FAIL" parameter passed into RPOLY.FOR has been eliminated.
// 3) RPOLY.FOR solves polynomials of degree up to 100, but does not explicitly
// state this limit.
//     rpoly_ak1 explicitly states this limit; uses the macro name MAXDEGREE to
//     specify this limit;
//     and does a check to ensure that the user input variable Degree is not
//     greater than MAXDEGREE
//     (if it is, an error message is output and rpoly_ak1 terminates). If a
//     user wishes to compute
//     roots of polynomials of degree greater than MAXDEGREE, using a macro name
//     like MAXDEGREE provides
//     the simplest way of offering this capability.
// 4) All "GO TO" statements have been eliminated.
//
// A small main program is included also, to provide an example of how to use
// rpoly_ak1. In this
// example, data is input from a file to eliminate the need for a user to type
// data in via
// the console.

#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <cfloat>

#include "mav_trajectory_generation/rpoly/rpoly_ak1.h"

namespace mav_trajectory_generation {

constexpr int kRpolyMaxDegree = 100;

// Wraps the call to rpoly_ak1.
void rpolyWrapper(double* coefficients_decreasing, int* degree,
                 double* roots_real, double* roots_imag);

int findLastNonZeroCoeff(const Eigen::VectorXd& coefficients) {
  int last_non_zero_coefficient = -1;

  // Find last non-zero coefficient:
  for (int i = coefficients.size() - 1; i != -1; i--) {
    if (std::abs(coefficients(i)) >= std::numeric_limits<double>::min()) {
      last_non_zero_coefficient = i;
      break;
    }
  }
  return last_non_zero_coefficient;
}

bool findRootsJenkinsTraub(const Eigen::VectorXd& coefficients_increasing,
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
  int degree = n_coefficients - 1;
  double* polynomial = new double[kRpolyMaxDegree + 1];
  double* roots_real = new double[kRpolyMaxDegree];
  double* roots_imag = new double[kRpolyMaxDegree];

  for (size_t i = 0; i < n_coefficients; i++) {
    polynomial[i] = coefficients_decreasing(i);
  }

  rpolyWrapper(polynomial, &degree, roots_real, roots_imag);
  if (degree > 0) {
    roots->resize(degree);
    for (int i = 0; i < degree; ++i) {
      (*roots)[i] = std::complex<double>(roots_real[i], roots_imag[i]);
    }
  }

  delete[] polynomial;
  delete[] roots_real;
  delete[] roots_imag;

  if (degree > 0) {
    return true;
  } else {
    return false;
  }
}

namespace rpoly_impl {

using namespace std;

#define MAXDEGREE 100
#define MDP1 MAXDEGREE + 1

void rpoly_ak1(double op[MDP1], int* Degree, double zeror[MAXDEGREE],
               double zeroi[MAXDEGREE]);
void Fxshfr_ak1(int L2, int* NZ, double sr, double bnd, double K[MDP1], int N,
                double p[MDP1], int NN, double qp[MDP1], double* lzi,
                double* lzr, double* szi, double* szr);
void QuadSD_ak1(int NN, double u, double v, double p[MDP1], double q[MDP1],
                double* a, double* b);
int calcSC_ak1(int N, double a, double b, double* a1, double* a3, double* a7,
               double* c, double* d, double* e, double* f, double* g, double* h,
               double K[MDP1], double u, double v, double qk[MDP1]);
void nextK_ak1(int N, int tFlag, double a, double b, double a1, double* a3,
               double* a7, double K[MDP1], double qk[MDP1], double qp[MDP1]);
void newest_ak1(int tFlag, double* uu, double* vv, double a, double a1,
                double a3, double a7, double b, double c, double d, double f,
                double g, double h, double u, double v, double K[MDP1], int N,
                double p[MDP1]);
void QuadIT_ak1(int N, int* NZ, double uu, double vv, double* szr, double* szi,
                double* lzr, double* lzi, double qp[MDP1], int NN, double* a,
                double* b, double p[MDP1], double qk[MDP1], double* a1,
                double* a3, double* a7, double* d, double* e, double* f,
                double* g, double* h, double K[MDP1]);
void RealIT_ak1(int* iFlag, int* NZ, double* sss, int N, double p[MDP1], int NN,
                double qp[MDP1], double* szr, double* szi, double K[MDP1],
                double qk[MDP1]);
void Quad_ak1(double a, double b1, double c, double* sr, double* si, double* lr,
              double* li);

void rpoly_ak1(double op[MDP1], int* Degree, double zeror[MAXDEGREE],
               double zeroi[MAXDEGREE]) {
  int i, j, jj, l, N, NM1, NN, NZ, zerok;

  double K[MDP1], p[MDP1], pt[MDP1], qp[MDP1], temp[MDP1];
  double bnd, df, dx, factor, ff, moduli_max, moduli_min, sc, x, xm;
  double aa, bb, cc, lzi, lzr, sr, szi, szr, t, xx, xxx, yy;

  const double RADFAC = 3.14159265358979323846 /
                        180;  // Degrees-to-radians conversion factor = pi/180
  const double lb2 = log(
      2.0);  // Dummy variable to avoid re-calculating this value in loop below
  const double lo = FLT_MIN / DBL_EPSILON;
  const double cosr = cos(94.0 * RADFAC);  // = -0.069756474
  const double sinr = sin(94.0 * RADFAC);  // = 0.99756405

  if ((*Degree) > MAXDEGREE) {
    cout << "\nThe entered Degree is greater than MAXDEGREE. Exiting rpoly. No "
            "further action taken.\n";
    *Degree = -1;
    return;
  }  // End ((*Degree) > MAXDEGREE)

  // Do a quick check to see if leading coefficient is 0
  if (op[0] != 0) {
    N = *Degree;
    xx = sqrt(0.5);  // = 0.70710678
    yy = -xx;

    // Remove zeros at the origin, if any
    j = 0;
    while (op[N] == 0) {
      zeror[j] = zeroi[j] = 0.0;
      N--;
      j++;
    }  // End while (op[N] == 0)

    NN = N + 1;

    // Make a copy of the coefficients
    for (i = 0; i < NN; i++) p[i] = op[i];

    while (N >= 1) {  // Main loop
      // Start the algorithm for one zero
      if (N <= 2) {
        // Calculate the final zero or pair of zeros
        if (N < 2) {
          zeror[(*Degree) - 1] = -(p[1] / p[0]);
          zeroi[(*Degree) - 1] = 0.0;
        }       // End if (N < 2)
        else {  // else N == 2
          Quad_ak1(p[0], p[1], p[2], &zeror[(*Degree) - 2],
                   &zeroi[(*Degree) - 2], &zeror[(*Degree) - 1],
                   &zeroi[(*Degree) - 1]);
        }  // End else N == 2
        break;
      }  // End if (N <= 2)

      // Find the largest and smallest moduli of the coefficients

      moduli_max = 0.0;
      moduli_min = FLT_MAX;

      for (i = 0; i < NN; i++) {
        x = fabs(p[i]);
        if (x > moduli_max) moduli_max = x;
        if ((x != 0) && (x < moduli_min)) moduli_min = x;
      }  // End for i

      // Scale if there are large or very small coefficients
      // Computes a scale factor to multiply the coefficients of the polynomial.
      // The scaling
      // is done to avoid overflow and to avoid undetected underflow interfering
      // with the
      // convergence criterion.
      // The factor is a power of the base.

      sc = lo / moduli_min;

      if (((sc <= 1.0) && (moduli_max >= 10)) ||
          ((sc > 1.0) && (FLT_MAX / sc >= moduli_max))) {
        sc = ((sc == 0) ? FLT_MIN : sc);
        l = (int)(log(sc) / lb2 + 0.5);
        factor = pow(2.0, l);
        if (factor != 1.0) {
          for (i = 0; i < NN; i++) p[i] *= factor;
        }  // End if (factor != 1.0)
      }    // End if (((sc <= 1.0) && (moduli_max >= 10)) || ((sc > 1.0) &&
           // (FLT_MAX/sc >= moduli_max)))

      // Compute lower bound on moduli of zeros

      for (i = 0; i < NN; i++) pt[i] = fabs(p[i]);
      pt[N] = -(pt[N]);

      NM1 = N - 1;

      // Compute upper estimate of bound

      x = exp((log(-pt[N]) - log(pt[0])) / (double)N);

      if (pt[NM1] != 0) {
        // If Newton step at the origin is better, use it
        xm = -pt[N] / pt[NM1];
        x = ((xm < x) ? xm : x);
      }  // End if (pt[NM1] != 0)

      // Chop the interval (0, x) until ff <= 0

      xm = x;
      do {
        x = xm;
        xm = 0.1 * x;
        ff = pt[0];
        for (i = 1; i < NN; i++) ff = ff * xm + pt[i];
      } while (ff > 0);  // End do-while loop

      dx = x;

      // Do Newton iteration until x converges to two decimal places

      while (fabs(dx / x) > 0.005) {
        df = ff = pt[0];
        for (i = 1; i < N; i++) {
          ff = x * ff + pt[i];
          df = x * df + ff;
        }  // End for i
        ff = x * ff + pt[N];
        dx = ff / df;
        x -= dx;
      }  // End while loop

      bnd = x;

      // Compute the derivative as the initial K polynomial and do 5 steps with
      // no shift

      for (i = 1; i < N; i++) K[i] = (double)(N - i) * p[i] / ((double)N);
      K[0] = p[0];

      aa = p[N];
      bb = p[NM1];
      zerok = ((K[NM1] == 0) ? 1 : 0);

      for (jj = 0; jj < 5; jj++) {
        cc = K[NM1];
        if (zerok) {
          // Use unscaled form of recurrence
          for (i = 0; i < NM1; i++) {
            j = NM1 - i;
            K[j] = K[j - 1];
          }  // End for i
          K[0] = 0;
          zerok = ((K[NM1] == 0) ? 1 : 0);
        }  // End if (zerok)

        else {  // else !zerok
          // Used scaled form of recurrence if value of K at 0 is nonzero
          t = -aa / cc;
          for (i = 0; i < NM1; i++) {
            j = NM1 - i;
            K[j] = t * K[j - 1] + p[j];
          }  // End for i
          K[0] = p[0];
          zerok = ((fabs(K[NM1]) <= fabs(bb) * DBL_EPSILON * 10.0) ? 1 : 0);
        }  // End else !zerok

      }  // End for jj

      // Save K for restarts with new shifts
      for (i = 0; i < N; i++) temp[i] = K[i];

      // Loop to select the quadratic corresponding to each new shift

      for (jj = 1; jj <= 20; jj++) {
        // Quadratic corresponds to a double shift to a non-real point and its
        // complex conjugate. The point has modulus BND and amplitude rotated
        // by 94 degrees from the previous shift.

        xxx = -(sinr * yy) + cosr * xx;
        yy = sinr * xx + cosr * yy;
        xx = xxx;
        sr = bnd * xx;

        // Second stage calculation, fixed quadratic

        Fxshfr_ak1(20 * jj, &NZ, sr, bnd, K, N, p, NN, qp, &lzi, &lzr, &szi,
                   &szr);

        if (NZ != 0) {
          // The second stage jumps directly to one of the third stage
          // iterations and
          // returns here if successful. Deflate the polynomial, store the zero
          // or
          // zeros, and return to the main algorithm.

          j = (*Degree) - N;
          zeror[j] = szr;
          zeroi[j] = szi;
          NN = NN - NZ;
          N = NN - 1;
          for (i = 0; i < NN; i++) p[i] = qp[i];
          if (NZ != 1) {
            zeror[j + 1] = lzr;
            zeroi[j + 1] = lzi;
          }  // End if (NZ != 1)
          break;
        }       // End if (NZ != 0)
        else {  // Else (NZ == 0)

          // If the iteration is unsuccessful, another quadratic is chosen after
          // restoring K
          for (i = 0; i < N; i++) K[i] = temp[i];
        }  // End else (NZ == 0)

      }  // End for jj

      // Return with failure if no convergence with 20 shifts

      if (jj > 20) {
        cout << "\nFailure. No convergence after 20 shifts. Program "
                "terminated.\n";
        *Degree -= N;
        break;
      }  // End if (jj > 20)

    }  // End while (N >= 1)

  }       // End if op[0] != 0
  else {  // else op[0] == 0
    cout << "\nThe leading coefficient is zero. No further action taken. "
            "Program terminated.\n";
    *Degree = 0;
  }  // End else op[0] == 0

  return;
}  // End rpoly_ak1

void Fxshfr_ak1(int L2, int* NZ, double sr, double bnd, double K[MDP1], int N,
                double p[MDP1], int NN, double qp[MDP1], double* lzi,
                double* lzr, double* szi, double* szr) {
  // Computes up to L2 fixed shift K-polynomials, testing for convergence in the
  // linear or
  // quadratic case. Initiates one of the variable shift iterations and returns
  // with the
  // number of zeros found.

  // L2 limit of fixed shift steps
  // NZ number of zeros found

  int fflag, i, iFlag, j, spass, stry, tFlag, vpass, vtry;
  double a, a1, a3, a7, b, betas, betav, c, d, e, f, g, h, oss, ots, otv, ovv,
      s, ss, ts, tss, tv, tvv, u, ui, v, vi, vv;
  double qk[MDP1], svk[MDP1];

  *NZ = 0;
  betav = betas = 0.25;
  u = -(2.0 * sr);
  oss = sr;
  ovv = v = bnd;

  // Evaluate polynomial by synthetic division
  QuadSD_ak1(NN, u, v, p, qp, &a, &b);

  tFlag =
      calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u, v, qk);

  for (j = 0; j < L2; j++) {
    // Calculate next K polynomial and estimate v
    nextK_ak1(N, tFlag, a, b, a1, &a3, &a7, K, qk, qp);
    tFlag =
        calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u, v, qk);
    newest_ak1(tFlag, &ui, &vi, a, a1, a3, a7, b, c, d, f, g, h, u, v, K, N, p);

    vv = vi;

    // Estimate s

    ss = ((K[N - 1] != 0.0) ? -(p[N] / K[N - 1]) : 0.0);

    ts = tv = 1.0;

    if ((j != 0) && (tFlag != 3)) {
      // Compute relative measures of convergence of s and v sequences

      tv = ((vv != 0.0) ? fabs((vv - ovv) / vv) : tv);
      ts = ((ss != 0.0) ? fabs((ss - oss) / ss) : ts);

      // If decreasing, multiply the two most recent convergence measures

      tvv = ((tv < otv) ? tv * otv : 1.0);
      tss = ((ts < ots) ? ts * ots : 1.0);

      // Compare with convergence criteria

      vpass = ((tvv < betav) ? 1 : 0);
      spass = ((tss < betas) ? 1 : 0);

      if ((spass) || (vpass)) {
        // At least one sequence has passed the convergence test.
        // Store variables before iterating

        for (i = 0; i < N; i++) svk[i] = K[i];

        s = ss;

        // Choose iteration according to the fastest converging sequence

        stry = vtry = 0;
        fflag = 1;

        do {
          iFlag = 1;  // Begin each loop by assuming RealIT will be called
                      // UNLESS iFlag changed below

          if ((fflag && ((fflag = 0) == 0)) &&
              ((spass) && (!vpass || (tss < tvv)))) {
            ;  // Do nothing. Provides a quick "short circuit".
          }    // End if (fflag)

          else {  // else !fflag
            QuadIT_ak1(N, NZ, ui, vi, szr, szi, lzr, lzi, qp, NN, &a, &b, p, qk,
                       &a1, &a3, &a7, &d, &e, &f, &g, &h, K);

            if ((*NZ) > 0) return;

            // Quadratic iteration has failed. Flag that it has been tried and
            // decrease the
            // convergence criterion

            vtry = 1;
            betav *= 0.25;

            // Try linear iteration if it has not been tried and the s sequence
            // is converging
            if (stry || (!spass)) {
              iFlag = 0;
            }  // End if (stry || (!spass))
            else {
              for (i = 0; i < N; i++) K[i] = svk[i];
            }  // End if (stry || !spass)

          }  // End else !fflag

          if (iFlag != 0) {
            RealIT_ak1(&iFlag, NZ, &s, N, p, NN, qp, szr, szi, K, qk);

            if ((*NZ) > 0) return;

            // Linear iteration has failed. Flag that it has been tried and
            // decrease the
            // convergence criterion

            stry = 1;
            betas *= 0.25;

            if (iFlag != 0) {
              // If linear iteration signals an almost double real zero, attempt
              // quadratic iteration

              ui = -(s + s);
              vi = s * s;
              continue;

            }  // End if (iFlag != 0)
          }    // End if (iFlag != 0)

          // Restore variables
          for (i = 0; i < N; i++) K[i] = svk[i];

          // Try quadratic iteration if it has not been tried and the v sequence
          // is converging

        } while (vpass && !vtry);  // End do-while loop

        // Re-compute qp and scalar values to continue the second stage

        QuadSD_ak1(NN, u, v, p, qp, &a, &b);
        tFlag = calcSC_ak1(N, a, b, &a1, &a3, &a7, &c, &d, &e, &f, &g, &h, K, u,
                           v, qk);

      }  // End if ((spass) || (vpass))

    }  // End if ((j != 0) && (tFlag != 3))

    ovv = vv;
    oss = ss;
    otv = tv;
    ots = ts;
  }  // End for j

  return;
}  // End Fxshfr_ak1

void QuadSD_ak1(int NN, double u, double v, double p[MDP1], double q[MDP1],
                double* a, double* b) {
  // Divides p by the quadratic 1, u, v placing the quotient in q and the
  // remainder in a, b

  int i;

  q[0] = *b = p[0];
  q[1] = *a = -((*b) * u) + p[1];

  for (i = 2; i < NN; i++) {
    q[i] = -((*a) * u + (*b) * v) + p[i];
    *b = (*a);
    *a = q[i];
  }  // End for i

  return;
}  // End QuadSD_ak1

int calcSC_ak1(int N, double a, double b, double* a1, double* a3, double* a7,
               double* c, double* d, double* e, double* f, double* g, double* h,
               double K[MDP1], double u, double v, double qk[MDP1]) {
  // This routine calculates scalar quantities used to compute the next K
  // polynomial and
  // new estimates of the quadratic coefficients.

  // calcSC - integer variable set here indicating how the calculations are
  // normalized
  //  to avoid overflow.

  int dumFlag = 3;  // TYPE = 3 indicates the quadratic is almost a factor of K

  // Synthetic division of K by the quadratic 1, u, v
  QuadSD_ak1(N, u, v, K, qk, c, d);

  if (fabs((*c)) <= (10.0 * DBL_EPSILON * fabs(K[N - 1]))) {
    if (fabs((*d)) <= (10.0 * DBL_EPSILON * fabs(K[N - 2]))) return dumFlag;
  }  // End if (fabs(c) <= (100.0*DBL_EPSILON*fabs(K[N - 1])))

  *h = v * b;
  if (fabs((*d)) >= fabs((*c))) {
    dumFlag = 2;  // TYPE = 2 indicates that all formulas are divided by d
    *e = a / (*d);
    *f = (*c) / (*d);
    *g = u * b;
    *a3 = (*e) * ((*g) + a) + (*h) * (b / (*d));
    *a1 = -a + (*f) * b;
    *a7 = (*h) + ((*f) + u) * a;
  }  // End if(fabs(d) >= fabs(c))
  else {
    dumFlag = 1;  // TYPE = 1 indicates that all formulas are divided by c;
    *e = a / (*c);
    *f = (*d) / (*c);
    *g = (*e) * u;
    *a3 = (*e) * a + ((*g) + (*h) / (*c)) * b;
    *a1 = -(a * ((*d) / (*c))) + b;
    *a7 = (*g) * (*d) + (*h) * (*f) + a;
  }  // End else

  return dumFlag;
}  // End calcSC_ak1

void nextK_ak1(int N, int tFlag, double a, double b, double a1, double* a3,
               double* a7, double K[MDP1], double qk[MDP1], double qp[MDP1]) {
  // Computes the next K polynomials using the scalars computed in calcSC_ak1

  int i;
  double temp;

  if (tFlag == 3) {  // Use unscaled form of the recurrence
    K[1] = K[0] = 0.0;

    for (i = 2; i < N; i++) K[i] = qk[i - 2];

    return;
  }  // End if (tFlag == 3)

  temp = ((tFlag == 1) ? b : a);

  if (fabs(a1) > (10.0 * DBL_EPSILON * fabs(temp))) {
    // Use scaled form of the recurrence

    (*a7) /= a1;
    (*a3) /= a1;
    K[0] = qp[0];
    K[1] = -((*a7) * qp[0]) + qp[1];

    for (i = 2; i < N; i++)
      K[i] = -((*a7) * qp[i - 1]) + (*a3) * qk[i - 2] + qp[i];

  }  // End if (fabs(a1) > (10.0*DBL_EPSILON*fabs(temp)))
  else {
    // If a1 is nearly zero, then use a special form of the recurrence

    K[0] = 0.0;
    K[1] = -(*a7) * qp[0];

    for (i = 2; i < N; i++) K[i] = -((*a7) * qp[i - 1]) + (*a3) * qk[i - 2];
  }  // End else

  return;

}  // End nextK_ak1

void newest_ak1(int tFlag, double* uu, double* vv, double a, double a1,
                double a3, double a7, double b, double c, double d, double f,
                double g, double h, double u, double v, double K[MDP1], int N,
                double p[MDP1]) {
  // Compute new estimates of the quadratic coefficients using the scalars
  // computed in calcSC_ak1

  double a4, a5, b1, b2, c1, c2, c3, c4, temp;

  (*vv) = (*uu) = 0.0;  // The quadratic is zeroed

  if (tFlag != 3) {
    if (tFlag != 2) {
      a4 = a + u * b + h * f;
      a5 = c + (u + v * f) * d;
    }       // End if (tFlag != 2)
    else {  // else tFlag == 2
      a4 = (a + g) * f + h;
      a5 = (f + u) * c + v * d;
    }  // End else tFlag == 2

    // Evaluate new quadratic coefficients

    b1 = -K[N - 1] / p[N];
    b2 = -(K[N - 2] + b1 * p[N - 1]) / p[N];
    c1 = v * b2 * a1;
    c2 = b1 * a7;
    c3 = b1 * b1 * a3;
    c4 = -(c2 + c3) + c1;
    temp = -c4 + a5 + b1 * a4;
    if (temp != 0.0) {
      *uu = -((u * (c3 + c2) + v * (b1 * a1 + b2 * a7)) / temp) + u;
      *vv = v * (1.0 + c4 / temp);
    }  // End if (temp != 0)

  }  // End if (tFlag != 3)

  return;
}  // End newest_ak1

void QuadIT_ak1(int N, int* NZ, double uu, double vv, double* szr, double* szi,
                double* lzr, double* lzi, double qp[MDP1], int NN, double* a,
                double* b, double p[MDP1], double qk[MDP1], double* a1,
                double* a3, double* a7, double* d, double* e, double* f,
                double* g, double* h, double K[MDP1]) {
  // Variable-shift K-polynomial iteration for a quadratic factor converges only
  // if the
  // zeros are equimodular or nearly so.

  int i, j = 0, tFlag, triedFlag = 0;
  double c, ee, mp, omp, relstp, t, u, ui, v, vi, zm;

  *NZ = 0;  // Number of zeros found
  u = uu;   // uu and vv are coefficients of the starting quadratic
  v = vv;

  do {
    Quad_ak1(1.0, u, v, szr, szi, lzr, lzi);

    // Return if roots of the quadratic are real and not close to multiple or
    // nearly
    // equal and of opposite sign.

    if (fabs(fabs(*szr) - fabs(*lzr)) > 0.01 * fabs(*lzr)) break;

    // Evaluate polynomial by quadratic synthetic division

    QuadSD_ak1(NN, u, v, p, qp, a, b);

    mp = fabs(-((*szr) * (*b)) + (*a)) + fabs((*szi) * (*b));

    // Compute a rigorous bound on the rounding error in evaluating p

    zm = sqrt(fabs(v));
    ee = 2.0 * fabs(qp[0]);
    t = -((*szr) * (*b));

    for (i = 1; i < N; i++) ee = ee * zm + fabs(qp[i]);

    ee = ee * zm + fabs((*a) + t);
    ee = (9.0 * ee + 2.0 * fabs(t) - 7.0 * (fabs((*a) + t) + zm * fabs((*b)))) *
         DBL_EPSILON;

    // Iteration has converged sufficiently if the polynomial value is less than
    // 20 times this bound

    if (mp <= 20.0 * ee) {
      *NZ = 2;
      break;
    }  // End if (mp <= 20.0*ee)

    j++;

    // Stop iteration after 20 steps
    if (j > 20) break;

    if (j >= 2) {
      if ((relstp <= 0.01) && (mp >= omp) && (!triedFlag)) {
        // A cluster appears to be stalling the convergence. Five fixed shift
        // steps are taken with a u, v close to the cluster.

        relstp = ((relstp < DBL_EPSILON) ? sqrt(DBL_EPSILON) : sqrt(relstp));

        u -= u * relstp;
        v += v * relstp;

        QuadSD_ak1(NN, u, v, p, qp, a, b);

        for (i = 0; i < 5; i++) {
          tFlag =
              calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
          nextK_ak1(N, tFlag, *a, *b, *a1, a3, a7, K, qk, qp);
        }  // End for i

        triedFlag = 1;
        j = 0;

      }  // End if ((relstp <= 0.01) && (mp >= omp) && (!triedFlag))

    }  // End if (j >= 2)

    omp = mp;

    // Calculate next K polynomial and new u and v

    tFlag = calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
    nextK_ak1(N, tFlag, *a, *b, *a1, a3, a7, K, qk, qp);
    tFlag = calcSC_ak1(N, *a, *b, a1, a3, a7, &c, d, e, f, g, h, K, u, v, qk);
    newest_ak1(tFlag, &ui, &vi, *a, *a1, *a3, *a7, *b, c, *d, *f, *g, *h, u, v,
               K, N, p);

    // If vi is zero, the iteration is not converging
    if (vi != 0) {
      relstp = fabs((-v + vi) / vi);
      u = ui;
      v = vi;
    }                 // End if (vi != 0)
  } while (vi != 0);  // End do-while loop

  return;

}  // End QuadIT_ak1

void RealIT_ak1(int* iFlag, int* NZ, double* sss, int N, double p[MDP1], int NN,
                double qp[MDP1], double* szr, double* szi, double K[MDP1],
                double qk[MDP1]) {
  // Variable-shift H-polynomial iteration for a real zero

  // sss  - starting iterate
  // NZ - number of zeros found
  // iFlag  - flag to indicate a pair of zeros near real axis

  int i, j = 0, nm1 = N - 1;
  double ee, kv, mp, ms, omp, pv, s, t;

  *iFlag = *NZ = 0;
  s = *sss;

  for (;;) {
    qp[0] = pv = p[0];

    // Evaluate p at s
    for (i = 1; i < NN; i++) qp[i] = pv = pv * s + p[i];

    mp = fabs(pv);

    // Compute a rigorous bound on the error in evaluating p

    ms = fabs(s);
    ee = 0.5 * fabs(qp[0]);
    for (i = 1; i < NN; i++) ee = ee * ms + fabs(qp[i]);

    // Iteration has converged sufficiently if the polynomial value is less than
    // 20 times this bound

    if (mp <= 20.0 * DBL_EPSILON * (2.0 * ee - mp)) {
      *NZ = 1;
      *szr = s;
      *szi = 0.0;
      break;
    }  // End if (mp <= 20.0*DBL_EPSILON*(2.0*ee - mp))

    j++;

    // Stop iteration after 10 steps

    if (j > 10) break;

    if (j >= 2) {
      if ((fabs(t) <= 0.001 * fabs(-t + s)) && (mp > omp)) {
        // A cluster of zeros near the real axis has been encountered;
        // Return with iFlag set to initiate a quadratic iteration

        *iFlag = 1;
        *sss = s;
        break;
      }  // End if ((fabs(t) <= 0.001*fabs(s - t)) && (mp > omp))

    }  // End if (j >= 2)

    // Return if the polynomial value has increased significantly

    omp = mp;

    // Compute t, the next polynomial and the new iterate
    qk[0] = kv = K[0];
    for (i = 1; i < N; i++) qk[i] = kv = kv * s + K[i];

    if (fabs(kv) > fabs(K[nm1]) * 10.0 * DBL_EPSILON) {
      // Use the scaled form of the recurrence if the value of K at s is
      // non-zero
      t = -(pv / kv);
      K[0] = qp[0];
      for (i = 1; i < N; i++) K[i] = t * qk[i - 1] + qp[i];
    }       // End if (fabs(kv) > fabs(K[nm1])*10.0*DBL_EPSILON)
    else {  // else (fabs(kv) <= fabs(K[nm1])*10.0*DBL_EPSILON)
      // Use unscaled form
      K[0] = 0.0;
      for (i = 1; i < N; i++) K[i] = qk[i - 1];
    }  // End else (fabs(kv) <= fabs(K[nm1])*10.0*DBL_EPSILON)

    kv = K[0];
    for (i = 1; i < N; i++) kv = kv * s + K[i];

    t = ((fabs(kv) > (fabs(K[nm1]) * 10.0 * DBL_EPSILON)) ? -(pv / kv) : 0.0);

    s += t;

  }  // End infinite for loop

  return;

}  // End RealIT_ak1

void Quad_ak1(double a, double b1, double c, double* sr, double* si, double* lr,
              double* li) {
  // Calculates the zeros of the quadratic a*Z^2 + b1*Z + c
  // The quadratic formula, modified to avoid overflow, is used to find the
  // larger zero if the
  // zeros are real and both zeros are complex. The smaller real zero is found
  // directly from
  // the product of the zeros c/a.

  double b, d, e;

  *sr = *si = *lr = *li = 0.0;

  if (a == 0) {
    *sr = ((b1 != 0) ? -(c / b1) : *sr);
    return;
  }  // End if (a == 0))

  if (c == 0) {
    *lr = -(b1 / a);
    return;
  }  // End if (c == 0)

  // Compute discriminant avoiding overflow

  b = b1 / 2.0;
  if (fabs(b) < fabs(c)) {
    e = ((c >= 0) ? a : -a);
    e = -e + b * (b / fabs(c));
    d = sqrt(fabs(e)) * sqrt(fabs(c));
  }       // End if (fabs(b) < fabs(c))
  else {  // Else (fabs(b) >= fabs(c))
    e = -((a / b) * (c / b)) + 1.0;
    d = sqrt(fabs(e)) * (fabs(b));
  }  // End else (fabs(b) >= fabs(c))

  if (e >= 0) {
    // Real zeros

    d = ((b >= 0) ? -d : d);
    *lr = (-b + d) / a;
    *sr = ((*lr != 0) ? (c / (*lr)) / a : *sr);
  }       // End if (e >= 0)
  else {  // Else (e < 0)
    // Complex conjugate zeros

    *lr = *sr = -(b / a);
    *si = fabs(d / a);
    *li = -(*si);
  }  // End else (e < 0)

  return;
}  // End Quad_ak1


}  // namespace rpoly_impl;

void rpolyWrapper(double* coefficients_decreasing, int* degree,
                 double* roots_real, double* roots_imag) {
  rpoly_impl::rpoly_ak1(coefficients_decreasing, degree, roots_real,
                        roots_imag);
}

}  // namespace mav_trajectory_generation
