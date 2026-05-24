#include "ExtraUtil.h"
#include "error.h"
#include "xCheck.h"
#include "Tests.h"
#include "varmapack.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

static int nearly_one(double x, double tol) { return fabs(x - 1.0) <= tol; }

static void checkP0(void) {
  double rho = varmapack_specrad(0, 4, 0);
  xCheck(almostSame(rho, 0));
  rho = varmapack_ma_specrad(0, 4, 0);
  xCheck(almostSame(rho, 0));
}

static void checkInvalidInput(void) {
  double A[] = {0};
  xCheck(isnan(varmapack_specrad(0, 1, 1)));
  xCheck(isnan(varmapack_specrad(A, 0, 1)));
  xCheck(isnan(varmapack_specrad(A, 1, -1)));
  xCheck(isnan(varmapack_specrad(A, 0, 0)));
  xCheck(isnan(varmapack_ma_specrad(0, 1, 1)));
  xCheck(isnan(varmapack_ma_specrad(A, 0, 1)));
  xCheck(isnan(varmapack_ma_specrad(A, 1, -1)));
  xCheck(isnan(varmapack_ma_specrad(A, 0, 0)));
}

static void checkScalarAR2ComplexRoots(void) {
  double A[] = {1, -0.5};
  double rho = varmapack_specrad(A, 1, 2);
  xCheck(fabs(rho - sqrt(0.5)) < 1e-14);
}

static void checkScalarUnstable(void) {
  double A[] = {1.25};
  double rho = varmapack_specrad(A, 1, 1);
  xCheck(almostSame(rho, 1.25));
  xCheck(rho > 1);
}

static void checkScalarMA(void) {
  double B1[] = {0.7};
  double B2[] = {-1.2};
  double Bma2[] = {0.4, -0.12};
  double rho = varmapack_ma_specrad(B1, 1, 1);
  xCheck(almostSame(rho, 0.7));
  rho = varmapack_ma_specrad(B2, 1, 1);
  xCheck(almostSame(rho, 1.2));
  rho = varmapack_ma_specrad(Bma2, 1, 2);
  xCheck(fabs(rho - 0.6) < 1e-14);
}

static void checkDiagonalMA(void) {
  int r = 2, q = 1;
  double B[] = {
    0.3, 0,
    0, -0.8
  };
  double rho = varmapack_ma_specrad(B, r, q);
  xCheck(almostSame(rho, 0.8));
}

static void check2x2(void) {
  int r = 2, p = 1;
  double A[] = {
    0.95, 0.40, 0.10, 0.20
  };
  double rho = varmapack_specrad(A, r, p);
  xCheck(nearly_one(rho, 1e-14));
  xCheck(rho <= 1 + 1e-14);
}

static void checkScaling(void) {
  int r = 2, p = 1;
  double A[] = {
    0.95 * 20, 0.40 * 20, 0.10 * 20, 0.20 * 20
  };
  double rho = varmapack_specrad(A, r, p);
  xCheck(almostSame(rho, 20));
  xCheck(rho > 19.999 && rho < 20.001);
}

static void checkZero(void) {
  int r = 3, p = 1;
  double A[9] = {0};
  double rho = varmapack_specrad(A, r, p);
  xCheck(almostSame(rho, 0));
}

static void checkAllNamedStationary(void) {
  const double strict_tol = 1e-12;
  int pmax = 0, qmax = 0, rmax = 0, maxcase = 0;
  varmapack_error error = varmapack_testcase(0, 0, 0, "max", &pmax, &qmax,
                                               &rmax, &maxcase, 0, 0);
  xCheck(!error);
  for (int icase = 1; icase <= maxcase; ++icase) {
    int p = 0, q = 0, r = 0;
    char name[32] = {0};
    error = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, 0);
    xCheck(!error);
    xCheck(r > 0);
    xCheck(p >= 0);
    xCheck(q >= 0);
    size_t nA = (size_t)r * (size_t)r * (size_t)(p > 0 ? p : 1);
    size_t nB = (size_t)r * (size_t)r * (size_t)(q > 0 ? q : 1);
    size_t nS = (size_t)r * (size_t)r;
    double *A = 0, *B = 0, *Sig = 0;
    xCheck(ALLOC(A, nA));
    xCheck(ALLOC(B, nB));
    xCheck(ALLOC(Sig, nS));
    error = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, 0);
    xCheck(!error);
    double rho = p == 0 ? 0 : varmapack_specrad(A, r, p);
    double maRho = q == 0 ? 0 : varmapack_ma_specrad(B, r, q);
    if (!(rho < 1 - strict_tol)) {
      fprintf(stderr, "[Testvarmapack_specrad] Nonstationary/borderline: case %d (%s), "
              "rho=%.15g (p=%d, q=%d, r=%d)\n", icase, name, rho, p, q, r);
    }
    xCheck(rho < 1 - strict_tol);
    xCheck(isfinite(maRho) && maRho >= 0);
    FREE(A);
    FREE(B);
    FREE(Sig);
  }
}

void Testvarmapack_specrad(void) {
  checkP0();
  checkInvalidInput();
  checkScalarAR2ComplexRoots();
  checkScalarUnstable();
  checkScalarMA();
  checkDiagonalMA();
  check2x2();
  checkScaling();
  checkZero();
  checkAllNamedStationary();
}
