#include "TestUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

static int nearly_one(double x, double tol) { return fabs(x - 1.0) <= tol; }

static void checkP0(void) {
  double rho = varmapack_specrad(0, 4, 0);
  checkVarmapackClean();
  xCheck(almostSame(rho, 0));
  rho = varmapack_ma_specrad(0, 4, 0);
  checkVarmapackClean();
  xCheck(almostSame(rho, 0));
}

static void checkInvalidInput(void) {
  double A[] = {0};
  checkVarmapackNaN(varmapack_specrad(0, 1, 1));
  checkVarmapackNaN(varmapack_specrad(A, 0, 1));
  checkVarmapackNaN(varmapack_specrad(A, 1, -1));
  checkVarmapackNaN(varmapack_specrad(A, 0, 0));
  checkVarmapackNaN(varmapack_ma_specrad(0, 1, 1));
  checkVarmapackNaN(varmapack_ma_specrad(A, 0, 1));
  checkVarmapackNaN(varmapack_ma_specrad(A, 1, -1));
  checkVarmapackNaN(varmapack_ma_specrad(A, 0, 0));
}

static void checkNonfiniteInput(void) {
  double bad[] = {NAN, INFINITY, -INFINITY};
  for (int i=0; i<LEN(bad); i++) {
    checkVarmapackNaN(varmapack_specrad(bad + i, 1, 1));
    checkVarmapackNaN(varmapack_ma_specrad(bad + i, 1, 1));
  }
}

static void checkOverflowingDimensions(void) {
  double A[] = {0};
  checkVarmapackNaN(varmapack_specrad(A, INT_MAX, 2));
  checkVarmapackNaN(varmapack_ma_specrad(A, INT_MAX, 2));
}

static void checkScalarAR2ComplexRoots(void) {
  double A[] = {1, -0.5};
  double rho = varmapack_specrad(A, 1, 2);
  checkVarmapackClean();
  xCheck(fabs(rho - sqrt(0.5)) < 1e-14);
}

static void checkScalarUnstable(void) {
  double A[] = {1.25};
  double rho = varmapack_specrad(A, 1, 1);
  checkVarmapackClean();
  xCheck(almostSame(rho, 1.25));
  xCheck(rho > 1);
}

static void checkScalarMA(void) {
  double B1[] = {0.7};
  double B2[] = {-1.2};
  double Bma2[] = {0.4, -0.12};
  double rho = varmapack_ma_specrad(B1, 1, 1);
  checkVarmapackClean();
  xCheck(almostSame(rho, 0.7));
  rho = varmapack_ma_specrad(B2, 1, 1);
  checkVarmapackClean();
  xCheck(almostSame(rho, 1.2));
  rho = varmapack_ma_specrad(Bma2, 1, 2);
  checkVarmapackClean();
  xCheck(fabs(rho - 0.6) < 1e-14);
}

static void checkDiagonalMA(void) {
  int r = 2, q = 1;
  double B[] = {
    0.3, 0,
    0, -0.8
  };
  double rho = varmapack_ma_specrad(B, r, q);
  checkVarmapackClean();
  xCheck(almostSame(rho, 0.8));
}

static void check2x2(void) {
  int r = 2, p = 1;
  double A[] = {
    0.95, 0.40, 0.10, 0.20
  };
  double rho = varmapack_specrad(A, r, p);
  checkVarmapackClean();
  xCheck(nearly_one(rho, 1e-14));
  xCheck(rho <= 1 + 1e-14);
}

static void checkScaling(void) {
  int r = 2, p = 1;
  double A[] = {
    0.95 * 20, 0.40 * 20, 0.10 * 20, 0.20 * 20
  };
  double rho = varmapack_specrad(A, r, p);
  checkVarmapackClean();
  xCheck(almostSame(rho, 20));
  xCheck(rho > 19.999 && rho < 20.001);
}

static void checkZero(void) {
  int r = 3, p = 1;
  double A[9] = {0};
  double rho = varmapack_specrad(A, r, p);
  checkVarmapackClean();
  xCheck(almostSame(rho, 0));
}

static void checkAllNamedStationary(void) {
  double strict_tol = 1e-12;
  int pmax = 0, qmax = 0, rmax = 0, maxcase = 0;
  bool ok = varmapack_testcase("max", &maxcase, 0, &pmax, &qmax,
                               &rmax, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  for (int icase=1; icase<=maxcase; icase++) {
    test_model model;
    loadIndexedTestModel(&model, icase);
    xCheck(model.r > 0);
    xCheck(model.p >= 0);
    xCheck(model.q >= 0);
    double rho = model.p == 0 ? 0 :
      varmapack_specrad(model.A, model.r, model.p);
    double maRho = model.q == 0 ? 0 :
      varmapack_ma_specrad(model.B, model.r, model.q);
    checkVarmapackClean();
    if (!(rho < 1 - strict_tol)) {
      fprintf(stderr, "[Testvarmapack_specrad] Nonstationary/borderline: "
              "case %d (%s), rho=%.15g (p=%d, q=%d, r=%d)\n", icase,
              model.name, rho, model.p, model.q, model.r);
    }
    xCheck(rho < 1 - strict_tol);
    xCheck(isfinite(maRho) && maRho >= 0);
    freeTestModel(&model);
  }
}

void Testvarmapack_specrad(void) {
  checkP0();
  checkInvalidInput();
  checkNonfiniteInput();
  checkOverflowingDimensions();
  checkScalarAR2ComplexRoots();
  checkScalarUnstable();
  checkScalarMA();
  checkDiagonalMA();
  check2x2();
  checkScaling();
  checkZero();
  checkAllNamedStationary();
}
