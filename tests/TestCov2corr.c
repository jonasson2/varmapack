#include <math.h>
#include "ExtraUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_copy(void) {
  double cov[] = {
    4, 2, 2, 9,
    2, 6, -3, 4.5,
    -4, 3, 12, -9
  };
  double original[] = {
    4, 2, 2, 9,
    2, 6, -3, 4.5,
    -4, 3, 12, -9
  };
  double expected[] = {
    1, 0.3333333333333333, 0.3333333333333333, 1,
    0.5, 1, -0.5, 0.5,
    -1, 0.5, 2, -1
  };
  double corr[12];
  bool ok = varmapack_cov2corr(cov, 2, 2, corr);
  checkVarmapackSuccess(ok);
  checkArrayTol(corr, expected, 12, 1e-15);
  checkArraySame(cov, original, 12);
}

static void check_in_place(void) {
  double cov[] = {
    4, 2, 2, 9,
    2, 6, -3, 4.5,
    -4, 3, 12, -9
  };
  double expected[] = {
    1, 0.3333333333333333, 0.3333333333333333, 1,
    0.5, 1, -0.5, 0.5,
    -1, 0.5, 2, -1
  };
  bool ok = varmapack_cov2corr(cov, 2, 2, cov);
  checkVarmapackSuccess(ok);
  checkArrayTol(cov, expected, 12, 1e-15);
}

static void check_corrected_autocovariance(void) {
  double X[] = {1, 2, 4};
  double cov[3], corr[3];
  double expected[] = {1, -0.035714285714285714, -1.4285714285714286};
  bool ok = varmapack_autocov("N", "C", 1, 3, X, 2, cov);
  checkVarmapackSuccess(ok);
  ok = varmapack_cov2corr(cov, 1, 2, corr);
  checkVarmapackSuccess(ok);
  checkArrayTol(corr, expected, 3, 1e-15);
  xCheck(corr[2] < -1);
}

static void check_invalid_input(void) {
  double cov[] = {1};
  double corr[] = {7};
  double zero[] = {0};
  double negative[] = {-1};
  double nonfinite[] = {NAN};
  bool ok;
  ok = varmapack_cov2corr(0, 1, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(cov, 1, 0, 0);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(cov, 0, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(cov, 1, -1, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(zero, 1, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(negative, 1, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(nonfinite, 1, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(cov, INT_MAX, 0, corr);
  checkVarmapackFailure(ok);
  ok = varmapack_cov2corr(cov, 1, INT_MAX, corr);
  checkVarmapackFailure(ok);
  xCheck(corr[0] == 7);
}

void TestCov2corr(void) {
  check_copy();
  check_in_place();
  check_corrected_autocovariance();
  check_invalid_input();
}
