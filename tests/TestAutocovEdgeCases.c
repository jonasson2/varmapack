#include <math.h>
#include "ExtraUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_scalar_maxlag(void) {
  double X[] = {1, 2, 4};
  double C[3];
  double ml[] = {14.0/9.0, -1.0/27.0, -20.0/27.0};
  double corrected[] = {14.0/9.0, -1.0/18.0, -20.0/9.0};
  varmapack_error error;
  error = varmapack_autocov("N", "ML", 1, 3, X, 2, C);
  xCheck(!error);
  checkArrayTol(C, ml, 3, 1e-12);
  error = varmapack_autocov("N", "C", 1, 3, X, 2, C);
  xCheck(!error);
  checkArrayTol(C, corrected, 3, 1e-12);
}

static void check_constant_series(void) {
  double X[] = {
    3, -1, 3, -1, 3, -1, 3, -1
  };
  double C[2*2*4];
  varmapack_error error;
  error = varmapack_autocov("N", "ML", 2, 4, X, 3, C);
  xCheck(!error);
  checkArrayZero(C, 16);
  error = varmapack_autocov("N", "C", 2, 4, X, 3, C);
  xCheck(!error);
  checkArrayZero(C, 16);
}

static void check_maxlag_zero(void) {
  double X[] = {
    1, 2, 3, 4, 5, 8
  };
  double C[4];
  double expected[] = {
    8.0/3.0, 4, 4, 56.0/9.0
  };
  varmapack_error error;
  error = varmapack_autocov("N", "ML", 2, 3, X, 0, C);
  xCheck(!error);
  checkArrayTol(C, expected, 4, 1e-12);
  error = varmapack_autocov("N", "C", 2, 3, X, 0, C);
  xCheck(!error);
  checkArrayTol(C, expected, 4, 1e-12);
}

static void check_single_observation(void) {
  double X[] = {2, -3};
  double C[4];
  varmapack_error error;
  error = varmapack_autocov("N", "ML", 2, 1, X, 0, C);
  xCheck(!error);
  checkArrayZero(C, 4);
  error = varmapack_autocov("N", "C", 2, 1, X, 0, C);
  xCheck(!error);
  checkArrayZero(C, 4);
}

void TestAutocovEdgeCases(void) {
  check_scalar_maxlag();
  check_constant_series();
  check_maxlag_zero();
  check_single_observation();
}
