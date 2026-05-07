#include <math.h>
#include <stdio.h>
#include "ExtraUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_acvf_scalar_ar1(void) {
  double A[] = {0.5};
  double Sig[] = {3};
  double Gamma[4];
  double expected[] = {4, 2, 1, 0.5};
  varmapack_error error = varmapack_acvf(A, 0, Sig, 1, 0, 1, Gamma, 3);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 4, 1e-12);
}

static void check_acvf_scalar_ma1(void) {
  double B[] = {0.5};
  double Sig[] = {2};
  double Gamma[3];
  double expected[] = {2.5, 1, 0};
  varmapack_error error = varmapack_acvf(0, B, Sig, 0, 1, 1, Gamma, 2);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 3, 1e-12);
}

static void check_acvf_scalar_white_noise(void) {
  double Sig[] = {1.25};
  double Gamma[4];
  double expected[] = {1.25, 0, 0, 0};
  varmapack_error error = varmapack_acvf(0, 0, Sig, 0, 0, 1, Gamma, 3);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 4, 1e-12);
}

static void check_acvf_scalar_arma11(void) {
  double A[] = {0.5};
  double B[] = {0.25};
  double Sig[] = {2};
  double Gamma[4];
  double expected[] = {3.5, 2.25, 1.125, 0.5625};
  varmapack_error error = varmapack_acvf(A, B, Sig, 1, 1, 1, Gamma, 3);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 4, 1e-12);
}

static void check_acvf_scalar_arma12(void) {
  double A[] = {0.5};
  double B[] = {0.25, -0.1};
  double Sig[] = {2};
  double Gamma[5];
  double expected[] = {
    3.3266666666666667,
    2.013333333333333,
    0.8066666666666668,
    0.4033333333333334,
    0.2016666666666667
  };
  varmapack_error error = varmapack_acvf(A, B, Sig, 1, 2, 1, Gamma, 4);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 5, 1e-12);
}

static void check_acvf_vector_var1_diagonal(void) {
  double A[] = {0.5, 0, 0, 0.25};
  double Sig[] = {3, 0, 0, 4};
  double Gamma[4*4];
  double expected[] = {
    4, 0, 0, 4.266666666666667,
    2, 0, 0, 1.066666666666667,
    1, 0, 0, 0.2666666666666667,
    0.5, 0, 0, 0.06666666666666667
  };
  varmapack_error error = varmapack_acvf(A, 0, Sig, 1, 0, 2, Gamma, 3);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 16, 1e-12);
}

static void check_acvf_vector_ma1(void) {
  double B[] = {0.2, 0.3, 0.1, -0.1};
  double Sig[] = {2, 0.5, 0.5, 1};
  double Gamma[4*3];
  double expected[] = {
    2.11, 0.615, 0.615, 1.16,
    0.45, 0.55, 0.2, 0.05,
    0, 0, 0, 0
  };
  varmapack_error error = varmapack_acvf(0, B, Sig, 0, 1, 2, Gamma, 2);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 12, 1e-12);
}

static void check_acvf_maxlag_equals_p(void) {
  double A[] = {0.5, -0.1};
  double Sig[] = {2};
  double Gamma[2];
  double expected[] = {2.5462962962962967, 1.1574074074074074};
  varmapack_error error = varmapack_acvf(A, 0, Sig, 2, 0, 1, Gamma, 2);
  xCheck(!error);
  checkArrayTol(Gamma, expected, 2, 1e-12);
}

static void check_acvf_nonstationary_failure(void) {
  double A[] = {1};
  double Sig[] = {1};
  double Gamma[1];
  varmapack_error error = varmapack_acvf(A, 0, Sig, 1, 0, 1, Gamma, 1);
  xCheck(error == VARMAPACK_NONSTATIONARY);
}

static void check_acvf_invalid_input(void) {
  double A[] = {0.5};
  double B[] = {0.25};
  double Sig[] = {1};
  double Gamma[2];
  varmapack_error error;
  error = varmapack_acvf(0, 0, Sig, 1, 0, 1, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, 0, Sig, 1, 1, 1, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, 0, 1, 1, 1, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, Sig, 1, 1, 1, 0, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, Sig, -1, 1, 1, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, Sig, 1, -1, 1, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, Sig, 1, 1, 0, Gamma, 1);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_acvf(A, B, Sig, 1, 1, 1, Gamma, 0);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
}

void TestAcvf(void) {
  check_acvf_scalar_ar1();
  check_acvf_scalar_ma1();
  check_acvf_scalar_white_noise();
  check_acvf_scalar_arma11();
  check_acvf_scalar_arma12();
  check_acvf_vector_var1_diagonal();
  check_acvf_vector_ma1();
  check_acvf_maxlag_equals_p();
  check_acvf_nonstationary_failure();
  check_acvf_invalid_input();
}
