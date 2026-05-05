#include <math.h>
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_acvf_scalar_ar1(void) {
  double A[] = {0.5};
  double Sig[] = {3};
  double Gamma[4];
  double expected[] = {4, 2, 1, 0.5};
  bool ok = varmapack_acvf(A, 0, Sig, 1, 0, 1, Gamma, 3);
  xCheck(ok);
  for (int i=0; i<4; i++) {
    xCheck(fabs(Gamma[i] - expected[i]) < 1e-12);
  }
}

static void check_acvf_scalar_ma1(void) {
  double B[] = {0.5};
  double Sig[] = {2};
  double Gamma[3];
  double expected[] = {2.5, 1, 0};
  bool ok = varmapack_acvf(0, B, Sig, 0, 1, 1, Gamma, 2);
  xCheck(ok);
  for (int i=0; i<3; i++) {
    xCheck(fabs(Gamma[i] - expected[i]) < 1e-12);
  }
}

static void check_acvf_scalar_white_noise(void) {
  double Sig[] = {1.25};
  double Gamma[4];
  double expected[] = {1.25, 0, 0, 0};
  bool ok = varmapack_acvf(0, 0, Sig, 0, 0, 1, Gamma, 3);
  xCheck(ok);
  for (int i=0; i<4; i++) {
    xCheck(fabs(Gamma[i] - expected[i]) < 1e-12);
  }
}

static void check_acvf_scalar_arma11(void) {
  double A[] = {0.5};
  double B[] = {0.25};
  double Sig[] = {2};
  double Gamma[4];
  double expected[] = {3.5, 2.25, 1.125, 0.5625};
  bool ok = varmapack_acvf(A, B, Sig, 1, 1, 1, Gamma, 3);
  xCheck(ok);
  for (int i=0; i<4; i++) {
    xCheck(fabs(Gamma[i] - expected[i]) < 1e-12);
  }
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
  bool ok = varmapack_acvf(A, B, Sig, 1, 2, 1, Gamma, 4);
  xCheck(ok);
  for (int i=0; i<5; i++) {
    xCheck(fabs(Gamma[i] - expected[i]) < 1e-12);
  }
}

void TestAcvf(void) {
  check_acvf_scalar_ar1();
  check_acvf_scalar_ma1();
  check_acvf_scalar_white_noise();
  check_acvf_scalar_arma11();
  check_acvf_scalar_arma12();
}
