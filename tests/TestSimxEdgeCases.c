#include <math.h>
#include "ExtraUtil.h"
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_white_noise_with_exog(void) {
  int p = 0, q = 0, s = 1, r = 1, n = 5, M = 2, h = 1;
  double C[] = {0.5};
  double Sig[] = {1};
  double z[] = {1, 2, 3, 4, 5};
  double X0[] = {2};
  double X[10], E[10];
  randompack_rng *rng = seededRng(100);
  varmapack_error error = varmapack_simx(0, 0, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X, E, rng);
  xCheck(!error);
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      xCheck(fabs(X[j*n + t] - E[j*n + t] - C[0]*z[t]) < 1e-12);
    }
  }
  randompack_free(rng);
}

static void check_arx_recursion(void) {
  int p = 1, q = 0, s = 1, r = 1, n = 6, M = 1, h = 2;
  double A[] = {0.6};
  double C[] = {0.4};
  double Sig[] = {1};
  double z[] = {1, -1, 2, -2, 3, -3};
  double X0[] = {0.5, -0.25};
  double X[6], E[6];
  randompack_rng *rng = seededRng(101);
  varmapack_error error = varmapack_simx(A, 0, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X, E, rng);
  xCheck(!error);
  xCheck(fabs(X[0] - X0[0]) < 1e-12);
  xCheck(fabs(X[1] - X0[1]) < 1e-12);
  for (int t=h; t<n; t++) {
    double pred = A[0]*X[t-1] + E[t] + C[0]*z[t];
    xCheck(fabs(X[t] - pred) < 1e-12);
  }
  randompack_free(rng);
}

static void check_no_return_shocks_reproducible(void) {
  int p = 1, q = 1, s = 2, r = 1, n = 8, M = 3, h = 3;
  double A[] = {0.4};
  double B[] = {0.2};
  double C[] = {0.3, -0.1};
  double Sig[] = {1};
  double z[] = {1, 0, -1, 2, -2, 1, 0, -1};
  double X0[] = {0.2, -0.1, 0.4};
  double X1[24], X2[24], E[24];
  randompack_rng *rng = seededRng(102);
  varmapack_error error = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X1, E, rng);
  xCheck(!error);
  randompack_free(rng);
  rng = seededRng(102);
  error = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1,
                         X2, 0, rng);
  xCheck(!error);
  checkArraySame(X1, X2, n*M);
  randompack_free(rng);
}

static void check_multipath_z_and_X0(void) {
  int p = 0, q = 0, s = 1, r = 1, n = 4, M = 2, h = 1;
  double C[] = {0.5};
  double Sig[] = {1};
  double z[] = {1, 2, 3, 4, -1, -2, -3, -4};
  double X0[] = {0.5, -0.5};
  double X[8], E[8];
  randompack_rng *rng = seededRng(104);
  varmapack_error error = varmapack_simx(0, 0, C, Sig, z, M, p, q, s, r, n, M,
                                         X0, h, M, X, E, rng);
  xCheck(!error);
  for (int j=0; j<M; j++) {
    xCheck(fabs(X[j*n] - X0[j]) < 1e-12);
    for (int t=h; t<n; t++) {
      xCheck(fabs(X[j*n + t] - E[j*n + t] - C[0]*z[j*n + t]) < 1e-12);
    }
  }
  randompack_free(rng);
}

static void check_invalid_input(void) {
  int p = 1, q = 1, s = 1, r = 1, n = 4, M = 2, h = 2;
  double A[] = {0.4}, B[] = {0.2}, C[] = {0.3}, Sig[] = {1}, z[] = {1, 2, 3, 4};
  double X0[] = {0, 0}, X[8], E[8];
  randompack_rng *rng = seededRng(103);
  xCheck(varmapack_simx(0, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, 0, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, 0, 1, p, q, s, r, n, M, X0, h, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, 0, h, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, 1, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 3, p, q, s, r, n, M, X0, h, 1, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 3, X, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, 0, E,
                        rng) ==
         VARMAPACK_INVALID_ARGUMENT);
  xCheck(varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E,
                        0) ==
         VARMAPACK_INVALID_ARGUMENT);
  randompack_free(rng);
}

void TestSimxEdgeCases(void) {
  check_white_noise_with_exog();
  check_arx_recursion();
  check_no_return_shocks_reproducible();
  check_multipath_z_and_X0();
  check_invalid_input();
}
