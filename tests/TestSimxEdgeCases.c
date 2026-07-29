#include <math.h>
#include <limits.h>
#include <string.h>
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
  bool ok = varmapack_simx(0, 0, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X, E, rng);
  checkVarmapackSuccess(ok);
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
  bool ok = varmapack_simx(A, 0, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  xCheck(fabs(X[0] - X0[0]) < 1e-12);
  xCheck(fabs(X[1] - X0[1]) < 1e-12);
  for (int t=h; t<n; t++) {
    double pred = A[0]*X[t-1] + E[t] + C[0]*z[t];
    xCheck(fabs(X[t] - pred) < 1e-12);
  }
  randompack_free(rng);
}

static void check_varma_without_exog(void) {
  int p = 2, q = 1, s = 0, r = 1, n = 7, M = 2, h = 3;
  double A[] = {0.4, -0.1};
  double B[] = {0.2};
  double Sig[] = {1};
  double X0[] = {0.5, -0.25, 0.1};
  double X[14], E[14];
  randompack_rng *rng = seededRng(105);
  bool ok = varmapack_simx(A, B, 0, Sig, 0, 1, p, q, s, r, n, M,
                           X0, h, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  for (int j=0; j<M; j++) {
    for (int t=h; t<n; t++) {
      int it = j*n + t;
      double prediction = E[it] + A[0]*X[it-1] + A[1]*X[it-2] + B[0]*E[it-1];
      xCheck(fabs(X[it] - prediction) < 1e-12);
    }
  }
  randompack_free(rng);
}

static void check_singular_sigma(void) {
  int p = 0, q = 0, s = 1, r = 2, n = 2, M = 1, h = 1;
  double C[] = {0, 0};
  double Sig[] = {1, 1, 1, 1};
  double z[] = {0, 0};
  double X0[] = {0, 0};
  double X[4];
  randompack_rng *rng = seededRng(107);
  bool ok = varmapack_simx(0, 0, C, Sig, z, 1, p, q, s, r, n, M,
                            X0, h, 1, X, 0, rng);
  checkVarmapackFailure(ok);
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
  bool ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M,
                                         X0, h, 1, X1, E, rng);
  checkVarmapackSuccess(ok);
  randompack_free(rng);
  rng = seededRng(102);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1,
                         X2, 0, rng);
  checkVarmapackSuccess(ok);
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
  bool ok = varmapack_simx(0, 0, C, Sig, z, M, p, q, s, r, n, M,
                                         X0, h, M, X, E, rng);
  checkVarmapackSuccess(ok);
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
  bool ok;
  ok = varmapack_simx(0, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, 0, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, 0, 1, p, q, s, r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, 0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, 1, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 3, p, q, s, r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 3, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, 0, Sig, 0, 3, p, q, 0, r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, 0, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(A, B, C, Sig, z, 1, p, q, s, r, n, M, X0, h, 1, X, E, 0);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

static void check_oversized_dimensions(void) {
  double Sig[] = {1}, X0[] = {0}, X[] = {0};
  randompack_rng *rng = seededRng(106);
  bool ok = varmapack_simx(0, 0, 0, Sig, 0, 1, 0, 0, 0, INT_MAX, 2, 1,
                           X0, 1, 1, X, 0, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_simx(0, Sig, 0, Sig, 0, 1, 0, INT_MAX, 0, 1, 1, 1,
                       X0, 1, 1, X, 0, rng);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

void TestSimxEdgeCases(void) {
  check_white_noise_with_exog();
  check_arx_recursion();
  check_varma_without_exog();
  check_singular_sigma();
  check_no_return_shocks_reproducible();
  check_multipath_z_and_X0();
  check_invalid_input();
  check_oversized_dimensions();
}
