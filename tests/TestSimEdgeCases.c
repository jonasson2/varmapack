#include <math.h>
#include <stdbool.h>
#include "ExtraUtil.h"
#include "Tests.h"
#include "VarmaUtilities.h"
#include "randompack.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_white_noise(void) {
  int p = 0, q = 0, r = 2, n = 4, M = 1;
  double Sig[] = {1, 0.25, 0.25, 2};
  double X[8], E[8];
  varmapack_error error;
  randompack_rng *rng = seededRng(7);
  error = varmapack_sim(0, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(!error);
  checkArraySame(X, E, r*n*M);
  randompack_free(rng);
}

static void check_minimal_length(void) {
  int p = 2, q = 1, r = 1, n = 2, M = 1;
  double A[] = {0.3, -0.1};
  double B[] = {0.2};
  double Sig[] = {1};
  double X[2], E[2];
  varmapack_error error;
  randompack_rng *rng = seededRng(11);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(!error);
  checkArrayFinite(X, n);
  checkArrayFinite(E, n);
  randompack_free(rng);
}

static void check_null_E(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {0.4};
  double Sig[] = {1};
  double X[5];
  varmapack_error error;
  randompack_rng *rng = seededRng(13);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X, 0);
  xCheck(!error);
  checkArrayFinite(X, n);
  randompack_free(rng);
}

static void check_nonstationary_without_X0(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {1.25};
  double Sig[] = {1};
  double X[5], E[5];
  varmapack_error error;
  randompack_rng *rng = seededRng(17);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_NONSTATIONARY);
  randompack_free(rng);
}

static void check_invalid_input(void) {
  int p = 1, q = 1, r = 1, n = 2, M = 1;
  double A[] = {0.5};
  double B[] = {0.2};
  double Sig[] = {1};
  double X0[] = {0};
  double X[2], E[2];
  randompack_rng *rng = randompack_create(0);
  varmapack_error error;
  xCheck(rng != 0);
  error = varmapack_sim(0, B, Sig, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, 0, 0, p, q, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 0, rng, 0, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 0, 0, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, -1, q, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, -1, r, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, 0, n, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, 0, M, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, 0, 0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 1, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, X0, 0, rng, X, E);
  xCheck(error == VARMAPACK_INVALID_ARGUMENT);
  randompack_free(rng);
}

static void check_nonstationary_with_X0(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {1.25};
  double Sig[] = {1};
  double X0[] = {2};
  double X[5], E[5];
  varmapack_error error;
  randompack_rng *rng = seededRng(19);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, X0, 1, rng, X, E);
  xCheck(!error);
  xCheck(fabs(X[0] - X0[0]) < 1e-15);
  for (int t=1; t<n; t++) {
    xCheck(fabs(X[t] - (A[0]*X[t-1] + E[t])) < 1e-12);
  }
  randompack_free(rng);
}

static void check_nonstationary_arma_with_X0(void) {
  int p = 1, q = 1, r = 1, n = 6, M = 1;
  double A[] = {1.25};
  double B[] = {0.5};
  double Sig[] = {1};
  double X0[] = {2};
  double X[6], E[6];
  varmapack_error error;
  randompack_rng *rng = seededRng(21);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, X0, 1, rng, X, E);
  xCheck(!error);
  xCheck(fabs(X[0] - X0[0]) < 1e-15);
  for (int t=1; t<n; t++) {
    double expected = A[0]*X[t-1] + E[t] + B[0]*E[t-1];
    xCheck(fabs(X[t] - expected) < 1e-12);
  }
  randompack_free(rng);
}

static void check_reproducibility(void) {
  int p = 1, q = 1, r = 1, n = 6, M = 2;
  double A[] = {0.4};
  double B[] = {0.2};
  double Sig[] = {1};
  double X1[12], E1[12], X2[12], E2[12];
  varmapack_error error;
  randompack_rng *rng = seededRng(23);
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 0, rng, X1, E1);
  xCheck(!error);
  xCheck(randompack_seed(23, 0, 0, rng));
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, 0, 0, rng, X2, E2);
  xCheck(!error);
  checkArraySame(X1, X2, n*M);
  checkArraySame(E1, E2, n*M);
  randompack_free(rng);
}

static void check_rng_state_advances(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {0.4};
  double Sig[] = {1};
  double X1[5], E1[5], X2[5], E2[5];
  varmapack_error error;
  int different = 0;
  randompack_rng *rng = seededRng(29);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X1, E1);
  xCheck(!error);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, 0, 0, rng, X2, E2);
  xCheck(!error);
  for (int i=0; i<n; i++) {
    if (fabs(X1[i] - X2[i]) > 1e-15 || fabs(E1[i] - E2[i]) > 1e-15) different = 1;
  }
  xCheck(different);
  randompack_free(rng);
}

static void check_nX0_longer_than_h(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {0.5};
  double Sig[] = {1};
  double X0[] = {2, 3, 4};
  double X[5], E[5];
  varmapack_error error;
  randompack_rng *rng = seededRng(31);
  error = varmapack_sim(A, 0, Sig, 0, p, q, r, n, M, X0, 3, rng, X, E);
  xCheck(!error);
  checkArraySame(X, X0, 3);
  for (int t=3; t<n; t++) {
    xCheck(fabs(X[t] - (A[0]*X[t-1] + E[t])) < 1e-12);
  }
  randompack_free(rng);
}

void TestSimEdgeCases(void) {
  check_white_noise();
  check_minimal_length();
  check_null_E();
  check_nonstationary_without_X0();
  check_nonstationary_with_X0();
  check_nonstationary_arma_with_X0();
  check_reproducibility();
  check_rng_state_advances();
  check_nX0_longer_than_h();
  check_invalid_input();
}
