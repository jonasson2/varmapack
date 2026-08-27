#include <math.h>
#include <stdbool.h>
#include "TestUtil.h"
#include "Tests.h"
#include "VarmaUtilities.h"
#include "randompack.h"
#include "varmapack.h"

static void check_white_noise(void) {
  int p = 0, q = 0, r = 2, n = 4, M = 1;
  double Sig[] = {1, 0.25, 0.25, 2};
  double X[8], E[8];
  bool ok;
  randompack_rng *rng = seededRng(7);
  ok = varmapack_sim(0, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X, E, r*n*M);
  randompack_free(rng);
}

static void check_time_dependent_mean(void) {
  int p = 0, q = 0, r = 2, n = 5, M = 1;
  double Sig[] = {1, 0.25, 0.25, 2};
  double mu[] = {
    10, 20, 11, 21, 12, 22
  };
  double X[10], E[10];
  bool ok;
  randompack_rng *rng = seededRng(8);
  ok = varmapack_sim(0, 0, Sig, mu, 3, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  for (int t=0; t<n; t++) {
    int imu = t < 3 ? t : 2;
    for (int i=0; i<r; i++) {
      xCheck(fabs(X[t*r+i] - E[t*r+i] - mu[imu*r+i]) < 1e-12);
    }
  }
  randompack_free(rng);
}

static void check_minimal_length(void) {
  int p = 2, q = 1, r = 1, n = 2, M = 1;
  double A[] = {0.3, -0.1};
  double B[] = {0.2};
  double Sig[] = {1};
  double X[2], E[2];
  bool ok;
  randompack_rng *rng = seededRng(11);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArrayFinite(X, n);
  checkArrayFinite(E, n);
  randompack_free(rng);
}

static void check_null_E(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {0.4};
  double Sig[] = {1};
  double X[5];
  bool ok;
  randompack_rng *rng = seededRng(13);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, 0, rng);
  checkVarmapackSuccess(ok);
  checkArrayFinite(X, n);
  randompack_free(rng);
}

static void check_null_E_reproducible_with_ma(void) {
  int p = 1, q = 2, r = 2, n = 9, M = 1;
  double A[] = {
    0.2, 0.05, 0.1, 0.3
  };
  double B[] = {
    0.1, -0.04, 0.03, 0.2, -0.05, 0.02, 0.06, -0.08
  };
  double Sig[] = {
    1, 0.2, 0.2, 1.5
  };
  double X1[18], X2[18];
  bool ok;
  randompack_rng *rng = seededRng(41);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X1, 0, rng);
  checkVarmapackSuccess(ok);
  reseedRng(rng, 41);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X2, 0, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X1, X2, r*n*M);
  randompack_free(rng);
}

static void check_null_E_multiple_replicates(void) {
  int p = 1, q = 2, r = 2, n = 9, M = 3;
  double A[] = {
    0.2, 0.05, 0.1, 0.3
  };
  double B[] = {
    0.1, -0.04, 0.03, 0.2, -0.05, 0.02, 0.06, -0.08
  };
  double Sig[] = {
    1, 0.2, 0.2, 1.5
  };
  double X1[54], X2[54], E[54];
  bool ok;
  randompack_rng *rng = seededRng(43);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X1, E, rng);
  checkVarmapackSuccess(ok);
  reseedRng(rng, 43);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X2, 0, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X1, X2, r*n*M);
  checkArrayFinite(E, r*n*M);
  randompack_free(rng);
}

static void check_example_2(void) {
  int n = 5;
  double expected[] = {
    1.3923, 2.0876, -1.1619, -1.9181, -1.7641,
    1.2999, 0.1571, -3.5218, 0.6110, -2.7306
  };
  double X[10];
  bool ok;
  test_model model;
  randompack_rng *rng = seededRng(123);
  loadNamedTestModel(&model, "smallAR1");
  ok = varmapack_sim(model.A, model.B, model.Sig, 0, 0, model.p, model.q,
                     model.r, n, 1, 0, 0, 1, X, 0, rng);
  checkVarmapackSuccess(ok);
  checkArrayTol(X, expected, 10, 5e-5);
  xCheck(fabs(varmapack_specrad(model.A, model.r, model.p) - 0.2) < 5e-4);
  xCheck(varmapack_ma_specrad(model.B, model.r, model.q) == 0);
  freeTestModel(&model);
  randompack_free(rng);
}

static void check_conditional_ma_recursion(void) {
  int p = 1, q = 3, r = 1, n = 8, M = 1;
  double A[] = {0.4};
  double B[] = {0.4, -0.2, 0.1};
  double Sig[] = {1};
  double X0[] = {2, 3, 4};
  double X[8], E[8];
  bool ok;
  randompack_rng *rng = seededRng(37);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, X0, 3, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  for (int t=3; t<n; t++) {
    double expected = A[0]*X[t-1] + E[t];
    expected += B[0]*E[t-1] + B[1]*E[t-2] + B[2]*E[t-3];
    xCheck(fabs(X[t] - expected) < 1e-12);
  }
  randompack_free(rng);
}

static void check_singular_sigma(void) {
  int p = 1, q = 1, r = 2, n = 7, M = 2;
  double A[] = {
    0.2, 0.05, -0.1, 0.25
  };
  double B[] = {
    0.15, -0.03, 0.04, 0.12
  };
  double Sig[] = {
    1, 1, 1, 1
  };
  double X[28], E[28];
  bool ok;
  randompack_rng *rng = seededRng(47);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArrayFinite(X, r*n*M);
  checkArrayFinite(E, r*n*M);
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      int i = j*r*n + t*r;
      xCheck(fabs(E[i] - E[i+1]) < 1e-12);
    }
    for (int t=1; t<n; t++) {
      int i = j*r*n + t*r;
      int im1 = j*r*n + (t-1)*r;
      double xPred0 = E[i] + A[0]*X[im1] + A[2]*X[im1+1] + B[0]*E[im1]
                      + B[2]*E[im1+1];
      double xPred1 = E[i+1] + A[1]*X[im1] + A[3]*X[im1+1] + B[1]*E[im1]
                      + B[3]*E[im1+1];
      xCheck(fabs(X[i] - xPred0) < 1e-12);
      xCheck(fabs(X[i+1] - xPred1) < 1e-12);
    }
  }
  randompack_free(rng);
}

static void check_singular_startup_covariance(void) {
  int p = 1, q = 0, r = 2, n = 5, M = 1;
  double A[] = {
    0.2, 0, 0, 0.2
  };
  double Sig[] = {
    1, 1, 1, 1
  };
  double X0[] = {0, 0};
  double X[10], E[10];
  bool ok;
  randompack_rng *rng = seededRng(53);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArrayFinite(X, r*n*M);
  checkArrayFinite(E, r*n*M);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, X0, 1, 1, X, E, rng);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

static void check_singular_conditional_covariance(void) {
  int p = 1, q = 0, r = 2, n = 5, M = 3;
  double A[] = {
    0.1, 0.1, 0.1, 0.1
  };
  double Sig[] = {
    2, 0, 0, 2
  };
  double X[30], E[30];
  bool ok;
  randompack_rng *rng = seededRng(59);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  for (int j=0; j<M; j++) {
    int offset = j*r*n;
    double u0 = X[offset] - E[offset];
    double u1 = X[offset + 1] - E[offset + 1];
    xCheck(fabs(u0 - u1) < 1e-12);
  }
  randompack_free(rng);
}

static void check_nonstationary_without_X0(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {1.25};
  double Sig[] = {1};
  double X[5], E[5];
  bool ok;
  randompack_rng *rng = seededRng(17);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

static void check_indefinite_sigma(void) {
  int p = 1, q = 0, r = 1, n = 3, M = 1;
  double A[] = {0.5};
  double Sig[] = {-1};
  double X[3], E[3];
  bool ok;
  randompack_rng *rng = seededRng(18);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

static void check_invalid_input(void) {
  int p = 1, q = 1, r = 1, n = 2, M = 1;
  double A[] = {0.5};
  double B[] = {0.2};
  double Sig[] = {1};
  double X0[] = {0};
  double X[2], E[2];
  randompack_rng *rng = seededRng(44);
  bool ok;
  ok = varmapack_sim(0, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, 0, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, 0, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X, E, 0);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, -1, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, -1, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, 0, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, 0, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, 0, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 1, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, X0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(0, 0, Sig, 0, 0, 0, 0, 1, 1, 1, X0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, 0, 1, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, Sig, n + 1, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, B, Sig, Sig, 0, p, q, r, n, M, 0, 0, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  randompack_free(rng);
}

static void check_problem_too_large(void) {
  double A[] = {0};
  double Sig[] = {1};
  double X[1], E[1];
  bool ok;
  randompack_rng *rng = seededRng(51);
  ok = varmapack_sim(0, 0, Sig, 0, 0, 0, 0, 50000, 50000, 1, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  ok = varmapack_sim(A, 0, Sig, 0, 0, 50000, 0, 1, 50000, 1, 0, 0, 1, X, E, rng);
  checkVarmapackFailure(ok);
  randompack_free(rng);
}

static void check_nonstationary_with_X0(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {1.25};
  double Sig[] = {1};
  double X0[] = {2};
  double X[5], E[5];
  bool ok;
  randompack_rng *rng = seededRng(19);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, X0, 1, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  xCheck(fabs(X[0] - X0[0]) < 1e-15);
  for (int t=1; t<n; t++) {
    xCheck(fabs(X[t] - (A[0]*X[t-1] + E[t])) < 1e-12);
  }
  randompack_free(rng);
}

static void check_multipath_X0(void) {
  int p = 1, q = 0, r = 1, n = 3, M = 2;
  double A[] = {1.2};
  double Sig[] = {1};
  double X0[] = {2, 4};
  double X[6], E[6];
  randompack_rng *rng = seededRng(45);
  bool ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, X0, 1,
                                        M, X, E, rng);
  checkVarmapackSuccess(ok);
  xCheck(fabs(X[0] - X0[0]) < 1e-15);
  xCheck(fabs(X[n] - X0[1]) < 1e-15);
  randompack_free(rng);
}

static void check_nonstationary_arma_with_X0(void) {
  int p = 1, q = 1, r = 1, n = 6, M = 1;
  double A[] = {1.25};
  double B[] = {0.5};
  double Sig[] = {1};
  double X0[] = {2, 3, 4};
  double X[6], E[6];
  bool ok;
  randompack_rng *rng = seededRng(21);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, X0, 3, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X, X0, 3);
  for (int t=1; t<n; t++) {
    xCheck(fabs(X[t] - (A[0]*X[t - 1] + E[t] + B[0]*E[t - 1])) < 1e-12);
  }
  randompack_free(rng);
}

static void check_nonstationary_arma_with_presample_shocks(void) {
  int p = 1, q = 3, r = 1, n = 6, M = 1;
  double A[] = {1.25};
  double B[] = {0.3, -0.1, 0.2};
  double Sig[] = {1};
  double X0[] = {2, 3, 4, 5};
  double X[6], E[6];
  bool ok;
  randompack_rng *rng = seededRng(22);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, X0, 4, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X, X0, 4);
  for (int t=3; t<n; t++) {
    double expected = A[0]*X[t - 1] + E[t] + B[0]*E[t - 1];
    expected += B[1]*E[t - 2] + B[2]*E[t - 3];
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
  bool ok;
  randompack_rng *rng = seededRng(23);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X1, E1, rng);
  checkVarmapackSuccess(ok);
  reseedRng(rng, 23);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X2, E2, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X1, X2, n*M);
  checkArraySame(E1, E2, n*M);
  randompack_free(rng);
}

static void check_rng_state_advances(void) {
  int p = 1, q = 0, r = 1, n = 5, M = 1;
  double A[] = {0.4};
  double Sig[] = {1};
  double X1[5], E1[5], X2[5], E2[5];
  bool ok;
  int different = 0;
  randompack_rng *rng = seededRng(29);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X1, E1, rng);
  checkVarmapackSuccess(ok);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, 0, 0, 1, X2, E2, rng);
  checkVarmapackSuccess(ok);
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
  bool ok;
  randompack_rng *rng = seededRng(31);
  ok = varmapack_sim(A, 0, Sig, 0, 0, p, q, r, n, M, X0, 3, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X, X0, 3);
  for (int t=3; t<n; t++) {
    xCheck(fabs(X[t] - (A[0]*X[t-1] + E[t])) < 1e-12);
  }
  randompack_free(rng);
}

static void check_nX0_longer_than_h_with_ma(void) {
  int p = 1, q = 1, r = 1, n = 6, M = 1;
  double A[] = {0.4};
  double B[] = {0.25};
  double Sig[] = {1};
  double X0[] = {2, 3, 4};
  double X[6], E[6];
  bool ok;
  randompack_rng *rng = seededRng(33);
  ok = varmapack_sim(A, B, Sig, 0, 0, p, q, r, n, M, X0, 3, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  checkArraySame(X, X0, 3);
  for (int t=3; t<n; t++) {
    double expected = A[0]*X[t-1] + E[t] + B[0]*E[t-1];
    xCheck(fabs(X[t] - expected) < 1e-12);
  }
  randompack_free(rng);
}

void TestSimEdgeCases(void) {
  check_white_noise();
  check_time_dependent_mean();
  check_minimal_length();
  check_null_E();
  check_null_E_reproducible_with_ma();
  check_null_E_multiple_replicates();
  check_example_2();
  check_conditional_ma_recursion();
  check_singular_sigma();
  check_singular_startup_covariance();
  check_singular_conditional_covariance();
  check_nonstationary_without_X0();
  check_indefinite_sigma();
  check_nonstationary_with_X0();
  check_multipath_X0();
  check_nonstationary_arma_with_X0();
  check_nonstationary_arma_with_presample_shocks();
  check_reproducibility();
  check_rng_state_advances();
  check_nX0_longer_than_h();
  check_nX0_longer_than_h_with_ma();
  check_invalid_input();
  check_problem_too_large();
}
