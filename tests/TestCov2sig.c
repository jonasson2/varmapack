#include <math.h>
#include <string.h>
#include "TestUtil.h"
#include "Tests.h"
#include "Lyapunov.h"
#include "error.h"
#include "varmapack.h"

static void check_recovery(double A[], double B[], double cov[], int p, int q, int r,
                            double expected[], double tol) {
  double *Sig = 0, *original = 0, *recovered = 0;
  int rr = r*r;
  ASSERT(ALLOC(Sig, rr) && ALLOC(original, rr) && ALLOC(recovered, rr),
         "allocation failed");
  memcpy(original, cov, rr*sizeof(*cov));
  bool ok = varmapack_cov2sig(A, B, cov, p, q, r, Sig);
  checkVarmapackSuccess(ok);
  if (ok) {
    checkArrayTol(Sig, expected, rr, tol);
    checkArraySame(cov, original, rr);
    ok = varmapack_acvf(A, B, Sig, p, q, r, recovered, 0);
    checkVarmapackSuccess(ok);
    if (ok) {
      double scale = 1;
      for (int i=0; i<rr; i++) scale = fmax(scale, fabs(cov[i]));
      checkArrayTol(recovered, cov, rr, tol*scale);
    }
  }
  FREE(recovered); FREE(original); FREE(Sig);
}

static void check_scalar_cases(void) {
  double A[] = {0.5, -0.1}, B[] = {0.25, -0.1};
  double V[] = {4}, Sig[] = {3};
  check_recovery(A, 0, V, 1, 0, 1, Sig, 1e-12);
  V[0] = 1.25;
  Sig[0] = 1.25;
  check_recovery(0, 0, V, 0, 0, 1, Sig, 1e-12);
  V[0] = 2.125;
  Sig[0] = 2;
  check_recovery(0, B, V, 0, 1, 1, Sig, 1e-12);
  V[0] = 3.5;
  check_recovery(A, B, V, 1, 1, 1, Sig, 1e-12);
  V[0] = 3.3266666666666667;
  check_recovery(A, B, V, 1, 2, 1, Sig, 1e-12);
  V[0] = 2.5462962962962967;
  check_recovery(A, 0, V, 2, 0, 1, Sig, 1e-12);
  double persistence[] = {-0.999, 0.99, 0.999, 1-1e-10};
  for (int i=0; i<LEN(persistence); i++) {
    A[0] = persistence[i];
    V[0] = Sig[0]/(1-A[0]*A[0]);
    check_recovery(A, 0, V, 1, 0, 1, Sig, 1e-10);
  }
}

static void check_vector_var1(void) {
  double A[] = {0.4, -0.1, 0.2, 0.3};
  double original[] = {0.4, -0.1, 0.2, 0.3};
  double V[] = {2, 0.3, 0.3, 1};
  double Sig[4];
  // Independent closed form: Sigma = V - A V A^T.
  for (int j=0; j<2; j++) {
    for (int i=0; i<2; i++) {
      Sig[i+2*j] = V[i+2*j];
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) Sig[i+2*j] -= A[i+2*k]*V[k+2*l]*A[j+2*l];
      }
    }
  }
  check_recovery(A, 0, V, 1, 0, 2, Sig, 1e-12);
  checkArraySame(A, original, 4);
}

static void check_named_round_trips(void) {
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "max";
  int count = 0, p = 0, q = 0, r = 0;
  bool ok = varmapack_testcase(name, &count, 0, &p, &q, &r, 0, 0, 0, 0);
  checkVarmapackSuccess(ok);
  for (int i=1; i<=count; i++) {
    test_model model;
    double *cov = 0;
    loadIndexedTestModel(&model, i);
    xCheckAddMsg(model.name);
    ASSERT(ALLOC(cov, model.r*model.r), "allocation failed");
    ok = varmapack_acvf(model.A, model.B, model.Sig, model.p, model.q, model.r,
                       cov, 0);
    checkVarmapackSuccess(ok);
    if (ok) {
      check_recovery(model.A, model.B, cov, model.p, model.q, model.r,
                      model.Sig, 1e-9);
    }
    FREE(cov);
    freeTestModel(&model);
  }
  xCheckAddMsg("");
}

static void check_impulse_recovery(int p, int q, int r) {
  double A[27] = {0}, B[27] = {0}, Sig[9] = {0}, cov[9] = {0};
  double Psi[301*9] = {0};
  int rr = r*r;
  for (int i=0; i<rr*p; i++) A[i] = 0.3*sin(i+1)/(r*p);
  for (int i=0; i<rr*q; i++) B[i] = 0.4*cos(i+1)/(r*q);
  for (int j=0; j<r; j++) {
    for (int i=0; i<r; i++) Sig[i+j*r] = i == j ? i+1 : 0.2;
    Psi[j+j*r] = 1;
  }
  // Direct impulse-response sum, independent of both covariance solvers.
  for (int t=0; t<=300; t++) {
    double *P = Psi+rr*t;
    if (t > 0) {
      for (int j=0; j<r; j++) {
        for (int i=0; i<r; i++) {
          if (t <= q) P[i+j*r] = B[i+j*r+(t-1)*rr];
          for (int lag=1; lag<=p && lag<=t; lag++) {
            for (int k=0; k<r; k++) {
              P[i+j*r] += A[i+k*r+(lag-1)*rr]*Psi[k+j*r+(t-lag)*rr];
            }
          }
        }
      }
    }
    for (int j=0; j<r; j++) {
      for (int i=0; i<r; i++) {
        for (int k=0; k<r; k++) {
          for (int l=0; l<r; l++) cov[i+j*r] += P[i+k*r]*Sig[k+l*r]*P[j+l*r];
        }
      }
    }
  }
  check_recovery(A, B, cov, p, q, r, Sig, 1e-11);
}

static void check_lyapunov_target(void) {
  int r = 16, rr = r*r, p = 2, q = 1;
  double A[512] = {0}, B[256] = {0}, Sig[256], S[768], C[512], G[512];
  for (int j=0; j<r; j++) {
    A[j+j*r] = 0.4;
    A[j+j*r+rr] = 0.1;
    A[(j+1)%r+j*r] = 0.03;
    B[j+j*r] = 0.2;
    B[(j+1)%r+j*r] = -0.1;
    for (int i=0; i<r; i++) Sig[i+j*r] = i == j ? 1 : 0.02;
  }
  bool ok = LyapunovFactorizeSolve(A, B, Sig, p, q, r, S, C, G);
  xCheck(ok);
  if (ok) check_recovery(A, B, S, p, q, r, Sig, 1e-10);
}

static void check_storage_and_scale(void) {
  double A[] = {0.5, 0, 0, 0.25};
  double B[] = {0.2, 0.3, 0.1, -0.1};
  double Bcopy[] = {0.2, 0.3, 0.1, -0.1};
  double V[] = {2.11, 0.615, 0.615, 1.16};
  double expected[] = {2, 0.5, 0.5, 1};
  double Sig[4];
  check_recovery(0, B, V, 0, 1, 2, expected, 1e-12);
  checkArraySame(B, Bcopy, 4);
  V[2] = NAN;
  bool ok = varmapack_cov2sig(0, B, V, 0, 1, 2, Sig);
  checkVarmapackSuccess(ok);
  if (ok) checkArrayTol(Sig, expected, 4, 1e-12);
  xCheck(isnan(V[2]));
  ok = varmapack_cov2sig(0, B, V, 0, 1, 2, V);
  checkVarmapackSuccess(ok);
  if (ok) checkArrayTol(V, expected, 4, 1e-12);
  double factors[] = {1e-200, 1, 1e200};
  for (int k=0; k<LEN(factors); k++) {
    V[0] = 4*factors[k];
    V[1] = 0.5*factors[k];
    V[2] = NAN;
    V[3] = 2*factors[k];
    ok = varmapack_cov2sig(A, 0, V, 1, 0, 2, Sig);
    checkVarmapackSuccess(ok);
    if (ok) {
      for (int i=0; i<4; i++) Sig[i] /= factors[k];
      double want[] = {3, 0.4375, 0.4375, 1.875};
      checkArrayTol(Sig, want, 4, 1e-12);
    }
  }
}

static void check_semidefinite(void) {
  double A[] = {0.5, 0, 0, 0.5};
  double V[] = {4, 4, 4, 4};
  double expected[] = {3, 3, 3, 3};
  check_recovery(A, 0, V, 1, 0, 2, expected, 1e-12);
  V[0] = V[1] = V[2] = V[3] = 0;
  expected[0] = expected[1] = expected[2] = expected[3] = 0;
  check_recovery(A, 0, V, 1, 0, 2, expected, 1e-12);
  V[0] = expected[0] = 1;
  check_recovery(0, 0, V, 0, 0, 2, expected, 1e-12);
  // A rank-one innovation covariance can produce a full-rank model covariance.
  double rotation[] = {0, 0.5, -0.5, 0};
  V[0] = 16.0/15;
  V[3] = 4.0/15;
  check_recovery(rotation, 0, V, 1, 0, 2, expected, 1e-12);
}

static void check_singular_and_incompatible(void) {
  double V[] = {1, 0, 0, 1};
  double B[] = {0, 1, -1, 0};
  double A[] = {0.5, 0, 2, 0.5};
  double Sig[] = {7, 8, 9, 10}, original[] = {7, 8, 9, 10};
  bool ok = varmapack_cov2sig(0, B, V, 0, 1, 2, Sig);
  checkVarmapackFailure(ok);
  xCheck(varmapack_last_error() && strstr(varmapack_last_error(), "singular"));
  checkArraySame(Sig, original, 4);
  B[1] = 1-1e-14;
  B[2] = -B[1];
  ok = varmapack_cov2sig(0, B, V, 0, 1, 2, Sig);
  checkVarmapackFailure(ok);
  xCheck(varmapack_last_error() && strstr(varmapack_last_error(), "ill-conditioned"));
  checkArraySame(Sig, original, 4);
  ok = varmapack_cov2sig(A, 0, V, 1, 0, 2, Sig);
  checkVarmapackFailure(ok);
  checkArraySame(Sig, original, 4);
  double copyV[] = {1, 0, 0, 1};
  ok = varmapack_cov2sig(A, 0, V, 1, 0, 2, V);
  checkVarmapackFailure(ok);
  checkArraySame(V, copyV, 4);
  B[0] = B[1] = B[2] = 0;
  B[3] = 2;
  V[1] = V[2] = 1;
  ok = varmapack_cov2sig(0, B, V, 0, 1, 2, Sig);
  checkVarmapackFailure(ok);
  checkArraySame(Sig, original, 4);
}

static void check_invalid_input(void) {
  double A[] = {0.5}, B[] = {0.25}, V[] = {1}, Sig[] = {7};
  checkVarmapackFailure(varmapack_cov2sig(0, B, V, 1, 1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, 0, V, 1, 1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, 0, 1, 1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, 0));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, -1, 1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, -1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 0, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, INT_MAX, 1, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, INT_MAX, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, INT_MAX, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 50000, 0, 1, Sig));
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1000, Sig));
  A[0] = 1;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  A[0] = -1.1;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  A[0] = NAN;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  A[0] = 0.5;
  B[0] = INFINITY;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  B[0] = 1e200;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  B[0] = 0.25;
  V[0] = NAN;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  V[0] = -1e-200;
  checkVarmapackFailure(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
  double indefinite[] = {1e-200, 2e-200, 2e-200, 1e-200};
  double output[4] = {7, 7, 7, 7};
  checkVarmapackFailure(varmapack_cov2sig(0, 0, indefinite, 0, 0, 2, output));
  indefinite[1] = INFINITY;
  checkVarmapackFailure(varmapack_cov2sig(0, 0, indefinite, 0, 0, 2, output));
  double extremes[] = {1e200, 0, 0, 1e-200};
  double recovered[4];
  checkVarmapackSuccess(varmapack_cov2sig(0, 0, extremes, 0, 0, 2, recovered));
  checkArraySame(recovered, extremes, 4);
  double V2[] = {1,0,0,.25}, B2[] = {0,1,0,0};
  checkVarmapackFailure(varmapack_cov2sig(0, B2, V2, 0, 1, 2, output));
  V2[0] = 1e16; B2[1] = 1e-8;
  checkVarmapackFailure(varmapack_cov2sig(0, B2, V2, 0, 1, 2, output));
  V2[3] = 1.25;
  checkVarmapackSuccess(varmapack_cov2sig(0, B2, V2, 0, 1, 2, recovered));
  xCheck(fabs(recovered[3]-.25) < 1e-14);
  V[0] = nextafter(0, 1);
  A[0] = 0.9;
  checkVarmapackFailure(varmapack_cov2sig(A, 0, V, 1, 0, 1, Sig));
  test_check_group unchanged = {0};
  for (int i=0; i<4; i++) groupCheck(&unchanged, output[i] == 7);
  checkGroup(&unchanged, "failed calls leave the output unchanged");
  xCheck(Sig[0] == 7);
  V[0] = 1;
  A[0] = 0.5;
  checkVarmapackSuccess(varmapack_cov2sig(A, B, V, 1, 1, 1, Sig));
}

void TestCov2sig(void) {
  check_scalar_cases();
  check_vector_var1();
  check_named_round_trips();
  for (int r=1; r<=3; r++) {
    for (int p=0; p<=3; p++) {
      for (int q=0; q<=3; q++) check_impulse_recovery(p, q, r);
    }
  }
  check_lyapunov_target();
  check_storage_and_scale();
  check_semidefinite();
  check_singular_and_incompatible();
  check_invalid_input();
}
