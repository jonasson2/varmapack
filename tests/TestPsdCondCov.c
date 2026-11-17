#include <math.h>
#include <stdbool.h>
#include "BlasGateway.h"
#include "TestUtil.h"
#include "Tests.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "randompack.h"
#include "varmapack.h"

static void conditional_moments(double A[], double B[], double Sig[], double X0[], int p,
                                int q, int r, int h, double R[], double e[]) {
  int rh = r*h, r2 = r*r, info;
  double *C = 0, *G = 0, *S = 0, *SS = 0, *CC = 0, *wrk = 0;
  ASSERT(ALLOC(C, r2*(q+1)), "allocation failed");
  ASSERT(ALLOC(G, r2*(q+1)), "allocation failed");
  ASSERT(ALLOC(S, r2*(p+1)), "allocation failed");
  ASSERT(ALLOC(SS, rh*rh), "allocation failed");
  ASSERT(ALLOC(CC, rh*rh), "allocation failed");
  ASSERT(ALLOC(wrk, rh), "allocation failed");
  FindC(A, B, Sig, p, q, r, C);
  FindG(B, C, q, r, G);
  bool ok = varmapack_acvf(A, B, Sig, p, q, r, S, p);
  checkVarmapackSuccess(ok);
  xCheck(SBuild("All", S, A, G, p, q, r, h, SS));
  potrf("Low", rh, SS, rh, &info);
  xCheck(info == 0);
  CCBuild(A, C, p, q, r, h, CC);
  trsm("Left", "Low", "NT", "NotUD", rh, rh, 1, SS, rh, CC, rh);
  copy(rh, X0, 1, wrk, 1);
  trsv("Lo", "NT", "NotUD", rh, SS, rh, wrk, 1);
  gemv("T", rh, rh, 1, CC, rh, wrk, 1, 0, e, 1);
  setzero(rh*rh, R);
  for (int j=0; j<h; j++) {
    lacpy("Low", r, r, Sig, r, R + j*r*(rh + 1), rh);
  }
  syrk("Low", "T", rh, rh, -1, CC, rh, 1, R, rh);
  copylowertoupper(rh, R, rh);
  FREE(wrk);
  FREE(CC);
  FREE(SS);
  FREE(S);
  FREE(G);
  FREE(C);
}

static void check_case_support(int icase, double support_tol) {
  test_model model;
  int h, rh, n = 5, M = 1, info, nulls = 0;
  double *X0 = 0, *R = 0, *e = 0, *Eig = 0, *lam = 0;
  double *work = 0, *X = 0, *E = 0;
  loadIndexedTestModel(&model, icase);
  h = imax(model.p, model.q);
  rh = model.r*h;
  ASSERT(ALLOC(X0, rh), "allocation failed");
  ASSERT(ALLOC(R, rh*rh), "allocation failed");
  ASSERT(ALLOC(e, rh), "allocation failed");
  ASSERT(ALLOC(Eig, rh*rh), "allocation failed");
  ASSERT(ALLOC(lam, rh), "allocation failed");
  ASSERT(ALLOC(work, 3*rh), "allocation failed");
  ASSERT(ALLOC(X, model.r*n*M), "allocation failed");
  ASSERT(ALLOC(E, model.r*n*M), "allocation failed");
  for (int i=0; i<rh; i++) X0[i] = (i+1)/10.0;
  conditional_moments(model.A, model.B, model.Sig, X0, model.p, model.q,
                      model.r, h, R, e);
  lacpy("All", rh, rh, R, rh, Eig, rh);
  syev("V", "Low", rh, Eig, rh, lam, work, 3*rh, &info);
  xCheck(info == 0);
  test_check_group psd = {0};
  for (int i=0; i<rh; i++) {
    if (lam[i] < 1e-10) nulls++;
    groupCheck(&psd, lam[i] > -1e-12);
  }
  checkGroup(&psd, "conditional covariance is positive semidefinite");
  xCheck(nulls > 0 && nulls < rh);
  randompack_rng *rng = seededRng(42);
  bool ok = varmapack_sim(model.A, model.B, model.Sig, 0, 0, model.p,
                          model.q, model.r, n, M, X0, h, 1, X, E, rng);
  checkVarmapackSuccess(ok);
  test_check_group nullspace = {0}, support = {0};
  for (int k=0; k<nulls; k++) {
    double nr = 0;
    for (int i=0; i<rh; i++) {
      nr += fabs(dot(rh, R + i, rh, Eig + k*rh, 1));
    }
    groupCheck(&nullspace, nr < support_tol);
  }
  checkGroup(&nullspace, "null eigenvectors annihilate the covariance matrix");
  for (int j=0; j<M; j++) {
    double *Ej = E + j*model.r*n;
    for (int k=0; k<nulls; k++) {
      double d = dot(rh, Eig + k*rh, 1, Ej, 1);
      d -= dot(rh, Eig + k*rh, 1, e, 1);
      groupCheck(&support, fabs(d) < support_tol);
    }
  }
  checkGroup(&support, "draws lie in the conditional covariance support");
  randompack_free(rng);
  FREE(E);
  FREE(X);
  FREE(work);
  FREE(lam);
  FREE(Eig);
  FREE(e);
  FREE(R);
  FREE(X0);
  freeTestModel(&model);
}

void TestPsdCondCov(void) {
  // Case 4 has a numerical null direction that varies modestly by platform.
  check_case_support(4, 1e-6);
  check_case_support(7, 1e-10);
  // Case 14 has a structural null direction from its second MA block.
  check_case_support(14, 1e-6);
}
