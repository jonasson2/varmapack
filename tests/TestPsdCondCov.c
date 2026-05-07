#include <math.h>
#include <stdbool.h>
#include "BlasGateway.h"
#include "ExtraUtil.h"
#include "Tests.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "randompack.h"
#include "varmapack.h"
#include "xCheck.h"

static void SExtend(double A[], double G[], double S[], double Scol[], int p, int q,
                    int r, int n) {
  int iScol = n*r;
  for (int j=0; j<p+1; j++) {
    lacpy("All", r, r, S + j*r*r, r, Scol + j*r, iScol);
  }
  for (int j=p+1; j<n; j++) {
    double *Scolj = Scol + j*r;
    if (j <= q) {
      lacpy("All", r, r, G + j*r*r, r, Scolj, iScol);
    }
    for (int i=0; i<p; i++) {
      gemm("N", "N", r, r, r, 1, A + i*r*r, r, Scolj - (i+1)*r, iScol, 1,
           Scolj, iScol);
    }
  }
}

static bool SBuild(double S[], double A[], double G[], int p, int q, int r, int n,
                   double SS[]) {
  if (n == 0) return true;
  int m = imax(p+1, n);
  double *Scol = 0;
  if (!ALLOC(Scol, r*m*r)) return false;
  SExtend(A, G, S, Scol, p, q, r, m);
  for (int j=0; j<n; j++) {
    double *SSj = SS + j*r*n*r + j*r;
    double *SSi = SSj + r*n*r;
    lacpy("All", r*(n-j), r, Scol, r*m, SSj, r*n);
    if (j < n-1) copytranspose(r*(n-j-1), r, Scol + r, r*m, SSi, r*n);
  }
  FREE(Scol);
  return true;
}

static void CCBuild(double A[], double C[], int p, int q, int r, int n, double CC[]) {
  setzero(r*n*r*n, CC);
  for (int j=0; j<=q && j<n; j++) {
    lacpy("All", r, r, C + j*r*r, r, CC + j*r, r*n);
  }
  for (int j=q+1; j<n; j++) {
    double *CCj = CC + j*r;
    for (int i=1; i<=j && i<=p; i++) {
      gemm("N", "N", r, r, r, 1, A + (i-1)*r*r, r, CC + (j-i)*r, r*n, 1,
           CCj, r*n);
    }
  }
  for (int j=1; j<n; j++) {
    lacpy("All", r*(n-j), r, CC, r*n, CC + j*r*n*r + j*r, r*n);
  }
}

static void conditional_moments(double A[], double B[], double Sig[], double X0[], int p,
                                int q, int r, int h, double R[], double e[]) {
  int rh = r*h, r2 = r*r, info;
  double *C = 0, *G = 0, *S = 0, *SS = 0, *CC = 0, *wrk = 0, *work = 0;
  xCheck(ALLOC(C, r2*(q+1)));
  xCheck(ALLOC(G, r2*(q+1)));
  xCheck(ALLOC(S, r2*(p+1)));
  xCheck(ALLOC(SS, rh*rh));
  xCheck(ALLOC(CC, rh*rh));
  xCheck(ALLOC(wrk, rh));
  FindCG(A, B, Sig, p, q, r, C, G);
  varmapack_error error = varmapack_acvf(A, B, Sig, p, q, r, S, p);
  xCheck(!error);
  xCheck(SBuild(S, A, G, p, q, r, h, SS));
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
  FREE(work);
  FREE(wrk);
  FREE(CC);
  FREE(SS);
  FREE(S);
  FREE(G);
  FREE(C);
}

static void check_case7_support(void) {
  int p = 0, q = 0, r = 0, icase = 7, h, rh, n = 5, M = 1, info, nulls = 0;
  char name[32] = "";
  double *A = 0, *B = 0, *Sig = 0, *X0 = 0, *R = 0, *e = 0, *Eig = 0, *lam = 0;
  double *work = 0, *X = 0, *E = 0;
  varmapack_error error = varmapack_testcase(0, 0, 0, name, &p, &q, &r,
                                               &icase, 0, 0);
  xCheck(!error);
  h = imax(p, q);
  rh = r*h;
  xCheck(ALLOC(A, r*r*(p > 0 ? p : 1)));
  xCheck(ALLOC(B, r*r*q));
  xCheck(ALLOC(Sig, r*r));
  xCheck(ALLOC(X0, rh));
  xCheck(ALLOC(R, rh*rh));
  xCheck(ALLOC(e, rh));
  xCheck(ALLOC(Eig, rh*rh));
  xCheck(ALLOC(lam, rh));
  xCheck(ALLOC(work, 3*rh));
  xCheck(ALLOC(X, r*n*M));
  xCheck(ALLOC(E, r*n*M));
  error = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, 0);
  xCheck(!error);
  for (int i=0; i<rh; i++) X0[i] = (i + 1)/10.0;
  conditional_moments(A, B, Sig, X0, p, q, r, h, R, e);
  lacpy("All", rh, rh, R, rh, Eig, rh);
  syev("V", "Low", rh, Eig, rh, lam, work, 3*rh, &info);
  xCheck(info == 0);
  for (int i=0; i<rh; i++) {
    if (lam[i] < 1e-10) nulls++;
    xCheck(lam[i] > -1e-12);
  }
  xCheck(nulls > 0 && nulls < rh);
  randompack_rng *rng = randompack_create(0);
  xCheck(randompack_seed(42, 0, 0, rng));
  error = varmapack_sim(A, B, Sig, 0, p, q, r, n, M, X0, h, rng, X, E);
  xCheck(!error);
  for (int k=0; k<nulls; k++) {
    double nr = 0;
    for (int i=0; i<rh; i++) {
      nr += fabs(dot(rh, R + i, rh, Eig + k*rh, 1));
    }
    xCheck(nr < 1e-10);
  }
  for (int j=0; j<M; j++) {
    double *Ej = E + j*r*n;
    for (int k=0; k<nulls; k++) {
      double d = dot(rh, Eig + k*rh, 1, Ej, 1);
      d -= dot(rh, Eig + k*rh, 1, e, 1);
      xCheck(fabs(d) < 1e-10);
    }
  }
  randompack_free(rng);
  FREE(E);
  FREE(X);
  FREE(work);
  FREE(lam);
  FREE(Eig);
  FREE(e);
  FREE(R);
  FREE(X0);
  FREE(Sig);
  FREE(B);
  FREE(A);
}

void TestPsdCondCov(void) {
  check_case7_support();
}
