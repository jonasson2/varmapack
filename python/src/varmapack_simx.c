// varmapack_simx — simulate VARMAX time series with fixed exogenous data
//
// See varmapack.h for parameter descriptions and storage conventions.

#include "BlasGateway.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"
#include "randompack.h"
#include "varmapack.h"
#include "varmapack_config.h"

static void buildH(double B[], int q, int r, int h, int t0, int firstShockTime,
                   double H[]);
static void startupResidual(double A[], double C[], double z[], double X0[], int p,
                            int s, int r, int h, int t0, double residual[]);
static varmapack_error drawStartupShocks(double A[], double B[], double C[], double Sig[],
                                         double z[], int Mz, double X0[], int MX0,
                                         double Eall[], int p, int q, int s, int r, int n,
                                         int h, int t0, int firstShockTime,
                                         int firstActiveShock, int ldE, int M,
                                         randompack_rng *rng);
static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int r, int n, int M, int h,
                        int firstShockTime, randompack_rng *rng);

varmapack_error varmapack_simx(double A[], double B[], double C[], double Sig[],
                               double z[], int Mz, int p, int q, int s, int r,
                               int n, int M, double X0[], int h, int MX0, double X[],
                               double E[], randompack_rng *rng)
{
  varmapack_error error = VARMAPACK_OK;
  double *Eall = 0, *Aflp = 0, *Bflp = 0;
  int zlag, t0, firstActiveShock, firstShockTime, nE, ldE, rn;
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || (s > 0 && C == 0) ||
      (s > 0 && z == 0) || Sig == 0 || X0 == 0 || X == 0 || rng == 0) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if (p < 0 || q < 0 || s < 0 || r <= 0 || n <= 0 || M <= 0 ||
      (s > 0 && Mz != 1 && Mz != M) || (MX0 != 1 && MX0 != M)) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  zlag = s > 0 ? s - 1 : 0;
  t0 = imax(p, zlag);
  if (h <= t0 || n < h) return VARMAPACK_INVALID_ARGUMENT;
  firstActiveShock = t0 - q;
  firstShockTime = imin(firstActiveShock, 0);
  nE = n - firstShockTime;
  ldE = r*nE;
  rn = r*n;
  if (!ALLOC(Eall, ldE*M)) goto alloc_fail;
  if (p > 0 && !ALLOC(Aflp, r*r*p)) goto alloc_fail;
  if (q > 0 && !ALLOC(Bflp, r*r*q)) goto alloc_fail;
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  for (int j=0; j<M; j++) {
    double *X0j = X0 + (MX0 == 1 ? 0 : j*r*h);
    copy(r*h, X0j, 1, X + j*rn, 1);
  }
  error = drawStartupShocks(A, B, C, Sig, z, Mz, X0, MX0, Eall, p, q, s, r, n, h, t0,
                            firstShockTime, firstActiveShock, ldE, M, rng);
  if (error) goto fail;
  if (!forwardSimx(Aflp, Bflp, C, Sig, z, Mz, Eall, ldE, X, p, q, s, r, n, M, h,
                   firstShockTime, rng)) {
    goto alloc_fail;
  }
  if (E != 0) {
    for (int j=0; j<M; j++) {
      lacpy("All", rn, 1, Eall + j*ldE - firstShockTime*r, ldE, E + j*rn, rn);
    }
  }
  FREE(Bflp); FREE(Aflp); FREE(Eall);
  return VARMAPACK_OK;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
fail:
  FREE(Bflp); FREE(Aflp); FREE(Eall);
  return error;
}

static void buildH(double B[], int q, int r, int h, int t0, int firstShockTime,
                   double H[]) {
  int m = h - t0;
  int nE0 = h - firstShockTime;
  int rm = r*m;
  laset("All", rm, r*nE0, 0, 0, H, rm);
  for (int t=t0; t<h; t++) {
    int row = r*(t - t0);
    for (int ell=firstShockTime; ell<h; ell++) {
      int lag = t - ell;
      int col = r*(ell - firstShockTime);
      double *Hij = H + col*rm + row;
      if (lag == 0) {
        for (int k=0; k<r; k++) Hij[k + k*rm] = 1;
      }
      else if (1 <= lag && lag <= q) {
        lacpy("All", r, r, B + (lag - 1)*r*r, r, Hij, rm);
      }
    }
  }
}

static void startupResidual(double A[], double C[], double z[], double X0[], int p,
                            int s, int r, int h, int t0, double residual[]) {
  for (int t=t0; t<h; t++) {
    double *rt = residual + r*(t - t0);
    copy(r, X0 + r*t, 1, rt, 1);
    for (int i=1; i<=p; i++) {
      gemv("NoT", r, r, -1, A + (i-1)*r*r, r, X0 + r*(t-i), 1, 1, rt, 1);
    }
    for (int i=1; i<=s; i++) {
      axpy(r, -z[t-i+1], C + (i-1)*r, 1, rt, 1);
    }
  }
}

static varmapack_error drawStartupShocks(double A[], double B[], double C[], double Sig[],
                                         double z[], int Mz, double X0[], int MX0,
                                         double Eall[], int p, int q, int s, int r, int n,
                                         int h, int t0, int firstShockTime,
                                         int firstActiveShock, int ldE, int M,
                                         randompack_rng *rng) {
  int info;
  int m = h - t0;
  int nPrefix = imax(firstActiveShock, 0);
  int activeOffset = r*(firstActiveShock - firstShockTime);
  int nActive = h - firstActiveShock;
  int rm = r*m;
  int re = r*nActive;
  varmapack_error error = VARMAPACK_OK;
  double *H = 0, *HD = 0, *Wlag = 0, *W = 0, *K = 0, *R = 0, *work = 0;
  double *residual = 0, *Ehat = 0;
  for (int j=0; j<M; j++) {
    if (nPrefix > 0 &&
        !randompack_mvn("T", 0, Sig, r, nPrefix, Eall + j*ldE, r, 0, rng)) {
      goto alloc_fail;
    }
  }
  if (q == 0) {
    if (!ALLOC(residual, rm)) goto alloc_fail;
    for (int j=0; j<M; j++) {
      double *zj = z + (Mz == 1 ? 0 : j*n);
      double *X0j = X0 + (MX0 == 1 ? 0 : j*r*h);
      startupResidual(A, C, zj, X0j, p, s, r, h, t0, residual);
      copy(rm, residual, 1, Eall + j*ldE + activeOffset, 1);
    }
    FREE(residual);
    return VARMAPACK_OK;
  }
  if (!ALLOC(H, rm*re)) goto alloc_fail;
  if (!ALLOC(HD, rm*re)) goto alloc_fail;
  if (!ALLOC(Wlag, r*r*(q+1))) goto alloc_fail;
  if (!ALLOC(W, rm*rm)) goto alloc_fail;
  if (!ALLOC(K, rm*re)) goto alloc_fail;
  if (!ALLOC(R, re*re)) goto alloc_fail;
  if (!ALLOC(work, imax(rm, re)*r)) goto alloc_fail;
  if (!ALLOC(residual, rm)) goto alloc_fail;
  if (!ALLOC(Ehat, re)) goto alloc_fail;
  buildH(B, q, r, h, t0, firstActiveShock, H);
  lacpy("All", rm, re, H, rm, HD, rm);
  postmultiplySigmaPrime(HD, rm, rm, nActive, Sig, r, work);
  if (!FindW(B, Sig, q, r, Wlag)) goto alloc_fail;
  WBuild(Wlag, q, r, m, W);
  potrf("Low", rm, W, rm, &info);
  if (info > 0) { error = VARMAPACK_SINGULAR; goto fail; }
  if (info < 0) { error = VARMAPACK_INTERNAL; goto fail; }
  lacpy("All", rm, re, HD, rm, K, rm);
  trsm("Left", "Low", "NoT", "NonUnit", rm, re, 1, W, rm, K, rm);
  trsm("Left", "Low", "Trans", "NonUnit", rm, re, 1, W, rm, K, rm);
  laset("All", re, re, 0, 0, R, re);
  for (int j=0; j<nActive; j++) lacpy("All", r, r, Sig, r, R + j*r*(re + 1), re);
  gemm("Trans", "NoT", re, re, rm, -1, HD, rm, K, rm, 1, R, re);
  for (int j=0; j<M; j++) {
    double *zj = z + (Mz == 1 ? 0 : j*n);
    double *X0j = X0 + (MX0 == 1 ? 0 : j*r*h);
    startupResidual(A, C, zj, X0j, p, s, r, h, t0, residual);
    trsv("Low", "NoT", "NonUnit", rm, W, rm, residual, 1);
    trsv("Low", "Trans", "NonUnit", rm, W, rm, residual, 1);
    gemv("Trans", rm, re, 1, HD, rm, residual, 1, 0, Ehat, 1);
    if (!randompack_mvn("T", Ehat, R, re, 1, Eall + j*ldE + activeOffset, ldE,
                        0, rng)) {
      goto alloc_fail;
    }
  }
  goto fail;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
fail:
  FREE(Ehat); FREE(residual); FREE(work); FREE(R); FREE(K); FREE(W); FREE(Wlag);
  FREE(HD); FREE(H);
  return error;
}

static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int r, int n, int M, int h,
                        int firstShockTime, randompack_rng *rng) {
  int rn = r*n;
  if (n > h) {
    for (int j=0; j<M; j++) {
      double *Eh = Eall + j*ldE + r*(h - firstShockTime);
      if (!randompack_mvn("T", 0, Sig, r, n-h, Eh, r, 0, rng)) return false;
    }
  }
  for (int t=h; t<n; t++) {
    int iX = r*t;
    int iA = r*(t - p);
    int iB = r*(t - q - firstShockTime);
    lacpy("All", r, M, Eall + r*(t - firstShockTime), ldE, X + iX, rn);
    if (p > 0) {
      gemm("NoT", "NoT", r, M, r*p, 1, Aflp, r, X + iA, rn, 1, X + iX, rn);
    }
    if (q > 0) {
      gemm("NoT", "NoT", r, M, r*q, 1, Bflp, r, Eall + iB, ldE, 1, X + iX, rn);
    }
    for (int i=1; i<=s; i++) {
      for (int j=0; j<M; j++) {
        double *zj = z + (Mz == 1 ? 0 : j*n);
        axpy(r, zj[t-i+1], C + (i-1)*r, 1, X + j*rn + iX, 1);
      }
    }
  }
  return true;
}
