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
static bool drawStartupShocks(double A[], double B[], double C[], double Sig[],
                              double z[], int Mz, double X0[], int MX0, double Eall[],
                              int p, int q, int s, int r, int n, int h, int t0,
                              int firstShockTime, int firstActiveShock, int ldE,
                              int M, randompack_rng *rng);
static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int r, int n, int M, int h, int rn,
                        int firstShockTime, randompack_rng *rng);
static bool validSimxSizes(int p, int q, int s, int r, int n, int M, int h,
                           int *t0, int *firstActiveShock, int *firstShockTime,
                           int *rh, int *rn, int *ldE, size_t *eCount,
                           size_t *aCount, size_t *bCount);
static bool requirePositiveDefinite(double Sig[], int r, size_t count);

bool varmapack_simx(double A[], double B[], double C[], double Sig[], double z[],
                    int Mz, int p, int q, int s, int r, int n, int M, double X0[],
                    int h, int MX0, double X[], double E[], randompack_rng *rng)
{
  double *Eall = 0, *Aflp = 0, *Bflp = 0;
  int t0, firstActiveShock, firstShockTime, rh, ldE, rn;
  size_t eCount, aCount, bCount, sigCount;
  clear_error();
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || (s > 0 && C == 0) ||
      (s > 0 && z == 0) || Sig == 0 || X0 == 0 || X == 0 || rng == 0) {
    return fail_error("invalid argument");
  }
  if (p < 0 || q < 0 || s < 0 || r <= 0 || n <= 0 || M <= 0 ||
      (Mz != 1 && Mz != M) || (MX0 != 1 && MX0 != M)) {
    return fail_error("invalid argument");
  }
  if (h <= imax(p, s > 0 ? s - 1 : 0) || n < h) {
    return fail_error("invalid argument");
  }
  if (!validSimxSizes(p, q, s, r, n, M, h, &t0, &firstActiveShock,
                      &firstShockTime, &rh, &rn, &ldE, &eCount, &aCount, &bCount)) {
    return fail_error("problem size too large");
  }
  if (!sizeProduct((size_t)r, (size_t)r, &sigCount))
    return fail_error("problem size too large");
  if (!requirePositiveDefinite(Sig, r, sigCount)) return false;
  if (!ALLOC(Eall, eCount)) goto alloc_fail;
  if (p > 0 && !ALLOC(Aflp, aCount)) goto alloc_fail;
  if (q > 0 && !ALLOC(Bflp, bCount)) goto alloc_fail;
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  for (int j=0; j<M; j++) {
    double *X0j = X0 + (MX0 == 1 ? 0 : (size_t)j*rh);
    copy(rh, X0j, 1, X + (size_t)j*rn, 1);
  }
  if (!drawStartupShocks(A, B, C, Sig, z, Mz, X0, MX0, Eall, p, q, s, r, n, h, t0,
                         firstShockTime, firstActiveShock, ldE, M, rng)) {
    goto fail;
  }
  if (!forwardSimx(Aflp, Bflp, C, Sig, z, Mz, Eall, ldE, X, p, q, s, r, n, M, h,
                   rn, firstShockTime, rng)) {
    goto fail;
  }
  if (E != 0) {
    for (int j=0; j<M; j++) {
      lacpy("All", rn, 1, Eall + (size_t)j*ldE + (size_t)(-firstShockTime)*r, ldE,
            E + (size_t)j*rn, rn);
    }
  }
  FREE(Bflp); FREE(Aflp); FREE(Eall);
  return true;
alloc_fail:
  fail_error("allocation failed");
fail:
  FREE(Bflp); FREE(Aflp); FREE(Eall);
  return false;
}

static bool validSimxSizes(int p, int q, int s, int r, int n, int M, int h,
                           int *t0, int *firstActiveShock, int *firstShockTime,
                           int *rh, int *rn, int *ldE, size_t *eCount,
                           size_t *aCount, size_t *bCount) {
  int r2, rp, rq, m, nActive, rm, re, nE;
  size_t hCount, wlagCount, wCount, rCount, workCount;
  *t0 = imax(p, s > 0 ? s - 1 : 0);
  if (q == INT_MAX) return false;
  *firstActiveShock = *t0 - q;
  *firstShockTime = imin(*firstActiveShock, 0);
  if (*firstShockTime < n - INT_MAX) return false;
  nE = n - *firstShockTime;
  m = h - *t0;
  nActive = h - *firstActiveShock;
  if (!intProduct(r, h, rh) || !intProduct(r, n, rn) || !intProduct(r, nE, ldE) ||
      !intProduct(r, r, &r2) || !intProduct(r, p, &rp) || !intProduct(r, q, &rq) ||
      !intProduct(r, m, &rm) || !intProduct(r, nActive, &re)) return false;
  if (!sizeProduct((size_t)*ldE, (size_t)M, eCount) ||
      !sizeProduct((size_t)r2, (size_t)p, aCount) ||
      !sizeProduct((size_t)r2, (size_t)q, bCount)) return false;
  if (q > 0 &&
      (!sizeProduct((size_t)rm, (size_t)re, &hCount) ||
       !sizeProduct((size_t)r2, (size_t)(q + 1), &wlagCount) ||
       !sizeProduct((size_t)rm, (size_t)rm, &wCount) ||
       !sizeProduct((size_t)re, (size_t)re, &rCount) ||
       !sizeProduct((size_t)imax(rm, re), (size_t)r, &workCount))) return false;
  return true;
}

static bool requirePositiveDefinite(double Sig[], int r, size_t count) {
  double *L = 0;
  bool triangular;
  if (!ALLOC(L, count)) return fail_error("allocation failed");
  if (!psdFactor(Sig, r, L, &triangular)) {
    FREE(L);
    return false;
  }
  FREE(L);
  if (!triangular) {
    return fail_error(
      "cannot condition on X0: shock covariance is not positive definite");
  }
  return true;
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
      double *Hij = H + (size_t)col*rm + row;
      if (lag == 0) {
        for (int k=0; k<r; k++) Hij[k + (size_t)k*rm] = 1;
      }
      else if (1 <= lag && lag <= q) {
        lacpy("All", r, r, B + (size_t)(lag - 1)*r*r, r, Hij, rm);
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
      gemv("NoT", r, r, -1, A + (size_t)(i-1)*r*r, r, X0 + r*(t-i), 1, 1, rt, 1);
    }
    for (int i=1; i<=s; i++) {
      axpy(r, -z[t-i+1], C + (size_t)(i-1)*r, 1, rt, 1);
    }
  }
}

static bool drawStartupShocks(double A[], double B[], double C[], double Sig[],
                              double z[], int Mz, double X0[], int MX0, double Eall[],
                              int p, int q, int s, int r, int n, int h, int t0,
                              int firstShockTime, int firstActiveShock, int ldE,
                              int M, randompack_rng *rng) {
  int info;
  int m = h - t0;
  int nPrefix = imax(firstActiveShock, 0);
  int activeOffset = r*(firstActiveShock - firstShockTime);
  int nActive = h - firstActiveShock;
  int rm = r*m;
  int re = r*nActive;
  bool ok = true;
  double *H = 0, *HD = 0, *Wlag = 0, *W = 0, *K = 0, *R = 0, *work = 0;
  double *residual = 0, *Ehat = 0;
  for (int j=0; j<M; j++) {
    if (nPrefix > 0 &&
        !randompack_mvn("T", 0, Sig, r, nPrefix, Eall + (size_t)j*ldE, r, 0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto done;
    }
  }
  if (q == 0) {
    if (!ALLOC(residual, rm)) goto alloc_fail;
    for (int j=0; j<M; j++) {
      double *zj = s > 0 ? z + (Mz == 1 ? 0 : (size_t)j*n) : 0;
      double *X0j = X0 + (MX0 == 1 ? 0 : (size_t)j*r*h);
      startupResidual(A, C, zj, X0j, p, s, r, h, t0, residual);
      copy(rm, residual, 1, Eall + (size_t)j*ldE + activeOffset, 1);
    }
    FREE(residual);
    return true;
  }
  size_t hCount, wlagCount, wCount, rCount, workCount;
  if (!sizeProduct((size_t)rm, (size_t)re, &hCount) ||
      !sizeProduct((size_t)r, (size_t)r, &wlagCount) ||
      !sizeProduct(wlagCount, (size_t)(q + 1), &wlagCount) ||
      !sizeProduct((size_t)rm, (size_t)rm, &wCount) ||
      !sizeProduct((size_t)re, (size_t)re, &rCount) ||
      !sizeProduct((size_t)imax(rm, re), (size_t)r, &workCount)) {
    return fail_error("problem size too large");
  }
  if (!ALLOC(H, hCount)) goto alloc_fail;
  if (!ALLOC(HD, hCount)) goto alloc_fail;
  if (!ALLOC(Wlag, wlagCount)) goto alloc_fail;
  if (!ALLOC(W, wCount)) goto alloc_fail;
  if (!ALLOC(K, hCount)) goto alloc_fail;
  if (!ALLOC(R, rCount)) goto alloc_fail;
  if (!ALLOC(work, workCount)) goto alloc_fail;
  if (!ALLOC(residual, rm)) goto alloc_fail;
  if (!ALLOC(Ehat, re)) goto alloc_fail;
  buildH(B, q, r, h, t0, firstActiveShock, H);
  lacpy("All", rm, re, H, rm, HD, rm);
  postmultiplySigmaPrime(HD, rm, rm, nActive, Sig, r, work);
  if (!FindW(B, Sig, q, r, Wlag)) {
    ok = varmapack_last_error() ? false : fail_error("allocation failed");
    goto done;
  }
  WBuild(Wlag, q, r, m, W);
  potrf("Low", rm, W, rm, &info);
  if (info > 0) { ok = fail_error("singular matrix"); goto done; }
  if (info < 0) { ok = fail_error("internal error"); goto done; }
  lacpy("All", rm, re, HD, rm, K, rm);
  trsm("Left", "Low", "NoT", "NonUnit", rm, re, 1, W, rm, K, rm);
  trsm("Left", "Low", "Trans", "NonUnit", rm, re, 1, W, rm, K, rm);
  laset("All", re, re, 0, 0, R, re);
  for (int j=0; j<nActive; j++)
    lacpy("All", r, r, Sig, r, R + (size_t)j*r*((size_t)re + 1), re);
  gemm("Trans", "NoT", re, re, rm, -1, HD, rm, K, rm, 1, R, re);
  for (int j=0; j<M; j++) {
    double *zj = s > 0 ? z + (Mz == 1 ? 0 : (size_t)j*n) : 0;
    double *X0j = X0 + (MX0 == 1 ? 0 : (size_t)j*r*h);
    startupResidual(A, C, zj, X0j, p, s, r, h, t0, residual);
    trsv("Low", "NoT", "NonUnit", rm, W, rm, residual, 1);
    trsv("Low", "Trans", "NonUnit", rm, W, rm, residual, 1);
    gemv("Trans", rm, re, 1, HD, rm, residual, 1, 0, Ehat, 1);
    if (!randompack_mvn("T", Ehat, R, re, 1,
                        Eall + (size_t)j*ldE + activeOffset, ldE,
                        0, rng)) {
      ok = fail_error(randompack_last_error(rng));
      goto done;
    }
  }
  goto done;
alloc_fail:
  ok = fail_error("allocation failed");
done:
  FREE(Ehat); FREE(residual); FREE(work); FREE(R); FREE(K); FREE(W); FREE(Wlag);
  FREE(HD); FREE(H);
  return ok;
}

static bool forwardSimx(double Aflp[], double Bflp[], double C[], double Sig[],
                        double z[], int Mz, double Eall[], int ldE, double X[], int p,
                        int q, int s, int r, int n, int M, int h, int rn,
                        int firstShockTime, randompack_rng *rng) {
  if (n > h) {
    for (int t=h; t<n; t++) {
      int iE = r*(t - firstShockTime);
      if (!randompack_mvn("T", 0, Sig, r, M, Eall + iE, ldE, 0, rng)) {
        return fail_error(randompack_last_error(rng));
      }
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
        double *zj = z + (Mz == 1 ? 0 : (size_t)j*n);
        axpy(r, zj[t-i+1], C + (size_t)(i-1)*r, 1,
             X + (size_t)j*rn + iX, 1);
      }
    }
  }
  return true;
}
