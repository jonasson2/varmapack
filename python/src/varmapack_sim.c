// varmapack_sim — simulate burn-in-free AR, VAR, ARMA or VARMA time series
//
// See varmapack.h for parameter descriptions and storage conventions
//
// GENERAL VARMA TIME SERIES MODEL:
//   The following is called a VARMA(p,q) model:
//
//             x(t) - mu(t) = A1·(x(t-1) - mu(t-1)) + ... +
//                            Ap·(x(t-p) - mu(t-p)) + y(t)
//   where:
//             y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
//
//   x(t) is r-dimensional and eps(t) is r-variate normal with mean 0 and
//   covariance Sig. The eps(t) are called shocks and they are uncorrelated in
//   time (i.e. eps(t) and eps(s) are independent if t and s differ). VARMA
//   stands for vector-autoregressive-moving-average.
//
// SPECIAL CASES, ARMA, VAR and AR:
//   If r=1 the model above is an ARMA model, if q=0 and r>1 it is a VAR model
//   and if q=0 and r=1 it is an AR model. Thus VAR and AR models may be
//   written:
//            x(t) - mu(t) = A1·(x(t-1) - mu(t-1)) + ... +
//                           Ap·(x(t-p) - mu(t-p)) + eps(t)
//
//   where x(t) and eps(t) are scalar for an AR series and vectors for a VAR
//   series.
//
// ABOUT THE PROGRAM:
//   The program is based on the Matlab program varmapack_sim described in
//   references [1], [2] and [3] (the report [3] contains some details not found
//   in the published papers), but uses the C Randompack interface for random
//   number generation and may use either VYW or Lyapunov setup internally.
//
// RANDOM NUMBER GENERATION
//   The random number generator used for the shocks is specified by the rng
//   parameter. Random number generation is provided by the external Randompack
//   library. See GitHub.com/jonasson2/randompack, and [4] for details.
//
// MEANS
//   If mu is zero or nmu is zero, the time series has zero mean. Otherwise mu
//   gives r by nmu time series means: E[x(t)] = mu(:,t) for t < nmu, and
//   E[x(t)] = mu(:,nmu-1) for t >= nmu.
//
// STARTING VALUES:
//   When starting values are specified in X0, each generated series starts with
//   X0 and subsequent terms are simulated conditionally on X0. Without X0, the
//   first h x(t) and eps(t) values, h = max(p,q), are drawn from their
//   stationary joint distribution, so there is no initial burn-in segment to
//   throw away. With X0 and a stationary model, h = max(p,q,nX0), and the first
//   h eps(t) values are drawn from the conditional distribution of the startup
//   shocks given the supplied x values. With X0 and a nonstationary pure AR
//   model, no latent startup shocks are needed; future shocks are drawn and the
//   recurrence is run forward from the supplied X0.
//
// MATRIX STORAGE SCHEME
//   All matrices are stored in column-major order, so that for example the
//   elements of A1, A2... are stored in the order A1(1,1),
//   A1(2,1),...,A1(r,r),..., Ap(r,r).
//
// STATIONARITY
//   Stationary models can be simulated with or without X0. Nonstationary pure
//   AR models can be simulated with X0, because no latent startup shocks are
//   needed for forward recursion. Models with MA terms and specified X0 must be
//   stationary.
//
// FAILURES
//   If the series is nonstationary without X0, or another recoverable error is
//   detected, the function returns a nonzero varmapack_error code.
//
// REFERENCES:
//   [1] K Jonasson and SE Ferrando 2008. Evaluating exact VARMA likelihood
//       and its gradient when data are incomplete. ACM Trans. Math. Softw.
//       35, No 1.
//   [2] K Jonasson 2008. Algorithm 878: Exact VARMA likelihood and its
//       gradient for complete and incomplete data with Matlab. ACM Trans.
//       Math. Softw. 35, No 1.
//   [3] K Jonasson and SE Ferrando 2006. Efficient likelihood evaluation for
//       VARMA processes with missing values. Report VHI-01-2006, Engineering
//       Research Institute, University of Iceland.
//   [4] K Jonasson 2026. Randompack: Cross-Platform Reproducible Random Number
//       Generation and Distribution Sampling, arXiv:2605.05099. Submitted to
//       ACM Trans. Math. Softw.

#include <stdbool.h>
#include "error.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "randompack.h"
#include "varmapack.h"

static void addMean(double X[], double mu[], int nmu, int r, int n, int M);
static void subtractMean(double X[], double mu[], int nmu, int r, int n, int M, int ldX);
static varmapack_error buildStartupCovar(double A[], double B[], double Sig[], int p,
                                         int q, int r, int h, double C[], double SS[]);
static varmapack_error startUnconditional(double A[], double B[], double Sig[],
                                          double SS[], double E[], int ldE,
                                          double X[], int p, int q, int r, int h,
                                          int n, int M, randompack_rng *rng);
static varmapack_error startFromX0(double A[], double Sig[], double X0[], double mu[],
                                   int nmu, double C[], double SS[], double E[],
                                   int ldE, bool returnE, bool stationary, double X[],
                                   int p, int q, int r, int h, int n, int M, int MX0,
                                   randompack_rng *rng);
static bool tailSimulate(double Aflp[], double Bflp[], double Sig[], double E[], bool
                         rollingE, int ldE, double X[], int p, int q, int r,
                         int n, int M, int h, randompack_rng *rng);

varmapack_error varmapack_sim(double A[], double B[], double Sig[], double mu[],
                              int nmu, int p, int q, int r, int n, int M, double X0[],
                              int nX0, int MX0, double X[], double E[],
                              randompack_rng *rng)
{
  varmapack_error error = VARMAPACK_OK;
  double *C = 0, *SS = 0;
  double *Aflp = 0, *Bflp = 0;
  bool returnE = E != 0;
  bool X0Given = X0 != 0;
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || Sig == 0 || X == 0 || rng == 0) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if (p < 0 || q < 0 || r <= 0 || n <= 0 || M <= 0 || n < imax(p, q)) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if (nmu < 0 || nmu > n || (mu == 0 && nmu > 0)) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if ((!X0Given && nX0 != 0) || (X0Given && (nX0 < imax(p, q) || nX0 > n)) ||
      (X0Given && MX0 != 1 && MX0 != M)) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  int h = imax(imax(p,q), nX0);
  int rn = r*n;  // Total number of observations
  int rh = r*h;
  int ldE = returnE ? rn : r*(h + 1);
  bool stationary = p == 0 || varmapack_specrad(A, r, p) < 1;
  if (!X0Given && !stationary)
    return VARMAPACK_NONSTATIONARY;
  if (X0Given && !stationary && q > 0)
    return VARMAPACK_NONSTATIONARY_MA;
  if (!X0Given && p == 0 && q == 0) { // White noise can be drawn directly into X.
    if (!randompack_mvn("T", 0, Sig, r, n*M, X, r, 0, rng)) goto alloc_fail;
    if (returnE) lacpy("All", rn, M, X, rn, E, rn); // and copy it to E
    addMean(X, mu, nmu, r, n, M);
    return VARMAPACK_OK;
  }
  if (!returnE)
    if (!ALLOC(E, ldE*M)) goto alloc_fail;
  if (!X0Given || stationary) { // compute C and SS
    if (!ALLOC(C, r*r*(q+1))) goto alloc_fail;
    if (rh > 0)
      if (!ALLOC(SS, rh*rh)) goto alloc_fail;
    error = buildStartupCovar(A, B, Sig, p, q, r, h, C, SS);
    if (error) goto fail;
  }
  if (X0Given) {
    error = startFromX0(A, Sig, X0, mu, nmu, C, SS, E, ldE, returnE, stationary, X, p, q,
                        r, h, n, M, MX0, rng);
    if (error) goto fail;
  }
  else { // Unconditional exact startup: draw X[:h] given E[:h].
    error = startUnconditional(A, B, Sig, SS, E, ldE, X, p, q, r, h, n, M, rng);
    if (error) goto fail;
  }
  FREE(C); FREE(SS);
  if (p > 0 && !ALLOC(Aflp, r*r*p)) goto alloc_fail;
  if (q > 0 && !ALLOC(Bflp, r*r*q)) goto alloc_fail;
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  bool rollingE = !returnE;
  if (!tailSimulate(Aflp, Bflp, Sig, E, rollingE, ldE, X, p, q, r, n, M, h, rng))
    goto alloc_fail;
  addMean(X, mu, nmu, r, n, M);
  FREE(Bflp); FREE(Aflp);
  if (!returnE) FREE(E);
  return VARMAPACK_OK;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
fail:
  FREE(Bflp); FREE(Aflp); FREE(SS); FREE(C);
  if (!returnE) FREE(E);
  return error;
}

static void addMean(double X[], double mu[], int nmu, int r, int n, int M) {
  int rn = r*n;
  if (mu == 0 || nmu == 0) return;
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      int imu = t < nmu ? t : nmu - 1;
      axpy(r, 1.0, mu + imu*r, 1, X + j*rn + t*r, 1);
    }
  }
}

static void subtractMean(double X[], double mu[], int nmu, int r, int n, int M, int ldX) {
  if (mu == 0 || nmu == 0) return;
  for (int j=0; j<M; j++) {
    for (int t=0; t<n; t++) {
      int imu = t < nmu ? t : nmu - 1;
      axpy(r, -1.0, mu + imu*r, 1, X + j*ldX + t*r, 1);
    }
  }
}

static varmapack_error buildStartupCovar(double A[], double B[], double Sig[], int p,
                                         int q, int r, int h, double C[], double SS[]) {
  double *G = 0;
  double *S = 0;
  if (!ALLOC(G, r*r*(q+1))) goto alloc_fail;
  if (!ALLOC(S, r*r*(p+1))) goto alloc_fail;
  if (!FindS(A, B, Sig, p, q, r, S, C, G)) {
    FREE(S); FREE(G);
    return VARMAPACK_INTERNAL;
  }
  if (!SBuild("Low", S, A, G, p, q, r, h, SS)) goto alloc_fail;
  FREE(S); FREE(G);
  return VARMAPACK_OK;
alloc_fail:
  FREE(S); FREE(G);
  return VARMAPACK_ALLOCATION;
}

static varmapack_error startUnconditional(double A[], double B[], double Sig[], double
                                          SS[], double E[], int ldE, double X[], int p,
                                          int q, int r, int h, int n, int M,
                                          randompack_rng *rng) {
  // Simulate M replicates X[:m] and E[:m] using joint distribution of (X,E)
  int rh = r*h;
  int rn = r*n;
  varmapack_error error = VARMAPACK_OK;
  double *Wrk = 0;
  double *Psi = 0;
  double *PsiHat = 0;
  double *R = 0;
  if (!ALLOC(Wrk, rh*M)) goto alloc_fail;
  if (!ALLOC(Psi, rh*rh)) goto alloc_fail;
  if (!ALLOC(PsiHat, rh*rh)) goto alloc_fail;
  if (!ALLOC(R, rh*rh)) goto alloc_fail;
  for (int j=0; j<M; j++) {
    if (!randompack_mvn("T", 0, Sig, r, h, E + j*ldE, r, 0, rng)) {
      goto alloc_fail;
    }
  }
  FindPsi(A, B, Psi, p, q, r);
  FindPsiHat(Psi, PsiHat, Sig, r, h);
  lacpy("Low", rh, rh, SS, rh, R, rh);
  syrk("Low", "NoT", rh, rh, -1.0, PsiHat, rh, 1.0, R, rh);
  if (!randompack_mvn("T", 0, R, rh, M, Wrk, rh, 0, rng)) {
    goto alloc_fail;
  }
  lacpy("All", rh, M, Wrk, rh, X, rn);
  gemm("NoT", "NoT", rh, M, rh, 1.0, Psi, rh, E, ldE, 1.0, X, rn);
  goto done;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
done:
  FREE(R); FREE(PsiHat); FREE(Psi); FREE(Wrk);
  return error;
}

static varmapack_error startFromX0(double A[], double Sig[], double X0[], double mu[],
                                   int nmu, double C[], double SS[], double E[], int ldE,
                                   bool returnE, bool stationary, double X[], int p, int
                                   q, int r, int h, int n, int M, int MX0,
                                   randompack_rng *rng) {
  // Simulate M replicates from supplied X0.
  int info;
  int rh = r*h;
  int rn = r*n;
  varmapack_error error = VARMAPACK_OK;
  double *CC = 0;
  double *x0bar = 0;
  double *e = 0;
  double *wrk = 0;
  double *R = 0;
  double *Chat;
  for (int j=0; j<M; j++) {
    double *X0j = X0 + (MX0 == 1 ? 0 : j*rh);
    copy(rh, X0j, 1, X + j*rn, 1);
    subtractMean(X + j*rn, mu, nmu, r, h, 1, rn);
  }
  if (!stationary) {
    if (returnE) {
      for (int j=0; j<M; j++) setzero(rh, E + j*ldE);
    }
    FREE(x0bar);
    return VARMAPACK_OK;
  }
  if (!ALLOC(CC, rh*rh)) goto alloc_fail;
  if (!ALLOC(x0bar, rh)) goto alloc_fail;
  if (!ALLOC(e, rh)) goto alloc_fail;
  if (!ALLOC(wrk, rh)) goto alloc_fail;
  if (!ALLOC(R, rh*rh)) goto alloc_fail;
  potrf("Low", rh, SS, rh, &info);
  if (info > 0) { error = VARMAPACK_NOT_POSITIVE_SEMIDEFINITE; goto fail; }
  if (info < 0) { error = VARMAPACK_INTERNAL; goto fail; }
  CCBuild(A, C, p, q, r, h, CC);
  Chat = CC;
  trsm("Left", "Low", "NT", "NotUD", rh, rh, 1.0, SS, rh, Chat, rh);
  for (int j=0; j<h; j++) {
    lacpy("Low", r, r, Sig, r, R + j*r*(rh + 1), rh);
  }
  syrk("L", "T", rh, rh, -1.0, Chat, rh, 1.0, R, rh);
  for (int j=0; j<M; j++) {
    copy(rh, X + j*rn, 1, x0bar, 1);
    copy(rh, x0bar, 1, wrk, 1);
    trsv("Lo", "NT", "NotUD", rh, SS, rh, wrk, 1);
    gemv("T", rh, rh, 1.0, Chat, rh, wrk, 1, 0.0, e, 1);
    if (!randompack_mvn("T", e, R, rh, 1, E + j*ldE, ldE, 0, rng)) {
      goto alloc_fail;
    }
  }
  FREE(R); FREE(wrk); FREE(e); FREE(x0bar); FREE(CC);
  return VARMAPACK_OK;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
fail:
  FREE(R); FREE(wrk); FREE(e); FREE(x0bar); FREE(CC);
  return error;
}

static bool tailSimulate(double Aflp[], double Bflp[], double Sig[], double E[], bool
                         rollingE, int ldE, double X[], int p, int q, int r,
                         int n, int M, int h, randompack_rng *rng) {
  // Simulate X[h:] with forward computation from X[:h] and E[:h]
  int rn = r*n;
  int rh = r*h;
  if (!rollingE) {
    for (int j=0; j<M; j++) {
      if (!randompack_mvn("T", 0, Sig, r, n - h, E + j*ldE + rh, r, 0, rng)) {
        return false;
      }
    }
    lacpy("All", (n-h)*r, M, E + rh, rn, X + rh, rn);
  }
  for (int t=h; t<n; t++) {
    int iX = r*t;
    int iA = r*(t - p);
    int iB = r*(t - q);
    if (rollingE) {
      int iE = r*(t%(h + 1));
      if (!randompack_mvn("T", 0, Sig, r, M, E + iE, ldE, 0, rng)) {
        return false;
      }
      lacpy("All", r, M, E + iE, ldE, X + iX, rn);
    }
    if (p > 0) // Apply AR
      gemm("NoT", "NoT", r, M, r*p, 1.0, Aflp, r, X + iA, rn, 1.0, X + iX, rn);
    if (q > 0 && rollingE) { // rolling MA to avoid allocating full E
      int firstSlot = (t - q)%(h + 1);
      int firstCount = h + 1 - firstSlot;
      if (firstCount > q) firstCount = q;
      gemm("NoT", "NoT", r, M, r*firstCount, 1.0, Bflp, r,
           E + r*firstSlot, ldE, 1.0, X + iX, rn);
      if (firstCount < q) {
        int secondCount = q - firstCount;
        gemm("NoT", "NoT", r, M, r*secondCount, 1.0,
             Bflp + firstCount*r*r, r, E, ldE, 1.0, X + iX, rn);
      }
    }
    else if (q > 0) // Using full E
      gemm("NoT", "NoT", r, M, r*q, 1.0, Bflp, r, E + iB, rn, 1.0, X + iX, rn);
  }
  return true;
}
