// varmapack_sim — simulate spin-up FREE AR, VAR, ARMA or VARMA time series
//
// See varmapack.h for parameter descriptions and storage conventions
//
// GENERAL VARMA TIME SERIES MODEL:
//   The following is called a VARMA(p,q) model:
//
//             x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)
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
//            x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + eps(t)
//
//   where x(t) and eps(t) are scalar for an AR series and vectors for a VAR
//   series.
//
// ABOUT THE PROGRAM:
//   The program is a translation of the Matlab program varmapack_sim described in
//   references [1], [2] and [3] (the report [3] contains some details not found
//   in the published papers).
//
// RANDOM NUMBER GENERATION
//   The random number generator used for the shocks is specified by the rng
//   parameter. Random number generation is provided by the external randompack
//   library. See randompack.h for details.
// 
// STARTING VALUES:
//   When starting values for the simulation are specified in X0, each generated
//   series will start with X0, and subsequent terms are simulated using X0 as
//   input. An important feature is that the simulated series have an accurate
//   distribution from term one (they are spin-up FREE). Thus there is no need
//   to throw a way an initial segent of the simulated series. If X0 is not
//   used, the first h x(t) and eps(t) values (with h = max(p,q)) are drawn from
//   the correct joint distribution of these variables, and if X0 is used the
//   first h eps(t) values are drawn from the correct conditional distribution
//   of eps(1)...eps(h)|x(1)...x(h).
//
// MATRIX STORAGE SCHEME
//   All matrices are stored in Fortran fashion (i.e. column major ordering), so
//   that for example the elements of A1, A2... are stored in the order A1(1,1),
//   A1(2,1),...,A1(r,r),..., Ap(r,r).
//
// STATIONARITY
//   If X0 is not specified, the parameter matrices Ai and Bi must represent a
//   stationary model, but if X0 is specified the series is allowed to be
//   non-stationary.
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

#include <stdbool.h>
#include "error.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "randompack.h"
#include "varmapack.h"
#include "VYW.h"
#include "printX.h"

static void CCBuild(double A[], double C[], int p, int q, int r, int n, double CC[]);

varmapack_error varmapack_sim(double A[], double B[], double Sig[], double mu[],
                              int p, int q, int r, int n, int M, double X0[],
                              int nX0, randompack_rng *rng, double X[],
                              double E[])
{  
  int info;
  varmapack_error error = VARMAPACK_OK;
  double *C = 0, *G = 0, *S = 0, *SS = 0, *R = 0;
  double *Wrk = 0, *Psi = 0, *PsiHat = 0;
  double *CC = 0, *x0bar = 0, *e = 0, *wrk = 0;
  double *Aflp = 0, *Bflp = 0;
  bool Ealloc = E==0;
  if ((p > 0 && A == 0) || (q > 0 && B == 0) || Sig == 0 || X == 0 || rng == 0) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if (p < 0 || q < 0 || r <= 0 || n <= 0 || M <= 0 || n < imax(p, q)) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  if ((X0 == 0 && nX0 != 0) || (X0 != 0 && (nX0 < imax(p, q) || nX0 > n))) {
    return VARMAPACK_INVALID_ARGUMENT;
  }
  int h = imax(imax(p,q), nX0);
  int rn = r*n;  // Total number of observations
  int rh = r*h;  // Observation count in starting segment, order of SS, CC and EE
  bool stationary = p == 0 || varmapack_specrad(A, r, p) < 1;
  if (X0 == 0 && !stationary) return VARMAPACK_NONSTATIONARY;
  if (Ealloc && !ALLOC(E, rn*M)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (!randompack_mvn("T", 0, Sig, r, n*M, E, r, 0, rng)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (p == 0 && q == 0) {
    lacpy("All", rn, M, E, rn, X, rn);
    if (mu != 0) {
      for (int j=0; j<M; j++) {
        for (int t=0; t<n; t++) {
          axpy(r, 1.0, mu, 1, X + j*rn + t*r, 1);
        }
      }
    }
    if (Ealloc) FREE(E);
    return VARMAPACK_OK;
  }
  if (X0 != 0 && (!stationary || nX0 > imax(p, q))) {
    if (p > 0 && !ALLOC(Aflp, r*r*p)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (q > 0 && !ALLOC(Bflp, r*r*q)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    flipmat(A, Aflp, r, p);
    flipmat(B, Bflp, r, q);
    for (int j=0; j<M; j++) {
      copy(rh, X0, 1, X + j*rn, 1);
      if (mu != 0) {
        for (int t=0; t<h; t++) axpy(r, -1.0, mu, 1, X + j*rn + t*r, 1);
      }
    }
    lacpy("All", (n-h)*r, M, E + rh, rn, X + rh, rn);
    for (int t=h; t<n; t++) {
      int iX = r*t;
      int iA = r*(t - p);
      int iB = r*(t - q);
      if (p > 0)
        gemm("NoT", "NoT", r, M, r*p, 1.0, Aflp, r, X + iA, rn, 1.0,
             X + iX, rn);
      if (q > 0)
        gemm("NoT", "NoT", r, M, r*q, 1.0, Bflp, r, E + iB, rn, 1.0,
             X + iX, rn);
    }
    if (mu != 0) {
      for (int j=0; j<M; j++) {
        for (int t=0; t<n; t++) {
          axpy(r, 1.0, mu, 1, X + j*rn + t*r, 1);
        }
      }
    }
    FREE(Bflp);
    FREE(Aflp);
    if (Ealloc) FREE(E);
    return VARMAPACK_OK;
  }
  // SOLVE VECTOR-YULE-WALKER EQUATIONS FOR COVARIANCE OF X
  if (!ALLOC(C, r*r*(q+1))) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (!ALLOC(G, r*r*(q+1))) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (!ALLOC(S, r*r*(p+1))) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (!VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) {
    error = VARMAPACK_INTERNAL;
    goto fail;
  }
  printM("S", S, r, r*(p+1));
  printM("G", G, r, r*(q+1));
  if (rh > 0 && !ALLOC(SS, rh*rh)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (!SBuild("Low", S, A, G, p, q, r, h, SS)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  printM("SS", SS, rh, rh);
  FREE(S); FREE(G);
  if (rh > 0 && !ALLOC(R, rh*rh)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  printM("E", E, rn, M);
  if (X0 == 0) {  // Start series from scratch
    if (rh > 0 && !ALLOC(Wrk, rh*M)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (rh > 0 && !ALLOC(Psi, rh*rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (rh > 0 && !ALLOC(PsiHat, rh*rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    FindPsi(A, B, Psi, p, q, r);
    FindPsiHat(Psi, PsiHat, Sig, r, h);
    //printM("PsiHat", PsiHat, rh, rh);
    lacpy("Low", rh, rh, SS, rh, R, rh);
    syrk("Low", "NoT", rh, rh, -1.0, PsiHat, rh, 1.0, R, rh);
    // printM("R", R, rh, rh);
    if (!randompack_mvn("T", 0, R, rh, M, Wrk, rh, 0, rng)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    printM("Wrk", Wrk, r, h*M);
    lacpy("All", rh, M, Wrk, rh, X, rn);    // copy to X1
    gemm("NoT", "NoT", rh, M, rh, 1.0, Psi, // X1 := Psi*E(1:h) + X1
         rh, E, rn, 1.0, X, rn);
    FREE(PsiHat); FREE(Psi); FREE(Wrk);
  }
  else { // initialize series with X0
    double *Chat, *LS;
    if (!ALLOC(CC, rh*rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (!ALLOC(x0bar, rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (!ALLOC(e, rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    if (!ALLOC(wrk, rh)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    LS = SS;
    printI("rh", rh);
    potrf("Low", rh, LS, rh, &info); // Cholesky factorize SS
    if (info != 0) {
      error = info > 0 ? VARMAPACK_NOT_POSITIVE_SEMIDEFINITE : VARMAPACK_INTERNAL;
      goto fail;
    }
    CCBuild(A, C, p, q, r, h, CC);
    Chat = CC; // Chat = LS\CC
    trsm("Left", "Low", "NT", "NotUD", rh, rh, 1.0, LS, rh, Chat, rh);
    copy(rh, X0, 1, x0bar, 1);
    if (mu != 0) for (int i=0; i<rh; i+=r) axpy(r, -1.0, mu, 1, x0bar + i, 1);
    copy(rh, x0bar, 1, wrk, 1);
    printM("x0bar", x0bar, 1, rh);
    
    trsv("Lo", "NT", "NotUD", rh, LS, rh, wrk, 1);
    gemv("T", rh, rh, 1.0, Chat, rh, wrk, 1, 0.0, e, 1);
    for (int j=0; j<h; j++) { // Put Sig on diagonal blocks of R
      lacpy("Low", r, r, Sig, r, R + j*r*(rh + 1), rh);
    }
    syrk("L", "T", rh, rh, -1.0, Chat, rh, 1.0, R, rh);
    printMT("e", e, rh, 1);
    printM("R", R, rh, rh);
    printM("E0-fyrir", E, rh, M);
    if (!randompack_mvn("T", e, R, rh, M, E, rn, 0, rng)) {
      error = VARMAPACK_ALLOCATION;
      goto fail;
    }
    printM("E0-eftir", E, rh, M);
    for (int j=0; j<M; j++) {
      copy(rh, x0bar, 1, X + j*rn, 1);
    }
    printMP("X1", X, rh, M, X, rn);
    FREE(wrk); FREE(e); FREE(x0bar); FREE(CC); 
  }
  FREE(C);
  FREE(R); FREE(SS); 
  printMT("E", E, rn, M);
  lacpy("All", (n-h)*r, M, E + rh, rn, X + rh, rn);
  if (p > 0 && !ALLOC(Aflp, r*r*p)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  if (q > 0 && !ALLOC(Bflp, r*r*q)) {
    error = VARMAPACK_ALLOCATION;
    goto fail;
  }
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  for (int t=h; t<n; t++) {    
    int
      iX = r*t,
      iA = r*(t - p),
      iB = r*(t - q);
    // print4I("t, iX, iA, iB", t, iX, iA, iB);
    if (p > 0)
      gemm("NoT", "NoT", r, M, r*p, 1.0, Aflp, r, X + iA, rn, 1.0, X + iX, rn);
    if (q > 0)
      gemm("NoT", "NoT", r, M, r*q, 1.0, Bflp, r, E + iB, rn, 1.0, X + iX, rn);
    printMT("X", X, rn, M);
  }
  if (mu != 0) {
    for (int j=0; j<M; j++) {
      for (int t=0; t<n; t++) {
        axpy(r, 1.0, mu, 1, X + j*rn + t*r, 1);
      }
    }
  }
  FREE(Bflp); FREE(Aflp);
  if (Ealloc) FREE(E);
  return VARMAPACK_OK;
fail:
  FREE(Bflp); FREE(Aflp);
  FREE(wrk); FREE(e); FREE(x0bar); FREE(CC);
  FREE(PsiHat); FREE(Psi); FREE(Wrk);
  FREE(R); FREE(SS);
  FREE(S); FREE(G); FREE(C);
  if (Ealloc) FREE(E);
  return error;
}

static void CCBuild( // Build covariance between terms and shocks of VARMA time series
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double C[],  // in   r × r·(q+1), = [C0...Cq], Ci = cov(x(t),eps(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n must be <= max(p+q)+1)
  double CC[]) // out  r·n × r·n, cov(x1'...xn',eps1'...epsn')
{
  //  DESCRIPTION: The time series is as shown in SBuild (e.g.). The CC-matrix
  //  is:
  //                   C0  0 ...... 0 
  //                   C1 C0  0 ... 0
  //                   :  C1  .     :
  //                   :         .  0
  //                   C(n-1).. C1 C0
  //
  //  For q < j <= p, Cj is found with the recurrence relation:
  //
  //      Cj = A1*C(j-1) + A2*C(j-2) + ........ + Aj*C0   (*)
  //
  int i,j;
  double *CCj;
  xAssert(n <= imax(p,q) + 1);
  // Fill CC with zeros:
  setzero(r*n*r*n, CC);
  // Copy C0..Cq to first q+1 block-rows in first block-column:
  for (j=0; j<=q && j<n; j++) lacpy("All", r, r, C+j*r*r, r, CC + j*r, r*n);
  // Calculate rest of first block-column using (*):
  for (j=q+1; j<n; j++) {
    CCj = CC + j*r;
    for (i=1; i<=j && i<=p; i++)  // (using that CC was set to 0)
      gemm("N", "N", r, r, r, 1.0, A+(i-1)*r*r, r, CC+(j-i)*r, r*n,1.0,CCj,r*n);
  }
  // Copy head of first block-column to other columns
  for (j=1; j<n; j++) {
    CCj = CC + j*r*n*r + j*r;
    lacpy("All", r*(n-j), r, CC, r*n, CCj, r*n);
  }
}
