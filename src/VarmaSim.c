// VarmaSim — simulate spin-up freem AR, VAR, ARMA or VARMA time series
//
// See VarmaSim.h for parameter descriptions and storage conventions
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
//   The program is a translation of the Matlab program varma_sim described in
//   references [1], [2] and [3] (the report [3] contains some details not found
//   in the published papers).
//
// RANDOM NUMBER GENERATION
//   The random number generator used for the shocks is specified by the rng
//   parameter. Several different generators are possible: The built-in R
//   generator (when called from R), an xorshift128+ generator implemented in
//   the package, and a Park-Miller generator allowing generation of identical
//   shock series across environments (R, Matlab, and C). See RandomNumbers.h
//   for details.
// 
// STARTING VALUES:
//   When starting values for the simulation are specified in X0, each generated
//   series will start with X0, and subsequent terms are simulated using X0 as
//   input. An important feature is that the simulated series have an accurate
//   distribution from term one (they are spin-up freem). Thus there is no need
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
//   If the series is nonstationary or Sig is not positive semidefinite then the
//   function returns (quietly) with ok = 0 (false). In exceptional cases,
//   (which should only be possible when X0 is specified and the series is
//   nonstationary) needed matrices may be singular, and then the function exits
//   by calling one of the function defined in xAssert.h, which are expected to
//   print an error message and terminate the program.
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
#include "allocate.h"
#include "BlasGateway.h"
#include "VYW.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "RandomNumbers.h"
#include "VarmaSim.h"
#include "xAssert.h"

void VarmaSim(double A[], double B[], double Sig[], double mu[], int p, int q,
              int r, int n, int M, double X0[], int nX0, RandRng *rng, double X[],
              double E[], bool *ok)
{
  int h = max(p,q), k = (X0 == 0 ? h : nX0), rk = r*k;
  int rn = r*n;  // Total number of observations
  int rh = r*h;  // Observation count in starting segment, order of SS, CC and EE
  bool Ealloc = E==0;
  xAssert(p>=0 && q>=0 && r>0 && M>0);
  xAssertMessage(n > h, "Illegal parameter in VarmaSim, n must be > max(p,q)");
  xAssertMessage(nX0 == 0 || (h <= nX0 && nX0 <= n),
                 "Illegal parameter in VarmaSim, max(p,q) ≤ nX0 ≤ n");
  *ok = true;
  if (Ealloc) allocate(E, rn*M);
                 
  // SOLVE VECTOR-YULE-WALKER EQUATIONS FOR COVARIANCE OF X
  double *vywFactors;
  int *piv, info;
  int nVYW = (p == 0) ? 0 : r*r*p - r*(r-1)/2; // order of matrix of VYW-equations
  allocate(vywFactors, nVYW * nVYW);
  allocate(piv, nVYW);
  VYWFactorize(A, vywFactors, piv, p, r, &info);
  if (info != 0) {
    freem(piv); freem(vywFactors);
    xErrorExit("VarmaSim: Singular Yule-Walker equations, unable to continue");
  }
  double *C, *G, *S;
  allocate(C, r*r*(q+1));
  allocate(G, r*r*(q+1));
  allocate(S, r*r*(p+1));
  FindCG(A, B, Sig, p, q, r, C, G);
  VYWSolve(A, vywFactors, S, G, 1, q+1, piv, p, r);
  double *SS;
  allocate(SS, rk*rk);
  SBuild("Low", S, A, G, p, q, r, k, SS);
  freem(S); freem(G); freem(C);
  freem(piv); freem(vywFactors);
  double *R;
  allocate(R, rk*rk);
  // GENERATE INITIAL SEGMENT OF X
  if (X0 == 0) {  // Start series from scratch
    double *Wrk, *Psi, *PsiHat;
    allocate(Wrk, rh*M);
    allocate(Psi, rh*rh);
    allocate(PsiHat, rh*rh);
    FindPsi(A, B, Psi, p, q, r);
    FindPsiHat(Psi, PsiHat, Sig, r, h);
    lacpy("Low", rh, rh, SS, rh, R, rh);
    syrk("Low", "NoT", rh, rh, -1.0, PsiHat, rh, 1.0, R, rh);
    RandNM("T", 0, Sig, r, h*M, E, 0, rng); // first h shocks
    RandNM("T", 0, R, rh, M, Wrk, 0, rng);  // draw Wrk from N(0, R)
    lacpy("All", rh, M, Wrk, rh, X, rh);    // copy to X1
    gemm("NoT", "NoT", rh, M, rh, 1.0, Psi, // X1 := Psi*E(1:h) + X1
         rh, E, rn, 1.0, X, rn);
    freem(PsiHat); freem(Psi); freem(Wrk);
  }
  else { // initialize series with X0
    double *CC, *Chat, *LS, *x0bar, *e;
    allocate(CC, rk*rk);
    allocate(x0bar, rk);
    allocate(e, rk*M);
    SBuild("Low", S, A, G, p, q, r, k, SS);
    LS = SS;
    potrf("Low", rk, LS, rk, &info); // Cholesky factorize SS
    xAssert(info == 0);
    CCBuild(A, C, p, q, r, k, CC);
    Chat = CC; // Chat = LS\CC
    trsm("Left", "Low", "NT", "NotUD", rk, rk, 1.0, LS, rk, Chat, rk);
    copy(rk, X0, 1, x0bar, 1);
    for (int i=0; i<rk; i+=r) axpy(r, -1.0, mu, 1, x0bar, 1);
    trsv("Lo", "NT", "NotUD", rk, LS, rk, x0bar, 1);
    gemv("T", rk, rk, 1.0, Chat, rk, x0bar, 1, 0.0, e, 1);
    for (int j=0; j<k; j++) { // Put Sig on diagonal blocks of R
      lacpy("Low", r, r, Sig, r, R + j*r*(rk + 1), rk);
    }
    syrk("L", "T", rk, rk, 1.0, Chat, rk, 1.0, R, rk);
    RandNM("T", e, R, rk, M, E, 0, rng); // first k shocks
    freem(e); freem(x0bar); freem(CC); 
  }
  freem(R); freem(SS); 
  RandNM("T", 0, Sig, r, (n-k)*M, E + rk, 0, rng); // remaining shocks
  copy((n-k)*r*M, E + rk, 1, X + rk, 1);
  double *Aflp, *Bflp;
  allocate(Aflp, r*p);
  allocate(Bflp, r*q);
  flipmat(A, Aflp, r, p);
  flipmat(B, Bflp, r, q);
  for (int t=k; t<n; t++) {    
    int
      iX = r*t,
      iA = r*(t - p),
      iB = r*(t - q);
    gemm("NoT", "NoT", r, M, r*p, 1.0, Aflp, r, X + iA, rn, 1.0, X + iX, rn);
    gemm("NoT", "NoT", r, M, r*q, 1.0, Bflp, r, X + iB, rn, 1.0, X + iX, rn);
  }
  freem(Bflp); freem(Aflp);
  if (Ealloc) freem(E);
}
