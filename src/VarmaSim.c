// VarmaSim — simulate spin-up free AR, VAR, ARMA or VARMA time series
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
//   distribution from term one (they are spin-up free). Thus there is no need
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

#include "allocate.h"
#include "BlasGateway.h"
#include "VYW.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"
#include "RandomNumbers.h"
#include "VarmaSim.h"
#include "xAssert.h"

void VarmaSim(double A[], double B[], double Sig[], double mu[], int p, int q,
              int r, int n, int M, double X0[], RandRng *rng, double X[],
              double eps[], int *ok)
{
  int *piv, info, i, k, t, CHOLF, I, Ip, Iq;
  double *vywFactors, *C, *G, *W, *S, *SS, *CC, *EE, *Xk, *mEps, *invSStimesCC;
  double del, *LSS, *Eps, *Aflp=0, *Bflp=0, *L, *XX;
  int h = max(p,q);
  int N = r*n;  // Total number of observations
  int g = r*h;  // Observation count in starting segment, order of SS, CC and EE
  int nVYW = r*r*p - r*(r-1)/2;            // order of matrix of VYW-equations
  xAssert(p>=0 && q>=0 && r>0 && M>0);
  xAssertMessage(n > h, "Illegal parameter in VarmaSim, n must be > max(p,q)");
  
  // SOLVE VECTOR-YULE-WALKER EQUATIONS FOR COVARIANCE OF X
  if (p == 0) nVYW = 0;
  allocate(vywFactors, nVYW * nVYW);
  allocate(piv, nVYW);
  VYWFactorize(A, vywFactors, piv, p, r, &info);
  if (info != 0) {
    freem(piv); freem(vywFactors);
    xErrorExit("VarmaSim: Singular Yule-Walker equations, unable to continue");
  }
  allocate(C, r*r*(q+1));
  allocate(G, r*r*(q+1));
  allocate(W, r*r*(q+1));
  allocate(S, r*r*(p+1));
  FindCGW(A, B, Sig, p, q, r, C, G, W);
  VYWSolve(A, vywFactors, S, G, 1, q+1, piv, p, r);
  freem(piv); freem(vywFactors);
  
  // COVARIANCES OF START SEGMENTS: SS = COV(X), CC = COV(X,Eps), EE = COV(Eps)
  allocate(SS, g*g);
  allocate(CC, g*g);
  allocate(EE, g*g);
  SBuild("Low", S, A, G, p, q, r, h, SS);
  CCBuild(A, C, p, q, r, h, CC);
  // EE = block diagonal matrix with all blocks = Sig
  for (i=0; i<h; i++) lacpy("All", r, r, Sig, r, EE + i*r*g + i*r, g);
  freem(S); freem(W); freem(G); freem(C);
  printSetNdec(15);
  
  // GENERATE INITIAL SEGMENT OF X
  allocate(LSS, g*g); // For (lower) Cholesky factor of SS
  allocate(XX, N*M);
  if (X0) { // initialize series with X0
    for (k=0; k<M; k++) {
      Xk = X + k*N;
      lacpy("All", r, h, X0, r, Xk, r);
      if (mu) for (t=0; t<h; t++) axpy(r, -1.0, mu, 1, Xk + t*r, 1);
    }
    copy(g*g, SS, 1, LSS, 1);
    potrf("Low", g, LSS, g, &info); // Try to Cholesky factorize SS
    CHOLF = (info == 0);
  }
  else {  // Start series from scratch
    RandNM(0, SS, M, g, XX, LSS, &del, rng, ok); // draw initial part from
    printM("SS", SS, g, g);
    printM("LSS", LSS, g, g);
    copytranspose(M, g, XX, M, X, N);            // correct distribution
    if (!*ok) {                                
      freem(XX); freem(LSS); freem(EE); freem(CC); freem(SS);
      return;
    }
    CHOLF = 1;
  }
  freem(SS);

  // COPY INITIAL PART OF X TO mEps (LATER TO STORE EXPECTED Eps VALUES GIVEN X)
  allocate(mEps, g*M);
  for (k=0; k<M; k++) copy(g, X + k*N, 1, mEps + k*g, 1);

  // DETERMINE COVARIANCE MATRIX OF INITIAL SHOCKS GIVEN Xinitial
  if (g > 0 && CHOLF) { // we have the Cholesky factor of S in LSS
    trsm("Left", "Low", "NT", "NotUD", g, g, 1.0, LSS, g, CC, g); // CC:=LSS\CC
    syrk("Low", "T", g, g, -1.0, CC, g, 1.0, EE, g); // EE := EE-CC'*inv(SS)*CC
    // mEps := LSS\Xinitial:
    trsm("Left", "Low", "NT", "NotUD", g, M, 1.0, LSS, g, mEps, g);
    printM("mEps", mEps, g, M);
    printM("EE", EE, g, g);
  }
  else if (g > 0) {  // Need to use Gaussian elimination
    copy(g*g, SS, 1, LSS, 1); // Use LSS to store factors
    allocate(piv, g);
    getrf(g, g, LSS, g, piv, &info);
    if (info != 0) {
      freem(piv); freem(mEps); freem(LSS); freem(EE); freem(CC);
      xErrorExit("VarmaSim: Error, unable to finish, singular SS matrix");
    }
    allocate(invSStimesCC, g*g);
    copy(g*g, CC, 1, invSStimesCC, 1);
    getrs("NoT", g, g, LSS, g, piv, invSStimesCC, g, &info);
    xAssert(info == 0); // only nonzero if an argument is wrong
    // EE <-- EE-CC'*inv(SS)*CC:
    gemm("T", "NT", g, g, g, -1.0, CC, g, invSStimesCC, g, 1.0, EE, g);
    getrs("NoT", g, M, LSS, g, piv, mEps, g, &info);  // mEps:= SS\Xinitial
    xAssert(info == 0); // only nonzero if an argument is wrong
    freem(invSStimesCC); freem(piv);
  }
  freem(LSS);
  
  // GENERATE INITIAL SHOCKS:
  if (!eps) { allocate(Eps, N*M); }
  else Eps = eps;
  RandNM(0, EE, M, g, XX, 0, &del, rng, ok);
  printM("EE", EE, g, g);
  copytranspose(M, g, XX, M, Eps, N);
  freem(EE);
  if (!*ok) {
    if (!eps) freem(Eps);
    freem(mEps); freem(CC);
    xErrorExit("VarmaSim: Error, singular covariance matrix of initial shocks");
  }
  // Either CC contains LSS\CC and mEps has LSS\Xinitial (if CHOLF) or CC is
  // still CC and mEps has SS\Xinitial. In both cases the following dgemm call
  // adds CC'*(SS\Xinitial) to Eps, to give the initial shocks.
  if (g) gemm("T", "NT", g, M, g, 1.0, CC, g, mEps, g, 1.0, Eps, N);
  freem(mEps); freem(CC);
  printM("Eps", Eps, M, g);

  // NOW GENERATE THE REST OF THE SHOCKS
  allocate(L,r*r);
  RandNM(0, Sig, (n-h)*M, r, XX, 0, &del, rng, ok);
  for (i=0; i<M; i++)
    copytranspose(n-h, r, XX + i*(n-h), M*(n-h), Eps + g + i*N, r);
  freem(XX);
  freem(L);

  // NOW GENERATE THE REST OF THE SERIES:
  if (p > 0) {
    allocate(Aflp, r*r*p);  // [Ap...A2 A1]
    for (i=0; i<p; i++) lacpy("All", r, r, A + i*r*r, r, Aflp + (p-1-i)*r*r, r);
  }
  if (q > 0) {
    allocate(Bflp, r*r*q);  // [Bq...B2 B1]
    for (i=0; i<q; i++) lacpy("All", r, r, B + i*r*r, r, Bflp + (q-1-i)*r*r, r);
  }
  for (k=0; k<M; k++) copy(N - g, Eps + k*N + g, 1, X + k*N + g, 1);
  for (I=g; I<N; I+=r) {
    Ip = I - p*r;
    Iq = I - q*r;
    gemm("N", "N", r, M, r*p, 1.0, Aflp, r, X + Ip, N, 1.0, X + I, N);
    gemm("N", "N", r, M, r*q, 1.0, Bflp, r, Eps + Iq, N, 1.0, X + I, N);
  }
  printM("X", X, r, n);
  if (mu)
    for (k=0; k<M; k++) // Add mu to each x(t)
      for (t=0; t<n; t++)
        axpy(r, 1.0, mu, 1, X + r*(t + k*n), 1);
  if (p > 0)
    freem(Aflp);
  if (q > 0)
    freem(Bflp);
  if (!eps) freem(Eps);
}
