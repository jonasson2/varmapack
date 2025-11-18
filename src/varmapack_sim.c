// varmapack_sim — simulate spin-up freem AR, VAR, ARMA or VARMA time series
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
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "RandomNumbers.h"
#include "varmapack.h"
#include "xAssert.h"
#include "varmapack_VYW.h"

static void SBuild( char *uplo, double S[], double A[], double G[], int p, int q, int r, int n,
             double SS[]);
static void CCBuild( double A[], double C[], int p, int q, int r, int n, double CC[]);

void varmapack_sim(double A[], double B[], double Sig[], double mu[], int p, int q,
              int r, int n, int M, double X0[], int nX0, RandRng *rng, double X[],
              double E[], bool *ok)
{
  int h = imax(p,q), k = (X0 == 0 ? h : nX0), rk = r*k;
  int rn = r*n;  // Total number of observations
  int rh = r*h;  // Observation count in starting segment, order of SS, CC and EE
  int info;
  bool Ealloc = E==0;
  xAssert(p>=0 && q>=0 && r>0 && M>0);
  xAssertMessage(n > h, "Illegal parameter in varmapack_sim, n must be > max(p,q)");
  xAssertMessage(nX0 == 0 || (h <= nX0 && nX0 <= n),
                 "Illegal parameter in varmapack_sim, max(p,q) ≤ nX0 ≤ n");
  *ok = true;
  if (Ealloc) allocate(E, rn*M);
                 
  // SOLVE VECTOR-YULE-WALKER EQUATIONS FOR COVARIANCE OF X
  double *C, *G, *S;
  allocate(C, r*r*(q+1));
  allocate(G, r*r*(q+1));
  allocate(S, r*r*(p+1));
  if (!varmapack_VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G)) {
    freem(S); freem(G); freem(C);
    xErrorExit("varmapack_sim: Singular Yule-Walker equations, unable to continue");
  }
  double *SS;
  allocate(SS, rk*rk);
  SBuild("Low", S, A, G, p, q, r, k, SS);
  freem(S);
  freem(G);
  double *R;
  allocate(R, rk*rk);
  // GENERATE INITIAL SEGMENT OF X
  if (X0 == 0) {  // Start series from scratch
    double *Wrk, *Psi, *PsiHat;
    allocate(Wrk, rh*M);
    allocate(Psi, rh*rh);
    allocate(PsiHat, rh*rh);
    varmapack_FindPsi(A, B, Psi, p, q, r);
    varmapack_FindPsiHat(Psi, PsiHat, Sig, r, h);
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
  freem(C);
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

static void SExtend ( // Extend Sj matrices to include S(p+1)...S(n-1)
  double A[],    // in   r × r × p, autoregressive parameter matrices
  double G[],    // in   r × r × (q+1), G0, G1,...,Gq
  double S[],    // in   r × r × (p+1), S0, S1,...,Sp
  double Scol[], // out  r·n × r, column of Si-matrices S0, S1,..., S(n-1)
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of x(t)
  int n)         // in   length of time series
{
  int iScol = n*r, i, j;
  double *Scolj, *Ai, *Scoli;
  for (j=0; j<p+1; j++) {
    lacpy("All", r, r, S + j*r*r, r, Scol + j*r, iScol);
  }
  for (j=p+1; j<n; j++) {
    Scolj = Scol + j*r;
    if (j<=q) {
      lacpy("All", r, r, G + j*r*r, r, Scolj, iScol);
    }
    for (i=0; i<p; i++) {
      Ai = A+i*r*r;
      Scoli = Scolj - (i+1)*r;
      gemm("N", "N", r, r, r, 1.0, Ai, r, Scoli, iScol, 1.0, Scolj, iScol);
    }
  }
}

static void SBuild( // Build covariance matrix of all the values of a VARMA time series
  char *uplo,  // in   Create all of SS (when "A") or lower part only (when "L")
  double S[],  // in   r × r·(p+1), = [S0...Sp], Si = Cov(x(t),x(t-i))
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double G[],  // in   r × r·(q+1), = [G0...Gq], Gi = cov(y(t),x(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n can be any nonnegative value, even 0)
  double SS[]) // out  r·n × r·n, covariance of [x1'...xn']'
{
  //  DESCRIPTION: The time series is given by
  //
  //                  x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)
  //  where
  //                  y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
  //
  //  and x(t), y(t) and eps(t) are r-dimensional with eps(t) N(0,Sig). S and G
  //  can (for example) have been obtained with varmapack_FindCG. The SS matrix is:
  //
  //                   S0  S1' S2'...Sn-1'
  //                   S1  S0  S1'...Sn-2'
  //                   S2               :
  //                   :                :
  //                   Sn-1 ...... S1  S0
  //
  //  and the Sj are found with the recurrence relation:
  //
  //      Sj = A1*S(j-1) + A2*S(j-2) + ........ + Ap*S(j-p) + Gj
  //
  //  with Gj = 0 for j > q.
  //
  // NOTE: This function may (of course) be used to determine the covariance 
  // matrix of a segment of a timeseries by specifying n smaller then the total
  // length of the series (for example an initial segment, but since the series
  // is stationary all segments of the same length have the same covariance).
  double *Scol, *SSj, *SSi;
  int j, m;
  if (n==0) return;
  m = imax(p+1,n);
  allocate(Scol, (r*m)*r);
  SExtend(A, G, S, Scol, p, q, r, m);
  for (j=0; j<n; j++) {
    SSj = SS + j*r*n*r + j*r;
    if (uplo[0] == 'A') {
      SSi = SSj + r*n*r;
      lacpy("All", r*(n-j), r, Scol, r*m, SSj, r*n);
      if (j < n-1) copytranspose(r*(n-j-1), r, Scol+r, r*m, SSi, r*n);
    }
    else {
      lacpy("Low", r*(n-j), r, Scol, r*m, SSj, r*n);
    }
  }
  freem(Scol);
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
