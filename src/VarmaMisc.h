// Declaration of miscellaneous functions called by VarmaLoglik & VarLoglik
#ifndef VARMAMISC_H
#define VARMAMISC_H

#include "mds.h"
void atbac( // Set C := A'·B·A + C and update derivative for symmetric B and C
  double A[],    // in      m×n matrix
  double B[],    // in      m×m symmetric matrix, only lower triangle used
  double C[],    // in/out  n×n symmetric matrix, only lower triangle used
  int m,         // in      row count of A, order of B
  int n,         // in      column count of A, order of C
  double Ad[],   // in      m×n×nPar, derivative of A (or null if nPar=0)
  double Bd[],   // in      m×m×nPar, derivative of B (or null)
  double Cd[],   // in/out  n×n×nPar, derivative of C (or null)
  int nPar);     // in      number of paramters to differentiate w.r.t.
  // Primary tester: test_atba_c.m

void ChngVar ( // Postmultiply gradient with variable-change Jacobian
  double **Fd,  // in/out  mF×nPar×Q in, mF×(nPar+nJ-mJ)×Q out, derivatives of F
  double J[],   // in      mJ×nJ, Jacobian of variable change
  int mF,       // in      row count of i-th plane of Fd
  int Q,        // in      number of planes in Fd
  int mJ,       // in      row count of J
  int nJ,       // in      column count of J
  int nPar);    // in      number of original parameters

void CMultiply(// Multiply with transpose of expanded C matrix
  double C[],  // in   r×r×(q+1), C matrices from FindCGW
  double X[],  // in   (r·n-M)-vector to be multiplied.
  double Y[],  // out  r·n-vector, returns product
  int q,       // in   number of moving average terms in time series
  int r,       // in   dimension of each x(t) and hence each Ci matrix
  int n,       // in   length of time series
  int miss[]); // in   r·n-vector, miss(j)=1 if x(j) is missing, otherwise 0
  // Primary tester: test_C_multiply

void FindAcol ( // Make column of Aj matrices (and associated derivatives)
  double A[],     // in   r x r×p, autoregressive parameter matrices
  double Acol[],  // out  r·p×r, returns [A1; A2;... Ap]
  int p,          // in   number of autoregressive terms
  int r,          // in   dimension of each xt
  double Acold[], // out  r×r×p×nPar, derivatives of Acol; null if nPar=0
  int nPar);      // in   number of parameters, nPar=0 for no derivatives
  // Primary tester: test_find_lambda_om

void FindCGW ( // Calculate the Gi and Wi matrices for VARMALLC and VARMALLM
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[],   // out  r×r×(q+1) with C0, C1...Cq where Ci = cov(x(t),eps(t-i))
  double G[],   // out  r×r×(q+1) with G0, G1...Gq where Gi = cov(y(t),x(t-i))
  double W[]);  // out  r×r×(q+1) with W0, W1...Wq where Wi = cov(y(t),y(t-i))
  // Primary tester: test_find_CGW.m

void FindCGWDeriv ( // Calculate Gi and Wi matrices along with derivatives
  double A[],     // in   r×r×p autoregressive parameter matrices
  double B[],     // in   r×r×q moving average parameter matrices
  double Sig[],   // in   r×r covariance of shocks
  int p,          // in   number of AR terms
  int q,          // in   number of MA terms
  int r,          // in   dimension of each x(t)
  struct mds **C, // out  Ci and their derivatives in MDS-format (see mds.h)
  struct mds **G, // out  Gi and their derivatives in MDS-format (see mds.h)
  struct mds **W);// out  Wi and their derivatives in MDS-format (see mds.h)
  // Primary tester: test_find_CGW_deriv.m

void FindGneg ( // Make Gj matrices for negative j
  double A[],     // in   r x r×p, autoregressive parameter matrices
  double B[],     // in   r×r×q, moving average parameter matrices
  double C[],     // in   r×r×(q+1), C matrices from FindCGW
  double Gneg[],  // out  r·(n-p-1)×r, column of Gj for negative j
  int p,          // in   number of autoregressive terms
  int q,          // in   number of moving average terms
  int r,          // in   dimension of each xt
  int n,          // in   length of time series
  double Cd[],    // in   r×r×(q+1)×nPar, derivat. of C, or null if nPar=0
  double Gnegd[], // out  r·(n-p-1)×r×nPar, derivat. of Gneg; null if nPar=0
  int nPar);      // in   number of parameters, nPar=0 for no derivatives
  // Primary tester: test_find_Vhat.m

void FindLambdaOm ( // Determine matrix Lam_om and optionally its derivative
  double Acol[],  // in   r·p×r, column of autoregressive parameter matrices
  double Laom[],  // out  mLaom×M, returns Lambda(obs,miss)
  int mLaom,      // in   first dimension of Laom and Laomd
  int miss[],     // in   r×n, miss(i,t)=1 if x(i,t) is missing, 0 if observed
  int ko[],       // in   ko[t]=number of observed values before time t+1, t<=n
  int km[],       // in   km[t]=number of missing values before time t+1, t<=n
  int p,          // in   number of autoregressive terms
  int r,          // in   dimension of x(t)
  int n,          // in   length of time series
  double Acold[], // in   r·p x r×nPar, derivatives of Acol, or null if nPar=0
  double Laomd[], // out  mLaom×M×nPar, derivat. of Laom, or null if nPar=0
  int nPar);      // in   no. of params, set to zero to suppress derivative calc
  // Primary tester: test_find_lambda_om
  // Also tested by: test_profile.m, test_varma_llm.m

int findmaxnz ( // Return maximum non-zero row-number of V or Lam_om
  int ko[], // in   ko[t]=number of observed values before time t+1,t<=n
  int km[], // in   km[t]=number of missing values before time t+1,t<=n
  int p,    // in   number of autoregressive terms
  int g,    // in   g=q for V and g=p for Lam_om
  int n);   // in   length of series
  // Primary tester: test_misc.c

void FindMissingInfo( // Construct arrays with info on which values are missing
  int miss[], // in   r×n, miss(i,j) = 1 if x(i,t) is missing, otherewise 0
  int ko[],   // out  n-vec, ko[t] = number of observed values before time t+1
  int km[],   // out  n-vec, km[t] = number of missing values before time t+1
  int r,      // in   dimension of x(t)
  int n);     // in   length of time series
  // Primary tester: test_misc.c

void FindRes ( // Find maximum likelihood estimate of residuals for VarmaLoglikC
  double A[],    // in   r×r×p, autoregressive parameter matrices
  double C[],    // in   r×r×(q+1), C matrices from FindCGW
  double Lu[],   // in   r·h×r·h, upper left part of L-factor of Omega
  double Ll[],   // in   (n-h)·r×(q+1)·r, lower part of L-factor of Omega
  double what[], // in   n·r-vector
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of xt
  int n,         // in   length of series
  double eps[]); // out  r×n, the estimated residuals
  // **** No tester yet ****

void FindResMiss (
  double A[],    // in   r×r×p, autoregressive parameter matrices
  double C[],    // in   r×r×(q+1), C matrices from FindCGW
  double Lu[],   // in   mSu×mSu, upper left part of L-factor of Omega
  double Ll[],   // in   mOlow×(q+1)·r, lower part of L-factor of Omega
  double what[], // in   N-vector, either what (for M=0) or wohat (if M>0)
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of xt
  int n,         // in   length of series
  int ko[],      // in   ko[t]=number of observed values before time t+1,t<=n
  int miss[],    // in   r×n, miss(i,t)=1 if X(i,t) is missing, 0 if observed
  double LR[],   // in   M×M, LR-matrix
  double LQ[],   // in   M×M, LQ-matrix
  double K[],    // in   M×M, K-matrix
  double u[],    // in   M-vector
  double v[],    // in   M-vector
  double Laomh[],// in   mL×M, Lambda_om_hat
  double Vhat[], // in   mV×M
  int mL,        // in   first dimension of Laomh
  int mV,        // in   first dimension of Vhat
  double LSm[],  // in   M×M, Cholesky factor of Sm-matrix
  double mu[],   // in   r-vector, mean of x(t) or null for zero mean
  double eps[],  // out  r×n, the estimated residuals
  double xm[]);  // out  mM-vector, estimated missing values, null to not find

void FindResMissAr ( // Maximum likelihood estimate of residuals for VARLoglik
  double A[],    // in   r×r×p, autoregressive parameter matrices
  double Sig[],  // in   r×r, covariance of eps
  double Lu[],   // in   r·h×r·h, upper left part of L-factor of Omega
  double LSig[], // in   r×r×npat, Cholesky factors of lower blocks of Omega
  double z[],    // in   N-vector with wohat
  int p,         // in   number of autoregressive terms
  int r,         // in   dimension of xt
  int n,         // in   length of series
  int ko[],      // in   ko[t] = nmbr of observed values before time t+1, t<=n
  int miss[],    // in   r×n, miss(i,t)=1 if X(i,t) is missing, 0 if observed
  double LR[],   // in   M×M, LR-matrix
  double LQ[],   // in   M×M, LQ-matrix
  double K[],    // in   M×M, K-matrix
  double u[],    // in   M-vector
  double v[],    // in   M-vector
  double Laomh[],// in   mL×M, Lambda_om_hat
  double Vhat[], // in   mV×M
  int mL,        // in   first dimension of Laomh
  int mV,        // in   first dimension of Vhat
  double Sm[],   // in   M×M, Sm-matrix
  int jpat[],    // in   jpat[i] = index of missing pattern in ith lower block
  double mu[],   // in   r-vector, mean of x(t) or null for zero mean
  double eps[],  // out  r×n, the estimated residuals (or null for only miss v.)
  double xm[]);  // out  M, the estimated missing values (or null for only res.)

void FindResMissMVN ( // Estimate residuals and missing values for MVNLoglik
  double Sig[],  // in   r×r, covariance of eps
  double LSig[], // in   r×r×npat, Cholesky factors of lower blocks of Omega
  double w[],    // in   N-vector with wo
  int r,         // in   dimension of xt
  int n,         // in   length of series
  int ko[],      // in   ko[t] = number of observed values before t+1, t<=n
  int miss[],    // in   r×n, miss(i,t)=1 if X(i,t) is missing, 0 if observed
  int jpat[],    // in   jpat[i] = index of missing pattern in ith block
  double mu[],   // in   r-vector, mean of x(t) or null for zero mean
  double eps[],  // out  r×n, the estimated residuals (or null for only miss v.)
  double xm[]);  // out  M, the estimated missing values (or null for only res.)
  // The determined estimated values are maximum likelihood estimates.

void FindSm ( // Determine matrix Sm and optionally its derivative
  double Scol[], // in   r·n×r, column of Si-matrices S0, S1,..., S(n-1)
  double Sm[],   // out  M×M, missing rows and columns of S matrix
  int miss[],    // in   r×n, miss(i,t) = 1 if x(i,t) is missing, 0 if observed
  int r,         // in   dimension of x(t)
  int n,         // in   length of time series
  double Scold[],// in   r·n×r×nPar, derivative of Scol (or null)
  double Smd[],  // out  M×M×nPar, derivatives of Sm
  int nPar);     // in   n of params, set to zero to suppress derivative calc
  // Primary tester: test_find_Sm

void FindV ( // Find the matrix V and optionally its derivative
  double G[],    // in   r×r×(q+1), Gi-matrices
  double Gneg[], // in   r·(n-p-1)×r, (column of) Gi-matrices for negative i
  double Scol[], // in   r·n×r, column of Si-matrices S0, S1,..., S(n-1)
  double V[],    // out  mV×M, M is total number of missing values
  int mV,        // in   number of rows in V and Vd (mV <= N)
  int miss[],    // in   r×n, miss[i,t] = true if the (i,t)-value is missing
  int km[],      // in   km[t] = number of missing values before time t+1, t<=n
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of x(t)
  int n,         // in   length of time series
  double Gd[],   // in   r×r×nPar×(q+1), derivat. of G (or null if nPar=0)
  double Gnegd[],// in   r·(n-p-1)×r×nPar, derivative of Gneg (or null)
  double Scold[],// in   r·n×r×nPar, derivative of Scol (or null)
  double Vd[],   // out  mV×M×nPar, returns derivative of V
  int nPar);     // in   n of params, set to zero to suppress derivative calc
  // Primary testers: test_find_Vhat.m, test_profile.m

int IsStationary ( // Return 1 if model is stationary
  double A[],  // in  r×r×p, autoregressive parameter matrices
  double Sig[],// in  r×r, covariance matrix of shocks
  double LU[], // in  N×N with N = r*r*p-r*(r-1)/2, factors from vyw_factorize
  int piv[],   // in  N-vector, pivots from vyw_factorize
  int p,       // in  number of autoregressive terms
  int r);      // in  dimension of each x(t)
  // Tested by test_varma_llm.m, test_varma_llm_deriv.m +++ better testing needed

void LambdaMultiply(// Multiply with Lambda or Lambda(obs,obs)
  double A[],  // in   r×r×p, autoregressive parameter matrices
  double X[],  // in   N×nc matrix to be multiplied, returns product
  int nc,      // in   column count for X and Y
  int p,       // in   number of autoregressive terms, last dimension of A
  int r,       // in   dimension of x(t), and hence each Ai matrix
  int n,       // in   length of time series (Lambda is n·r×n·r)
  int miss[],  // in   r×n, miss(i,j) = 1 if x(i,t) is missing, otherewise 0
  int nPar,    // in   number of parameters, set to 0 to suppress derivatives
  double Xd[]);// out  mY×nc×nPar, returns derivatives of Y (see NOTE B)
  // Primary tester: test_lambda_multiply.m

void LambdaTMultiply(// Multiply with transpose of Lambda or Lambda(obs,obs)
  double A[],  // in      r×r×p, autoregressive parameter matrices
  double X[],  // in/out  mX-vector to be multiplied. Also returns product
  int p,       // in      number of autoregressive terms, last dimension of A
  int r,       // in      dimension of x(t), and hence each Ai matrix
  int n,       // in      length of time series (Lambda is n·r×n·r)
  int miss[]); // in      r×n, miss(i,j) = 1 if x(i,t) is missing, otherewise 0
  // Primary tester: test_lambda_multiply.m

double LogdetL(  // returns log(det(L))
  double L[],        // in   an m×m lower triangular matrix
  int iL,            // in   leading defining dimension of L and Ld
  int m,             // in   dimension of used part of L and Ld
  int nPar,          // in   number of parameters (set to 0 for no derivatives)
  double Ld[],       // in   derivatives of L w.r.t. nPar parameters, m×m×nPar
  int iLd2,          // in   second defining dimension of Ld
  double logdetd[]); // out  nPar-vector, derivatives of log det(L)
  // Tested by: test_omega.m, test_omega_ar.m, test_omega_deriv.m

void SExtend ( // Extend Sj matrices to include S(p+1)...S(n-1)
  double A[],    // in   r×r×p, autoregressive parameter matrices
  double G[],    // in   r×r×(q+1), G0, G1,...,Gq
  double S[],    // in   r×r×(p+1), S0, S1,...,Sp
  double Scol[], // out  r·n×r, column of Si-matrices S0, S1,..., S(n-1)
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of x(t)
  int n,         // in   length of time series
  double Gd[],   // in   r×r×nPar×(q+1), derivat. of G (or null if nPar=0)
  double Sd[],   // in   r×r×nPar×(p+1), derivat. of S (or null if nPar=0)
  double Scold[],// out  r·n×r×nPar, derivative of Scol (or null)
  int nPar);     // in   n of params, set to zero to suppress derivative calc
  // It is assumed that Scol and Scold are zero on entry.
  // Primary testers: test_find_Sm.m, test_find_Vhat.m

void SBuild( // Build covariance matrix of all the values of a VARMA time series
  char *uplo,  // in   Create all of SS (when "A") or lower part only (when "L")
  double S[],  // in   r×r·(p+1), = [S0...Sp], Si = Cov(x(t),x(t-i))
  double A[],  // in   r×r×p, A=[A1...Ap], autoregressive parameter matrices
  double G[],  // in   r×r·(q+1), = [G0...Gq], Gi = cov(y(t),x(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n can be any nonnegative value, even 0)
  double SS[]);// out  r·n×r·n, covariance of [x1'...xn']'

void CCBuild( // Build covariance between terms and shocks of VARMA time series
  double A[],  // in   r×r×p, A=[A1...Ap], autoregressive parameter matrices
  double C[],  // in   r×r·(q+1), = [C0...Cq], Ci = cov(x(t),eps(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n must be <= max(p+q)+1)
  double CC[]);// out  r·n×r·n, cov(x1'...xn',eps1'...epsn')

void unique01 ( // Find unique misssing patterns
  int miss[], // in,  r×n, miss[i,t] = 1 if i-th value at time t is missing
  int r,      // in,  row count of miss
  int n,      // in,  column count of miss, defining dimension of ipat and jpat
  int ipat[], // out  *npat-vector, ipat[j]=t if miss[.,t] has pattern number j
  int jpat[], // out  n-vector, jpat[t] = pattern number of miss[.,t]
  int *npat); // out  number of unique missing patterns
  // Primary tester: test_unique01.m (via t_unique01.c)

void XmuDeriv( // Derivatives of X minus mu.
  int r,         // in   dimension of each x(i)
  int n,         // in   length of time series
  int nPar,      // in   N. of parameters to diff. w.r.t., r^2·(p+q)+r·(r+3)/2
  int miss[],    // in   r×n, true if corresponding X-value is missing
  double xmud[]);// out  nO×nPar, derivatives of observed components of x - mu
  // Primary tester: test_misc.c

#endif
