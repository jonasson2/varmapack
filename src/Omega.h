// Omega.h:  Declaration of functions associated with the matrix Omega
#ifndef OMEGA_H
#define OMEGA_H

#include "BlasGateway.h"
void OmegaBuild ( // Build Omega = cov(w) with w' = (w1'...wn')
  double Su[],   // out  h·r×h·r, upper left part of Omega, where h = max(p,q)
  double Olow[], // out  (n-h)r×(q+1)r, block-diagonals of lower part of Omega
  double S[],    // in   r×r·(p+1), = [S0 S1...Sp] (lagged) covariances of xt
  double G[],    // in   r×r·(q+1), = [G0...Gq], Gi = cov(y(t),y(t-i))
  double W[],    // in   r×r·(q+1), = [W0...Wq], Wi = cov(y(t),x(t-i))
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of xt
  int n);        // in   length of series

void OmegaBuildDeriv ( // Find derivative of Omega
  double Sud[],   // out  r·h×r·h×nPar, derivatives of upper left part of Omega
  double Olowd[], // out  r(n-h)×r(q+1)×nPar, derivatives of lower part of Omega
  double Sd[],    // in   r×r×nPar×(p+1), derivatives of S0, S1...Sp
  double Gd[],    // in   r×r×nPar×(q+1), derivatives of G0, G1...Gq
  double Wd[],    // in   r×r×nPar×(q+1), derivatives of W0, W1...Wq
  int p,          // in   number of autoregressive terms
  int q,          // in   number of moving average terms
  int r,          // in   dimension of xt
  int nPar,       // in   total number of parameters
  int n);         // in   length of series

void OmegaRemoveMiss ( // Find Omega_o (Remove rows/columns of missing values)
  double Su[],   // in/out  Upper left part of Omega, in: r·h×r·h, out: mSu×mSu
  double Olow[], // in/out  Omega lower, in: (n-h)·r×(q+1)·r, out: mOlow×(q+1)·r
  double Sud[],  // in/out  In:r·h×r·h×nPar, out:mSu×mSu×nPar, derivatives of Su
  double Olowd[],// in/out  Derivatives of Olow, dimension: (same as Olow)×nPar 
  int h,         // in      Number of block rows in Su
  int q,         // in      Number of moving average terms in time series model
  int r,         // in      Dimension of each x(t)
  int n,         // in      Length of time series
  int nPar,      // in      Number of parameters differentiated w.r.t.
  int miss[],    // in      r×n, miss(i,t) = 1 if x(i,t) is missing, otherwise 0
  int ko[]);     // out     (n+1)-vec, ko[t] = n. of obs. values before time t+1

void OmegaFactor (
  double Su[],    // in/out  upper left part of Omega, dimension mSu×mSu
  double Olow[],  // in/out  mOlow×nOlow, block-diagonals of lower part of Omega
  int p,          // in      number of autoregressive terms
  int q,          // in      number of moving average terms
  int n,          // in      length of time series
  int ko[],       // in      ko[t] = N of observed values before time t+1, t<=n
  int *info);     // out     0 if ok, otherwise k for first nonpositive Ltt

void OmegaFactorDeriv ( // Derivative of Cholesky-factor of Omega or Omega_o
  double Sud[],   // in/out  derivatives of Su on entry, of Lu on return
  double Olowd[], // in/out  derivatives of Olow on entry, of Ll on return
  double Lu[],    // in      Cholesky factor of Su (from OmegaFactor)
  double Ll[],    // in      Cholesky factor of Olow (from OmegaFactor)
  int p,          // in      number of autoregressive terms
  int q,          // in      number of moving average terms
  int r,          // in      dimension of each time t observation, x(t)
  int n,          // in      length of time series
  int ko[],       // in      ko[t] = N of observed values before time t+1, t<=n
  int nPar);      // in      number of parameters to differentiate w.r.t.

void OmegaForward ( // Forward substitution with Cholesky-factor of Omega, L
  char *Ytransp,// in      String beginning with "T" if op(Y)=Y'; "N" if op(Y)=Y
  double Lu[],  // in      upper left part of L-factor of Omega
  double Ll[],  // in      lower part of L (block diagonals and subdiagonals)
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving average terms
  int ko[],     // in      ko[t]=number of observed values before time t+1, t<=n
  int n,        // in      length of time series
  double Y[],   // in/out  right hand side / solution
  int iY,       // in      leading dimension of Y
  int nY,       // in      number of columns in op(Y)
  int tmin,     // in      >0 used for partial solution (see note below)
  int tmax);    // in      =n-1 unless partial solution wanted

void OmegaForwardDeriv ( // Derivative of solution to L·op(X) = op(Y)
  char *Ytransp,// in      String beginning with "T" if op(Y)=Y'; "N" if op(Y)=Y
  double Lu[],  // in      upper left part of L-factor of Omega
  double Ll[],  // in      lower part of L (block diagonals and subdiagonals)
  double Lud[], // in      derivatives of Lu
  double Lld[], // in      derivatives of Ll
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving average terms
  int r,        // in      dimension of time series observation, x(k)
  int ko[],     // in      ko[k]=number of observed values before time k+1, k<=n
  int n,        // in      length of time series
  double X[],   // in      the solution from OmegaForward
  double Yd[],  // in/out  derivatives of Y on entry, of X on return
  int i1,       // in      leading dimension of X, 1st dimension of Yd
  int i2,       // in      2nd dimension of (X, Y and) Yd
  int nY,       // in      number of columns in op(Y)
  int nPar,     // in      number of parameters to diff. w.r.k., 3rd dim. of Yd
  int tmin,     // in      >0 used for partial solution (see note below)
  int tmax);    // in      =n-1 unless partial solution wanted

void OmegaBackSub ( // Back substitution with Cholesky-factor of Omega, L
  double Lu[],  // in      e×e, upper left part of L-factor of Omega
  double Ll[],  // in      (N-e)×(q+1)r, L-lower-part (block diags & subdiags)
  double Lud[], // in      e×e×nPar, derivative of Lu (or null if nPar=0)
  double Lld[], // in      (N-e)×(q+1)r×nPar, der. of Ll (or null if nPar=0)
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving average terms
  int r,        // in      length of each xt
  int ko[],     // in      ko[t]=number of observed values before time t+1, t<=n
  int n,        // in      length of time series
  int nY,       // in      number of columns of Y
  int tmax,     // in      index of last nonzero block-row of Y
  double Y[],   // in/out  k×nY, right hand side / solution
  double Yd[],  // in/out  k×nY×nPar, derivative of Y (or null if nPar=0)
  int iY,       // in      first declared dimension of Y and Yd
  int jY,       // in      second declared dimension of Yd
  int nPar);    // in      number of parameters or 0 for no derivative

void OmegaDiag ( // finds diagonal of Cholesky factor of Omega + derivatives
  double Lu[],      // in   upper left part of L-factor of Omega (square lower)
  double Ll[],      // in   lower part of L (block diagonals and subdiagonals)
  double Lud[],     // in   derivatives of Lu
  double Lld[],     // in   derivatives of Ll
  int p,            // in   number of autoregressive terms
  int q,            // in   number of moving average terms
  int r,            // in   dimension of time series observation, x(t)
  int n,            // in   length of time series
  int ko[],         // in   ko[t]=number of observed values before time t+1,t<=n
  int nPar,         // in   number of parameters (set to 0 for no derivatives)
  double Ldiag[],   // out  diagonal of L
  double Ldiagd[]); // out  derivatives Ldiag

double OmegaLogdet ( // returns log(det(Omega)) or log(det(Omega_o))
  double Lu[],      // in   upper left part of L-factor of Omega (square lower)
  double Ll[],      // in   lower part of L (block diagonals and subdiagonals)
  double Lud[],     // in   derivatives of Lu
  double Lld[],     // in   derivatives of Ll
  int p,            // in   number of autoregressive terms
  int q,            // in   number of moving average terms
  int r,            // in   dimension of time series observation, x(t)
  int n,            // in   length of time series
  int ko[],         // in   ko[t]=number of observed values before time t+1,t<=n
  int nPar,         // in   number of parameters (set to 0 for no derivatives)
  double logdetd[]);// out  gradient of log(det(Omega))

void OmegaLtl ( // LTL-factorize Omega, or, for missing values, Omega_o
  double Lu[],  // in/out  e×e, upper left part of Omega / L
  double Ll[],  // in/out  mLl×nLl, block-diagonals of lower part of Omega / L
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving average terms
  int r,        // in      dimension of each time t observation, x(t)
  int n,        // in      length of time series
  int ko[],     // in      ko[t] = N of observed values before time t+1, t<=n
  int *info,    // out     0 if ok, otherwise t for first nonpositive Ltt
  double Lud[], // in/out  e×e×nPar, derivatives of Lu, or null if nPar=0
  double Lld[], // in/out  mLl×nLl×nPar, derivatives of Ll or null if nPar=0
  int nPar);    // in      number of paramters to differentiate w.r.t.

#endif
