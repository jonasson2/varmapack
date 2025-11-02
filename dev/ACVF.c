// ACVF — compute the theoretical autocovariance function (ACVF) of the stationary
// time-series process {x_t} defined by a VARMA(p,q) model, using the vector Yule–Walker
// (VYW) equations.
//
// Inputs:
//   A     : r×r×p array of autoregressive coefficient matrices
//   B     : r×r×q array of moving-average coefficient matrices
//   Sig   : r×r innovation (shock) covariance matrix
//   p,q,r : model orders and system dimension
//
// Output:
//   S   : r×r×(maxlag+1) array containing the first (maxlag+1) matrices of the
//         theoretical autocovariance function (ACVF) of the process {x_t},
//         S[:,:,k] = Cov(x_t, x_{t−k}),  k = 0,…,maxlag
//   ok  : 1 on success, 0 on failure (non-
//
// Returns true on success, false on failure.

#include <stdbool.h>
#include "ACVF.h"
#include "BlasGateway.h"
#include "allocate.h"
#include "VYW.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"

bool ACVF(double A[], double B[], double Sig[], int p, int q, int r, double Gamma[],
              int maxlag) {
  int *piv, info, nVYW;
  double *vywFactors, *C, *G, *W, *S;
  S = Gamma;
  xAssert(p >= 0 && q >= 0 && r > 0);
  nVYW = (p == 0) ? 0 : (r*r*p - r*(r-1)/2);
  allocate(vywFactors, nVYW*nVYW);
  allocate(piv, nVYW);
  if (p > 0) {
    VYWFactorize(A, vywFactors, piv, p, r, &info);
    if (info != 0) {
      freem(piv);
      freem(vywFactors);
      return false;
    }
  }
  allocate(C, r*r*(q+1));
  allocate(G, r*r*(q+1));
  allocate(W, r*r*(q+1));
  FindCGW(A, B, Sig, p, q, r, C, G, W);
  if (p > 0) {
    VYWSolve(A, vywFactors, S, G, 1, q+1, piv, p, r);
  }
  else {
    copy(r*r, G, 1, S, 1);
  }
  for (int j=p+1; j<=maxlag; j++) {
    double *Sj = S + (j - 1) * r*r;
    setzero(r*r, Sj);
    // AR-recursion, see method article, last line of eq. (10)
    for (int i = 1; i <= p; ++i) {
      double *Ai  = A + (i - 1) * r*r;
      double *Sji = S + (j - i - 1) * r*r;
      gemm("N", "N", r, r, r, 1.0, Ai, r, Sji, r, 1.0, Sj, r);
    }
  }
  freem(W);
  freem(G);
  freem(C);
  freem(piv);
  freem(vywFactors);
  return true;
}
