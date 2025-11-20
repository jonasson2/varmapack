// varmapack_acvf — compute the theoretical autocovariance function (varmapack_acvf) of the stationary
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
//         theoretical autocovariance function (varmapack_acvf) of the process {x_t},
//         S[:,:,k] = Cov(x_t, x_{t−k}),  k = 0,…,maxlag
//   ok  : 1 on success, 0 on failure (non-
//
// Returns true on success, false on failure.

#include <stdbool.h>
#include "varmapack.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "varmapack_VYW.h"

bool varmapack_acvf(double A[], double B[], double Sig[], int p, int q, int r, double Gamma[],
              int maxlag) {
  double *S = Gamma;
  xAssert(p >= 0 && q >= 0 && r > 0);
  if (!vpack_VYWFactorizeSolve(A, B, Sig, p, q, r, S, 0, 0)) {
    return false;
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
  return true;
}
