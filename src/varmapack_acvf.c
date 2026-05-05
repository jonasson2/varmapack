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
  double *G = 0;
  int r2 = r*r;
  xAssert(p >= 0 && q >= 0 && r > 0);
  xAssert(maxlag >= p);
  if (!ALLOC(G, r2*(q+1))) {
    return false;
  }
  if (!vpack_VYWFactorizeSolve(A, B, Sig, p, q, r, S, 0, G)) {
    FREE(G);
    return false;
  }
  if (p == 0) {
    int qcopy = imin(q, maxlag);
    copy(r2*(qcopy+1), G, 1, S, 1);
    if (qcopy < maxlag) setzero(r2*(maxlag - qcopy), S + r2*(qcopy+1));
    FREE(G);
    return true;
  }
  for (int j=p+1; j<=maxlag; j++) {
    double *Sj = S + j*r2;
    if (j <= q) {
      copy(r2, G + j*r2, 1, Sj, 1);
    }
    else {
      setzero(r2, Sj);
    }
    // AR-recursion, see method article, last line of eq. (10)
    for (int i = 1; i <= p; ++i) {
      double *Ai  = A + (i - 1)*r2;
      double *Sji = S + (j - i)*r2;
      gemm("N", "N", r, r, r, 1.0, Ai, r, Sji, r, 1.0, Sj, r);
    }
  }
  FREE(G);
  return true;
}
