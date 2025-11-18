// varmapack_covar — covariance utilities for VARMA models
//
// varmapack_acvf computes the theoretical autocovariance function (varmapack_acvf) of the stationary
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
  if (!varmapack_VYWFactorizeSolve(A, B, Sig, p, q, r, S, 0, 0)) {
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

void varmapack_autocov(const char *transp, const char *norm, int r, int n,
                       double X[], int maxlag, double C[]) {
  xAssert(transp != 0 && norm != 0);
  bool ML = (norm[0] == 'M' || norm[0] == 'm');
  bool CORRECTED = (norm[0] == 'C' || norm[0] == 'c');
  xAssert(ML || CORRECTED);
  bool TRANSP = !(transp[0] == 'N' || transp[0] == 'n');
  xAssert((!TRANSP && (transp[0] == 'N' || transp[0] == 'n')) || (TRANSP && (transp[0] == 'T' || transp[0] == 't')));
  xAssert(maxlag <= n-1);
  double *Xc;
  double *mu;
  allocate(Xc, r*n);
  allocate(mu, r);
  copy(r*n, X, 1, Xc, 1);
  if (TRANSP) {
    meanmat("N", n, r, Xc, n, mu);
    for (int i = 0; i < n; i++) {
      double *row = Xc + i;
      axpy(r, -1.0, mu, 1, row, n);
    }
  }
  else {
    meanmat("T", r, n, Xc, r, mu);
    for (int j = 0; j < n; j++) {
      double *col = Xc + j*r;
      axpy(r, -1.0, mu, 1, col, 1);
    }
  }
  double *Y = Xc;
  for (int k = 0; k <= maxlag; k++) {
    double fctr = ML ? 1.0/n : 1.0/(n - k);
    double *Ck = C + r*r*k;
    if (TRANSP) {
      double *Z = Xc + k;
      gemm("T", "N", r, r, n-k, fctr, Y, n, Z, n, 0.0, Ck, r);
    }
    else {
      double *Z = Xc + r*k;
      gemm("N", "T", r, r, n-k, fctr, Y, r, Z, r, 0.0, Ck, r);
    }
  }
  freem(mu);
  freem(Xc);
}
