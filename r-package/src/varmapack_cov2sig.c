// Recover innovation covariance from lag-zero covariance using inverse VYW.
// Derivation: ../output/pdf/cov2sig.tex and cov2sig.pdf.

#include <math.h>
#include "varmapack.h"
#include "VYW.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"

bool varmapack_cov2sig(double A[], double B[], double cov[], int p, int q, int r,
                       double Sig[]) {
  double *V = 0, *candidate = 0, *factor = 0, *S = 0, *C = 0, *G = 0;
  int rr, nS, nG, square, worksize;
  double scale = 0;
  bool ok = false;
  clear_error();
  if (p < 0 || q < 0 || r <= 0 || cov == 0 || Sig == 0 ||
      (p > 0 && A == 0) || (q > 0 && B == 0)) {
    return fail_error("invalid argument");
  }
  if (p == INT_MAX || q == INT_MAX || !intProduct(r, r, &rr) ||
      !intProduct(rr, p+1, &nS) || !intProduct(rr, q+1, &nG) ||
      (p > 0 && !intProduct(rr, rr, &square))) {
    return fail_error("problem size too large");
  }
  int N = r*(r+1)/2 + rr*imax(p-1, 0);
  if (!intProduct(N, N, &square) || !intProduct(4, N, &worksize)) {
    return fail_error("problem size too large");
  }
  for (int i=0; i<nG-rr; i++) {
    if (!isfinite(B[i])) return fail_error("nonfinite coefficient");
  }
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      if (!isfinite(cov[i+j*r])) return fail_error("nonfinite covariance");
      scale = fmax(scale, fabs(cov[i+j*r]));
    }
  }
  double rho = varmapack_specrad(A, r, p);
  if (!isfinite(rho)) {
    return varmapack_last_error() ? false : fail_error("nonfinite spectral radius");
  }
  if (rho >= 1) return fail_error("nonstationary model");
  if (!ALLOC(V, rr) || !ALLOC(candidate, rr) || !ALLOC(factor, rr) ||
      !ALLOC(S, nS) || !ALLOC(C, nG) || !ALLOC(G, nG)) {
    fail_error("allocation failed");
    goto done;
  }
  // Normalize so PSD tolerances do not depend on the units of the covariance.
  if (scale == 0) scale = 1;
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double value = cov[i+j*r]/scale;
      if (value == 0 && cov[i+j*r] != 0) {
        fail_error("covariance scale range is too large");
        goto done;
      }
      V[i+j*r] = V[j+i*r] = value;
    }
  }
  if (!psdFactor(V, r, factor, 0)) goto done;
  if (!VYWInverse(A, B, V, p, q, r, candidate)) goto done;
  if (!psdFactor(candidate, r, factor, 0)) goto done;
  if (!VYWFactorizeSolve(A, B, candidate, p, q, r, S, C, G)) {
    fail_error("forward covariance check failed");
    goto done;
  }
  double tol = 1000*N*lamch("E");
  for (int i=0; i<rr; i++) {
    if (!isfinite(S[i]) || fabs(S[i]-V[i]) > tol) {
      fail_error("forward covariance check failed");
      goto done;
    }
    double value = candidate[i]*scale;
    if (!isfinite(value) || (value == 0 && candidate[i] != 0)) {
      fail_error("innovation covariance is outside the representable range");
      goto done;
    }
    candidate[i] = value;
  }
  copy(rr, candidate, 1, Sig, 1);
  ok = true;
done:
  FREE(G); FREE(C); FREE(S); FREE(factor); FREE(candidate); FREE(V);
  return ok;
}
