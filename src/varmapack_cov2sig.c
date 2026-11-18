// Recover innovation covariance from lag-zero covariance using inverse VYW.
// Derivation: "Recovering innovation covariance" in python/docs/math-v0.2-draft.md.

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
  double *As = 0, *Bs = 0;
  int *ex = 0;
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
    }
  }
  if (!ALLOC(V, rr) || !ALLOC(candidate, rr) || !ALLOC(factor, rr) ||
      !ALLOC(ex, r) || (p && !ALLOC(As, nS-rr)) || (q && !ALLOC(Bs, nG-rr)) ||
      !ALLOC(S, nS) || !ALLOC(C, nG) || !ALLOC(G, nG)) {
    fail_error("allocation failed");
    goto done;
  }
  // Solve for D^-1 X with power-of-two scales close to each marginal SD.
  for (int i=0; i<r; i++) {
    double v = cov[i+i*r];
    if (v < 0) {
      fail_error("covariance is not positive semidefinite: negative marginal variance");
      goto done;
    }
    ex[i] = v == 0 ? 0 : ilogb(v)/2;
  }
  for (int k=0; k<imax(p,q); k++) {
    for (int j=0; j<r; j++) {
      for (int i=0; i<r; i++) {
        int u = k*rr+j*r+i, shift = ex[j]-ex[i];
        if (k < p) As[u] = scalbn(A[u], shift);
        if (k < q) Bs[u] = scalbn(B[u], shift);
        if ((k < p && (!isfinite(As[u]) || (A[u] && !As[u]))) ||
            (k < q && (!isfinite(Bs[u]) || (B[u] && !Bs[u])))) {
          fail_error("scaled coefficient is not representable");
          goto done;
        }
      }
    }
  }
  double rho = varmapack_specrad(As, r, p);
  if (!isfinite(rho) || rho >= 1) {
    if (!varmapack_last_error()) fail_error("nonstationary model");
    goto done;
  }
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      double value = scalbn(cov[i+j*r], -ex[i]-ex[j]);
      if (cov[i+j*r] != 0 && (cov[i+i*r] == 0 || cov[j+j*r] == 0)) {
        fail_error("covariance is not positive semidefinite: nonzero zero-variance row");
        goto done;
      }
      if (!isfinite(value) || (value == 0 && cov[i+j*r] != 0)) {
        fail_error("covariance scale range is too large");
        goto done;
      }
      V[i+j*r] = V[j+i*r] = value;
    }
  }
  if (!psdFactor(V, r, factor, 0)) goto done;
  if (!VYWInverse(As, Bs, V, p, q, r, candidate)) goto done;
  bool triangular;
  if (!psdFactor(candidate, r, factor, &triangular)) goto done;
  if (!triangular) {
    // Return the PSD matrix actually accepted by the rank-revealing factorization.
    syrk("Low", "NoT", r, r, 1, factor, r, 0, candidate, r);
    copylowertoupper(r, candidate, r);
  }
  if (!VYWFactorizeSolve(As, Bs, candidate, p, q, r, S, C, G)) {
    fail_error("forward covariance check failed");
    goto done;
  }
  double tol = 1000*N*lamch("E");
  for (int i=0; i<rr; i++) {
    if (!isfinite(S[i]) || fabs(S[i]-V[i]) > tol) {
      fail_error("forward covariance check failed");
      goto done;
    }
    double value = scalbn(candidate[i], ex[i%r]+ex[i/r]);
    if (!isfinite(value) || (value == 0 && candidate[i] != 0)) {
      fail_error("innovation covariance is outside the representable range");
      goto done;
    }
    candidate[i] = value;
  }
  copy(rr, candidate, 1, Sig, 1);
  ok = true;
done:
  FREE(Bs); FREE(As); FREE(ex);
  FREE(G); FREE(C); FREE(S); FREE(factor); FREE(candidate); FREE(V);
  return ok;
}
