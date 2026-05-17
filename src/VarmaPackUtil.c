#include "BlasGateway.h"
#include "error.h"
#include "printX.h"
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "varmapack_config.h"
#include "VarmaUtilities.h"

HIDDEN void FindCG ( // Calculate the Ci and Gi matrices for VARMASIM
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[],   // out  r×r×(q+1) with C0, C1...Cq where Ci = cov(x(t),eps(t-i))
  double G[])   // out  r×r×(q+1) with G0, G1...Gq where Gi = cov(y(t),x(t-i))
{
  int i, j, rr = r*r;
  double *Bj, *Cj, *Gj, *Ai, *Bi, *Cimj, *Cjmi;
  //
  // CALCULATE C0,... Cq:
  copy(rr, Sig, 1, C, 1);
  for (j=1; j<=q; j++) {
    Cj = C + j*rr;
    Bj = B + (j-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1.0, Bj, r, Sig, r, 0.0, Cj, r);
    for (i=1; i<=p && i<=j; i++) {
      Cjmi = C + (j-i)*rr;
      Ai = A + (i-1)*rr;
      gemm("NoT", "NoT", r, r, r, 1.0, Ai, r, Cjmi, r, 1.0, Cj, r);
    }
  }
  // NOW TURN ATTENTION TO THE Gj-MATRICES:
  setzero(rr*(q+1), G);
  for (j=0; j<=q; j++) {
    Gj = G + j*rr;
    Cj = C + j*rr;
    for (i=j; i<=q; i++) {
      Cimj = C + (i-j)*rr;
      if (i==0) {
        addmat("All", r, r, Cimj, r, Gj, r);
      }
      else {
        Bi = B + (i-1)*rr;
        gemm("NoT", "T", r, r, r, 1.0, Bi, r, Cimj, r, 1.0, Gj, r);
      }
    }
  }
}

HIDDEN void FindPsi(double *A, double *B, double *Psi, int p, int q, int r) {
  // Prepare
  int h = imax(p, q), rr = r*r, hr = h*r, i, j, k, l;
  double *Aflp = 0, *Psi_jj;
  if (h == 0) return;
  if (p > 0) xAssert(ALLOC(Aflp, p*rr));
  flipmat(A, Aflp, r, p);
  laset("All", hr, hr, 0.0, 0.0, Psi, hr);
  
  // Compute first block column
  setEye(r, Psi, hr);
  for (i=1; i<h; i++) {
    l = i*r;
    k = p*r - i*r;
    if (i <= p) {
      gemm("NoT", "NoT", r, r, l, 1.0, Aflp + k*r, r, Psi, hr, 1.0, Psi + i*r, hr);
    }
    if (i <= q) {
      addmat("All", r, r, B + (i - 1)*rr, r, Psi + i*r, hr);
    }
  }
  // Copy to remaining lower block columns
  for (j=1; j<h; j++) {
    Psi_jj = Psi + j*hr*r + j*r;
    lacpy("All", hr - j*r, r, Psi, hr, Psi_jj, hr);
  }
  FREE(Aflp);
}

HIDDEN void FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h) {
  double *LSig, *Psi_hat_kk, *W;
  int hr = h*r, rr = r*r, k, nrow;
  bool triangular;
  xAssert(ALLOC(LSig, rr));
  xAssert(ALLOC(W, hr*r));
  xAssert(psdFactor(Sig, r, LSig, &triangular) == VARMAPACK_OK);
  lacpy("Low", hr, hr, Psi, hr, Psi_hat, hr);
  for (k=0; k<h; k++) {
    Psi_hat_kk = Psi_hat + k*h*rr + k*r;
    nrow = hr - k*r;
    if (triangular) {
      trmm("Right", "Low", "NoT", "NonUnit", nrow, r, 1, LSig, r, Psi_hat_kk, hr);
    }
    else {
      lacpy("All", nrow, r, Psi_hat_kk, hr, W, nrow);
      gemm("NoT", "NoT", nrow, r, r, 1, W, nrow, LSig, r, 0, Psi_hat_kk, hr);
    }
  }
  FREE(W);
  FREE(LSig);
}

HIDDEN varmapack_error psdFactor(double Sig[], int r, double L[], bool *triangular) {
  int info, rank, rr = r*r;
  int *piv = 0;
  double *C = 0;
  double *work = 0;
  double *recon = 0;
  double tol;
  varmapack_error error = VARMAPACK_OK;
  if (triangular) *triangular = false;
  if (!ALLOC(C, rr)) goto alloc_fail;
  if (!ALLOC(piv, r)) goto alloc_fail;
  if (!ALLOC(work, 2*r)) goto alloc_fail;
  if (!ALLOC(recon, rr)) goto alloc_fail;
  lacpy("All", r, r, Sig, r, C, r);
  potrf("Low", r, C, r, &info);
  if (info < 0) { error = VARMAPACK_INTERNAL; goto fail; }
  if (info == 0) {
    setzero(rr, L);
    lacpy("Low", r, r, C, r, L, r);
    if (triangular) *triangular = true;
    goto fail;
  }
  lacpy("All", r, r, Sig, r, C, r);
  pstrf("Low", r, C, r, piv, &rank, -1, work, &info);
  if (info < 0) { error = VARMAPACK_INTERNAL; goto fail; }
  setzero(rr, L);
  for (int k=0; k<rank; k++) {
    for (int i=k; i<r; i++) {
      L[piv[i] + k*r] = C[i + k*r];
    }
  }
  setzero(rr, recon);
  if (rank > 0) syrk("Low", "NoT", r, rank, 1, L, r, 0, recon, r);
  copylowertoupper(r, recon, r);
  tol = 100*r*lamch("E");
  if (relabsdiff(Sig, recon, rr) > tol) {
    error = VARMAPACK_NOT_POSITIVE_SEMIDEFINITE;
  }
fail:
  FREE(recon);
  FREE(work);
  FREE(piv);
  FREE(C);
  return error;
alloc_fail:
  error = VARMAPACK_ALLOCATION;
  goto fail;
}
