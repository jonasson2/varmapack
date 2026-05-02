#include "BlasGateway.h"
#include "allocate.h"
#include "printX.h"
#include "VarmaUtilities.h"

void vpack_FindCG ( // Calculate the Ci and Gi matrices for VARMASIM
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

void vpack_FindPsi(double *A, double *B, double *Psi, int p, int q, int r) {
  // Prepare
  int h = imax(p, q), rr = r*r, hr = h*r, i, j, k, l;
  double *Aflp, *Psi_jj;
  if (h == 0) return;
  allocate(Aflp, p*rr);
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

void vpack_FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h) {
  double *LSig, *Psi_hat_kk;
  int info, hr = h*r, rr = r*r, k, nrow;
  allocate(LSig, rr);
  lacpy("Low", r, r, Sig, r, LSig, r);
  potrf("Low", r, LSig, r, &info);
  xAssert(info == 0);
  lacpy("Low", hr, hr, Psi, hr, Psi_hat, hr);
  for (k=0; k<h; k++) {
    Psi_hat_kk = Psi_hat + k*h*rr + k*r;
    nrow = hr - k*r;
    trmm("Right", "Low", "NoT", "NonUnit", nrow, r, 1.0, LSig, r, Psi_hat_kk, hr);
  }
  FREE(LSig);
}
