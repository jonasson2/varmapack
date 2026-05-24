#include "BlasGateway.h"
#include "Lyapunov.h"
#include "error.h"
#include "printX.h"
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "varmapack_config.h"
#include "VarmaUtilities.h"
#include "VYW.h"

static int slicotCutoff(int p, int q);
static void SExtend(double A[], double G[], double S[], double Scol[], int p, int q,
                    int r, int n);

HIDDEN void FindC( // Calculate Ci = cov(x(t), eps(t-i))
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[])   // out  r×r×(q+1) with C0, C1...Cq
{
  int i, j, rr = r*r;
  double *Bj, *Cj, *Ai, *Cjmi;
  copy(rr, Sig, 1, C, 1);
  for (j=1; j<=q; j++) {
    Cj = C + j*rr;
    Bj = B + (j-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1, Bj, r, Sig, r, 0, Cj, r);
    for (i=1; i<=p && i<=j; i++) {
      Cjmi = C + (j-i)*rr;
      Ai = A + (i-1)*rr;
      gemm("NoT", "NoT", r, r, r, 1, Ai, r, Cjmi, r, 1, Cj, r);
    }
  }
}

HIDDEN void FindG( // Calculate Gi = cov(y(t), x(t-i))
  double B[],   // in   r×r×q, moving average parameter matrices
  double C[],   // in   r×r×(q+1), Ci = cov(x(t), eps(t-i))
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double G[])   // out  r×r×(q+1) with G0, G1...Gq
{
  int i, j, rr = r*r;
  double *Bi, *Cimj, *Gj;
  setzero(rr*(q+1), G);
  for (j=0; j<=q; j++) {
    Gj = G + j*rr;
    for (i=j; i<=q; i++) {
      Cimj = C + (i-j)*rr;
      if (i==0) {
        addmat("All", r, r, Cimj, r, Gj, r);
      }
      else {
        Bi = B + (i-1)*rr;
        gemm("NoT", "T", r, r, r, 1, Bi, r, Cimj, r, 1, Gj, r);
      }
    }
  }
}

HIDDEN bool FindS(double A[], double B[], double Sig[], int p, int q, int r,
                  double S[], double C[], double G[]) {
  int rc = slicotCutoff(p, q);
  if (rc > 0 && r >= rc) return LyapunovFactorizeSolve(A, B, Sig, p, q, r, S, C, G);
  return VYWFactorizeSolve(A, B, Sig, p, q, r, S, C, G);
}

HIDDEN bool SBuild(char *uplo, double S[], double A[], double G[], int p, int q, int r,
                   int n, double SS[]) {
  double *Scol = 0;
  double *SSj, *SSi;
  int j, m;
  if (n == 0) return true;
  m = imax(p+1, n);
  if (!ALLOC(Scol, (r*m)*r)) return false;
  SExtend(A, G, S, Scol, p, q, r, m);
  for (j = 0; j < n; j++) {
    SSj = SS + j*r*n*r + j*r;
    if (uplo[0] == 'A') {
      SSi = SSj + r*n*r;
      lacpy("All", r*(n-j), r, Scol, r*m, SSj, r*n);
      if (j < n-1) copytranspose(r*(n-j-1), r, Scol+r, r*m, SSi, r*n);
    }
    else {
      lacpy("Low", r*(n-j), r, Scol, r*m, SSj, r*n);
    }
  }
  FREE(Scol);
  return true;
}

HIDDEN void CCBuild( // Build covariance between terms and shocks
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double C[],  // in   r × r·(q+1), [C0...Cq], Ci = cov(x(t),eps(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series
  double CC[]) // out  r·n × r·n, cov(x1'...xn',eps1'...epsn')
{
  int i, j;
  double *CCj;
  setzero(r*n*r*n, CC);
  for (j=0; j<=q && j<n; j++) {
    lacpy("All", r, r, C+j*r*r, r, CC + j*r, r*n);
  }
  for (j=q+1; j<n; j++) {
    CCj = CC + j*r;
    for (i=1; i<=j && i<=p; i++) {
      gemm("N", "N", r, r, r, 1, A+(i-1)*r*r, r, CC+(j-i)*r, r*n, 1, CCj, r*n);
    }
  }
  for (j=1; j<n; j++) {
    CCj = CC + j*r*n*r + j*r;
    lacpy("All", r*(n-j), r, CC, r*n, CCj, r*n);
  }
}

static int slicotCutoff(int p, int q) {
  int pc, qc;
  // Average break-even r for choosing SLICOT over VYW; 0 means always use VYW.
  // The numbers are averages across Mac, XEON, Ubuntu; OpenBLAS, Accelerate, MKL.
  int rc[7][8] = {
    {16, 26, 0, 0, 0, 0, 0, 0}, {10, 13, 18, 21, 24, 24, 26, 28},
    {9, 13, 15, 17, 19, 21, 22, 23}, {11, 13, 14, 15, 16, 18, 18, 19},
    {11, 13, 14, 15, 15, 16, 17, 18}, {11, 13, 13, 14, 14, 15, 15, 16},
    {12, 13, 13, 13, 14, 14, 14, 15}
  };
  if (p <= 0 || q < 0) return 0;
  pc = p < 7 ? p : 7;
  qc = q < 7 ? q : 7;
  return rc[pc-1][qc];
}

static void SExtend(double A[], double G[], double S[], double Scol[], int p, int q,
                    int r, int n) {
  int iScol = n*r, i, j;
  double *Scolj, *Ai, *Scoli;
  for (j = 0; j < p+1; j++) {
    lacpy("All", r, r, S + j*r*r, r, Scol + j*r, iScol);
  }
  for (j = p+1; j < n; j++) {
    Scolj = Scol + j*r;
    if (j <= q) {
      lacpy("All", r, r, G + j*r*r, r, Scolj, iScol);
    }
    for (i = 0; i < p; i++) {
      Ai = A+i*r*r;
      Scoli = Scolj - (i+1)*r;
      gemm("N", "N", r, r, r, 1, Ai, r, Scoli, iScol, 1, Scolj, iScol);
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
