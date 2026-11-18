// Inverse VYW: see "Recovering innovation covariance" in python/docs/math-v0.2-draft.md.

#include "VYW.h"

#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "error.h"

static bool vyw_build_matrix(double A[], double F[], int p, int r);
static void vyw_rhs(double A[], double Y[], int nY, int p, int r, double g[]);
static bool vyw_solve(double A[], double LU[], double S[], double Y[], int nrhs, int nY,
                      int piv[], int p, int r);
static void kronecker(int n, double alpha, double A[], double B[], double C[]);
static void kronI(int n, double alpha, double A[], double C[]);

HIDDEN bool VYWFactorizeSolve(double A[], double B[], double Sig[], int p, int q, int r,
                        double S[], double C[], double G[])
{
  xAssert((A != 0 || p == 0) && Sig != 0 && S != 0);
  xAssert(p >= 0 && q >= 0 && r > 0);
  xAssert(B != 0 || q == 0);
  xAssert(C != 0 && G != 0);
  int rr = r*r;
  double *vywFactors = 0;
  int *piv = 0;
  FindC(A, B, Sig, p, q, r, C);
  FindG(B, C, q, r, G);
  if (p == 0) {
    copy(rr, G, 1, S, 1);
    return true;
  }
  int nVYW = rr*p - r*(r-1)/2;
  if (!ALLOC(vywFactors, nVYW*nVYW)) goto fail;
  if (!ALLOC(piv, nVYW)) goto fail;
  int info;
  if (!vyw_build_matrix(A, vywFactors, p, r)) goto fail;
  getrf(nVYW, nVYW, vywFactors, nVYW, piv, &info);
  if (info == 0) {
    if (!vyw_solve(A, vywFactors, S, G, 1, q+1, piv, p, r)) goto fail;
  }
  FREE(piv);
  FREE(vywFactors);
  return info == 0;
fail:
  FREE(piv);
  FREE(vywFactors);
  return false;
}

HIDDEN bool VYWInverse(double A[], double B[], double cov[], int p, int q, int r,
                        double Sig[]) {
  int rr = r*r, m = r*(r+1)/2, N = m + rr*imax(p-1, 0);
  double *F = 0, *rhs = 0, *g = 0, *E = 0, *C = 0, *G = 0, *work = 0;
  int *piv = 0, *iwork = 0;
  int info, col = 0;
  double anorm = 0, rcond;
  bool ok = false;
  if (!ALLOC(F, N*N) || !ALLOC(rhs, N) || !ALLOC(g, imax(N, rr)) ||
      !ALLOC(E, rr) || !ALLOC(C, rr*(q+1)) || !ALLOC(G, rr*(q+1)) ||
      !ALLOC(piv, N) || !ALLOC(work, 4*N) || !ALLOC(iwork, N)) goto alloc_fail;
  if (p > 0) {
    if (!vyw_build_matrix(A, F, p, r)) goto alloc_fail;
    for (int j=0; j<r; j++) {
      for (int i=j; i<r; i++) g[col++] = cov[i+j*r]*(i == j ? 0.5 : 1);
    }
    // Exchange the known Gamma_0 and unknown Sigma: [R, -F1] u = F0 v.
    gemv("N", N, m, 1, F, N, g, 1, 0, rhs, 1);
    scal(N*(N-m), -1, F+N*m, 1);
  }
  else {
    for (int j=0; j<r; j++) {
      for (int i=j; i<r; i++) rhs[col++] = cov[i+j*r];
    }
  }
  col = 0;
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) {
      E[i+j*r] = E[j+i*r] = 1;
      FindC(A, B, E, p, q, r, C);
      FindG(B, C, q, r, G);
      vyw_rhs(A, G, q+1, p, r, g);
      copy(N, g, 1, F+N*col++, 1);
      E[i+j*r] = E[j+i*r] = 0;
    }
  }
  for (int j=0; j<N; j++) {
    double sum = 0;
    for (int i=0; i<N; i++) {
      if (!isfinite(F[i+j*N])) goto nonfinite;
      sum += fabs(F[i+j*N]);
    }
    anorm = fmax(anorm, sum);
    if (!isfinite(rhs[j])) goto nonfinite;
  }
  if (!isfinite(anorm)) goto nonfinite;
  getrf(N, N, F, N, piv, &info);
  if (info > 0) {
    fail_error("singular inverse covariance system");
    goto done;
  }
  if (info < 0) goto internal_fail;
  gecon(N, F, anorm, &rcond, work, iwork, &info);
  if (info != 0) goto internal_fail;
  if (!isfinite(rcond) || rcond <= 100*N*lamch("E")) {
    fail_error("ill-conditioned inverse covariance system");
    goto done;
  }
  getrs("N", N, 1, F, N, piv, rhs, N, &info);
  if (info != 0) goto internal_fail;
  for (int i=0; i<N; i++) {
    if (!isfinite(rhs[i])) goto nonfinite;
  }
  col = 0;
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) Sig[i+j*r] = Sig[j+i*r] = rhs[col++];
  }
  ok = true;
  goto done;
alloc_fail:
  fail_error("allocation failed");
  goto done;
nonfinite:
  fail_error("nonfinite inverse covariance system");
  goto done;
internal_fail:
  fail_error("internal error");
done:
  FREE(iwork); FREE(work); FREE(piv); FREE(G); FREE(C); FREE(E);
  FREE(g); FREE(rhs); FREE(F);
  return ok;
}

static bool vyw_build_matrix(double A[], double F[], int p, int r) {
  int rr = r*r;
  int N = rr*p - r*(r-1)/2;
  int N0 = r*(r+1)/2;
  int Nc = rr*(p-1);
  int Nr = rr;
  double *F0r, *F0c, *F00, *F11, *Ap, *Fii, *Apmi, *F0ri, *F00jj, *FKK;
  double *FKJ, *Ai, *Aj, *FiimjKK, *FijmiKJ, *F0ciKK, *F0ciKJ, *F00KK1;
  int i, j, k, i1, j1, I, J;
  if (p == 0) return true;
  F0r = F0c = F00 = 0;
  if (Nc > 0 && !ALLOC(F0r, Nr*Nc)) goto fail;
  if (Nc > 0 && !ALLOC(F0c, Nc*Nr)) goto fail;
  if (!ALLOC(F00, Nr*Nr)) goto fail;
  F11 = p > 1 ? F + N0*(1 + N) : F;
  setzero(N*N, F);
  Ap = A + rr*(p-1);
  Fii = F11;
  for (i = 0; i < p-1; i++) {
    for (j = 0; j < rr; j++) {
      *Fii = 1.0;
      Fii += (N+1);
    }
    Ai = A + rr*i;
    Apmi = A + rr*(p-i-2);
    F0ri = F0r + Nr*rr*i;
    kronI(r, -1.0, Ai, F0ri);
    kronecker(r, -1.0, Ap, Apmi, F0ri);
  }
  for (j = 0; j < rr; j++) {
    F00jj = F00 + (1 + Nr)*j;
    *F00jj = 1.0;
  }
  kronecker(r, -1.0, Ap, Ap, F00);
  for (k = 0; k < r; k++) {
    if (p > 1) {
      FKK = F11 + r*k + N*r*k;
      FKJ = F11 + r*k + N*k;
      for (i = 0; i < p-1; i++) {
        for (j = 0; j < i; j++) {
          Aj = A + rr*j;
          FiimjKK = FKK + rr*i + N*rr*(i-j-1);
          subtractmat(r, r, Aj, r, FiimjKK, N);
        }
        for (j = i+1; j < p; j++) {
          Aj = A + rr*j;
          FijmiKJ = FKJ + rr*i + N*rr*(j-i-1);
          subtractmat(r, r, Aj, r, FijmiKJ, N*r);
        }
        Ai = A + rr*i;
        F0ciKK = F0c + rr*i + r*k + Nc*r*k;
        F0ciKJ = F0c + rr*i + r*k + Nc*k;
        subtractmat(r, r, Ai, r, F0ciKK, Nc);
        subtractmat(r, r, Ai, r, F0ciKJ, Nc*r);
      }
    }
    for (int k1 = 0; k1 < r; k1++) {
      F00KK1 = F00 + r*k + Nr*r*k1;
      ger(r, r, -1.0, Ap + r*k1, 1, Ap + k, r, F00KK1, Nr);
    }
    *(F00 + (1 + Nr)*(r*k + k)) += 1;
  }
  j1 = 0;
  for (j = 0; j < r; j++) {
    J = r*j + j;
    i1 = 0;
    for (i = 0; i < r; i++) {
      I = i*r + i;
      lacpy("All", r-i, r-j, F00+I+Nr*J, Nr, F+i1+N*j1, N);
      i1 += r-i;
    }
    if (Nc > 0) {
      lacpy("All", r-j, Nc, F0r+J, Nr, F+N0*N+j1, N);
      lacpy("All", Nc, r-j, F0c+Nc*J, Nc, F+N0+N*j1, N);
    }
    j1 += r-j;
  }
  FREE(F0r);
  FREE(F0c);
  FREE(F00);
  return true;
fail:
  FREE(F0r);
  FREE(F0c);
  FREE(F00);
  return false;
}

// g has room for max(r*r, m + max(p-1,0)*r*r); diagonal entries are not halved.
static void vyw_rhs(double A[], double Y[], int nY, int p, int r, double g[]) {
  int rr = r*r, m = r*(r+1)/2, N = m + rr*imax(p-1, 0);
  setzero(imax(N, rr), g);
  copytranspose(r, r, Y, r, g, r);
  if (p > 0 && nY > p) {
    gemm("N", "T", r, r, r, 1, Y+rr*p, r, A+rr*(p-1), r, 1, g, r);
  }
  int k = 0;
  for (int j=0; j<r; j++) {
    for (int i=j; i<r; i++) g[k++] = g[i+j*r];
  }
  setzero(N-m, g+m);
  for (int j=1; j<nY && j<p; j++) copy(rr, Y+rr*j, 1, g+m+(j-1)*rr, 1);
}

static bool vyw_solve(double A[], double LU[], double S[], double Y[], int nrhs, int nY,
                      int piv[], int p, int r)
{
  xAssert(nrhs == 1);
  int rr = r*r;
  double *g, *g1;
  int mg = r*r*(p+1);
  int N = r*r*p - r*(r-1)/2;
  if (!ALLOC(g, mg*nrhs)) return false;
  g1 = g + r*(r+1)/2;
  vyw_rhs(A, Y, nY, p, r, g);
  if (p > 0) {
    int info;
    getrs("NoT", N, nrhs, LU, N, piv, g, mg, &info);
    xAssert(info == 0);
  }
  if (p > 0) {
    for (int k = 0; k < nrhs; k++) {
      double *Sk = S + rr*k;
      double *gptr = g + mg*k;
      for (int i = 0; i < r; i++) {
        double *Sii0k = Sk + (1+r)*i;
        copy(r-i, gptr, 1, Sii0k, 1);
        copy(r-i-1, gptr+1, 1, Sii0k+r, r);
        *Sii0k *= 2;
        gptr += r-i;
      }
    }
  }
  for (int j = 1; j < p; j++) {
    double *gj = g1 + (j-1)*rr;
    lacpy("All", rr, nrhs, gj, mg, S + rr*nrhs*j, rr);
  }
  // In alignment with the original vauto Matlab code and TOMS article, Sp is
  // computed here, although it could also be computed by SExtend.
  double *Sp = S + rr*nrhs*p;
  if (nY > p) {
    lacpy("All", rr, nrhs, Y+rr*nrhs*p, rr, Sp, rr);
  }
  else {
    setzero(rr*nrhs, Sp);
  }
  for (int k = 0; k < nrhs; k++) {
    double *Skp = Sp + rr*k;
    for (int i = 1; i <= p; i++) {
      double *Skpmi = S + rr*k + rr*nrhs*(p-i);
      double *Ai = A + rr*(i-1);
      gemm("NoT", "NoT", r, r, r, 1.0, Ai, r, Skpmi, r, 1.0, Skp, r);
    }
  }
  FREE(g);
  return true;
}

static void kronecker(int n, double alpha, double A[], double B[], double C[])
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double *Aij = A + i + n*j;
      double *Cij = C + n*i + n*n*n*j;
      for (int k = 0; k < n; k++) {
        double *Bk = B + n*k;
        double *Cijk = Cij + n*n*k;
        axpy(n, alpha*(*Aij), Bk, 1, Cijk, 1);
      }
    }
  }
}

static void kronI(int n, double alpha, double A[], double C[])
{
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double *Aij = A + i + n*j;
      double *Cij = C + n*i + n*n*n*j;
      double aa = alpha*(*Aij);
      for (int k = 0; k < n; k++) {
        double *Cijkk = Cij + (1 + n*n)*k;
        *Cijkk += aa;
      }
    }
  }
}
