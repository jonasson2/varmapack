#include "varmapack_VYW.h"

#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "allocate.h"
#include "error.h"

static void vyw_factorize(double A[], double LU[], int piv[], int p, int r, int *info);
static void vyw_solve(double A[], double LU[], double S[], double Y[], int nrhs, int nY,
                      int piv[], int p, int r);
static void kronecker(int n, double alpha, double A[], double B[], double C[]);
static void kronI(int n, double alpha, double A[], double C[]);

bool vpack_VYWFactorizeSolve(double A[], double B[], double Sig[],
                        int p, int q, int r,
                        double S[], double C[], double G[])
{
  xAssert(A != 0 && Sig != 0 && S != 0);
  xAssert(p >= 0 && q >= 0 && r > 0);
  xAssert(B != 0 || q == 0);

  int rr = r*r;
  double *Cbuf = C;
  double *Gbuf = G;
  bool ownC = false;
  bool ownG = false;
  if (Cbuf == 0) {
    allocate(Cbuf, rr*(q+1));
    ownC = true;
  }
  if (Gbuf == 0) {
    allocate(Gbuf, rr*(q+1));
    ownG = true;
  }

  vpack_FindCG(A, B, Sig, p, q, r, Cbuf, Gbuf);

  if (p == 0) {
    copy(rr, Gbuf, 1, S, 1);
    if (ownG) freem(Gbuf);
    if (ownC) freem(Cbuf);
    return true;
  }

  int nVYW = rr*p - r*(r-1)/2;
  double *vywFactors;
  int *piv;
  allocate(vywFactors, nVYW*nVYW);
  allocate(piv, nVYW);

  int info;
  vyw_factorize(A, vywFactors, piv, p, r, &info);
  if (info == 0) {
    vyw_solve(A, vywFactors, S, Gbuf, 1, q+1, piv, p, r);
  }

  freem(piv);
  freem(vywFactors);
  if (ownG) freem(Gbuf);
  if (ownC) freem(Cbuf);
  return info == 0;
}

static void vyw_factorize(double A[], double LU[], int piv[], int p, int r, int *info)
{
  int rr = r*r;
  int N = rr*p - r*(r-1)/2;
  int N0 = r*(r+1)/2;
  int Nc = rr*(p-1);
  int Nr = rr;
  double *F, *F0r, *F0c, *F00, *F11, *Ap, *Fii, *Apmi, *F0ri, *F00jj, *FKK;
  double *FKJ, *Ai, *Aj, *FiimjKK, *FijmiKJ, *F0ciKK, *F0ciKJ, *F00KK1;
  int i, j, k, i1, j1, I, J;
  *info = 0;
  if (p == 0) return;
  F = LU;
  allocate(F0r, Nr*Nc);
  allocate(F0c, Nc*Nr);
  allocate(F00, Nr*Nr);
  F11 = F + N0*(1 + N);
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
    lacpy("All", r-j, Nc, F0r+J, Nr, F+N0*N+j1, N);
    lacpy("All", Nc, r-j, F0c+Nc*J, Nc, F+N0+N*j1, N);
    j1 += r-j;
  }
  getrf(N, N, F, N, piv, info);
  freem(F0r);
  freem(F0c);
  freem(F00);
}

static void vyw_solve(double A[], double LU[], double S[], double Y[], int nrhs, int nY,
                      int piv[], int p, int r)
{
  xAssert(nrhs == 1);
  int rr = r*r;
  double *g, *g1, *Yp;
  int mg = r*r*(p+1);
  int N = r*r*p - r*(r-1)/2;
  allocate(g, mg*nrhs);
  g1 = g + r*(r+1)/2;
  Yp = Y + rr*nrhs*p;
  for (int k = 0; k < nrhs; k++) {
    double *g0k = g + mg*k;
    copytranspose(r, r, Y+rr*k, r, g0k, r);
    if (nY > p && p > 0) {
      double *Ykp = Yp + rr*k;
      double *Ap = A + rr*(p-1);
      gemm("NoT", "T", r, r, r, 1.0, Ykp, r, Ap, r, 1.0, g0k, r);
    }
  }
  for (int k = 0; k < nrhs; k++) {
    double *gptr = g + mg*k;
    for (int i = 0; i < r; i++) {
      double *gii0k = g + (1+r)*i + mg*k;
      for (int j = 0; j < r-i; j++) *gptr++ = *gii0k++;
    }
    for (int i = 0; i < r*(r-1)/2; i++) *gptr++ = 0;
  }
  for (int j = 1; j < nY && j < p; j++) {
    double *gj = g1 + (j-1)*rr;
    lacpy("All", rr, nrhs, Y + rr*nrhs*j, rr, gj, mg);
  }
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
  double *Sp = S + rr*nrhs*p;
  if (nY > p) lacpy("All", rr, nrhs, Yp, rr, Sp, rr);
  for (int k = 0; k < nrhs; k++) {
    double *Skp = Sp + rr*k;
    for (int i = 1; i <= p; i++) {
      double *Skpmi = S + rr*k + rr*nrhs*(p-i);
      double *Ai = A + rr*(i-1);
      gemm("NoT", "NoT", r, r, r, 1.0, Ai, r, Skpmi, r, 1.0, Skp, r);
    }
  }
  freem(g);
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
