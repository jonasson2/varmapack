#include "allocate.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "VYW.h"
static void kronecker(int n, double alpha, double A[], double B[], double C[]);
static void kronI(int n, double alpha, double A[], double C[]);

void VYWFactorize ( // Set up and factorize matrix of Yule-Walker equations
  double A[],   // in   r × r × p, autoregressive parameter matrices
  double LU[],  // out  LU factors of the N×N matrix F where N = r*r*p-r*(r-1)/2
  int piv[],    // out  N-vector with pivoting information
  int p,        // in   number of autoregressive terms
  int r,        // in   dimension of each x(t)
  int *info)    // out  nonzero if factorization fails due to singularity
{
  int rr = r*r;
  int N = rr*p-r*(r-1)/2; // leading dimension of F
  int N0 = r*(r+1)/2;     // rows in lower triangle part of F00 and F0r
  int Nc = rr*(p-1);      // leading dimension of F0c
  int Nr = rr;            // leading dimension of F0r and F00
  double *F, *F0r, *F0c, *F00, *F11, *Ap, *Fii, *Apmi, *F0ri, *F00jj, *FKK;
  double *FKJ, *Ai, *Aj, *FiimjKK, *FijmiKJ, *F0ciKK, *F0ciKJ, *F00KK1;
  int i, j, k, i1, j1, k1, I, J;
  *info = 0;
  if (p==0) return;
  F = LU; // F and its factorization share memory
  allocate(F0r, Nr*Nc);
  allocate(F0c, Nc*Nr); // set to zero
  allocate(F00, Nr*Nr);
  F11 = F + N0*(1 + N); // location of F11 in F matrix
  // SET F TO THE IDENTITY MATRIX AND INITIALIZE F0r
  setzero(N*N, F);
  Ap = A + rr*(p-1);
  Fii = F11;
  for (i=0; i<p-1; i++) {
    for (j=0; j<rr; j++) {
      *Fii = 1.0;
      Fii += (N+1);
    }
    Ai = A + rr*i;
    Apmi = A + rr*(p-i-2);
    F0ri = F0r + Nr*rr*i;
    kronI(r, -1.0, Ai, F0ri);
    kronecker(r, -1.0, Ap, Apmi, F0ri);
  }
  for (j=0; j<rr; j++) { // INITIALIZE F00
    F00jj = F00 + (1 + Nr)*j;
    *F00jj = 1.0;
  }
  kronecker(r, -1.0, Ap, Ap, F00);
  // MAIN LOOP, OVER SUBMATRICES (K,J), BLOCK-ROWS (i) AND BLOCK-COLUMNS (j)
  for (k=0; k<r; k++) {
    FKK = F11 + r*k + N*r*k;
    FKJ = F11 + r*k + N*k;
    for (i=0; i<p-1; i++) {
      for (j=0; j<i; j++) { // SUBTRACT AHAT MATRICES FROM F
        Aj = A + rr*j;
        FiimjKK = FKK + rr*i + N*rr*(i-j-1);
        subtractmat(r, r, Aj, r, FiimjKK, N);
      }
      for (j=i+1; j<p; j++) { // SUBTRACT ATILDA MATRICES FROM F
        Aj = A + rr*j;
        FijmiKJ = FKJ + rr*i + N*rr*(j-i-1);
        subtractmat(r, r, Aj, r, FijmiKJ, N*r);
      }
      Ai = A + rr*i;
      F0ciKK = F0c + rr*i + r*k + Nc*r*k;
      F0ciKJ = F0c + rr*i + r*k + Nc*k;
      subtractmat(r, r, Ai, r, F0ciKK, Nc);   // SUBTRACT AHATi FROM COLUMN 0
      subtractmat(r, r, Ai, r, F0ciKJ, Nc*r); // SUBTRACT ATILDAi FROM COLUMN 0
    }
    // SUBTRACT ATILDAp·AHATp FROM F00 & ADD TO DIAGONAL OF F00
    for (k1=0; k1<r; k1++) {
      F00KK1 = F00 + r*k + Nr*r*k1;
      ger(r, r, -1.0, Ap + r*k1, 1, Ap + k, r, F00KK1, Nr);
    }
    *(F00 + (1 + Nr)*(r*k + k)) += 1;
  }
  // COPY ROWS CORRESPONDING TO LOWER TRIANGLE FROM F0r, F0c AND F00 TO F
  j1 = 0;
  for (j=0; j<r; j++) {
    J = r*j + j;
    i1 = 0;
    for (i=0; i<r; i++) {
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

static void kronecker(int n, double alpha, double A[], double B[], double C[]) { 
  // Add alpha times kronecker product of n×n matrices A and B to n·n×n·n mat. C 
  int i, j, k, nn = n*n;
  double *Aij, *Cij, *Cijk, *Bk;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Aij = A + i + n*j;
      Cij = C + n*i + nn*n*j;
      for (k=0; k<n; k++) {
        Bk = B + n*k;
        Cijk = Cij + nn*k;
        axpy(n, alpha*(*Aij), Bk, 1, Cijk, 1);
      }
    }
  }
}

static void kronI(int n, double alpha, double A[], double C[]) { 
  // Add alpha times the kronecker product of the n×n matrix A and the n×n 
  // identity matrix to the n·n × n·n matrix C 
  int i, j, k, nn = n*n;
  double *Aij, *Cij, *Cijkk, aa;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      Aij = A + i + n*j;
      Cij = C + n*i + nn*n*j;
      aa = alpha*(*Aij);
      for (k=0; k<n; k++) {
        Cijkk = Cij + (1+nn)*k;
        *Cijkk += aa;
      }
    }
  }
}
