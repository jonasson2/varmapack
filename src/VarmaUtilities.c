// Utilities for VarmaLoglik
#include "BlasGateway.h"
#include "xAssert.h"
#include "allocate.h"
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "VarmaUtilities.h"

void copytranspose(int m, int n, double A[], int ldA, double B[], int ldB){
  // Set B to the transpose of the m×n matrix A. Leading dimensions are in ldA
  // and ldB. Matrices must not overlap.
  int j;
  double *A1j = A, *Bj1 = B;
  for (j=0; j<n; j++) {
    copy(m, A1j, 1, Bj1, ldB);
    A1j += ldA;
    Bj1++;
  }  
}

void setzero(int n, double *x) {   // set x[0]...x[n-1] to 0.0
  memset(x, 0, n*sizeof(double));
}

void subtractmat(int m, int n, double A[], int ldA, double B[], int ldB) {
  // Subtract m×n matrix A from m×n matrix B
  int i;
  double *Ai, *Bi;
  for (i=0; i<n; i++) {
    Ai = A + ldA*i;
    Bi = B + ldB*i;    
    axpy(m, -1.0, Ai, 1, Bi, 1);
  }
}

void subtracttranspose(int m, int n, double A[], int ldA, double B[], int ldB) {
  // Subtract transpose of n × m matrix A from m × n matrix B
  int i;
  for (i=0; i<n; i++) {
    axpy(m, -1.0, A + i, ldA, B + ldB*i, 1);
  }
}

void addmat
(char *uplo, int m, int n, double A[], int ldA, double B[], int ldB){
  // With uplo beginning with "A", add m×n matrix A to m×n matrix B.
  // If uplo begins with "U" or "L" only the upper or lower triangles are added.
  int i;
  double *Ai, *Bi;
  switch (uplo[0]) {
    case('A'): case('a'):
      if (ldA == m && ldB == m) {
        axpy(m*n, 1.0, A, 1, B, 1);
      }
      else {
        for (i=0; i<n; i++) {
          Ai = A + ldA*i;
          Bi = B + ldB*i;
          axpy(m, 1.0, Ai, 1, Bi, 1);
        }
      }
      break;
    case('L'): case('l'):
      for (i=0; i<n; i++) {
        Ai = A + i + ldA*i;
        Bi = B + i + ldB*i;
        axpy(m-i, 1.0, Ai, 1, Bi, 1);
      }
      break;
    case('U'): case('u'):      
      for (i=0; i<n; i++) {
        Ai = A + ldA*i;
        Bi = B + ldB*i;
        axpy(i+1, 1.0, Ai, 1, Bi, 1);
      };
      break;
    default: xAssert(0);
  }
}

void meanmat(char *transp, int m, int n, double A[], int ldA, double mu[]) {
  // Find means of columns (if *transp begins with 'N', mu is an n-vector) or
  // rows (if it begins with 'T', mu is an m-vector) of the m by n matrix A.
  double *Ai;
  int i, j;
  switch (transp[0]) {
    case('N'): case('n'):
      setzero(n, mu);
      for (i=0; i<n; i++) {
        Ai = A + i*ldA;
        for (j=0; j<m; j++) mu[i] += Ai[j];
        mu[i] /= m;
      }
      break;
    case('T'): case('t'):
      setzero(m, mu);
      for (i=0; i<n; i++) {
        Ai = A + i*ldA;
        for (j=0; j<m; j++) mu[j] += Ai[j];
      }
      scal(m, 1.0/n, mu, 1);
      break;
  }
}

void aplusat(int n, double A[], int ldA) {
  // Add A' to A and return lower triangle of result in lower triangle of A
  int i;
  for (i=1; i<n; i++) {
    axpy(n-i, 1.0, A + (i-1) + i*ldA, ldA, A + i + (i-1)*ldA, 1);
  }
  scal(n, 2.0, A, ldA+1); // Scale diagonal
}

void copylowertoupper(int n, double A[], int ldA) {
  // Set strictly upper triangle of n × n matrix A to strictly lower triangle of 
  // its transpose.
  int i;
  for (i=1; i<n; i++) {
    copy(n-i, A + i + (i-1)*ldA, 1, A + (i-1) + i*ldA, ldA);
  }
}

bool anynan(int n, double A[]) { // return true if A[i] != A[i] for any i
  int i;
  for (i=0; i<n; i++) if (A[i] != A[i]) return true;
  return false;
}

bool islowermat(int n, double A[]) { // true iff A is lower triangular
  int i, j;
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      if (A[i + j*n] != 0) return false;
  return true;
}

void makeposdef(int n, double A[], double del0) { // make sure A is positive
  // definite, begin by adding del*I and increase del if needed
  int i, info = 1, imx;
  double *Aii, *B, del;
  allocate(B, n*n);
  imx = iamax(n, A, n+1);
  del = max(1.0, fabs(A[imx*(n+1)]))*max(del0, 1.0e-15);
  while (true) {
    lacpy("Low", n, n, A, n, B, n);
    potrf("Low", n, B, n, &info);  // Cholesky factorize B
    if (info == 0) break;
    for (i=0, Aii=A; i<n; i++, Aii += (n+1)) *Aii += del;
    del *= 10.0;
  }
  freem(B);
}
