#ifndef VARMAUTILITIES_H
#define VARMAUTILITIES_H

// Following will on some platforms define min and max, so including it here
// prevents macro redefinition warning message
#include <stdlib.h>

// Utilities for VarmaLoglik

void copytranspose(int m, int n, double A[], int ldA, double B[], int ldB);
  // Set B to the transpose of the m×n matrix A. Leading dimensions are in ldA
  // and ldB.

void setzero(int n, double *x);
  // set x[0]...x[n-1] to 0.0

void subtractmat(int m, int n, double A[], int ldA, double B[], int ldB);
  // Subtract m×n matrix A from m×n matrix B

void addmat(char *uplo, int m, int n, double A[], int ldA, double B[], int ldB);
  // With uplo beginning with "A", add m×n matrix A to m×n matrix B.
  // If uplo begins with "U" or "L" only the upper or lower triangles are added.

void subtracttranspose(int m, int n, double A[], int ldA, double B[], int ldB);
  // Subtract transpose of n×m matrix A from m×n matrix B

void meanmat(char *transp, int m, int n, double A[], int ldA, double mu[]);
  // Find means of columns (if *transp begins with 'N', mu is an n-vector) or
  // rows (if it begins with 'T', mu is an m-vector) of the m by n matrix A.

void aplusat(int n, double A[], int ldA);
  // Add A' to A and return lower triangle of result in lower triangle of A

void copylowertoupper(int n, double A[], int ldA);
  // Set strictly upper triangle of n×n matrix A to strictly lower triangle of 
  // its transpose.

int anynan(int n, double A[]); // return true if A[i] != A[i] for any i

void makeposdef(int n, double A[], double del);
  // make sure A is positive definite

#ifndef min
#define min(x, y) ( (x) < (y) ? (x) : (y) )
#endif
#ifndef max
#define max(x, y) ( (x) > (y) ? (x) : (y) )
#endif

#endif
