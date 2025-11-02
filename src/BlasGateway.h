#ifndef BLASGATEWAY_H
#define BLASGATEWAY_H
#include "BlasUpper.h"

void axpy(int n, double alpha, double x[], int incx, double y[], int incy);

void copy(int n, double x[], int incx, double y[], int incy);

double dot(int n, double x[], int incx, double y[], int incy);

double dnrm2(int n, double x[], int incx);

void geev(char *jobvl, char *jobvr, int n, double a[], int lda, double wr[], double wi[],
          double vl[], int ldvl, double vr[], int ldvr, double work[], int lwork, int *info);

void gemm(char *transa, char *transb, int m, int n, int k, double alpha,
          double a[], int lda, double b[], int ldb, double beta, double c[],
          int ldc);

void gemv(char *transa, int m, int n, double alpha, double a[], int lda,
          double x[], int incx, double beta, double y[], int incy);

void ger(int m, int n, double alpha, double x[], int incx, double y[], int incy,
         double a[], int lda);

void getrf(int m, int n, double a[], int lda, int ipiv[], int *info);

void getrs(char *transa, int n, int nrhs, double a[], int lda, int ipiv[],
           double b[], int ldb, int *info);

int iamax(int n, double dx[], int incx); // -->idamax, returns 0 for 1st elem.

void lacpy(char *uplo, int m, int n, double a[], int lda, double b[], int ldb);

void potrf(char *uplo, int n, double a[], int lda, int *info);

void scal(int m, double alpha, double *x, int incx);

void symm (char *side, char *uplo, int m, int n, double alpha, double a[], 
           int lda, double b[], int ldb, double beta, double c[], int ldc);

void symv(char *uplo, int n, double alpha, double a[], int lda, double x[],
          int incx, double beta, double y[], int incy);

void syr(char *uplo, int n, double alpha, double x[], int incx, double a[],
         int lda);

void syr2k(char *uplo, char *trans, int m, int n, double alpha, double a[], 
           int lda, double b[], int ldb, double beta, double c[], int ldc);

void syrk(char *uplo, char *trans, int m, int n, double alpha, double a[], 
          int lda, double beta, double c[], int ldc);

void trmm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb);

void trmv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx);

void trsm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb);

void trsv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx);
#endif
