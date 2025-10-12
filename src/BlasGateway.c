// Gateway to Fortran Blas functions used by the VarmaLoglik / VarLoglik package
#include "BlasGateway.h"
#include "blasref.h"

// If STRPAIR is defined, character string length arguments come immediately
// after the string (instead of all together at the end of the argument list).
// This is the default with the ifort compiler with the /iface:cvf flag and the
// Nag f95 compiler with the -f77 flag. These flags are used (and needed) when
// linking against the NAG library on Microsoft Windows, to obtain the stdcall
// calling convention (the -mrtd flag of g95, that also effects stdcall calling,
// does not change the placement of the length arguments).
//
// If STRLEN is defined the Blas are assumed to require character string length
// arguments. Normally they are absent, and that applies to how the CLapack
// reference Blas are implemented (the ones made with f2c, with names such as
// daxpy_, not the ones having names such as cblas_daxpy, as implemented in the
// gsl), and also to OpenBLAS, Accelerate, and MKL.


void axpy(int n, double alpha, double x[], int incx, double y[], int incy) {
  daxpy_(&n, &alpha, x, &incx, y, &incy);
}

void copy(int n, double x[], int incx, double y[], int incy) {
  dcopy_(&n, x, &incx, y, &incy);
}

double dot(int n, double x[], int incx, double y[], int incy) {
  return ddot_(&n, x, &incx, y, &incy);
}

double nrm2(int n, double x[], int incx) {
  return dnrm2_(&n, x, &incx);
}

void gemm(char *transa, char *transb, int m, int n, int k, double alpha,
          double a[], int lda, double b[], int ldb, double beta, double c[], int ldc) {
#ifdef STRPAIR
  dgemm_(transa, 1, transb, 1, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#elif defined(STRLEN)
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#else
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}

void gemv(char *transa, int m, int n, double alpha, double a[], int lda,
          double x[], int incx, double beta, double y[], int incy) {
#ifdef STRPAIR
  dgemv_(transa, 1, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#elif defined(STRLEN)
  dgemv_(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
#else
  dgemv_(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
}

void ger(int m, int n, double alpha, double x[], int incx, double y[], int incy,
         double a[], int lda) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

void getrf(int m, int n, double a[], int lda, int ipiv[], int *info) {
  dgetrf_(&m, &n, a, &lda, ipiv, info);
}

void getrs(char *transa, int n, int nrhs, double a[], int lda, int ipiv[],
           double b[], int ldb, int *info) {
#ifdef STRPAIR
  dgetrs_(transa, 1, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
#elif defined(STRLEN)
  dgetrs_(transa, &n, &nrhs, a, &lda, ipiv, b, &ldb, info, 1);
#else
  dgetrs_(transa, &n, &nrhs, a, &lda, ipiv, b, &ldb, info);
#endif
}

int iamax(int n, double dx[], int incx) { // NOTE: Returns 0 for 1st elem. etc.
  return idamax_(&n, dx, &incx) - 1;
}

void lacpy(char *uplo, int m, int n, double a[], int lda, double b[], int ldb) {
  // Ideally this function should call dlacpy from lapack, but because the acml
  // library does not include dlacpy, a function based on the one from clapack is
  // included here instead.
  int i, j;
  if (*uplo == 'U' || *uplo == 'u') {
    for (j = 0; j < n; j++) {
      for (i = 0; i <= j && i < m; i++) b[i + j * ldb] = a[i + j * lda];
    }
  } else if (*uplo == 'L' || *uplo == 'l') {
    for (j = 0; j < n; j++) {
      for (i = j; i < m; i++) b[i + j * ldb] = a[i + j * lda];
    }
  } else {
    for (j = 0; j < n; j++) {
      for (i = 0; i < m; i++) b[i + j * ldb] = a[i + j * lda];
    }
  }
}

void potrf(char *uplo, int n, double a[], int lda, int *info) {
#ifdef STRPAIR
  dpotrf_(uplo, 1, &n, a, &lda, info);
#elif defined(STRLEN)
  dpotrf_(uplo, &n, a, &lda, info, 1);
#else
  dpotrf_(uplo, &n, a, &lda, info);
#endif
}

void scal(int m, double alpha, double *x, int incx) {
  dscal_(&m, &alpha, x, &incx);
}

void symm(char *side, char *uplo, int m, int n, double alpha, double a[], 
          int lda, double b[], int ldb, double beta, double c[], int ldc) {
#ifdef STRPAIR
  dsymm_(side, 1, uplo, 1, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#elif defined(STRLEN)
  dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#else
  dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}

void symv(char *uplo, int n, double alpha, double a[], int lda, double x[],
          int incx, double beta, double y[], int incy) {
#ifdef STRPAIR
  dsymv_(uplo, 1, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#elif defined(STRLEN)
  dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
#else
  dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
#endif
}

void syr(char *uplo, int n, double alpha, double x[], int incx, double a[], int lda) {
#ifdef STRPAIR
  dsyr_(uplo, 1, &n, &alpha, x, &incx, a, &lda);
#elif defined(STRLEN)
  dsyr_(uplo, &n, &alpha, x, &incx, a, &lda, 1);
#else
  dsyr_(uplo, &n, &alpha, x, &incx, a, &lda);
#endif
}

void syr2k(char *uplo, char *trans, int n, int k, double alpha, double a[],
           int lda, double b[], int ldb, double beta, double c[], int ldc) {
#ifdef STRPAIR
  dsyr2k_(uplo, 1, trans, 1, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#elif defined(STRLEN)
  dsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
#else
  dsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
#endif
}

void syrk(char *uplo, char *trans, int m, int n, double alpha, double a[], 
          int lda, double beta, double c[], int ldc) {
#ifdef STRPAIR
  dsyrk_(uplo, 1, trans, 1, &m, &n, &alpha, a, &lda, &beta, c, &ldc);
#elif defined(STRLEN)
  dsyrk_(uplo, trans, &m, &n, &alpha, a, &lda, &beta, c, &ldc, 1, 1);
#else
  dsyrk_(uplo, trans, &m, &n, &alpha, a, &lda, &beta, c, &ldc);
#endif
}

void trmm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
#ifdef STRPAIR
  dtrmm_(side, 1, uplo, 1, transa, 1, diag, 1, &m, &n, &alpha, a, &lda, b, &ldb);
#elif defined(STRLEN)
  dtrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
#else
  dtrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}    

void trmv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
#ifdef STRPAIR
  dtrmv_(uplo, 1, transa, 1, diag, 1, &n, a, &lda, x, &incx);
#elif defined(STRLEN)
  dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
#else
  dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx);
#endif
}

void trsm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
#ifdef STRPAIR
  dtrsm_(side, 1, uplo, 1, transa, 1, diag, 1, &m, &n, &alpha, a, &lda, b, &ldb);
#elif defined(STRLEN)
  dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
#else
  dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
#endif
}

void trsv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
#ifdef STRPAIR
  dtrsv_(uplo, 1, transa, 1, diag, 1, &n, a, &lda, x, &incx);
#elif defined(STRLEN)
  dtrsv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
#else
  dtrsv_(uplo, transa, diag, &n, a, &lda, x, &incx);
#endif
}
