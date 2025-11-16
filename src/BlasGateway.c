// Gateway to Fortran Blas functions used by the VarmaLoglik / VarLoglik package
#include "BlasGateway.h"
#include "blasref.h"
#include "printX.h"

// Note that lapack_dpstrf_ is used instead of dpstrf_ because the latter routine is
// faulty in Accelerate; lapack_dpstrf.f with Netlib's official code must be compiled and
// linked against.

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
  dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

void gemv(char *transa, int m, int n, double alpha, double a[], int lda,
          double x[], int incx, double beta, double y[], int incy) {
  dgemv_(transa, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
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
  dgetrs_(transa, &n, &nrhs, a, &lda, ipiv, b, &ldb, info, 1);
}

void geev(char *jobvl, char *jobvr, int n, double a[], int lda, double wr[], double wi[],
          double vl[], int ldvl, double vr[], int ldvr, double work[], int lwork, int
          *info) {
  dgeev_(jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, info, 1, 1);
}

int iamax(int n, double dx[], int incx) { // NOTE: Returns 0 for 1st elem. etc.
  return idamax_(&n, dx, &incx) - 1;
}

void lacpy(char *uplo, int m, int n, double a[], int lda, double b[], int ldb) {
  dlacpy_(uplo, &m, &n, a, &lda, b, &ldb, 1);
}

void laset(char *uplo, int m, int n, double alpha, double beta, double a[], int lda) {
  dlaset_(uplo, &m, &n, &alpha, &beta, a, &lda, 1);
}

void potrf(char *uplo, int n, double a[], int lda, int *info) {
  dpotrf_(uplo, &n, a, &lda, info, 1);
}

void posv(char *uplo, int n, int nrhs, double a[], int lda, double b[], int ldb,
	  int *info) {
  dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, info, 1);
}

void pstrf(char *uplo, int n, double a[], int lda, int piv[], int *rank, double tol,
	   double work[], int *info) {
  printMsg("calling lapack_dpstrf");
  lapack_dpstrf_(uplo, &n, a, &lda, piv, rank, &tol, work, info, 1);
  for (int i=0; i<n; i++) piv[i]--;
}

void scal(int m, double alpha, double *x, int incx) {
  dscal_(&m, &alpha, x, &incx);
}

void symm(char *side, char *uplo, int m, int n, double alpha, double a[], 
          int lda, double b[], int ldb, double beta, double c[], int ldc) {
  dsymm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

void symv(char *uplo, int n, double alpha, double a[], int lda, double x[],
          int incx, double beta, double y[], int incy) {
  dsymv_(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy, 1);
}

void syr(char *uplo, int n, double alpha, double x[], int incx, double a[], int lda) {
  dsyr_(uplo, &n, &alpha, x, &incx, a, &lda, 1);
}

void syr2k(char *uplo, char *trans, int n, int k, double alpha, double a[],
           int lda, double b[], int ldb, double beta, double c[], int ldc) {
  dsyr2k_(uplo, trans, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc, 1, 1);
}

void syrk(char *uplo, char *trans, int m, int n, double alpha, double a[], 
          int lda, double beta, double c[], int ldc) {
  dsyrk_(uplo, trans, &m, &n, &alpha, a, &lda, &beta, c, &ldc, 1, 1);
}

void trmm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
  dtrmm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
}    

void trmv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
  dtrmv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
}

void trsm(char *side, char *uplo, char *transa, char *diag, int m, int n,
          double alpha, double a[], int lda, double b[], int ldb) {
  dtrsm_(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb, 1, 1, 1, 1);
}

void trsv(char *uplo, char *transa, char *diag, int n, double a[], int lda, 
          double x[], int incx) {
  dtrsv_(uplo, transa, diag, &n, a, &lda, x, &incx, 1, 1, 1);
}
