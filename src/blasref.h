// Include file declaring the reference blas functions
#ifndef BLASREF_H
#define BLASREF_H
#include "BlasUpper.h"

void daxpy_(int *n, double *alpha, double x[], int *incx, double y[], int *incy);

void dcopy_(int *n, double x[], int *incx, double y[], int *incy);

double ddot_(int *n, double x[], int *incx, double y[], int *incy);

#ifdef STRPAIR
void dgemm_(char *transa, int lentransa, char *transb, int lentransb, int *m, int *n,
            int *k, double *alpha, double a[], int *lda, double b[], int *ldb,
            double *beta, double c[], int *ldc);
#elif defined(STRLEN)
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
            double a[], int *lda, double b[], int *ldb, double *beta, double c[], 
            int *ldc, int lentransa, int lentransb);
#else
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha,
            double a[], int *lda, double b[], int *ldb, double *beta, double c[], 
            int *ldc);
#endif

#ifdef STRPAIR
void dgemv_(char *transa, int lentransa, int *m, int *n, double *alpha, double a[],
            int *lda, double x[], int *incx, double *beta, double y[], int *incy);
#elif defined(STRLEN)
void dgemv_(char *transa, int *m, int *n, double *alpha, double a[], int *lda,
            double x[], int *incx, double *beta, double y[], int *incy, int lentransa);
#else
void dgemv_(char *transa, int *m, int *n, double *alpha, double a[], int *lda,
            double x[], int *incx, double *beta, double y[], int *incy);
#endif

void dger_(int *m, int *n, double *alpha, double x[], int *incx, double y[], 
           int *incy, double a[], int *lda);

void dgetrf_(int *m, int *n, double a[], int *lda, int ipiv[], int *info);

#ifdef STRPAIR
void dgetrs_(char *transa, int lentransa, int *n, int *nrhs, double a[], int *lda,
             int ipiv[], double b[], int *ldb, int *info);
#elif defined(STRLEN)
void dgetrs_(char *transa, int *n, int *nrhs, double a[], int *lda, int ipiv[],
             double b[], int *ldb, int *info, int lentransa);
#else
void dgetrs_(char *transa, int *n, int *nrhs, double a[], int *lda, int ipiv[],
             double b[], int *ldb, int *info);
#endif

double dnrm2_(int *n, double x[], int *incx);

#ifdef STRPAIR
void dpotrf_(char *uplo, int lenuplo, int *n, double a[], int *lda, int *info);
#elif defined(STRLEN)
void dpotrf_(char *uplo, int *n, double a[], int *lda, int *info, int lenuplo);
#else
void dpotrf_(char *uplo, int *n, double a[], int *lda, int *info);
#endif

void dscal_(int *m, double *alpha, double *x, int *incx);

#ifdef STRPAIR
void dspr_(char *uplo, int lenuplo, int *n, double *alpha, double x[], int *incx,
           double ap[]);
#elif defined(STRLEN)
void dspr_(char *uplo, int *n, double *alpha, double x[], int *incx, double ap[],
           int lenuplo);
#else
void dspr_(char *uplo, int *n, double *alpha, double x[], int *incx, double ap[]);
#endif

#ifdef STRPAIR
void dspr2_(char *uplo, int lenuplo, int *n, double *alpha, double x[], int *incx,
            double y[], int *incy, double ap[]);
#elif defined(STRLEN)
void dspr2_(char *uplo, int *n, double *alpha, double x[], int *incx, double y[],
            int *incy, double ap[], int lenuplo);
#else
void dspr2_(char *uplo, int *n, double *alpha, double x[], int *incx, double y[],
            int *incy, double ap[]);
#endif

#ifdef STRPAIR
void dsymm_(char *side, int lenside, char *uplo, int lenuplo, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb, double *beta,
            double c[], int *ldc);
#elif defined(STRLEN)
void dsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double a[], 
            int *lda, double b[], int *ldb, double *beta, double c[], int *ldc,
            int lenside, int lenuplo);
#else
void dsymm_(char *side, char *uplo, int *m, int *n, double *alpha, double a[], 
            int *lda, double b[], int *ldb, double *beta, double c[], int *ldc);
#endif

#ifdef STRPAIR
void dsymv_(char *uplo, int lenuplo, int *n, double *alpha, double a[], int *lda,
            double x[], int *incx, double *beta, double y[], int *incy);
#elif defined(STRLEN)
void dsymv_(char *uplo, int *n, double *alpha, double a[], int *lda, double x[],
            int *incx, double *beta, double y[], int *incy, int lenuplo);
#else
void dsymv_(char *uplo, int *n, double *alpha, double a[], int *lda, double x[],
            int *incx, double *beta, double y[], int *incy);
#endif

#ifdef STRPAIR
void dsyr_(char *uplo, int lenuplo, int *n, double *alpha, double x[], int *incx,
          double a[], int *lda);
#elif defined(STRLEN)
void dsyr_(char *uplo, int *n, double *alpha, double x[], int *incx, double a[],
          int *lda, int lenuplo);
#else
void dsyr_(char *uplo, int *n, double *alpha, double x[], int *incx, double a[],
          int *lda);
#endif

#ifdef STRPAIR
void dsyr2k_(char *uplo, int lenuplo, char *trans, int lentrans, int *m, int *n,
             double *alpha, double a[], int *lda, double b[], int *ldb, double *beta,
             double c[], int *ldc);
#elif defined(STRLEN)
void dsyr2k_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], 
             int *lda, double b[], int *ldb, double *beta, double c[], int *ldc, 
             int lenuplo, int lentrans);
#else
void dsyr2k_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], 
             int *lda, double b[], int *ldb, double *beta, double c[], int *ldc);
#endif

#ifdef STRPAIR
void dsyrk_(char *uplo, int lenuplo, char *trans, int lentrans, int *m, int *n,
            double *alpha, double a[], int *lda, double *beta, double c[], int *ldc);
#elif defined(STRLEN)
void dsyrk_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], 
            int *lda, double *beta, double c[], int *ldc, int lenuplo, int lentrans);
#else
void dsyrk_(char *uplo, char *trans, int *m, int *n, double *alpha, double a[], 
            int *lda, double *beta, double c[], int *ldc);
#endif

#ifdef STRPAIR
void dtrmm_(char *side, int lenside, char *uplo, int lenuplo, char *transa,
            int lentransa, char *diag, int lendiag, int *m, int *n, double *alpha,
            double a[], int *lda, double b[], int *ldb);
#elif defined(STRLEN)
void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb, int lenside, 
            int lenuplo, int lentransa, int lendiag);
#else
void dtrmm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb);
#endif

#ifdef STRPAIR
void dtrmv_(char *uplo, int lenuplo, char *trans, int lentrans, char *diag,
            int lendiag, int *n, double A[], int *lda, double x[], int *incx);
#elif defined(STRLEN)
void dtrmv_(char *uplo, char *trans, char *diag, int *n, double A[], int *lda,
            double x[], int *incx, int lenuplo, int lentrans, int lendiag);
#else
void dtrmv_(char *uplo, char *trans, char *diag, int *n, double A[], int *lda,
            double x[], int *incx);
#endif

#ifdef STRPAIR
void dtrsm_(char *side, int lenside, char *uplo, int lenuplo, char *transa,
            int lentransa, char *diag, int lendiag, int *m, int *n, double *alpha,
            double a[], int *lda, double b[], int *ldb);
#elif defined(STRLEN)
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb, int lenside, 
            int lenuplo, int lentransa, int lendiag);
#else
void dtrsm_(char *side, char *uplo, char *transa, char *diag, int *m, int *n,
            double *alpha, double a[], int *lda, double b[], int *ldb);
#endif

#ifdef STRPAIR
void dtrsv_(char *uplo, int lenuplo, char *transa, int lentransa, char *diag,
            int lendiag, int *n, double a[], int *lda, double x[], int *incx);
#elif defined(STRLEN)
void dtrsv_(char *uplo, char *transa, char *diag, int *n, double a[], int *lda, 
            double x[], int *incx, int lenuplo, int lentransa, int lendiag);
#else
void dtrsv_(char *uplo, char *transa, char *diag, int *n, double a[], int *lda, 
            double x[], int *incx);
#endif

int idamax_(int *n, double dx[], int *incx);

#endif // BLASREF_H
