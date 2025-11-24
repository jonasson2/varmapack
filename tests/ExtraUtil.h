// ExtraUtils.h — Utilities used by example/demo and test programs (means, (co)variances,
// and approximate-equality helpers)

#ifndef EXTRAUTIL_H
#define EXTRAUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>
double mean(const double *x, int n);               // Mean of vector
double var(const double *x, int n, double mu);     // Unbiased variance of vector
int almostSame(double a, double b);                // are a and b are almost equal?
int almostEqual(double a[], double b[], int n);    // is rel.diff. beween a and b < 5e-14?
int almostAllSame(double a[], int n);              // is max diff. among x-elements < 5e-14?
int almostZero(double a[], int n);                 // is a ≈ 0?
bool TestMatlabMatrix(char *filename, double *A, int m, int n);
double condnum(double *A, int n);                  // norm-2 condition number of sym. PD matrix

void cov(char *transp, int m, int n, double X[], double C[]);
// C := covariance between columns of op(X). X is an m by n matrix. If transp begins with
// N, op(X) = X, and C is n by n, but if it begins with T, op(X) = X^T and C is m by m.

void SCbuild(double A[], double B[], double Sig[], int p, int q, int r, int n,
             double CC[], double SS[]);
// CC = cov(x,eps) and SS = cov(x,x) (s = t-i). Used by testing functions. See also
// CCBuild.c

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* EXTRAUTIL_H */
