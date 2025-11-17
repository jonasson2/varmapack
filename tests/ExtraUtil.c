// Utilities used by example programs and tests for varmapack_sim
#include "VarmaUtilities.h"
#include "VarmaPackUtil.h"
#include "allocate.h"
#include "ExtraUtil.h"
#include "BlasGateway.h"
#include "xAssert.h"
#include "VYW.h"
#include "Tests.h"
#include "varmapack.h"
#include <math.h>
#include <stdio.h>

void write_matrix(const char *fname, int m, int n, const double *A) {
  FILE *fp = fopen(fname, "w");
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      fprintf(fp, "%.17g ", A[i + j*n]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

double mean(const double *x, int n) {
  // Mean of vector
  if (n <= 0) return 0.0;
  double s = 0.0; for (int i=0; i<n; i++) s += x[i];
  return s/n;
}

double var(const double *x, int n, double mu) {
  // Unbiased variance of vector
  if (n <= 1) return 0.0;
  double s = 0.0;
  for (int i=0; i<n; i++) { double d = x[i] - mu; s += d*d; }
  return s/(n - 1);
}

double relabsdiff(double a[], double b[], int n) {
  // max(relative diff, absolute diff) where diff is difference between 
  // vectors a and b
  int ia, ib, ic;
  double *c, rmx, r;
  if (n == 0) return 0.0;
  allocate(c, n);
  copy(n, a, 1, c, 1);
  axpy(n, -1.0, b, 1, c, 1);
  ia = iamax(n, a, 1);
  ib = iamax(n, b, 1);
  ic = iamax(n, c, 1);
  rmx = fmax(fabs(a[ia]), fabs(b[ib]));
  rmx = fmax(1.0, rmx);
  r = fabs(c[ic])/rmx;
  freem(c);
  return r;
}

int almostSame(double a, double b) {
  // true if a and b are almost equal
  return relabsdiff(&a, &b, 1) < 5.0e-14;
}

int almostZero(double a[], int n) {
  // true if a ≈ 0
  int ia = iamax(n, a, 1);
  return fabs(a[ia]) < 5.0e-14;
}

int almostEqual(double a[], double b[], int n) {
  // true if relative difference between vectors a and b < 5e-14
  return relabsdiff(a, b, n) < 5.0e-14;
}

int almostAllSame(double a[], int n) {
  // true if max difference among elements of a is < 5e-14
  double minel, maxel, amax;
  int i;
  if (n == 0) return 1;
  minel = maxel = a[0];
  for (i=1; i<n; i++) {
    minel = fmin(minel, a[i]);
    maxel = fmax(maxel, a[i]);
  }
  amax = fmax(fabs(minel),fabs(maxel));
  return (maxel - minel) <= 5e-14*fmax(1,amax);
}

void mean3(double X[], int m, int n, int k, int idim, double mu[]) {
  // Calculate the mean of the m by n by k array X along the dimension idim.
  // The mean is returned in mu which is n by k, m by k or m by n for idim =
  // 1, 2 and 3 respectively.
  int i;
  switch(idim) {
  case(1):
    for (i=0; i<k; i++) meanmat("N", m, n, X + i*m*n, m, mu + i*n);
    break;
  case(2):
    for (i=0; i<k; i++) meanmat("T", m, n, X + i*m*n, m, mu + i*m);
    break;
  case(3):
    for (i=0; i<n; i++) meanmat("T", m, k, X + i*m, m*n, mu + i*m);
  }
}

void cov(char *transp, int m, int n, double X[], double C[]) {
  // C := covariance between columns of the m by n matrix op(X). If transp begins with N,
  // op(X) = X, but if it begins with T, op(X) = X^T. C is n by n.
  double *mu, *Xm;
  int i, tmp;
  if (transp[0] == 'T') { tmp = n; n = m; m = tmp; }
  setzero(n*n, C);
  if (m <= 1) return;
  allocate(mu, n);
  allocate(Xm, m*n);
  setzero(n, mu);
  if (transp[0] == 'T') copytranspose(n, m, X, n, Xm, m);
  else copy(m*n, X, 1, Xm, 1);
  for (i=0; i<m; i++) axpy(n, 1.0, Xm + i, m, mu, 1);
  scal(n, 1.0/m, mu, 1); // mu[i] = mean of X(:,i)
  for (i=0; i<m; i++) axpy(n, -1.0, mu, 1, Xm + i, m); // subtract mu from rows
  syrk("Low", "T", n, m, 1.0/(m-1), Xm, m, 0.0, C, n);
  copylowertoupper(n, C, n);
  freem(Xm); freem(mu);
}

void autocov(char *transp, char *norm, int r, int n, double X[], int maxlag, double C[]) {
  // Compute Ck = sample lag-k autocovariance matrix of X for k = 0,…,maxlag.
  //
  // When transp begins with "N", X is an r×n matrix whose column t contains the t-th
  // observation of an r-dimensional time series. When transp begins with "T", X is n×r,
  // with the t-th observation in row t.
  //
  // C is an r×r×(maxlag + 1) array returning Ck in its k-th subarray, C(:,:,k).
  //
  // When norm begins with "B" the estimation is *biased* (ML), normalizing all Ck by 1/n.
  // When norm begins with "U" it is *unbiased*, normalizing Ck by 1/(n−k−1).
  //
  // Mathematically:
  //
  //   Ck(i,j) = Cov(xi_t, xj_{t−k})
  //   Ck      = Cov(x_t, x_{t−k}) = Z·Y' / N
  //
  // where Y has the leading (n−k) columns of X and Z has the trailing (n−k) columns,
  // and N is the normalization factor: N = n (biased) or N = n−k−1 (unbiased).
  double fctr, *Y, *Z, *Ck;
  bool biased = (*norm == 'B');
  xAssert(*transp == 'N' || *transp == 'T');
  xAssert(biased || *norm == 'U');
  if (biased) xAssert(maxlag <= n-1);
  else        xAssert(maxlag <= n-2);
  Y = X;
  for (int k=0; k <= maxlag; k++) {
    fctr = biased ? 1.0/n : 1.0/(n - k - 1);
    Ck = C + r*r*k;
    if (*transp == 'N') {
      Z = X + r*k;
      gemm("N", "T", r, r, n-k, fctr, Z, r, Y, r, 0.0, Ck, r);
    }
    else {
      Z = X + k;
      gemm("T", "N", r, r, n-k, fctr, Z, n, Y, n, 0.0, Ck, r);
    }
  }
}
