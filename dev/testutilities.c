// Utilities used by tests of VarmaC
#include "VarmaSim.h"
#include "VarmaUtilities.h"
#include "allocate.h"
#include "testutilities.h"
#include "ExampleUtil.h"
#include "BlasGateway.h"
#include "xAssert.h"
#include <math.h>
#include <stdio.h>

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
  rmx = max(fabs(a[ia]), fabs(b[ib]));
  rmx = max(1.0, rmx);
  r = fabs(c[ic])/rmx;
  freem(c);
  return r;
}

int almostSame(double a, double b) {
  // true if a and b are almost equal
  return relabsdiff(&a, &b, 1) < 5.0e-14;
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
    minel = min(minel, a[i]);
    maxel = max(maxel, a[i]);
  }
  amax = max(fabs(minel),fabs(maxel));
  return (maxel - minel) <= 5e-14*max(1,amax);
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

double mean(int n, double x[]) { // Return mean of x
  double sum = 0.0;
  int i;
  if (n == 0) return 0.0;
  for (i=0; i<n; i++) sum += x[i];
  return sum/n;
}

void cov(char *transp, int m, int n, double X[], double C[]) {
  // C := covariance between columns of op(X). X is an m by n matrix.
  // If transp begins with N, op(X) = X, and C is n by n, but if
  // it begins with T, op(X) = X^T and C is m by m.
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

double difftest(fun f, int n, double x[], int ipar[], double par[]) {
  // Return max difference between numerical and analytic gradient
  double *xm, *xp, eps, fm, fp, dmax, *g, *gnum;
  int i;
  eps = sqrt(dot(n, x, 1, x, 1))*2.0e-6;
  allocate(xm, n);
  allocate(xp, n);
  allocate(g, n);
  allocate(gnum, n);
  copy(n, x, 1, xm, 1);
  copy(n, x, 1, xp, 1);
  f(n, x, g, ipar, par);
  dmax = 0;
  for (i=0; i<n; i++) {
    xp[i] += eps;
    xm[i] -= eps;
    fm = f(n, xm, 0, ipar, par);
    fp = f(n, xp, 0, ipar, par);
    gnum[i] = (fp - fm)/(2*eps);
    xm[i] = xp[i] = x[i];
  }
  dmax = relabsdiff(g, gnum, n);
  freem(gnum); freem(g); freem(xp); freem(xm);
  return dmax;
}

void lehmer(int n, double A[]) { // Return Lehmer matrix of order n
  int i, j;
  for (i=0; i<n; i++) {
    for (j=i; j<n; j++) A[j + i*n] = (double) (i+1) / (double) (j+1);
  }
  copylowertoupper(n, A, n);
}

void minij(int n, double A[]) { // Return matrix with A(i,j) = min(i,j)
  int i, j;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) A[i + j*n] = min(i+1, j+1);
  }
}

static void checkmeancov(void) {
  double X[] = {1.0, 3.0, 5.0,
                1.0, 4.0, 4.0};
  double Y[] = {1.0, 1.0,
                3.0, 4.0,
                5.0, 4.0};
  double C[4], D[4];
  cov("N", 3, 2, X, C);
  xAssert(almostSame(mean(3, X), 3.0));
  xAssert(almostSame(C[0], 4.0));
  xAssert(almostSame(C[1], 3.0));
  xAssert(almostSame(C[2], 3.0));
  xAssert(almostSame(C[3], 3.0));
  cov("T", 2, 3, Y, D);
  xAssert(almostEqual(C, D, 4));
}

static void checkmean3(void) {
  // Check correctness of mean3 with a simple 2 by 3 by 4 array.
  double X[2*3*4], mu[3*4];
  int i;
  for (i=0; i<2*3*4; i++) X[i] = i+1;
  mean3(X, 2, 3, 4, 1, mu);
  for (i=0; i<12; i++) xAssert(fabs(mu[i] - (1.5 + 2*i)) < 1e-14);
  mean3(X, 2, 3, 4, 2, mu);
  for (i=0; i<8; i+=2)
    xAssert(fabs(mu[i] - (3 + 3*i)) + fabs(mu[i+1] - (4 + 3*i)) < 1e-14);
  mean3(X, 2, 3, 4, 3, mu);
  for (i=0; i<6; i++) xAssert(fabs(mu[i] - (10 + i)) < 1e-14);
}

static double dtf(int n, double x[], double g[], int ipar[], double par[]) {
  // function for checkdifftest
  double a = par[0], fx;
  int i = ipar[0];
  (void) n; // suppress unused parameter warning
  fx = sin(x[0] + a)*x[1]*i;
  if (g) {
    g[0] = cos(x[0] + a)*x[1]*i;
    g[1] = sin(x[0] + a)*i;
  }
  return fx;
}

static void checkdifftest(void) {
  double par[1] = {2.5}, x[2] = {1.2, 1.7}, dmax;
  int ipar[1] = {3};
  dmax = difftest(dtf, 2, x, ipar, par);
  xAssert(dmax < 1.0e-8);
}

static void checkgallery(void) {
  // Quickly check lehmer and minij which are based on Matlab's gallery function
  double Lehmer[] = {1.0, 0.5, 1.0/3.0, 0.5, 1.0, 2.0/3.0, 1.0/3.0, 2.0/3.0, 1.0};
  double Minij[] = {1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 1.0, 2.0, 3.0}, A[9];  
  lehmer(3, A);
  xAssert(almostEqual(A, Lehmer, 9));
  minij(3, A);
  xAssert(almostEqual(A, Minij, 9));
}

void checktestutilities(void) {
  double a[] = {1.0, 2.0, 3.0}, b[] = {2.0, 2.0, 2.0 + 1e-15};
  xAssert(almostEqual(a, a, 3));
  xAssert(!almostEqual(a, b, 3));
  xAssert(almostAllSame(b, 3));
  xAssert(!almostAllSame(a, 3));
  xAssert(relabsdiff(a, a, 3) < 1e-15);
  xAssert(relabsdiff(a, b, 3) > 0.1);
  xAssert(almostSame(b[1], b[2]));
  xAssert(!almostSame(a[1], a[2]));
  checkdifftest();
  checkgallery();
  checkmean3();
  checkmeancov();
}
