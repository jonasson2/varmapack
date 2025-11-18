#include "ExtraUtil.h"
#include "allocate.h"
#include "xAssert.h"
#include "xCheck.h"
#include <math.h>
#include <stdio.h>
#include "Tests.h"

static void checkmeancov(void) {
  double X[] = {1.0, 3.0, 5.0,
                1.0, 4.0, 4.0};
  double Y[] = {1.0, 1.0,
                3.0, 4.0,
                5.0, 4.0};
  double C[4], D[4];
  cov("N", 3, 2, X, C);
  xCheck(almostSame(mean(X, 3), 3.0));
  xCheck(almostSame(C[0], 4.0));
  xCheck(almostSame(C[1], 3.0));
  xCheck(almostSame(C[2], 3.0));
  xCheck(almostSame(C[3], 3.0));
  cov("T", 2, 3, Y, D);
  xCheck(almostEqual(C, D, 4));
}

static void checkvar(void) {
  // Also checks mean
  double x[] = {1.0, 2.0, 3.0, 4.0};
  int n = 4;
  double mu = mean(x, n);
  xCheck(almostSame(mu, 2.5));
  double v = var(x, n, mu);
  xCheck(almostSame(v, 5.0/3.0));
  // all equal values (should be 0):
  double x2[] = {3.0, 3.0, 3.0, 3.0, 3.0};
  double mu2 = mean(x2, 5);
  xCheck(almostSame(mu2, 3.0));
  double v2 = var(x2, 5, mu2);
  xCheck(almostSame(v2, 0.0));
}

static void checkdiff(void) {
  double a[] = {1, 2, 3};
  double b[] = {1, 2 + 1e-14, 3 - 1e-14}; // small perturbations
  double c[] = {1, 2 + 1e-12, 3};         // larger perturbations
  int n = 3;
  xCheck(relabsdiff(a, b, n) < 1e-12);
  xCheck(almostEqual(a, a, n));
  xCheck(almostEqual(b, b, n));
  xCheck(almostEqual(a, b, n));
  xCheck(!almostEqual(a, c, n));
  double d[] = {5, 5 + 1e-14, 5 - 1e-14};
  xCheck(almostAllSame(d, 3));
}

void TestExtraUtil(void) {
  double a[] = {1.0, 2.0, 3.0}, b[] = {2.0, 2.0, 2.0 + 1e-15};
  xCheck(almostEqual(a, a, 3));
  xCheck(!almostEqual(a, b, 3));
  xCheck(almostAllSame(b, 3));
  xCheck(!almostAllSame(a, 3));
  xCheck(relabsdiff(a, a, 3) < 1e-15);
  xCheck(relabsdiff(a, b, 3) > 0.1);
  xCheck(almostSame(b[1], b[2]));
  xCheck(!almostSame(a[1], a[2]));
  checkdiff(); // more checks of the above
  checkmeancov();
  checkvar();
}
