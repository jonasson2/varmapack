#include <math.h>
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_autocov_example(void) {
  int r = 2, n = 5, maxlag = 2;
  double tol = 1e-12;
  double X[] = {
    1, 0,
    2, 2,
    3, 1,
    2, 2,
    1, 3
  };
  double XT[] = {
    1, 2, 3, 2, 1,
    0, 2, 1, 2, 3
  };
  double C[2*2*3];
  double ml[] = {
    0.5600, -0.0800, -0.0800,  1.0400,
    0.0320, -0.0560,  0.0640, -0.1120,
   -0.3760, -0.2720,  0.4480,  0.0560
  };
  double corrected[] = {
    0.5600, -0.0800, -0.0800,  1.0400,
    0.0400, -0.0700,  0.0800, -0.1400,
   -0.6266666666666667, -0.4533333333333333,
    0.7466666666666666,  0.09333333333333334
  };
  varmapack_autocov("N", "ML", r, n, X, maxlag, C);
  for (int i=0; i<r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - ml[i]) < tol);
  }
  varmapack_autocov("N", "C", r, n, X, maxlag, C);
  for (int i=0; i<r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - corrected[i]) < tol);
  }
  varmapack_autocov("T", "ML", r, n, XT, maxlag, C);
  for (int i=0; i<r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - ml[i]) < tol);
  }
  varmapack_autocov("T", "C", r, n, XT, maxlag, C);
  for (int i=0; i<r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - corrected[i]) < tol);
  }
}

void TestAutocov(void) {
  check_autocov_example();
}
