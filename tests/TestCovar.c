#include <math.h>
#include "Tests.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_autocov_example(void) {
#define R 2
#define N 5
  int r = R, n = N; 
  const int maxlag = 2;
  const double tol = 1e-12;
  // Columns are observations at times t = 1..n
  double X[R*N] = {
    1.0, 0.0,
    2.0, 2.0,
    3.0, 1.0,
    2.0, 2.0,
    1.0, 3.0
  };
  double C[ r*r*(maxlag+1) ];
  const double Cb[] = {
    0.5600, -0.0800, -0.0800,  1.0400,
    0.0320, -0.0560,  0.0640, -0.1120,
   -0.3760, -0.2720,  0.4480,  0.0560
  };
  const double Cu[] = {
    0.5600, -0.0800, -0.0800,  1.0400,
    0.0400, -0.0700,  0.0800, -0.1400,
   -0.6266666666666667, -0.4533333333333333,
    0.7466666666666666,  0.09333333333333334
  };

  varmapack_autocov("N", "ML", r, n, X, maxlag, C);
  for (int i = 0; i < r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - Cb[i]) < tol);
  }

  varmapack_autocov("N", "C", r, n, X, maxlag, C);
  for (int i = 0; i < r*r*(maxlag+1); i++) {
    xCheck(fabs(C[i] - Cu[i]) < tol);
  }
}

void TestCovar(void) {
  check_autocov_example();
}
