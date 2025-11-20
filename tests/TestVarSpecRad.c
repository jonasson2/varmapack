// Testvarmapack_specrad.c
#include "ExtraUtil.h"
#include "allocate.h"
#include "error.h"
#include "xCheck.h"
#include "Tests.h"

#include "varmapack.h"   // varmapack_testcase(...)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// Prototype under test
double varmapack_specrad(double A[], int r, int p);

// --- helpers ---------------------------------------------------------------
static int nearly_one(double x, double tol) { return fabs(x - 1.0) <= tol; }

// ---------------------------------------------------------------------------
//  check2x2 — boundary case (rho == 1)
// ---------------------------------------------------------------------------
static void check2x2(void) {
  const int r = 2, p = 1;

  // Column-major A(:,:,1):
  double A[] = {
    0.95, 0.40,  // col 1
    0.10, 0.20   // col 2
  };

  double rho = varmapack_specrad(A, r, p);
  xCheck(nearly_one(rho, 1e-14));
  xCheck(rho <= 1.0 + 1e-14);
}

// ---------------------------------------------------------------------------
//  checkScaling — rho scales linearly with scalar multiplication
// ---------------------------------------------------------------------------
static void checkScaling(void) {
  const int r = 2, p = 1;
  double A[] = {
    0.95 * 20, 0.40 * 20,
    0.10 * 20, 0.20 * 20
  };

  double rho = varmapack_specrad(A, r, p);
  xCheck(almostSame(rho, 20.0));
  xCheck(rho > 19.999 && rho < 20.001);
}

// ---------------------------------------------------------------------------
//  checkZero — zero AR => rho = 0
// ---------------------------------------------------------------------------
static void checkZero(void) {
  const int r = 3, p = 1;
  double A[9] = {0.0};
  double rho = varmapack_specrad(A, r, p);
  xCheck(almostSame(rho, 0.0));
}

// ---------------------------------------------------------------------------
//  checkAllNamedStationary
//  For each icase=1..12:
//   1) query dimensions with A=B=Sig=0 (OK by API)
//   2) allocate A (r*r*p), B (r*r*q), Sig (r*r) — ALL nonzero
//   3) fetch the testcase
//   4) assert varmapack_specrad(A,r,p) < 1 (strict stationarity of AR part)
// ---------------------------------------------------------------------------
static void checkAllNamedStationary(void) {
  const double strict_tol = 1e-12;

  for (int icase = 1; icase <= 12; ++icase) {
    int p = 0, q = 0, r = 0;
    char name[32] = {0};

    // Step 1: query dims (A=B=Sig=0 is required when any is 0)
    bool ok = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stderr);
    xCheck(ok);

    // Some named cases might be ARMA(p,q) with p >= 0. We can handle p==0.
    xCheck(r > 0);
    xCheck(p >= 0);
    xCheck(q >= 0);

    // Step 2: allocate all three since we’ll pass nonzero pointers
    size_t nA = (size_t)r * (size_t)r * (size_t)(p > 0 ? p : 1);
    size_t nB = (size_t)r * (size_t)r * (size_t)(q > 0 ? q : 1);
    size_t nS = (size_t)r * (size_t)r;

    double *A = 0, *B = 0, *Sig = 0;
    allocate(A, nA);
    allocate(B, nB);
    allocate(Sig, nS);

    // Step 3: fetch the data (A,B,Sig all nonzero per API)
    ok = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stderr);
    xCheck(ok);

    // Step 4: spectral radius of AR companion (treat p==0 as rho=0)
    double rho = (p == 0) ? 0.0 : varmapack_specrad(A, r, p);

    if (!(rho < 1.0 - strict_tol)) {
      fprintf(stderr,
              "[Testvarmapack_specrad] Nonstationary/borderline: case %d (%s), "
              "rho=%.15g (p=%d, q=%d, r=%d)\n",
              icase, name, rho, p, q, r);
    }
    xCheck(rho < 1.0 - strict_tol);

    freem(A);
    freem(B);
    freem(Sig);
  }
}

// ---------------------------------------------------------------------------
//  Entry point
// ---------------------------------------------------------------------------
void Testvarmapack_specrad(void) {
  check2x2();                // boundary (rho == 1)
  checkScaling();            // scaling sanity
  checkZero();               // degenerate
  checkAllNamedStationary(); // the 12 named varmapack_testcase models
}
