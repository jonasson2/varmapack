// TestCompareWithMatlab.c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "VarmaUtilities.h"
#include "ExtraUtil.h"
#include "allocate.h"
#include "error.h"
#include "xCheck.h"
#include "printX.h"
#include "Tests.h"
#include "varmapack.h"
#include "DebugUtil.h"
#include "FromMatlab.h"

static int getdim(FILE *fid, char *base, int k) {
  char name[32];
  double scalar = 0;
  bool ok;
  STRSETF(name, "%s%d", base, k);
  ok = DoubleFromMatlab(fid, name, &scalar);
  xCheck(ok);
  return scalar;
}

bool TestAgainstMatlab(void) {
  // Bring in MATLAB reference data
  int n, p, q, r, pk, qk, rk, h, ncases, nMs;
  double *cases, *Ms, *A, *B, *Sig, *x0, condR, tol;
  bool ok;
  char name[64] = {0};
  char comparefile[] = "matlabcompare.txt";
  FILE *fid = fopen(comparefile, "r");
  ASSERT(fid != 0, "File %s cannot be read", comparefile);
  ok = IntFromMatlab(fid, "#cases", &ncases);
  ok &= IntFromMatlab(fid, "#Ms", &nMs);
  ok &= IntFromMatlab(fid, "n", &n);
  xCheck(ok && ncases < 1000);
  ALLOC(cases, ncases);
  ALLOC(Ms, nMs);
  ok = MatrixFromMatlab(fid, "cases", cases, (int)ncases, 1);
  xCheck(ok);

  // Loop over cases
  for (int k = 0; k < ncases; k++) {
    int icase = cases[k];
    xCheck(0 <= icase && icase < 100);
    
    // Get dimensions from Matlab file
    p = getdim(fid, "p", k);
    q = getdim(fid, "q", k);
    r = getdim(fid, "r", k);
    h = imax(p, q);

    // Get condition number, set tolerance
    ok = DoubleFromMatlab(fid, "condR", &condR);
    xCheck(ok);
    tol = 5e-16;

    // Get starting vector
    ALLOC(x0, r*h);
    STRSETF(name, "x0%d", k);
    ok = MatrixFromMatlab(fid, name, x0, r, h);
    xCheck(ok);

    // Set additional message to testcase number
    printI("TESTCASE NUMBER", icase);
    char addmsg[32];
    snprintf(addmsg, 32, "varmapack_testcase number %2d", icase);
    xCheckAddMsg(addmsg);

    // Query testcase dimensions and check them
    ok = varmapack_testcase(0, 0, 0, name, &pk, &qk, &rk, &icase, 0, stderr);
    xCheck(ok);
    xCheck(pk == p);
    xCheck(qk == q);
    xCheck(rk == r);
    xCheck(h <= n);

    // Allocate and fetch actual testcase coefficients
    int r2 = r*r;
    ALLOC(A, r2*p);
    ALLOC(B, r2*q);
    ALLOC(Sig, r2);
    STRSET(name, "");
    ok = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stderr);
    xCheck(ok);

    // Compare A, B, Sig against MATLAB reference
    double diffA, diffB, diffSig;
    
    STRSETF(name, "A%d", k);   CompareWithMatlab(fid, name, A, r, r*p, &diffA);
    STRSETF(name, "B%d", k);   CompareWithMatlab(fid, name, B, r, r*q, &diffB);
    STRSETF(name, "Sig%d", k); CompareWithMatlab(fid, name, Sig, r, r, &diffSig);
    xCheck(diffA < tol && diffB < tol && diffSig < tol);

    // Simulate with varmapack_sim and compare the simulated series with Matlab's
    for (int j=0; j<(int)nMs; j++) {
      int M = Ms[j];
      double *X, *E, *X0, *E0;
      bool sim_ok;
      h = imax(p, q);
    
      ALLOC(X, r*n);
      ALLOC(X0, r*n*2);
      ALLOC(E, r*n);
      ALLOC(E0, r*n*2);
      randompack_rng *rng = randompack_create("PM", 42);
      printM("A", A, r, r*p);
      printM("B", B, r, r*q);
      printM("Sig", Sig, r, r);
      double *mu = 0;
      varmapack_sim(A, B, Sig, mu, p, q, r, n, M, 0, 0, rng, X, E, &sim_ok);
      varmapack_sim(A, B, Sig, mu, p, q, r, n, M, x0, h, rng, X0, E0, &sim_ok);
      printM("X0", X0, r, n*M);
      xCheck(sim_ok == 1);
      // Compare E and X
      
      double diffE, diffX, diffE0, diffX0;
      STRSETF(name, "E%d", k);  ok = CompareWithMatlab(fid, name, E, r*n, M, &diffE);
      STRSETF(name, "X%d", k);  ok |= CompareWithMatlab(fid, name, X, r*n, M, &diffX);
      STRSETF(name, "E0%d", k); ok |= CompareWithMatlab(fid, name, E0, r*n, M, &diffE0);
      STRSETF(name, "X0%d", k); ok |= CompareWithMatlab(fid, name, X0, r*n, M, &diffX0);
      ASSERT(ok, "error comparing with Matlab");
      printD(">>> diffX", diffX);
      printD(">>> diffE", diffE);
      printD(">>> diffX0", diffX0);
      printD(">>> diffE0", diffE0);
      printD(">>> tolerance", tol*condR);
      xCheck(diffE < tol*condR);
      xCheck(diffX < tol*condR);
      xCheck(diffE0 < tol*condR);
      xCheck(diffX0 < tol*condR);

      // Free memory
      FREE(E0);
      FREE(X0);
      FREE(E);
      FREE(X);
      randompack_free(rng);
    }
    FREE(A);
    FREE(B);
    FREE(Sig);
    FREE(x0);
    xCheckAddMsg("");
  }
  FREE(cases);
  fclose(fid);
  return true;
}
