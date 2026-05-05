#include <stdio.h>
#include "Tests.h"
#include "error.h"
#include "xCheck.h"
#include "varmapack.h"
#include "VarmaPackUtil.h"
#include "ExtraUtil.h"

void TestFindCG(void) {
  // Just a rudimentary test. Tested thoroughly by Testvarmapack_sim.
  int p, q, r, icase;
  char name[16] = "";
  double *A, *B, *Sig, *C, *G;
  bool ok;

  // First call: query dimensions for smallARMA1
  icase = 8;
  ok = varmapack_testcase(0, 0, 0, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  xCheck(p == 1);
  xCheck(q == 1);
  xCheck(r == 2);

  // Allocate and fetch testcase data
  xCheck(ALLOC(A, r*r*p));
  xCheck(ALLOC(B, r*r*q));
  xCheck(ALLOC(Sig, r*r));
  xCheck(ALLOC(C, r*r*(q + 1)));
  xCheck(ALLOC(G, r*r*(q + 1)));

  ok = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, stdout);
  xCheck(ok);
  vpack_FindCG(A, B, Sig, p, q, r, C, G);

  // Expected values from MATLAB for smallARMA1.
  // [ 2.0  1.0  1.4  1.3
  //   1.0  2.0  1.7  1.6 ]
  double Cexp[] = {2.0, 1.0, 1.0, 2.0, 1.4, 1.7, 1.3, 1.6};

  // [ 2.67  1.82  0.70  0.80
  //   1.67  2.82  0.70  0.80 ]
  double Gexp[] = {2.67, 1.67, 1.82, 2.82, 0.70, 0.70, 0.80, 0.80};

  xCheck(almostEqual(C, Cexp, r*r*(q + 1)));
  xCheck(almostEqual(G, Gexp, r*r*(q + 1)));

  FREE(A);
  FREE(B);
  FREE(Sig);
  FREE(C);
  FREE(G);
}
