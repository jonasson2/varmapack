#include "Tests.h"
#include "TestUtil.h"
#include "varmapack.h"
#include "VarmaPackUtil.h"

void TestFindCFindG(void) {
  // Just a rudimentary test. Tested thoroughly by Testvarmapack_sim.
  test_model model;
  double *C, *G, *W;
  loadIndexedTestModel(&model, 8);
  xCheck(model.p == 1);
  xCheck(model.q == 1);
  xCheck(model.r == 2);
  int n = model.r*model.r*(model.q+1);
  xCheck(ALLOC(C, n));
  xCheck(ALLOC(G, n));
  xCheck(ALLOC(W, n));
  FindC(model.A, model.B, model.Sig, model.p, model.q, model.r, C);
  FindG(model.B, C, model.q, model.r, G);
  // Expected values from MATLAB for smallARMA1.
  // [ 2.0  1.0  1.4  1.3
  //   1.0  2.0  1.7  1.6 ]
  double Cexp[] = {2, 1, 1, 2, 1.4, 1.7, 1.3, 1.6};
  // [ 2.67  1.82  0.70  0.80
  //   1.67  2.82  0.70  0.80 ]
  double Gexp[] = {2.67, 1.67, 1.82, 2.82, 0.70, 0.70, 0.80, 0.80};
  xCheck(almostEqual(C, Cexp, n));
  xCheck(almostEqual(G, Gexp, n));
  xCheck(FindW(model.B, model.Sig, model.q, model.r, W));
  // Lower-lag convention: Wi = cov(y(t), y(t-i)).
  // [ 2.38  1.38  0.70  0.80
  //   1.38  2.38  0.70  0.80 ]
  double Wexp[] = {2.38, 1.38, 1.38, 2.38, 0.70, 0.70, 0.80, 0.80};
  xCheck(almostEqual(W, Wexp, n));
  FREE(W);
  FREE(G);
  FREE(C);
  freeTestModel(&model);
}
