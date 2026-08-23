#include "TestUtil.h"
#include "Lyapunov.h"
#include "VYW.h"
#include "varmapack.h"

static void check_case(const char *name) {
  test_model model;
  double *Sv = 0;
  double *Cv = 0;
  double *Gv = 0;
  double *Sl = 0;
  double *Cl = 0;
  double *Gl = 0;
  loadNamedTestModel(&model, name);
  int rr = model.r*model.r;
  xCheck(ALLOC(Sv, rr*(model.p+1)));
  xCheck(ALLOC(Cv, rr*(model.q+1)));
  xCheck(ALLOC(Gv, rr*(model.q+1)));
  xCheck(ALLOC(Sl, rr*(model.p+1)));
  xCheck(ALLOC(Cl, rr*(model.q+1)));
  xCheck(ALLOC(Gl, rr*(model.q+1)));
  xCheck(VYWFactorizeSolve(model.A, model.B, model.Sig, model.p, model.q,
                           model.r, Sv, Cv, Gv));
  xCheck(LyapunovFactorizeSolve(model.A, model.B, model.Sig, model.p,
                                model.q, model.r, Sl, Cl, Gl));
  checkArrayTol(Sl, Sv, rr*(model.p+1), 1e-9);
  checkArrayTol(Cl, Cv, rr*(model.q+1), 1e-9);
  checkArrayTol(Gl, Gv, rr*(model.q+1), 1e-9);
  FREE(Gl);
  FREE(Cl);
  FREE(Sl);
  FREE(Gv);
  FREE(Cv);
  FREE(Sv);
  freeTestModel(&model);
}

void TestLyapunov(void) {
  check_case("tinyARMA");
  check_case("smallARMA2");
  check_case("mediumARMA1");
}
