#include "ExtraUtil.h"
#include "Lyapunov.h"
#include "VYW.h"
#include "error.h"
#include "randompack.h"
#include "varmapack.h"
#include "xCheck.h"

static void check_case(char name[]) {
  int p = 0;
  int q = 0;
  int r = 0;
  int icase = 0;
  double *A = 0;
  double *B = 0;
  double *Sig = 0;
  double *Sv = 0;
  double *Cv = 0;
  double *Gv = 0;
  double *Sl = 0;
  double *Cl = 0;
  double *Gl = 0;
  randompack_rng *rng = randompack_create(0);
  varmapack_error error = varmapack_testcase(0, 0, 0, name, &p, &q, &r,
                                               &icase, 0, rng);
  int rr = r*r;
  xCheck(!error);
  xCheck(ALLOC(A, rr*(p > 0 ? p : 1)));
  xCheck(ALLOC(B, rr*(q > 0 ? q : 1)));
  xCheck(ALLOC(Sig, rr));
  xCheck(ALLOC(Sv, rr*(p+1)));
  xCheck(ALLOC(Cv, rr*(q+1)));
  xCheck(ALLOC(Gv, rr*(q+1)));
  xCheck(ALLOC(Sl, rr*(p+1)));
  xCheck(ALLOC(Cl, rr*(q+1)));
  xCheck(ALLOC(Gl, rr*(q+1)));
  error = varmapack_testcase(A, B, Sig, name, &p, &q, &r, &icase, 0, rng);
  xCheck(!error);
  xCheck(VYWFactorizeSolve(A, B, Sig, p, q, r, Sv, Cv, Gv));
  xCheck(LyapunovFactorizeSolve(A, B, Sig, p, q, r, Sl, Cl, Gl));
  checkArrayTol(Sl, Sv, rr*(p+1), 1e-9);
  checkArrayTol(Cl, Cv, rr*(q+1), 1e-9);
  checkArrayTol(Gl, Gv, rr*(q+1), 1e-9);
  randompack_free(rng);
  FREE(Gl);
  FREE(Cl);
  FREE(Sl);
  FREE(Gv);
  FREE(Cv);
  FREE(Sv);
  FREE(Sig);
  FREE(B);
  FREE(A);
}

int main(void) {
  xCheckInit("Lyapunov");
  check_case("tinyARMA");
  check_case("smallARMA2");
  check_case("mediumARMA1");
  if (xCheckNFailures() > 0) return 1;
  return 0;
}
