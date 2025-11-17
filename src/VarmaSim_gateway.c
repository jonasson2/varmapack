#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>
#include "varmapack.h"
#include "RandomNumbers.h"
#include "printX.h"

SEXP varmapack_sim_gateway(SEXP A, SEXP B, SEXP Sig, SEXP mu, SEXP p, SEXP q, SEXP r,
                      SEXP n, SEXP M, SEXP X0, SEXP X, SEXP eps, SEXP ok,
                      SEXP park_miller, SEXP seed) {
  // R-to-C gateway function to be called by .Call in varmapack_sim.R
  double
    *Ap = REAL(A),
    *Bp = REAL(B),
    *Sigp = REAL(Sig),
    *mup = Rf_isNull(mu) ? 0 : REAL(mu),
    *X0p = Rf_isNull(X0) ? 0 : REAL(X0),
    *Xp = REAL(X),
    *epsp = REAL(eps);
  
  int
    pv = asInteger(p),
    qv = asInteger(q),
    rv = asInteger(r),
    nv = asInteger(n),
    Mv = asInteger(M),
    s = asInteger(seed),
    *okp = INTEGER(ok);

  bool PM = LOGICAL(park_miller)[0];
  RandRng *rng = rand_create();
  if (PM) {
    rand_settype(PARKMILLER, rng);
    rand_setPMseed(s, rng);
  }

  // printI("pv", pv);
  // printI("qv", qv);
  // printI("rv", rv);
  // printI("nv", nv);
  // printI("Mv", Mv);
  // printM("A", Ap, rv, rv*pv);
  // printM("B", Bp, rv, rv*qv);
  // printM("Sig", Sigp, rv, rv);
  // printV("mu", mup, rv);
  // printP("X0p", X0p);  
  // printM("X", Xp, rv, nv);
  // printM("eps", epsp, rv, nv);
  // printI("okp", *okp);
  varmapack_sim(Ap, Bp, Sigp, mup, pv, qv, rv, nv, Mv, X0p, rng, Xp, epsp, okp);
  rand_free(rng);  
  return R_NilValue; // .Call wants a return value.
}
