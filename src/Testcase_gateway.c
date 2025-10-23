#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdint.h>
#include "VarmaSim.h"
#include "VarmaMisc.h"
#include "RandomNumbers.h" // RandRng, rand_create, rand_settype, rand_setPMseed, rand_free

// Name table for named cases (order must match your C Testcase’s enumeration)
static const char *CASE_NAMES[] = {
  "tinyAR","tinyMA","tinyARMA",
  "smallAR","smallMA","smallARMA1","smallARMA2",
  "mediumAR","mediumARMA1","mediumARMA2","mediumMA",
  "largeAR"
};

static int name_to_icase(const char *s) {
  int N = (int)(sizeof(CASE_NAMES) / sizeof(CASE_NAMES[0]));
  for (int i = 0; i < N; ++i) if (strcmp(s, CASE_NAMES[i]) == 0) return i + 1; // 1-based
  return 0; // treat unknown as random (your R wrapper won’t pass unknowns)
}

static SEXP make_array3(int r, int k) {
  SEXP dims = PROTECT(allocVector(INTSXP, 3));
  INTEGER(dims)[0] = r;
  INTEGER(dims)[1] = r;
  INTEGER(dims)[2] = (k > 0 ? k : 0);
  SEXP arr = PROTECT(allocArray(REALSXP, dims));
  UNPROTECT(2);
  return arr;
}

// .Call("Testcase_gateway", name, p, q, r, park_miller, seed)
SEXP Testcase_gateway(SEXP name_, SEXP p_, SEXP q_, SEXP r_,
                      SEXP park_miller_, SEXP seed_) {
  // RNG setup (no checks)
  RandRng *rng = rand_create();
  if (LOGICAL(park_miller_)[0]) {
    rand_settype(PARKMILLER, rng);
    rand_setPMseed((uint32_t)INTEGER(seed_)[0], rng);
  } else {
    rand_settype(DEFAULT_RNG, rng);
  }

  int p = 0, q = 0, r = 0;
  int icase = 0;
  const char *in_name = NULL;

  if (!isNull(name_)) {
    in_name = CHAR(STRING_ELT(name_, 0));
    icase = name_to_icase(in_name);

    // pre-query dims for named case
    char tmpname[64] = {0};
    Testcase(NULL, NULL, NULL, tmpname, &p, &q, &r, icase, rng);
  } else {
    // random case: use provided dims
    p = INTEGER(p_)[0];
    q = INTEGER(q_)[0];
    r = INTEGER(r_)[0];
    icase = 0;
  }

  // allocate outputs
  SEXP A   = PROTECT(make_array3(r, p));
  SEXP B   = PROTECT(make_array3(r, q));
  SEXP Sig = PROTECT(allocMatrix(REALSXP, r, r));

  // fill via C
  char outname[64] = {0};
  Testcase((p > 0) ? REAL(A) : NULL,
           (q > 0) ? REAL(B) : NULL,
           REAL(Sig),
           outname,
           &p, &q, &r, icase, rng);

  // pack list(A,B,Sig,p,q,r,name)
  SEXP Rp    = PROTECT(ScalarInteger(p));
  SEXP Rq    = PROTECT(ScalarInteger(q));
  SEXP Rr    = PROTECT(ScalarInteger(r));
  SEXP Rname = PROTECT(mkString(icase > 0 ? (outname[0] ? outname : in_name) : "random"));

  SEXP out = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(out, 0, A);
  SET_VECTOR_ELT(out, 1, B);
  SET_VECTOR_ELT(out, 2, Sig);
  SET_VECTOR_ELT(out, 3, Rp);
  SET_VECTOR_ELT(out, 4, Rq);
  SET_VECTOR_ELT(out, 5, Rr);
  SET_VECTOR_ELT(out, 6, Rname);

  SEXP nms = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(nms, 0, mkChar("A"));
  SET_STRING_ELT(nms, 1, mkChar("B"));
  SET_STRING_ELT(nms, 2, mkChar("Sig"));
  SET_STRING_ELT(nms, 3, mkChar("p"));
  SET_STRING_ELT(nms, 4, mkChar("q"));
  SET_STRING_ELT(nms, 5, mkChar("r"));
  SET_STRING_ELT(nms, 6, mkChar("name"));
  setAttrib(out, R_NamesSymbol, nms);

  rand_free(rng);
  UNPROTECT(9);
  return out;
}
