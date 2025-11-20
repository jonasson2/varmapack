#include "randompack.h"

#ifdef USING_R

#include <R.h>
#include <Rinternals.h>

static void rp_rng_finalizer(SEXP ext) {
  randompack_rng *rng = (randompack_rng *) R_ExternalPtrAddr(ext);
  if (rng) {
    randompack_free(rng);
    R_ClearExternalPtr(ext);
  }
}

static randompack_rng *rp_get_rng(SEXP ext) {
  if (TYPEOF(ext) != EXTPTRSXP)
    Rf_error("rp: expected external pointer");
  randompack_rng *rng = (randompack_rng *) R_ExternalPtrAddr(ext);
  if (!rng)
    Rf_error("rp: RNG pointer is NULL (already freed?)");
  return rng;
}

SEXP rp_create(SEXP typeSEXP, SEXP seedSEXP) {
  const char *type = NULL;
  if (!Rf_isNull(typeSEXP)) {
    if (!Rf_isString(typeSEXP) || Rf_length(typeSEXP) != 1)
      Rf_error("rp_create: `type` must be NULL or a length-1 character vector");
    type = CHAR(STRING_ELT(typeSEXP, 0));
  }
  int seed = Rf_isNull(seedSEXP) ? 0 : Rf_asInteger(seedSEXP);
  randompack_rng *rng = randompack_create(type, seed);
  if (!rng)
    Rf_error("rp_create: failed to allocate RNG");
  SEXP ext = PROTECT(R_MakeExternalPtr(rng, R_NilValue, R_NilValue));
  R_RegisterCFinalizer(ext, rp_rng_finalizer);
  UNPROTECT(1);
  return ext;
}

SEXP rp_free(SEXP rng_xptr) {
  randompack_rng *rng = rp_get_rng(rng_xptr);
  randompack_free(rng);
  R_ClearExternalPtr(rng_xptr);
  return R_NilValue;
}

SEXP rp_u01(SEXP rng_xptr, SEXP nSEXP) {
  randompack_rng *rng = rp_get_rng(rng_xptr);
  int n = Rf_asInteger(nSEXP);
  if (n < 0)
    Rf_error("rp_u01: `n` must be non-negative");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  if (n > 0)
    randompack_u01(REAL(out), n, rng);
  UNPROTECT(1);
  return out;
}

#else

typedef int rp_gateway_requires_using_r;

#endif /* USING_R */
