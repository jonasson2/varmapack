#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP varmapack_sim_gateway(SEXP A, SEXP B, SEXP Sig, SEXP mu, SEXP p, SEXP q, SEXP r,
                      SEXP n, SEXP M, SEXP X0, SEXP X, SEXP eps, SEXP ok,
                      SEXP park_miller, SEXP seed);
SEXP OmegaBuild_gateway(SEXP S_, SEXP G_, SEXP W_, SEXP p_, SEXP q_, SEXP r_,
                        SEXP n_);
SEXP VYWFactorize_gateway(SEXP A_, SEXP p_, SEXP r_);
SEXP VYWSolve_gateway(SEXP A_, SEXP LU_, SEXP piv_, SEXP Y_, SEXP p_, SEXP r_);
SEXP vpack_FindCG_gateway(SEXP A_, SEXP B_, SEXP Sig_);
SEXP rp_create(SEXP type, SEXP seed);
SEXP rp_free(SEXP rng);
SEXP rp_u01(SEXP rng, SEXP n);
SEXP Testcase_gateway(SEXP name_, SEXP p_, SEXP q_, SEXP r_, SEXP park_miller_, SEXP seed_);

static const R_CallMethodDef CallEntries[] = {
  {"varmapack_sim_gateway",     (DL_FUNC) &varmapack_sim_gateway,     15},
  {"OmegaBuild_gateway",   (DL_FUNC) &OmegaBuild_gateway,   7},
  {"VYWFactorize_gateway", (DL_FUNC) &VYWFactorize_gateway, 3},
  {"VYWSolve_gateway",     (DL_FUNC) &VYWSolve_gateway,     6},
  {"vpack_FindCG_gateway",      (DL_FUNC) &vpack_FindCG_gateway,      3},
  {"rp_create",            (DL_FUNC) &rp_create,            2},
  {"rp_free",              (DL_FUNC) &rp_free,              1},
  {"rp_u01",               (DL_FUNC) &rp_u01,               2},
  {"Testcase_gateway",     (DL_FUNC) &Testcase_gateway,     6},
  {0, 0, 0}
};

void R_init_varmapack(DllInfo *dll) {
    R_registerRoutines(dll, 0, CallEntries, 0, 0);
    R_useDynamicSymbols(dll, FALSE);
}
