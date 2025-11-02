#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP VarmaSim_gateway(SEXP A, SEXP B, SEXP Sig, SEXP mu, SEXP p, SEXP q, SEXP r,
                      SEXP n, SEXP M, SEXP X0, SEXP X, SEXP eps, SEXP ok,
                      SEXP park_miller, SEXP seed);
SEXP OmegaBuild_gateway(SEXP S_, SEXP G_, SEXP W_, SEXP p_, SEXP q_, SEXP r_,
                        SEXP n_);
SEXP VYWFactorize_gateway(SEXP A_, SEXP p_, SEXP r_);
SEXP VYWSolve_gateway(SEXP A_, SEXP LU_, SEXP piv_, SEXP Y_, SEXP p_, SEXP r_);
SEXP FindCGW_gateway(SEXP A_, SEXP B_, SEXP Sig_);
SEXP Testcase_gateway(SEXP name_, SEXP p_, SEXP q_, SEXP r_, SEXP park_miller_, SEXP seed_);

static const R_CallMethodDef CallEntries[] = {
  {"VarmaSim_gateway",     (DL_FUNC) &VarmaSim_gateway,     15},
  {"OmegaBuild_gateway",   (DL_FUNC) &OmegaBuild_gateway,   7},
  {"VYWFactorize_gateway", (DL_FUNC) &VYWFactorize_gateway, 3},
  {"VYWSolve_gateway",     (DL_FUNC) &VYWSolve_gateway,     6},
  {"FindCGW_gateway",      (DL_FUNC) &FindCGW_gateway,      3},
  {"Testcase_gateway",     (DL_FUNC) &Testcase_gateway,     6},
  {NULL, NULL, 0}
};

void R_init_varmapack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

