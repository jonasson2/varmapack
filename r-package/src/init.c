#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP varmapack_sim_R(SEXP A, SEXP B, SEXP C, SEXP Sig, SEXP mu, SEXP X0, SEXP z,
                     SEXP n, SEXP M, SEXP return_shocks, SEXP rng);
SEXP varmapack_acvf_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag);
SEXP varmapack_psi_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag);
SEXP varmapack_irf_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag);
SEXP varmapack_specrad_R(SEXP A, SEXP Sig);
SEXP varmapack_ma_specrad_R(SEXP B, SEXP Sig);
SEXP varmapack_autocov_R(SEXP X, SEXP maxlag, SEXP norm);
SEXP varmapack_cov2corr_R(SEXP cov);
SEXP varmapack_testcase_R(SEXP name, SEXP index, SEXP rho, SEXP p, SEXP q, SEXP r,
                          SEXP rng);
SEXP varmapack_testcases_R(void);

static const R_CallMethodDef CallEntries[] = {
  {"varmapack_sim_R", (DL_FUNC)&varmapack_sim_R, 11},
  {"varmapack_acvf_R", (DL_FUNC)&varmapack_acvf_R, 4},
  {"varmapack_psi_R", (DL_FUNC)&varmapack_psi_R, 4},
  {"varmapack_irf_R", (DL_FUNC)&varmapack_irf_R, 4},
  {"varmapack_specrad_R", (DL_FUNC)&varmapack_specrad_R, 2},
  {"varmapack_ma_specrad_R", (DL_FUNC)&varmapack_ma_specrad_R, 2},
  {"varmapack_autocov_R", (DL_FUNC)&varmapack_autocov_R, 3},
  {"varmapack_cov2corr_R", (DL_FUNC)&varmapack_cov2corr_R, 1},
  {"varmapack_testcase_R", (DL_FUNC)&varmapack_testcase_R, 7},
  {"varmapack_testcases_R", (DL_FUNC)&varmapack_testcases_R, 0},
  {0, 0, 0}
};

void R_init_varmapack(DllInfo *dll) {
  R_registerRoutines(dll, 0, CallEntries, 0, 0);
  R_useDynamicSymbols(dll, FALSE);
}
