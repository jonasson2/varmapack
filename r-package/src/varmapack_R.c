#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <string.h>
#include "varmapack.h"
#include "slepian_internal.h"

static int get_ndim(SEXP x) {
  return Rf_length(Rf_getAttrib(x, R_DimSymbol));
}

static int get_dim(SEXP x, int i) {
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  return INTEGER(dim)[i];
}

static SEXP new_series_array(int r, int n, int M) {
  R_xlen_t len = (R_xlen_t)r*n*M;
  SEXP x = PROTECT(Rf_allocVector(REALSXP, len));
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, 3));
  INTEGER(dim)[0] = r;
  INTEGER(dim)[1] = n;
  INTEGER(dim)[2] = M;
  Rf_setAttrib(x, R_DimSymbol, dim);
  UNPROTECT(2);
  return x;
}

static int lag_count_R(SEXP lags, int r) {
  R_xlen_t n = Rf_xlength(lags);
  if (TYPEOF(lags) != REALSXP || get_ndim(lags) != 0 || n > INT_MAX)
    Rf_error("lags must be a real vector of supported length");
  if (n && (R_xlen_t)r*r > R_XLEN_T_MAX/n)
    Rf_error("problem size too large");
  return n;
}

static randompack_rng *get_rng(SEXP rng_, bool *own_rng) {
  randompack_rng *rng;
  *own_rng = false;
  if (Rf_isNull(rng_)) {
    rng = randompack_create(0);
    *own_rng = true;
  }
  else {
    rng = randompack_rng_from_R(rng_);
  }
  return rng;
}

static void free_rng(randompack_rng *rng, bool own_rng) {
  if (own_rng) randompack_free(rng);
}

static void check_ok(bool ok, int nprotected) {
  char *message;
  if (ok) return;
  message = varmapack_last_error();
  UNPROTECT(nprotected);
  Rf_error("%s", message ? message : "Varmapack operation failed");
}

static SEXP named_list2(SEXP x, SEXP y, char *xname, char *yname) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_VECTOR_ELT(out, 0, x);
  SET_VECTOR_ELT(out, 1, y);
  SET_STRING_ELT(names, 0, Rf_mkChar(xname));
  SET_STRING_ELT(names, 1, Rf_mkChar(yname));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

static int ou_dimension_R(SEXP parameter, SEXP cov, SEXP mu) {
  if (TYPEOF(parameter) != REALSXP || TYPEOF(cov) != REALSXP ||
      TYPEOF(mu) != REALSXP || get_ndim(cov) != 2)
    Rf_error("invalid OU model dimensions or types");
  int r = get_dim(cov, 0), ndim = get_ndim(parameter);
  bool vector = ndim == 0 && Rf_xlength(parameter) == r;
  bool matrix = ndim == 2 && get_dim(parameter, 0) == r &&
                get_dim(parameter, 1) == r;
  if (r < 1 || get_dim(cov, 1) != r || Rf_xlength(mu) != r || !(vector || matrix))
    Rf_error("invalid OU model dimensions or types");
  return r;
}

static varmapack_ou_model *ou_model_R(SEXP parameter, SEXP cov, SEXP mu, int r) {
  if (get_ndim(parameter) == 2)
    return varmapack_ou_create_drift(REAL(parameter), REAL(cov), REAL(mu), r);
  return varmapack_ou_create(REAL(parameter), REAL(cov), REAL(mu), r);
}

SEXP varmapack_ou_validate_R(SEXP parameter, SEXP cov, SEXP mu) {
  int r = ou_dimension_R(parameter, cov, mu);
  varmapack_ou_model *model = ou_model_R(parameter, cov, mu, r);
  check_ok(model != 0, 0);
  varmapack_ou_free(model);
  return R_NilValue;
}

SEXP varmapack_ou_sim_R(SEXP parameter, SEXP cov, SEXP mu, SEXP t,
                        SEXP M_, SEXP rng_) {
  int r = ou_dimension_R(parameter, cov, mu), n = Rf_length(t), M = Rf_asInteger(M_);
  bool own_rng;
  if (TYPEOF(t) != REALSXP || n < 1 || M < 1) Rf_error("invalid OU simulation");
  SEXP X = PROTECT(new_series_array(r, n, M));
  randompack_rng *rng = get_rng(rng_, &own_rng);
  if (!rng) Rf_error("could not obtain a Randompack generator");
  varmapack_ou_model *model = ou_model_R(parameter, cov, mu, r);
  bool ok = model && varmapack_ou_sim(model, REAL(t), n, M, REAL(X), rng);
  varmapack_ou_free(model);
  free_rng(rng, own_rng);
  check_ok(ok, 1);
  UNPROTECT(1);
  return X;
}

static int slepian_dimension_R(SEXP theta, SEXP cov, SEXP mu) {
  if (TYPEOF(theta) != REALSXP || TYPEOF(cov) != REALSXP ||
      TYPEOF(mu) != REALSXP || get_ndim(theta) != 2 || get_ndim(cov) != 2)
    Rf_error("invalid Slepian model dimensions or types");
  int r = get_dim(theta, 0);
  if (r < 1 || get_dim(theta, 1) != r || get_dim(cov, 0) != r ||
      get_dim(cov, 1) != r || Rf_xlength(mu) != r)
    Rf_error("invalid Slepian model dimensions or types");
  return r;
}

SEXP varmapack_ou_acvf_R(SEXP parameter, SEXP cov, SEXP mu, SEXP lags) {
  int r = ou_dimension_R(parameter, cov, mu), n = lag_count_R(lags, r);
  SEXP out = PROTECT(new_series_array(r, r, n));
  varmapack_ou_model *model = ou_model_R(parameter, cov, mu, r);
  bool ok = model && varmapack_ou_acvf(model, REAL(lags), n, REAL(out));
  varmapack_ou_free(model);
  check_ok(ok, 1);
  UNPROTECT(1);
  return out;
}

SEXP varmapack_ou_specrad_R(SEXP parameter, SEXP cov, SEXP mu, SEXP gap) {
  int r = ou_dimension_R(parameter, cov, mu);
  if (TYPEOF(gap) != REALSXP || Rf_xlength(gap) != 1)
    Rf_error("gap must be a real scalar");
  varmapack_ou_model *model = ou_model_R(parameter, cov, mu, r);
  check_ok(model != 0, 0);
  double rho = varmapack_ou_specrad(model, REAL(gap)[0]);
  varmapack_ou_free(model);
  check_ok(!isnan(rho), 0);
  return Rf_ScalarReal(rho);
}

SEXP varmapack_slepian_acvf_R(SEXP theta, SEXP factor, SEXP mu, SEXP cov,
                              SEXP lags) {
  int r = slepian_dimension_R(theta, factor, mu), n = lag_count_R(lags, r);
  slepian_dimension_R(theta, cov, mu);
  SEXP out = PROTECT(new_series_array(r, r, n));
  varmapack_slepian_model model = {r, REAL(theta), REAL(mu), REAL(factor), REAL(cov)};
  check_ok(varmapack_slepian_acvf(&model, REAL(lags), n, REAL(out)), 1);
  UNPROTECT(1);
  return out;
}

SEXP varmapack_slepian_factor_R(SEXP theta, SEXP cov, SEXP mu) {
  int r = slepian_dimension_R(theta, cov, mu);
  SEXP factor = PROTECT(Rf_allocMatrix(REALSXP, r, r));
  varmapack_slepian_model *model = varmapack_slepian_create(REAL(theta), REAL(cov),
                                                          REAL(mu), r);
  check_ok(model != 0, 0);
  memcpy(REAL(factor), model->factor, (size_t)r*r*sizeof(double));
  varmapack_slepian_free(model);
  UNPROTECT(1);
  return factor;
}

SEXP varmapack_slepian_sim_R(SEXP theta, SEXP factor, SEXP mu, SEXP t,
                             SEXP M_, SEXP rng_) {
  int r = slepian_dimension_R(theta, factor, mu);
  int n = Rf_length(t), M = Rf_asInteger(M_);
  bool own_rng;
  if (TYPEOF(t) != REALSXP || n < 1 || M < 1)
    Rf_error("invalid Slepian simulation");
  if ((R_xlen_t)r*n > R_XLEN_T_MAX/M) Rf_error("problem size too large");
  SEXP X = PROTECT(new_series_array(r, n, M));
  randompack_rng *rng = get_rng(rng_, &own_rng);
  if (!rng) Rf_error("could not obtain a Randompack generator");
  varmapack_slepian_model model = {r, REAL(theta), REAL(mu), REAL(factor), 0};
  bool ok = varmapack_slepian_sim(&model, REAL(t), n, M, REAL(X), rng);
  free_rng(rng, own_rng);
  check_ok(ok, 1);
  UNPROTECT(1);
  return X;
}

SEXP varmapack_gamma_times_R(SEXP n_, SEXP shape, SEXP mean_gap, SEXP rng_) {
  int n = Rf_asInteger(n_);
  bool own_rng;
  SEXP t = PROTECT(Rf_allocVector(REALSXP, n));
  randompack_rng *rng = get_rng(rng_, &own_rng);
  if (!rng) Rf_error("could not obtain a Randompack generator");
  bool ok = varmapack_gamma_times(n, Rf_asReal(shape), Rf_asReal(mean_gap),
                                   REAL(t), rng);
  free_rng(rng, own_rng);
  check_ok(ok, 1);
  UNPROTECT(1);
  return t;
}

static SEXP named_list3(SEXP x, SEXP y, SEXP z, char *xname, char *yname,
                        char *zname) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_VECTOR_ELT(out, 0, x);
  SET_VECTOR_ELT(out, 1, y);
  SET_VECTOR_ELT(out, 2, z);
  SET_STRING_ELT(names, 0, Rf_mkChar(xname));
  SET_STRING_ELT(names, 1, Rf_mkChar(yname));
  SET_STRING_ELT(names, 2, Rf_mkChar(zname));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

SEXP varmapack_sim_R(SEXP A, SEXP B, SEXP C, SEXP Sig, SEXP mu, SEXP X0, SEXP z,
                     SEXP n_, SEXP M_, SEXP return_shocks_, SEXP rng_) {
  int n = Rf_asInteger(n_);
  int M = Rf_asInteger(M_);
  int return_shocks = Rf_asLogical(return_shocks_);
  int r = get_dim(Sig, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  int nmu = Rf_isNull(mu) ? 0 : get_dim(mu, 1);
  int nX0 = Rf_isNull(X0) ? 0 : get_dim(X0, 1);
  int MX0 = Rf_isNull(X0) || get_ndim(X0) == 2 ? 1 : get_dim(X0, 2);
  int s = Rf_isNull(C) ? 0 : get_dim(C, 2);
  int d = Rf_isNull(C) ? 0 : get_dim(C, 1);
  int Mz = Rf_isNull(z) || get_ndim(z) == 2 ? 1 : get_dim(z, 2);
  int h = Rf_isNull(X0) ? 0 : get_dim(X0, 1);
  int nprotected = 0;
  bool ok;
  bool own_rng;
  randompack_rng *rng;
  SEXP X = PROTECT(new_series_array(r, n, M));
  SEXP E = R_NilValue;
  nprotected++;
  if (return_shocks) {
    E = PROTECT(new_series_array(r, n, M));
    nprotected++;
  }
  rng = get_rng(rng_, &own_rng);
  if (!rng) {
    UNPROTECT(nprotected);
    Rf_error("could not obtain a Randompack random number generator");
  }
  if (s) {
    ok = varmapack_simx(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                        REAL(C), REAL(Sig), REAL(z), Mz, p, q, s, d, r, n, M,
                        Rf_isNull(X0) ? 0 : REAL(X0), h, MX0, REAL(X),
                        return_shocks ? REAL(E) : 0, rng);
  }
  else {
    ok = varmapack_sim(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                        REAL(Sig), Rf_isNull(mu) ? 0 : REAL(mu), nmu, p, q, r, n, M,
                        Rf_isNull(X0) ? 0 : REAL(X0), nX0, MX0, REAL(X),
                        return_shocks ? REAL(E) : 0, rng);
  }
  free_rng(rng, own_rng);
  check_ok(ok, nprotected);
  if (!return_shocks) {
    UNPROTECT(nprotected);
    return X;
  }
  SEXP out = PROTECT(named_list2(X, E, "X", "E"));
  UNPROTECT(nprotected + 1);
  return out;
}

SEXP varmapack_acvf_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag_) {
  int r = get_dim(Sig, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  int maxlag = Rf_asInteger(maxlag_);
  SEXP Gamma = PROTECT(Rf_alloc3DArray(REALSXP, r, r, maxlag + 1));
  bool ok = varmapack_acvf(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                           REAL(Sig), p, q, r, REAL(Gamma), maxlag);
  check_ok(ok, 1);
  UNPROTECT(1);
  return Gamma;
}

SEXP varmapack_cov2sig_R(SEXP A, SEXP B, SEXP cov) {
  int r = get_dim(cov, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  SEXP Sig = PROTECT(Rf_allocMatrix(REALSXP, r, r));
  bool ok = varmapack_cov2sig(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                             REAL(cov), p, q, r, REAL(Sig));
  check_ok(ok, 1);
  UNPROTECT(1);
  return Sig;
}

SEXP varmapack_psi_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag_) {
  int r = get_dim(Sig, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  int maxlag = Rf_asInteger(maxlag_);
  SEXP Psi = PROTECT(Rf_alloc3DArray(REALSXP, r, r, maxlag + 1));
  bool ok = varmapack_psi(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                          p, q, r, maxlag, REAL(Psi));
  check_ok(ok, 1);
  UNPROTECT(1);
  return Psi;
}

SEXP varmapack_irf_R(SEXP A, SEXP B, SEXP Sig, SEXP maxlag_) {
  int r = get_dim(Sig, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  int maxlag = Rf_asInteger(maxlag_);
  SEXP Theta = PROTECT(Rf_alloc3DArray(REALSXP, r, r, maxlag + 1));
  bool ok = varmapack_irf(Rf_isNull(A) ? 0 : REAL(A), Rf_isNull(B) ? 0 : REAL(B),
                          REAL(Sig), p, q, r, maxlag, REAL(Theta));
  check_ok(ok, 1);
  UNPROTECT(1);
  return Theta;
}

SEXP varmapack_specrad_R(SEXP A, SEXP Sig) {
  int r = get_dim(Sig, 0);
  int p = Rf_isNull(A) ? 0 : get_dim(A, 2);
  double rho = varmapack_specrad(Rf_isNull(A) ? 0 : REAL(A), r, p);
  if (isnan(rho)) {
    char *message = varmapack_last_error();
    Rf_error("%s", message ? message : "Varmapack spectral-radius calculation failed");
  }
  return Rf_ScalarReal(rho);
}

SEXP varmapack_ma_specrad_R(SEXP B, SEXP Sig) {
  int r = get_dim(Sig, 0);
  int q = Rf_isNull(B) ? 0 : get_dim(B, 2);
  double rho = varmapack_ma_specrad(Rf_isNull(B) ? 0 : REAL(B), r, q);
  if (isnan(rho)) {
    char *message = varmapack_last_error();
    Rf_error("%s", message ? message : "Varmapack spectral-radius calculation failed");
  }
  return Rf_ScalarReal(rho);
}

SEXP varmapack_autocov_R(SEXP X, SEXP maxlag_, SEXP norm) {
  int r = get_dim(X, 0);
  int n = get_dim(X, 1);
  int maxlag = Rf_asInteger(maxlag_);
  const char *norm_string = CHAR(STRING_ELT(norm, 0));
  SEXP C = PROTECT(Rf_alloc3DArray(REALSXP, r, r, maxlag + 1));
  bool ok = varmapack_autocov("N", norm_string, r, n, REAL(X), maxlag, REAL(C));
  check_ok(ok, 1);
  UNPROTECT(1);
  return C;
}

SEXP varmapack_cov2corr_R(SEXP cov) {
  int ndim = get_ndim(cov);
  int r = get_dim(cov, 0);
  int maxlag = ndim == 2 ? 0 : get_dim(cov, 2) - 1;
  SEXP corr = PROTECT(Rf_duplicate(cov));
  bool ok = varmapack_cov2corr(REAL(cov), r, maxlag, REAL(corr));
  check_ok(ok, 1);
  UNPROTECT(1);
  return corr;
}

SEXP varmapack_testcase_R(SEXP name_, SEXP index_, SEXP rho_, SEXP p_, SEXP q_, SEXP r_,
                          SEXP rng_) {
  char name[VARMAPACK_TESTCASE_NAME_LEN];
  int index = Rf_asInteger(index_);
  int p = Rf_asInteger(p_);
  int q = Rf_asInteger(q_);
  int r = Rf_asInteger(r_);
  double rho = Rf_asReal(rho_);
  double dummy = 0;
  bool own_rng = false;
  bool ok;
  randompack_rng *rng = 0;
  SEXP A = R_NilValue;
  SEXP B = R_NilValue;
  SEXP Sig;
  if (strlen(CHAR(STRING_ELT(name_, 0))) >= VARMAPACK_TESTCASE_NAME_LEN) {
    Rf_error("testcase name is too long");
  }
  strcpy(name, CHAR(STRING_ELT(name_, 0)));
  if ((name[0] && strcmp(name, "rho") != 0) || index > 0) {
    ok = varmapack_testcase(name, &index, rho, &p, &q, &r, 0, 0, 0, 0);
    check_ok(ok, 0);
  }
  if (index == 0 && strcmp(name, "rho") != 0) {
    rng = get_rng(rng_, &own_rng);
    if (!rng) Rf_error("could not obtain a Randompack random number generator");
  }
  if (p) A = PROTECT(Rf_alloc3DArray(REALSXP, r, r, p));
  if (q) B = PROTECT(Rf_alloc3DArray(REALSXP, r, r, q));
  Sig = PROTECT(Rf_allocMatrix(REALSXP, r, r));
  ok = varmapack_testcase(name, &index, rho, &p, &q, &r,
                          Rf_isNull(A) ? &dummy : REAL(A),
                          Rf_isNull(B) ? &dummy : REAL(B), REAL(Sig), rng);
  free_rng(rng, own_rng);
  check_ok(ok, 1 + (p > 0) + (q > 0));
  SEXP out = PROTECT(named_list3(A, B, Sig, "A", "B", "Sig"));
  UNPROTECT(2 + (p > 0) + (q > 0));
  return out;
}

SEXP varmapack_testcases_R(void) {
  char name[VARMAPACK_TESTCASE_NAME_LEN] = "max";
  int index = 0;
  int p = 0;
  int q = 0;
  int r = 0;
  bool ok = varmapack_testcase(name, &index, 0, &p, &q, &r, 0, 0, 0, 0);
  check_ok(ok, 0);
  int ncase = index;
  SEXP indices = PROTECT(Rf_allocVector(INTSXP, ncase));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, ncase));
  SEXP ps = PROTECT(Rf_allocVector(INTSXP, ncase));
  SEXP qs = PROTECT(Rf_allocVector(INTSXP, ncase));
  SEXP rs = PROTECT(Rf_allocVector(INTSXP, ncase));
  for (int i=0; i<ncase; i++) {
    name[0] = '\0';
    index = i + 1;
    p = q = r = 0;
    ok = varmapack_testcase(name, &index, 0, &p, &q, &r, 0, 0, 0, 0);
    check_ok(ok, 5);
    INTEGER(indices)[i] = index;
    SET_STRING_ELT(names, i, Rf_mkChar(name));
    INTEGER(ps)[i] = p;
    INTEGER(qs)[i] = q;
    INTEGER(rs)[i] = r;
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SEXP outnames = PROTECT(Rf_allocVector(STRSXP, 5));
  SET_VECTOR_ELT(out, 0, indices);
  SET_VECTOR_ELT(out, 1, names);
  SET_VECTOR_ELT(out, 2, ps);
  SET_VECTOR_ELT(out, 3, qs);
  SET_VECTOR_ELT(out, 4, rs);
  SET_STRING_ELT(outnames, 0, Rf_mkChar("index"));
  SET_STRING_ELT(outnames, 1, Rf_mkChar("name"));
  SET_STRING_ELT(outnames, 2, Rf_mkChar("p"));
  SET_STRING_ELT(outnames, 3, Rf_mkChar("q"));
  SET_STRING_ELT(outnames, 4, Rf_mkChar("r"));
  Rf_setAttrib(out, R_NamesSymbol, outnames);
  UNPROTECT(7);
  return out;
}
