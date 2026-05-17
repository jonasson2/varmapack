#include <stdint.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include "mex.h"
#include "randompack.h"
#include "varmapack.h"

static randompack_rng *get_rng(const mxArray *arg);
static int get_int_scalar(const mxArray *arg, const char *name);
static double get_double_scalar(const mxArray *arg, const char *name);
static void get_string(const mxArray *arg, char text[64], const char *name);
static void get_count(int *count);
static void query_case(char name[64], int *p, int *q, int *r, int *icase);
static void create_outputs(mxArray *plhs[], int p, int q, int r);
static void set_extra_outputs(mxArray *plhs[], int nlhs, int nrhs, bool indexed,
                              char name[64], int p, int q, int r, int icase);
static void check_varmapack_error(varmapack_error error);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int p = 0, q = 0, r = 0, icase = 0;
  double rho = 0;
  char name[64] = "";
  bool indexed = false;
  randompack_rng *rng = 0;
  double Adummy = 0, Bdummy = 0, *Aout, *Bout;
  varmapack_error error;
  if (nrhs != 1 && nrhs != 3 && nrhs != 4)
    mexErrMsgIdAndTxt("varmapack:testcase:nrhs", "Wrong number of inputs");
  if (nrhs == 1 && mxIsChar(prhs[0])) {
    get_string(prhs[0], name, "name");
    if (strcmp(name, "count") == 0) {
      get_count(&icase);
      if (nlhs > 1)
        mexErrMsgIdAndTxt("varmapack:testcase:nlhs", "Too many outputs");
      if (nlhs == 1) plhs[0] = mxCreateDoubleScalar(icase);
      return;
    }
    query_case(name, &p, &q, &r, &icase);
  }
  else if (nrhs == 1) {
    indexed = true;
    icase = get_int_scalar(prhs[0], "i");
    query_case(name, &p, &q, &r, &icase);
  }
  else {
    p = get_int_scalar(prhs[0], "p");
    q = get_int_scalar(prhs[1], "q");
    r = get_int_scalar(prhs[2], "r");
    if (r <= 0)
      mexErrMsgIdAndTxt("varmapack:testcase:argument", "r must be positive");
    icase = nrhs == 4 ? 0 : -1;
    if (nrhs == 4 && mxIsUint64(prhs[3])) {
      rng = get_rng(prhs[3]);
    }
    else if (nrhs == 4) {
      rho = get_double_scalar(prhs[3], "rho");
      strcpy(name, "rho");
    }
  }
  if (nlhs < 3)
    mexErrMsgIdAndTxt("varmapack:testcase:nlhs", "At least 3 outputs required");
  if ((nrhs == 1 && nlhs > 7) || (nrhs != 1 && nlhs > 6))
    mexErrMsgIdAndTxt("varmapack:testcase:nlhs", "Too many outputs");
  create_outputs(plhs, p, q, r);
  Aout = mxGetPr(plhs[0]);
  Bout = mxGetPr(plhs[1]);
  if (Aout == 0) Aout = &Adummy;
  if (Bout == 0) Bout = &Bdummy;
  error = varmapack_testcase(Aout, Bout, mxGetPr(plhs[2]), name, &p, &q, &r,
                             &icase, rho, rng);
  check_varmapack_error(error);
  set_extra_outputs(plhs, nlhs, nrhs, indexed, name, p, q, r, icase);
}

static randompack_rng *get_rng(const mxArray *arg) {
  uint64_t value;
  if (!mxIsUint64(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:testcase:rng", "rng must be a uint64 scalar");
  value = *(uint64_t *)mxGetData(arg);
  if (value == 0)
    mexErrMsgIdAndTxt("varmapack:testcase:rng", "rng handle is zero");
  return (randompack_rng *)(uintptr_t)value;
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:testcase:argument", "%s must be a numeric scalar", name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:testcase:argument",
                      "%s must be a nonnegative integer", name);
  return (int)x;
}

static double get_double_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:testcase:argument", "%s must be a numeric scalar", name);
  x = mxGetScalar(arg);
  if (x < 0)
    mexErrMsgIdAndTxt("varmapack:testcase:argument", "%s must be nonnegative", name);
  return x;
}

static void get_string(const mxArray *arg, char text[64], const char *name) {
  if (!mxIsChar(arg))
    mexErrMsgIdAndTxt("varmapack:testcase:argument", "%s must be a string", name);
  if (mxGetString(arg, text, 64) != 0)
    mexErrMsgIdAndTxt("varmapack:testcase:argument", "%s is too long", name);
}

static void get_count(int *count) {
  int p = 0, q = 0, r = 0;
  char maxname[64] = "max";
  varmapack_error error;
  *count = 0;
  error = varmapack_testcase(0, 0, 0, maxname, &p, &q, &r, count, 0, 0);
  check_varmapack_error(error);
}

static void query_case(char name[64], int *p, int *q, int *r, int *icase) {
  varmapack_error error;
  error = varmapack_testcase(0, 0, 0, name, p, q, r, icase, 0, 0);
  check_varmapack_error(error);
}

static void create_outputs(mxArray *plhs[], int p, int q, int r) {
  plhs[0] = mxCreateDoubleMatrix(r, r*p, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(r, r*q, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(r, r, mxREAL);
}

static void set_extra_outputs(mxArray *plhs[], int nlhs, int nrhs, bool indexed,
                              char name[64], int p, int q, int r, int icase) {
  if (nlhs <= 3) return;
  if (nrhs == 1 && nlhs == 4) {
    plhs[3] = indexed ? mxCreateString(name) : mxCreateDoubleScalar(icase);
    return;
  }
  plhs[3] = mxCreateDoubleScalar(p);
  if (nlhs >= 5) plhs[4] = mxCreateDoubleScalar(q);
  if (nlhs >= 6) plhs[5] = mxCreateDoubleScalar(r);
  if (nlhs >= 7)
    plhs[6] = indexed ? mxCreateString(name) : mxCreateDoubleScalar(icase);
}

static void check_varmapack_error(varmapack_error error) {
  if (error == VARMAPACK_OK) return;
  switch (error) {
    case VARMAPACK_INVALID_ARGUMENT:
      mexErrMsgIdAndTxt("varmapack:testcase:invalidArgument", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_UNKNOWN_TESTCASE:
      mexErrMsgIdAndTxt("varmapack:testcase:unknownTestcase", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_ALLOCATION:
      mexErrMsgIdAndTxt("varmapack:testcase:allocation", "%s", varmapack_strerror(error));
      break;
    default:
      mexErrMsgIdAndTxt("varmapack:testcase:error", "%s", varmapack_strerror(error));
      break;
  }
}
