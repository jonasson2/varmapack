#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static int get_int_scalar(const mxArray *arg, const char *name);
static void check_varmapack_error(bool ok);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int s, d, r, n;
  bool ok;
  if (nrhs != 4)
    mexErrMsgIdAndTxt("varmapack:testcasex:nrhs", "Wrong number of input arguments");
  if (nlhs != 2)
    mexErrMsgIdAndTxt("varmapack:testcasex:nlhs", "Two output arguments are required");
  s = get_int_scalar(prhs[0], "s");
  d = get_int_scalar(prhs[1], "d");
  r = get_int_scalar(prhs[2], "r");
  n = get_int_scalar(prhs[3], "n");
  if (r == 0)
    mexErrMsgIdAndTxt("varmapack:testcasex:r", "r must be positive");
  plhs[0] = mxCreateDoubleMatrix(r, (mwSize)d*s, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(d, n, mxREAL);
  ok = varmapack_testcasex(s, d, r, n, mxGetPr(plhs[0]), mxGetPr(plhs[1]));
  check_varmapack_error(ok);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:testcasex:argument", "%s must be a numeric scalar",
                      name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:testcasex:argument", "%s must be a nonnegative integer",
                      name);
  return (int)x;
}

static void check_varmapack_error(bool ok) {
  char *message;
  if (ok) return;
  message = varmapack_last_error();
  mexErrMsgIdAndTxt("varmapack:testcasex:error", "%s",
                    message ? message : "varmapack error");
}
