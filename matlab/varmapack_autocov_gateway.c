#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static double *get_double_matrix(const mxArray *arg, const char *name);
static int get_int_scalar(const mxArray *arg, const char *name);
static void get_norm(const mxArray *arg, char norm[16]);
static void check_varmapack_error(varmapack_error error);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int r, n, maxlag;
  char norm[16];
  double *X, *C;
  mwSize dim[3];
  varmapack_error error;
  if (nrhs != 3)
    mexErrMsgIdAndTxt("varmapack:autocov:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:autocov:nlhs", "Wrong number of output arguments");
  X = get_double_matrix(prhs[0], "X");
  if (mxGetM(prhs[0]) > INT_MAX || mxGetN(prhs[0]) > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:autocov:X", "X is too large");
  r = (int)mxGetM(prhs[0]);
  n = (int)mxGetN(prhs[0]);
  maxlag = get_int_scalar(prhs[1], "maxlag");
  get_norm(prhs[2], norm);
  dim[0] = r;
  dim[1] = r;
  dim[2] = maxlag + 1;
  plhs[0] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
  C = mxGetPr(plhs[0]);
  error = varmapack_autocov("N", norm, r, n, X, maxlag, C);
  check_varmapack_error(error);
}

static double *get_double_matrix(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg) || mxGetNumberOfDimensions(arg) != 2)
    mexErrMsgIdAndTxt("varmapack:autocov:argument", "%s must be a real double matrix",
                      name);
  return mxGetPr(arg);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:autocov:argument", "%s must be a numeric scalar",
                      name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:autocov:argument",
                      "%s must be a nonnegative integer", name);
  return (int)x;
}

static void get_norm(const mxArray *arg, char norm[16]) {
  if (!mxIsChar(arg))
    mexErrMsgIdAndTxt("varmapack:autocov:norm", "norm must be a string");
  if (mxGetString(arg, norm, 16) != 0)
    mexErrMsgIdAndTxt("varmapack:autocov:norm", "norm string is too long");
  if (norm[0] == 'c' || norm[0] == 'C') {
    norm[0] = 'C';
    norm[1] = 0;
  }
  else if (norm[0] == 'm' || norm[0] == 'M') {
    norm[0] = 'M';
    norm[1] = 'L';
    norm[2] = 0;
  }
}

static void check_varmapack_error(varmapack_error error) {
  if (error == VARMAPACK_OK) return;
  switch (error) {
    case VARMAPACK_INVALID_ARGUMENT:
      mexErrMsgIdAndTxt("varmapack:autocov:invalidArgument", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_ALLOCATION:
      mexErrMsgIdAndTxt("varmapack:autocov:allocation", "%s",
                        varmapack_strerror(error));
      break;
    default:
      mexErrMsgIdAndTxt("varmapack:autocov:error", "%s", varmapack_strerror(error));
      break;
  }
}
