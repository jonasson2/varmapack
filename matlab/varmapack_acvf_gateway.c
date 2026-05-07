#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static double *get_double_array(const mxArray *arg, const char *name);
static int get_int_scalar(const mxArray *arg, const char *name);
static void check_dims(const mxArray *A, const mxArray *B, const mxArray *Sig,
                       int *p, int *q, int *r);
static void check_varmapack_error(varmapack_error error);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int p, q, r, maxlag;
  double *A, *B, *Sig, *Gamma;
  mwSize dim[3];
  varmapack_error error;
  if (nrhs != 4)
    mexErrMsgIdAndTxt("varmapack:acvf:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:acvf:nlhs", "Wrong number of output arguments");
  A = get_double_array(prhs[0], "A");
  B = get_double_array(prhs[1], "B");
  Sig = get_double_array(prhs[2], "Sig");
  maxlag = get_int_scalar(prhs[3], "maxlag");
  check_dims(prhs[0], prhs[1], prhs[2], &p, &q, &r);
  dim[0] = r;
  dim[1] = r;
  dim[2] = maxlag + 1;
  plhs[0] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
  Gamma = mxGetPr(plhs[0]);
  error = varmapack_acvf(A, B, Sig, p, q, r, Gamma, maxlag);
  check_varmapack_error(error);
}

static double *get_double_array(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg))
    mexErrMsgIdAndTxt("varmapack:acvf:argument", "%s must be a real double array",
                      name);
  return mxGetPr(arg);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:acvf:argument", "%s must be a numeric scalar",
                      name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:acvf:argument",
                      "%s must be a nonnegative integer", name);
  return (int)x;
}

static void check_dims(const mxArray *A, const mxArray *B, const mxArray *Sig,
                       int *p, int *q, int *r) {
  mwSize rSig = mxGetM(Sig);
  mwSize cSig = mxGetN(Sig);
  mwSize cA = mxGetN(A);
  mwSize cB = mxGetN(B);
  if (rSig == 0 || cSig != rSig)
    mexErrMsgIdAndTxt("varmapack:acvf:Sig", "Sig must be square");
  if (rSig > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:acvf:Sig", "Sig is too large");
  if (mxGetM(A) != rSig)
    mexErrMsgIdAndTxt("varmapack:acvf:A", "A must have r rows");
  if (mxGetM(B) != rSig)
    mexErrMsgIdAndTxt("varmapack:acvf:B", "B must have r rows");
  if (cA % rSig != 0)
    mexErrMsgIdAndTxt("varmapack:acvf:A", "size(A,2) must be a multiple of r");
  if (cB % rSig != 0)
    mexErrMsgIdAndTxt("varmapack:acvf:B", "size(B,2) must be a multiple of r");
  if (cA/rSig > INT_MAX || cB/rSig > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:acvf:argument", "A or B is too large");
  *r = (int)rSig;
  *p = (int)(cA/rSig);
  *q = (int)(cB/rSig);
}

static void check_varmapack_error(varmapack_error error) {
  if (error == VARMAPACK_OK) return;
  switch (error) {
    case VARMAPACK_INVALID_ARGUMENT:
      mexErrMsgIdAndTxt("varmapack:acvf:invalidArgument", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_ALLOCATION:
      mexErrMsgIdAndTxt("varmapack:acvf:allocation", "%s", varmapack_strerror(error));
      break;
    case VARMAPACK_NONSTATIONARY:
      mexErrMsgIdAndTxt("varmapack:acvf:nonstationary", "%s",
                        varmapack_strerror(error));
      break;
    default:
      mexErrMsgIdAndTxt("varmapack:acvf:error", "%s", varmapack_strerror(error));
      break;
  }
}
