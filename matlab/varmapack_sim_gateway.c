#include <stdint.h>
#include <limits.h>
#include "mex.h"
#include "randompack.h"
#include "varmapack.h"

static randompack_rng *get_rng(const mxArray *arg);
static double *get_double_array(const mxArray *arg, const char *name);
static int get_int_scalar(const mxArray *arg, const char *name);
static void check_dims(const mxArray *A, const mxArray *B, const mxArray *Sig,
                       int *p, int *q, int *r);
static void check_varmapack_error(varmapack_error error);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int p, q, r, n, M, nX0 = 0;
  varmapack_error error;
  double *A, *B, *Sig, *mu = 0, *X0 = 0, *X, *E;
  randompack_rng *rng;
  mxArray *Eout = 0;
  mwSize dim[3];
  if (nrhs != 8)
    mexErrMsgIdAndTxt("varmapack:sim:nrhs", "Wrong number of input arguments");
  if (nlhs > 2)
    mexErrMsgIdAndTxt("varmapack:sim:nlhs", "Wrong number of output arguments");
  A = get_double_array(prhs[0], "A");
  B = get_double_array(prhs[1], "B");
  Sig = get_double_array(prhs[2], "Sig");
  check_dims(prhs[0], prhs[1], prhs[2], &p, &q, &r);
  n = get_int_scalar(prhs[3], "n");
  if (!mxIsEmpty(prhs[4])) {
    mu = get_double_array(prhs[4], "mu");
    if ((int)mxGetNumberOfElements(prhs[4]) != r)
      mexErrMsgIdAndTxt("varmapack:sim:mu", "mu must be empty or length r");
  }
  M = get_int_scalar(prhs[5], "M");
  if (!mxIsEmpty(prhs[6])) {
    X0 = get_double_array(prhs[6], "X0");
    if ((int)mxGetM(prhs[6]) != r)
      mexErrMsgIdAndTxt("varmapack:sim:X0", "X0 must have r rows");
    nX0 = (int)mxGetN(prhs[6]);
  }
  rng = get_rng(prhs[7]);
  dim[0] = r;
  dim[1] = n;
  dim[2] = M;
  plhs[0] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
  if (nlhs >= 2) {
    plhs[1] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
    Eout = plhs[1];
  }
  X = mxGetPr(plhs[0]);
  E = Eout ? mxGetPr(Eout) : 0;
  error = varmapack_sim(A, B, Sig, mu, p, q, r, n, M, X0, nX0, rng, X, E);
  check_varmapack_error(error);
}

static randompack_rng *get_rng(const mxArray *arg) {
  uint64_t value;
  if (!mxIsUint64(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:sim:rng", "rng must be a uint64 scalar");
  value = *(uint64_t *)mxGetData(arg);
  if (value == 0)
    mexErrMsgIdAndTxt("varmapack:sim:rng", "rng handle is zero");
  return (randompack_rng *)(uintptr_t)value;
}

static double *get_double_array(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg))
    mexErrMsgIdAndTxt("varmapack:sim:argument", "%s must be a real double array",
                      name);
  return mxGetPr(arg);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:sim:argument", "%s must be a numeric scalar",
                      name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:sim:argument",
                      "%s must be a nonnegative integer", name);
  return (int)x;
}

static void check_dims(const mxArray *A, const mxArray *B, const mxArray *Sig,
                       int *p, int *q, int *r) {
  mwSize rSig, cSig, cA, cB;
  rSig = mxGetM(Sig);
  cSig = mxGetN(Sig);
  if (rSig == 0 || cSig != rSig)
    mexErrMsgIdAndTxt("varmapack:sim:Sig", "Sig must be square");
  if (rSig > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:sim:Sig", "Sig is too large");
  *r = (int)rSig;
  if (mxGetM(A) != rSig)
    mexErrMsgIdAndTxt("varmapack:sim:A", "A must have r rows");
  if (mxGetM(B) != rSig)
    mexErrMsgIdAndTxt("varmapack:sim:B", "B must have r rows");
  cA = mxGetN(A);
  cB = mxGetN(B);
  if (cA % rSig != 0)
    mexErrMsgIdAndTxt("varmapack:sim:A", "size(A,2) must be a multiple of r");
  if (cB % rSig != 0)
    mexErrMsgIdAndTxt("varmapack:sim:B", "size(B,2) must be a multiple of r");
  if (cA/rSig > INT_MAX || cB/rSig > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:sim:argument", "A or B is too large");
  *p = (int)(cA/rSig);
  *q = (int)(cB/rSig);
}

static void check_varmapack_error(varmapack_error error) {
  if (error == VARMAPACK_OK) return;
  switch (error) {
    case VARMAPACK_INVALID_ARGUMENT:
      mexErrMsgIdAndTxt("varmapack:sim:invalidArgument", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_ALLOCATION:
      mexErrMsgIdAndTxt("varmapack:sim:allocation", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_NONSTATIONARY:
      mexErrMsgIdAndTxt("varmapack:sim:nonstationary", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_NOT_POSITIVE_SEMIDEFINITE:
      mexErrMsgIdAndTxt("varmapack:sim:notPositiveSemidefinite", "%s",
                        varmapack_strerror(error));
      break;
    case VARMAPACK_INTERNAL:
      mexErrMsgIdAndTxt("varmapack:sim:internal", "%s",
                        varmapack_strerror(error));
      break;
    default:
      mexErrMsgIdAndTxt("varmapack:sim:error", "%s", varmapack_strerror(error));
      break;
  }
}
