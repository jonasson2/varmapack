#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static double *get_double_array(const mxArray *arg, const char *name);
static int get_int_scalar(const mxArray *arg, const char *name);
static void check_dims(const mxArray *A, const mxArray *B, int *p, int *q, int *r);
static void check_varmapack_error(bool ok);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int p, q, r, maxlag;
  double *A, *B, *Psi;
  mwSize dim[3];
  bool ok;
  if (nrhs != 3)
    mexErrMsgIdAndTxt("varmapack:psi:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:psi:nlhs", "Wrong number of output arguments");
  A = get_double_array(prhs[0], "A");
  B = get_double_array(prhs[1], "B");
  maxlag = get_int_scalar(prhs[2], "maxlag");
  check_dims(prhs[0], prhs[1], &p, &q, &r);
  dim[0] = r;
  dim[1] = r;
  dim[2] = maxlag + 1;
  plhs[0] = mxCreateNumericArray(3, dim, mxDOUBLE_CLASS, mxREAL);
  Psi = mxGetPr(plhs[0]);
  ok = varmapack_psi(A, B, p, q, r, maxlag, Psi);
  check_varmapack_error(ok);
}

static double *get_double_array(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg))
    mexErrMsgIdAndTxt("varmapack:psi:argument", "%s must be a real double array", name);
  return mxGetPr(arg);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:psi:argument", "%s must be a numeric scalar", name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:psi:argument", "%s must be a nonnegative integer", name);
  return (int)x;
}

static void check_dims(const mxArray *A, const mxArray *B, int *p, int *q, int *r) {
  mwSize mA = mxGetM(A), mB = mxGetM(B);
  mwSize nA = mxGetN(A), nB = mxGetN(B);
  mwSize m = mA > 0 ? mA : mB;
  if (m == 0) m = 1;
  if (m > INT_MAX || nA > INT_MAX || nB > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:psi:argument", "A or B is too large");
  if ((mA != 0 && mA != m) || (mB != 0 && mB != m))
    mexErrMsgIdAndTxt("varmapack:psi:argument", "A and B must have the same row count");
  if (nA % m != 0 || nB % m != 0)
    mexErrMsgIdAndTxt("varmapack:psi:argument", "A and B must have r columns per term");
  *r = (int)m;
  *p = (int)(nA/m);
  *q = (int)(nB/m);
}

static void check_varmapack_error(bool ok) {
  char *message;
  if (ok) return;
  message = varmapack_last_error();
  mexErrMsgIdAndTxt("varmapack:psi:error", "%s", message ? message : "varmapack error");
}
