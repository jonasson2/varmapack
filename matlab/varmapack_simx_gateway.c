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
static void check_varmapack_error(bool ok);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int p, q, s, d, r, n, M, h, Mz = 1, MX0 = 1;
  bool ok;
  double *A, *B, *C, *Sig, *z = 0, *X0, *X, *E;
  randompack_rng *rng;
  mxArray *Eout = 0;
  mwSize dim[3];
  if (nrhs != 10)
    mexErrMsgIdAndTxt("varmapack:simx:nrhs", "Wrong number of input arguments");
  if (nlhs > 2)
    mexErrMsgIdAndTxt("varmapack:simx:nlhs", "Wrong number of output arguments");
  A = get_double_array(prhs[0], "A");
  B = get_double_array(prhs[1], "B");
  C = get_double_array(prhs[2], "C");
  Sig = get_double_array(prhs[4], "Sig");
  check_dims(prhs[0], prhs[1], prhs[4], &p, &q, &r);
  n = get_int_scalar(prhs[5], "n");
  M = get_int_scalar(prhs[6], "M");
  h = get_int_scalar(prhs[8], "h");
  X0 = get_double_array(prhs[7], "X0");
  mwSize ndimsX0 = mxGetNumberOfDimensions(prhs[7]);
  const mwSize *dimsX0 = mxGetDimensions(prhs[7]);
  if (ndimsX0 > 3 || dimsX0[0] != r || dimsX0[1] != h) {
    mexErrMsgIdAndTxt("varmapack:simx:X0", "X0 must have shape r-by-h or r-by-h-by-M");
  }
  if (ndimsX0 == 3) {
    mwSize MX0size = dimsX0[2];
    if (MX0size > INT_MAX)
      mexErrMsgIdAndTxt("varmapack:simx:X0", "X0 is too large");
    MX0 = (int)MX0size;
    if (MX0 != 1 && MX0 != M)
      mexErrMsgIdAndTxt("varmapack:simx:X0", "size(X0,3) must be 1 or M");
  }
  if (mxGetNumberOfDimensions(prhs[2]) != 2 || mxGetM(prhs[2]) != r ||
      mxGetN(prhs[2]) > INT_MAX) {
    mexErrMsgIdAndTxt("varmapack:simx:C", "C must be an r-by-(d*s) matrix");
  }
  if (mxGetN(prhs[2]) == 0) {
    s = 0;
    d = 0;
  }
  else {
    mwSize ndimsZ = mxGetNumberOfDimensions(prhs[3]);
    const mwSize *dimsZ = mxGetDimensions(prhs[3]);
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
        (ndimsZ != 2 && ndimsZ != 3) || dimsZ[0] == 0 || dimsZ[1] < n ||
        dimsZ[0] > INT_MAX) {
      mexErrMsgIdAndTxt("varmapack:simx:z", "z must have shape d-by-n or d-by-n-by-M");
    }
    d = (int)dimsZ[0];
    if (mxGetN(prhs[2]) % d != 0)
      mexErrMsgIdAndTxt("varmapack:simx:C", "C must have d columns per term");
    s = (int)(mxGetN(prhs[2])/d);
    if (ndimsZ == 3 && dimsZ[2] > INT_MAX)
      mexErrMsgIdAndTxt("varmapack:simx:z", "z is too large");
    Mz = ndimsZ == 3 ? (int)dimsZ[2] : 1;
    if (Mz != 1 && Mz != M)
      mexErrMsgIdAndTxt("varmapack:simx:z", "size(z,3) must be 1 or M");
    z = mxGetPr(prhs[3]);
  }
  rng = get_rng(prhs[9]);
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
  ok = varmapack_simx(A, B, C, Sig, z, Mz, p, q, s, d, r, n, M, X0, h, MX0,
                       X, E, rng);
  check_varmapack_error(ok);
}

static randompack_rng *get_rng(const mxArray *arg) {
  uint64_t value;
  if (!mxIsUint64(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:simx:rng", "rng must be a uint64 scalar");
  value = *(uint64_t *)mxGetData(arg);
  if (value == 0)
    mexErrMsgIdAndTxt("varmapack:simx:rng", "rng handle is zero");
  return (randompack_rng *)(uintptr_t)value;
}

static double *get_double_array(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg))
    mexErrMsgIdAndTxt("varmapack:simx:argument", "%s must be a real double array", name);
  return mxGetPr(arg);
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  double x;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:simx:argument", "%s must be a numeric scalar", name);
  x = mxGetScalar(arg);
  if (x < 0 || x > INT_MAX || x != (int)x)
    mexErrMsgIdAndTxt("varmapack:simx:argument",
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
    mexErrMsgIdAndTxt("varmapack:simx:Sig", "Sig must be square");
  if (rSig > INT_MAX || cA > INT_MAX || cB > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:simx:argument", "A, B, or Sig is too large");
  if (mxGetM(A) != rSig || mxGetM(B) != rSig)
    mexErrMsgIdAndTxt("varmapack:simx:argument", "A and B must have r rows");
  if (cA % rSig != 0 || cB % rSig != 0)
    mexErrMsgIdAndTxt("varmapack:simx:argument", "A and B must have r columns per term");
  *r = (int)rSig;
  *p = (int)(cA/rSig);
  *q = (int)(cB/rSig);
}

static void check_varmapack_error(bool ok) {
  char *message;
  if (ok) return;
  message = varmapack_last_error();
  mexErrMsgIdAndTxt("varmapack:simx:error", "%s",
                    message ? message : "varmapack error");
}
