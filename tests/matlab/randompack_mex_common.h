#ifndef RANDOMPACK_MEX_COMMON_H
#define RANDOMPACK_MEX_COMMON_H

#include <stdint.h>
#include <limits.h>
#include <string.h>
#include "mex.h"
#include "randompack.h"

static randompack_rng *get_rng(const mxArray *arg) {
  if (!mxIsUint64(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("randompack:handle", "RNG handle must be a uint64 scalar");
  uint64_t value = *(uint64_t *)mxGetData(arg);
  if (value == 0)
    mexErrMsgIdAndTxt("randompack:handle", "RNG handle is zero");
  return (randompack_rng *)(uintptr_t)value;
}

static int get_int_scalar(const mxArray *arg, const char *name) {
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("randompack:argument", "%s must be a numeric scalar", name);
  double x = mxGetScalar(arg);
  if (x < (double)INT32_MIN || x > (double)INT32_MAX)
    mexErrMsgIdAndTxt("randompack:argument", "%s is outside int range", name);
  return (int)x;
}

static size_t get_size_scalar(const mxArray *arg, const char *name) {
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("randompack:argument", "%s must be a numeric scalar", name);
  double x = mxGetScalar(arg);
  if (x < 0)
    mexErrMsgIdAndTxt("randompack:argument", "%s must be nonnegative", name);
  return (size_t)x;
}

static uint64_t get_uint64_scalar(const mxArray *arg, const char *name) {
  if (mxIsUint64(arg) && !mxIsComplex(arg) && mxGetNumberOfElements(arg) == 1)
    return *(uint64_t *)mxGetData(arg);
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("randompack:argument", "%s must be a numeric scalar", name);
  double x = mxGetScalar(arg);
  if (x < 0)
    mexErrMsgIdAndTxt("randompack:argument", "%s must be nonnegative", name);
  return (uint64_t)x;
}

static double *get_real_double_matrix(const mxArray *arg, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg))
    mexErrMsgIdAndTxt("randompack:argument", "%s must be a real double array", name);
  return mxGetPr(arg);
}

static void check_nargs(int nlhs, int nrhs, int nlhs_max, int nrhs_need) {
  if (nrhs != nrhs_need)
    mexErrMsgIdAndTxt("randompack:nrhs", "Wrong number of input arguments");
  if (nlhs > nlhs_max)
    mexErrMsgIdAndTxt("randompack:nlhs", "Wrong number of output arguments");
}

static void fail_if_false(bool ok, randompack_rng *rng, const char *msg) {
  if (ok) return;
  const char *err = rng ? randompack_last_error(rng) : 0;
  mexErrMsgIdAndTxt("randompack:call", "%s%s%s", msg, err ? ": " : "",
                    err ? err : "");
}

#endif
