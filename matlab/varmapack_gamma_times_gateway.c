#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "varmapack.h"

static double scalar(const mxArray *arg) {
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:gamma_times:argument", "Expected a real scalar");
  return mxGetScalar(arg);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 4 || nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:gamma_times:arguments", "Wrong number of arguments");
  double n = scalar(prhs[0]), shape = scalar(prhs[1]), gap = scalar(prhs[2]);
  if (!isfinite(n) || n < 1 || n > INT_MAX || n != floor(n))
    mexErrMsgIdAndTxt("varmapack:gamma_times:n", "n must be a positive integer");
  if (!mxIsUint64(prhs[3]) || mxGetNumberOfElements(prhs[3]) != 1)
    mexErrMsgIdAndTxt("varmapack:gamma_times:rng", "Invalid RNG handle");
  randompack_rng *rng =
    (randompack_rng *)(uintptr_t)(*(uint64_t *)mxGetData(prhs[3]));
  mxArray *out = mxCreateDoubleMatrix(n, 1, mxREAL);
  if (!varmapack_gamma_times(n, shape, gap, mxGetPr(out), rng)) {
    mxDestroyArray(out);
    mexErrMsgIdAndTxt("varmapack:gamma_times:error", "%s", varmapack_last_error());
  }
  if (nlhs) plhs[0] = out;
  else mxDestroyArray(out);
}
