#include <math.h>
#include <stdint.h>
#include <string.h>
#include "mex.h"
#include "slepian_internal.h"

static double *real_array(const mxArray *arg) {
  if (!mxIsDouble(arg) || mxIsComplex(arg) || mxIsSparse(arg))
    mexErrMsgIdAndTxt("varmapack:slepian:argument", "Expected a full real double array");
  return mxGetPr(arg);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if ((nrhs != 3 && nrhs != 6) || nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:slepian:arguments", "Wrong number of arguments");
  double *theta = real_array(prhs[0]);
  double *cov = real_array(prhs[1]);
  double *mu = real_array(prhs[2]);
  mwSize r = mxGetM(prhs[0]);
  if (r == 0 || r > INT_MAX || mxGetNumberOfDimensions(prhs[0]) != 2 ||
      mxGetN(prhs[0]) != r || mxGetNumberOfDimensions(prhs[1]) != 2 ||
      mxGetM(prhs[1]) != r || mxGetN(prhs[1]) != r ||
      mxGetNumberOfElements(prhs[2]) != r)
    mexErrMsgIdAndTxt("varmapack:slepian:dimensions", "Invalid Slepian dimensions");
  mxArray *out = 0;
  double *t = 0;
  int n = 0, M = 0;
  randompack_rng *rng = 0;
  if (nrhs == 6) {
    t = real_array(prhs[3]);
    mwSize nt = mxGetNumberOfElements(prhs[3]);
    if (nt < 1 || nt > INT_MAX || !mxIsNumeric(prhs[4]) ||
        mxIsComplex(prhs[4]) || mxGetNumberOfElements(prhs[4]) != 1)
      mexErrMsgIdAndTxt("varmapack:slepian:dimensions", "Invalid times or nrep");
    double m = mxGetScalar(prhs[4]);
    if (!isfinite(m) || m < 1 || m > INT_MAX || m != floor(m))
      mexErrMsgIdAndTxt("varmapack:slepian:nrep", "nrep must be a positive integer");
    n = nt;
    M = m;
    if (!mxIsUint64(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1)
      mexErrMsgIdAndTxt("varmapack:slepian:rng", "Invalid RNG handle");
    rng = (randompack_rng *)(uintptr_t)(*(uint64_t *)mxGetData(prhs[5]));
    mwSize dims[] = {r, nt, M};
    out = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  }
  bool ok;
  if (out) {
    // In simulation mode the second matrix is the prepared innovation factor.
    varmapack_slepian_model model = {r, theta, mu, cov};
    ok = varmapack_slepian_sim(&model, t, n, M, mxGetPr(out), rng);
  }
  else {
    out = mxCreateDoubleMatrix(r, r, mxREAL);
    varmapack_slepian_model *model = varmapack_slepian_create(theta, cov, mu, r);
    ok = model != 0;
    if (ok) memcpy(mxGetPr(out), model->factor, r*r*sizeof(double));
    varmapack_slepian_free(model);
  }
  if (!ok) {
    if (out) mxDestroyArray(out);
    mexErrMsgIdAndTxt("varmapack:slepian:error", "%s", varmapack_last_error());
  }
  if (nlhs) plhs[0] = out;
  else mxDestroyArray(out);
}
