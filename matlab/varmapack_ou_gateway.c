#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "varmapack.h"

static double *real_array(const mxArray *arg) {
  if (!mxIsDouble(arg) || mxIsComplex(arg) || mxIsSparse(arg))
    mexErrMsgIdAndTxt("varmapack:ou:argument", "Expected a full real double array");
  return mxGetPr(arg);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if ((nrhs != 4 && nrhs != 7) || nlhs > 1 || !mxIsLogicalScalar(prhs[3]))
    mexErrMsgIdAndTxt("varmapack:ou:arguments", "Wrong number of arguments");
  double *parameter = real_array(prhs[0]);
  double *cov = real_array(prhs[1]);
  double *mu = real_array(prhs[2]);
  bool matrix_drift = mxIsLogicalScalarTrue(prhs[3]);
  mwSize r = mxGetM(prhs[1]);
  bool parameter_ok = matrix_drift ?
    mxGetNumberOfDimensions(prhs[0]) == 2 && mxGetM(prhs[0]) == r &&
      mxGetN(prhs[0]) == r :
    mxGetNumberOfElements(prhs[0]) == r;
  if (r == 0 || r > INT_MAX || mxGetNumberOfDimensions(prhs[1]) != 2 ||
      mxGetN(prhs[1]) != r || !parameter_ok ||
      mxGetNumberOfElements(prhs[2]) != r)
    mexErrMsgIdAndTxt("varmapack:ou:dimensions", "Invalid OU dimensions");
  mxArray *out = 0;
  double *t = 0;
  int n = 0, M = 0;
  randompack_rng *rng = 0;
  if (nrhs == 7) {
    t = real_array(prhs[4]);
    mwSize nt = mxGetNumberOfElements(prhs[4]);
    if (nt < 1 || nt > INT_MAX || !mxIsNumeric(prhs[5]) ||
        mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1)
      mexErrMsgIdAndTxt("varmapack:ou:dimensions", "Invalid times or nrep");
    double m = mxGetScalar(prhs[5]);
    if (!isfinite(m) || m < 1 || m > INT_MAX || m != floor(m))
      mexErrMsgIdAndTxt("varmapack:ou:nrep", "nrep must be a positive integer");
    n = nt;
    M = m;
    if (!mxIsUint64(prhs[6]) || mxGetNumberOfElements(prhs[6]) != 1)
      mexErrMsgIdAndTxt("varmapack:ou:rng", "Invalid RNG handle");
    rng = (randompack_rng *)(uintptr_t)(*(uint64_t *)mxGetData(prhs[6]));
    mwSize dims[] = {r, nt, M};
    out = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
  }
  varmapack_ou_model *model = matrix_drift ?
    varmapack_ou_create_drift(parameter, cov, mu, r) :
    varmapack_ou_create(parameter, cov, mu, r);
  bool ok = model != 0;
  if (ok && nrhs == 7) ok = varmapack_ou_sim(model, t, n, M, mxGetPr(out), rng);
  varmapack_ou_free(model);
  if (!ok) {
    if (out) mxDestroyArray(out);
    mexErrMsgIdAndTxt("varmapack:ou:error", "%s", varmapack_last_error());
  }
  if (!out) out = mxCreateDoubleMatrix(0, 0, mxREAL);
  if (nlhs) plhs[0] = out;
  else mxDestroyArray(out);
}
