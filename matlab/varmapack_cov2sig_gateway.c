#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static int coefficient_order(const mxArray *arg, mwSize r, const char *name) {
  if (!mxIsDouble(arg) || mxIsComplex(arg) || mxIsSparse(arg) ||
      mxGetNumberOfDimensions(arg) != 2 || mxGetM(arg) != r ||
      mxGetN(arg) % r != 0 || mxGetN(arg)/r > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:cov2sig:argument",
                      "%s must be a real full double r-by-(r*order) matrix", name);
  return mxGetN(arg)/r;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 3)
    mexErrMsgIdAndTxt("varmapack:cov2sig:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:cov2sig:nlhs", "Wrong number of output arguments");
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) ||
      mxGetNumberOfDimensions(prhs[2]) != 2 || mxGetM(prhs[2]) == 0 ||
      mxGetN(prhs[2]) != mxGetM(prhs[2]) || mxGetM(prhs[2]) > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:cov2sig:cov",
                      "cov must be a non-empty real full double square matrix");
  int r = mxGetM(prhs[2]);
  int p = coefficient_order(prhs[0], r, "A");
  int q = coefficient_order(prhs[1], r, "B");
  mxArray *Sig = mxCreateDoubleMatrix(r, r, mxREAL);
  bool ok = varmapack_cov2sig(mxGetPr(prhs[0]), mxGetPr(prhs[1]),
                             mxGetPr(prhs[2]), p, q, r, mxGetPr(Sig));
  if (!ok) {
    char *message = varmapack_last_error();
    mxDestroyArray(Sig);
    mexErrMsgIdAndTxt("varmapack:cov2sig:error", "%s",
                      message ? message : "varmapack error");
  }
  if (nlhs == 0) mxDestroyArray(Sig);
  else plhs[0] = Sig;
}
