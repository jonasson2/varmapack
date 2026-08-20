#include <limits.h>
#include "mex.h"
#include "varmapack.h"

static void check_varmapack_error(bool ok);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *cov, *corr;
  int ndim, r, maxlag;
  const mwSize *dim;
  bool ok;
  if (nrhs != 1)
    mexErrMsgIdAndTxt("varmapack:cov2corr:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:cov2corr:nlhs", "Wrong number of output arguments");
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]))
    mexErrMsgIdAndTxt("varmapack:cov2corr:cov", "cov must be a real double array");
  ndim = (int)mxGetNumberOfDimensions(prhs[0]);
  dim = mxGetDimensions(prhs[0]);
  if (ndim > 3 || dim[0] == 0 || dim[1] != dim[0] ||
      (ndim == 3 && dim[2] == 0))
    mexErrMsgIdAndTxt("varmapack:cov2corr:cov",
                      "cov must be an r-by-r matrix or r-by-r-by-lag array");
  if (dim[0] > INT_MAX || (ndim == 3 && dim[2] - 1 > INT_MAX))
    mexErrMsgIdAndTxt("varmapack:cov2corr:cov", "cov is too large");
  r = (int)dim[0];
  maxlag = ndim == 2 ? 0 : (int)dim[2] - 1;
  cov = mxGetPr(prhs[0]);
  plhs[0] = mxDuplicateArray(prhs[0]);
  corr = mxGetPr(plhs[0]);
  ok = varmapack_cov2corr(cov, r, maxlag, corr);
  check_varmapack_error(ok);
}

static void check_varmapack_error(bool ok) {
  char *message;
  if (ok) return;
  message = varmapack_last_error();
  mexErrMsgIdAndTxt("varmapack:cov2corr:error", "%s",
                    message ? message : "varmapack error");
}
