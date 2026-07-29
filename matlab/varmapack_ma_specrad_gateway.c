#include <limits.h>
#include <math.h>
#include "mex.h"
#include "varmapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int r, q;
  mwSize m, n;
  double rho;
  if (nrhs != 1)
    mexErrMsgIdAndTxt("varmapack:ma_specrad:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:ma_specrad:nlhs", "Wrong number of output arguments");
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgIdAndTxt("varmapack:ma_specrad:B", "B must be a real double matrix");
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if (m == 0) {
    plhs[0] = mxCreateDoubleScalar(0);
    return;
  }
  if (m > INT_MAX || n > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:ma_specrad:B", "B is too large");
  if (n % m != 0)
    mexErrMsgIdAndTxt("varmapack:ma_specrad:B",
                      "size(B,2) must be a multiple of size(B,1)");
  r = (int)m;
  q = (int)(n/m);
  rho = varmapack_ma_specrad(mxGetPr(prhs[0]), r, q);
  if (isnan(rho)) {
    char *message = varmapack_last_error();
    mexErrMsgIdAndTxt("varmapack:ma_specrad:error", "%s",
                      message ? message : "varmapack error");
  }
  plhs[0] = mxCreateDoubleScalar(rho);
}
