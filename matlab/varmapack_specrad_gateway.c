#include <math.h>
#include <limits.h>
#include "mex.h"
#include "varmapack.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  int r, p;
  mwSize m, n;
  double rho;
  if (nrhs != 1)
    mexErrMsgIdAndTxt("varmapack:specrad:nrhs", "Wrong number of input arguments");
  if (nlhs > 1)
    mexErrMsgIdAndTxt("varmapack:specrad:nlhs", "Wrong number of output arguments");
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
    mexErrMsgIdAndTxt("varmapack:specrad:A", "A must be a real double matrix");
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  if (m == 0) {
    plhs[0] = mxCreateDoubleScalar(0);
    return;
  }
  if (m > INT_MAX || n > INT_MAX)
    mexErrMsgIdAndTxt("varmapack:specrad:A", "A is too large");
  if (n % m != 0)
    mexErrMsgIdAndTxt("varmapack:specrad:A", "size(A,2) must be a multiple of size(A,1)");
  r = (int)m;
  p = (int)(n/m);
  rho = varmapack_specrad(mxGetPr(prhs[0]), r, p);
  if (isnan(rho))
    mexErrMsgIdAndTxt("varmapack:specrad:error", "invalid argument");
  plhs[0] = mxCreateDoubleScalar(rho);
}
