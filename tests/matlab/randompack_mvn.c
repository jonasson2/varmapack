#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  check_nargs(nlhs, nrhs, 1, 4);
  randompack_rng *rng = get_rng(prhs[0]);
  double *Sig = get_real_double_matrix(prhs[1], "Sig");
  mwSize r = mxGetM(prhs[1]);
  if (r == 0 || mxGetN(prhs[1]) != r)
    mexErrMsgIdAndTxt("randompack:mvn", "Sig must be square");
  size_t len = get_size_scalar(prhs[2], "len");
  double *mu = get_real_double_matrix(prhs[3], "mu");
  if (mxGetNumberOfElements(prhs[3]) != r)
    mexErrMsgIdAndTxt("randompack:mvn", "mu must have length size(Sig,1)");
  plhs[0] = mxCreateDoubleMatrix(r, (mwSize)len, mxREAL);
  fail_if_false(randompack_mvn("T", mu, Sig, (int)r, len, mxGetPr(plhs[0]),
                               (int)r, 0, rng), rng, "randompack_mvn failed");
}
