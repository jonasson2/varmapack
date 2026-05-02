#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  check_nargs(nlhs, nrhs, 1, 4);
  randompack_rng *rng = get_rng(prhs[0]);
  size_t len = get_size_scalar(prhs[1], "len");
  double mu = mxGetScalar(prhs[2]);
  double sigma = mxGetScalar(prhs[3]);
  plhs[0] = mxCreateDoubleMatrix((mwSize)len, 1, mxREAL);
  fail_if_false(randompack_normal(mxGetPr(plhs[0]), len, mu, sigma, rng), rng,
                "randompack_normal failed");
}
