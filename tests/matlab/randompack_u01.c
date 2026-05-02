#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  check_nargs(nlhs, nrhs, 1, 2);
  randompack_rng *rng = get_rng(prhs[0]);
  size_t len = get_size_scalar(prhs[1], "len");
  plhs[0] = mxCreateDoubleMatrix((mwSize)len, 1, mxREAL);
  fail_if_false(randompack_u01(mxGetPr(plhs[0]), len, rng), rng,
                "randompack_u01 failed");
}
