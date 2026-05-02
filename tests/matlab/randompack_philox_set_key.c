#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)plhs;
  check_nargs(nlhs, nrhs, 0, 2);
  randompack_rng *rng = get_rng(prhs[0]);
  if (!mxIsUint64(prhs[1]) || mxIsComplex(prhs[1]) ||
      mxGetNumberOfElements(prhs[1]) != 2)
    mexErrMsgIdAndTxt("randompack:philox_set_key", "key must be uint64 with length 2");
  uint64_t *key = (uint64_t *)mxGetData(prhs[1]);
  fail_if_false(randompack_philox_set_key(key, rng), rng,
                "randompack_philox_set_key failed");
}
