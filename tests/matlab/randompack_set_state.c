#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)plhs;
  check_nargs(nlhs, nrhs, 0, 2);
  randompack_rng *rng = get_rng(prhs[0]);
  int n = (int)mxGetNumberOfElements(prhs[1]);
  uint64_t *state = 0;
  if (!mxIsUint64(prhs[1]) || mxIsComplex(prhs[1]))
    mexErrMsgIdAndTxt("randompack:set_state", "state must be a uint64 array");
  state = (uint64_t *)mxGetData(prhs[1]);
  fail_if_false(randompack_set_state(state, n, rng), rng, "randompack_set_state failed");
}
