#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)plhs;
  if (nrhs != 2 && nrhs != 3)
    mexErrMsgIdAndTxt("randompack:nrhs", "Wrong number of input arguments");
  if (nlhs > 0)
    mexErrMsgIdAndTxt("randompack:nlhs", "Wrong number of output arguments");
  randompack_rng *rng = get_rng(prhs[0]);
  int seed = get_int_scalar(prhs[1], "seed");
  uint32_t *spawn_key = 0;
  int nkey = 0;
  if (nrhs == 3) {
    if (!mxIsUint32(prhs[2]) || mxIsComplex(prhs[2]))
      mexErrMsgIdAndTxt("randompack:seed", "spawn_key must be a uint32 array");
    spawn_key = (uint32_t *)mxGetData(prhs[2]);
    nkey = (int)mxGetNumberOfElements(prhs[2]);
  }
  fail_if_false(randompack_seed(seed, spawn_key, nkey, rng), rng,
                "randompack_seed failed");
}
