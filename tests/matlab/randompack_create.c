#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  check_nargs(nlhs, nrhs, 1, 1);
  char engine[64];
  if (mxGetString(prhs[0], engine, sizeof(engine)) != 0)
    mexErrMsgIdAndTxt("randompack:create", "Engine name is too long");
  randompack_rng *rng = randompack_create(engine);
  if (!rng)
    mexErrMsgIdAndTxt("randompack:create", "Could not create RNG");
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *(uint64_t *)mxGetData(plhs[0]) = (uint64_t)(uintptr_t)rng;
}
