#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nlhs > 1)
    mexErrMsgIdAndTxt("randompack:nlhs", "Wrong number of output arguments");
  if (nrhs > 1)
    mexErrMsgIdAndTxt("randompack:nrhs", "Wrong number of input arguments");
  char engine[64];
  char *enginep = 0;
  if (nrhs == 1) {
    if (mxGetString(prhs[0], engine, sizeof(engine)) != 0)
      mexErrMsgIdAndTxt("randompack:create", "Engine name is too long");
    if (engine[0] != 0)
      enginep = engine;
  }
  randompack_rng *rng = randompack_create(enginep);
  if (!rng)
    mexErrMsgIdAndTxt("randompack:create", "Could not create RNG");
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *(uint64_t *)mxGetData(plhs[0]) = (uint64_t)(uintptr_t)rng;
}
