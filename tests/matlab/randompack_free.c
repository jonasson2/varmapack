#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)plhs;
  check_nargs(nlhs, nrhs, 0, 1);
  randompack_free(get_rng(prhs[0]));
}
