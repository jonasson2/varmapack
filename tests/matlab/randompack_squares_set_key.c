#include "randompack_mex_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  (void)plhs;
  check_nargs(nlhs, nrhs, 0, 2);
  randompack_rng *rng = get_rng(prhs[0]);
  uint64_t key = get_uint64_scalar(prhs[1], "key");
  fail_if_false(randompack_squares_set_key(key, rng), rng,
                "randompack_squares_set_key failed");
}
