#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include "mex.h"
#include "randompack.h"

static randompack_rng *get_rng(const mxArray *arg);
static int get_seed(const mxArray *arg);
static void fail_rng(randompack_rng *rng, char *fallback);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  char command[16];
  char error[256];
  char *message;
  randompack_rng *rng;
  if (nrhs < 1 || !mxIsChar(prhs[0]) ||
      mxGetString(prhs[0], command, sizeof(command)) != 0)
    mexErrMsgIdAndTxt("varmapack:rng:command",
                      "First argument must be a command string");
  if (strcmp(command, "create") == 0) {
    if (nrhs != 1 || nlhs != 1)
      mexErrMsgIdAndTxt("varmapack:rng:create",
                        "create requires one output and no arguments");
    rng = randompack_create(0);
    if (!rng)
      mexErrMsgIdAndTxt("varmapack:rng:create", "Could not create Randompack RNG");
    message = randompack_last_error(rng);
    if (message) {
      snprintf(error, sizeof(error), "%s", message);
      randompack_free(rng);
      mexErrMsgIdAndTxt("varmapack:rng:create", "%s", error);
    }
    mexLock();
    plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *(uint64_t *)mxGetData(plhs[0]) = (uint64_t)(uintptr_t)rng;
    return;
  }
  if (strcmp(command, "seed") == 0) {
    if (nrhs != 3 || nlhs != 0)
      mexErrMsgIdAndTxt("varmapack:rng:seed", "seed requires an RNG and a seed");
    rng = get_rng(prhs[1]);
    if (!randompack_seed(get_seed(prhs[2]), 0, 0, rng))
      fail_rng(rng, "Could not seed Randompack RNG");
    return;
  }
  if (strcmp(command, "free") == 0) {
    if (nrhs != 2 || nlhs != 0)
      mexErrMsgIdAndTxt("varmapack:rng:free", "free requires an RNG and no output");
    randompack_free(get_rng(prhs[1]));
    mexUnlock();
    return;
  }
  mexErrMsgIdAndTxt("varmapack:rng:command", "Unknown command");
}

static randompack_rng *get_rng(const mxArray *arg) {
  uint64_t value;
  if (!mxIsUint64(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:rng:handle", "RNG handle must be a uint64 scalar");
  value = *(uint64_t *)mxGetData(arg);
  if (value == 0)
    mexErrMsgIdAndTxt("varmapack:rng:handle", "RNG handle is zero");
  return (randompack_rng *)(uintptr_t)value;
}

static int get_seed(const mxArray *arg) {
  double seed;
  if (!mxIsNumeric(arg) || mxIsComplex(arg) || mxGetNumberOfElements(arg) != 1)
    mexErrMsgIdAndTxt("varmapack:rng:seed", "Seed must be an integer scalar");
  seed = mxGetScalar(arg);
  if (seed < INT_MIN || seed > INT_MAX || seed != (int)seed)
    mexErrMsgIdAndTxt("varmapack:rng:seed", "Seed must be an integer in the C int range");
  return (int)seed;
}

static void fail_rng(randompack_rng *rng, char *fallback) {
  char *message = randompack_last_error(rng);
  mexErrMsgIdAndTxt("varmapack:rng:error", "%s", message ? message : fallback);
}
