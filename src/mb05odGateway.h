#ifndef MB05ODGATEWAY_H
#define MB05ODGATEWAY_H

#include "varmapack_config.h"

typedef struct {
  int mdig;  // Conservative estimate of accurate decimal digits in the 1-norm.
  int idig;  // Estimate of accurate digits at the routine's 95% confidence level.
  int iwarn; // 0: none; 1/2: possible/severe inaccuracy; 3: balancing not used.
} mb05od_accuracy;

// Compute E = exp(delta*A), using MB05OD with balancing and Pade order 9.
// Matrices are n*n, column-major. n=0 is allowed; otherwise A and E are required.
// Inputs must be finite, including delta*A. A is preserved unless E == A;
// other overlapping input/output storage is unsupported. E is unchanged on failure.
// Returns false for invalid input, allocation, SLICOT INFO, or nonfinite output.
// INFO=0 with IWARN!=0 still returns true: inspect accuracy before using the result.
// Optional accuracy is zeroed on entry and filled when SLICOT has been called.
// Errors use varmapack_last_error(). Workspaces are local to each call.
HIDDEN bool mb05odGateway(double A[], double delta, int n, double E[],
                          mb05od_accuracy *accuracy);

// Variants that construct a matrix already scaled for exponentiation can
// disable balancing explicitly.
HIDDEN bool mb05odExponential(double A[], double delta, int n, double E[],
                              bool balance, mb05od_accuracy *accuracy);

#endif
