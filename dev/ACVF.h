// ACVF — compute the theoretical autocovariance function (ACVF) of the stationary
// time-series process {x_t} defined by a VARMA(p,q) model, using the vector Yule–Walker
// (VYW) equations.
//
// Inputs:
//   A     : r×r×p array of autoregressive coefficient matrices
//   B     : r×r×q array of moving-average coefficient matrices
//   Sig   : r×r innovation (shock) covariance matrix
//   p,q,r : model orders and system dimension
//
// Output:
//   S   : r×r×(maxlag+1) array containing the first (maxlag+1) matrices of the
//         theoretical autocovariance function (ACVF) of the process {x_t},
//         S[:,:,k] = Cov(x_t, x_{t−k}),  k = 0,…,maxlag
//   ok  : 1 on success, 0 on failure (non-
//
// Returns true on success, false on failure.

#include <stdbool.h>

bool ACVF(double A[], double B[], double Sig[], int p, int q, int r, double Gamma[],
	  int maxlag);

