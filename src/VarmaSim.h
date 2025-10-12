#include "random.h"
void VarmaSim( //  Generate simulated AR, VAR, ARMA or VARMA time series
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  double mu[],  // in   r-vector, mean of x(t); null for zero-mean
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  int n,        // in   length of each generated series
  int M,        // in   number of replicates to generate
  double X0[],  // in   r×h with h = max(p,q), series starter, or null. If it
  //            //      is not null then each generated series starts with X0.
  rand_rng *rng,// in   random number generator
  double X[],   // out  r×n×M, the generated series
  double eps[], // out  r×n×M, the shock series (or null to not return them)
  int *ok);     // out  1 if the series is stationary or X0 is non-null, else 0
  // See comments in VarmaSim.c
void Testcase ( // Make test data for VARMA simulation or likelihood calculation
  double A[],     // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],     // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],   // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,     // out     string w. defined length >= 12, name of testcase
  int *p,         // in/out  Number of autoregressive terms
  int *q,         // in/out  Number of moving avg. terms
  int *r,         // in/out  dimension of each x(t)
  int icase,      // in      index of named testcase to create or 0 to use p,q,r
  rand_rng *rng); // int     Random number generator
  // Notes: When icase = 0 a random model with dimensions p, q, r is created.
  // When 1 <= icase <= 12 and A, B and Sig are null, p, q and r return
  // dimensions of one of 12 predefined testcases, and when they are non-null
  // they return the data for these cases. See more detailes in Testcase.c.
