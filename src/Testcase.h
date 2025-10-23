// Make test data for VARMA simulation (or likelihood calculation)
#ifndef TESTCASE_H
#define TESTCASE_H
#ifdef __cplusplus
extern "C" {
#endif

#include "RandomNumbers.h"

bool Testcase (
  double A[],     // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],     // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],   // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,     // out     string w. defined length >= 12, name of testcase
  int *p,         // in/out  Number of autoregressive terms
  int *q,         // in/out  Number of moving avg. terms
  int *r,         // in/out  dimension of each x(t)
  int *icase,     // in      index of named testcase to create or 0 to use p,q,r
  RandRng *rng,   // int     Random number generator
  FILE *errfp);  
  // Notes: When icase = 0 a random model with dimensions p, q, r is created.
  // When 1 <= icase <= 12 and A, B and Sig are null, p, q and r return
  // dimensions of one of 12 predefined testcases, and when they are non-null
  // they return the data for these cases. See more detailes in Testcase.c.

#ifdef __cplusplus
}
#endif
#endif
