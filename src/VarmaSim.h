// VarmaSim — simulate spin-up free AR, VAR, ARMA or VARMA time series
//
//   Model:
//       x(t) - μ = A1(x(t-1) - μ) + … + Ap(x(t-p) - μ) + y(t),
//       y(t) = ε(t) + B1ε(t-1) + … + Bqε(t-q),
//       ε(t) ~ N(0, Σ), independent over time.
//
//   Special cases:
//       r=1, q>0   → ARMA(p,q)
//       r>1, q=0   → VAR(p)
//       r=1, q=0   → AR(p)
//
//   Matrices use Fortran (column-major) storage.  When X0 is NULL, the
//   model must be stationary; otherwise non-stationary models are allowed.
//   Each simulated series has the correct distribution from term one
//   (it is “spin-up free”).
//
//   See VarmaSim.c for full discussion, algorithm notes, and references
//   and parameter descriptions below for further details

#ifndef VARMASIM_H
#define VARMASIM_H
#ifdef __cplusplus
extern "C" {
#endif

#include "RandomNumbers.h"
#include <stdbool.h>

void VarmaSim(
  double A[],   // in      r×r×p, autoregressive parameter matrices
  double B[],   // in      r×r×q, moving average parameter matrices
  double Sig[], // in      r×r, covariance of the shock terms ε(t)
  double mu[],  // in      r-vector, mean of x(t); NULL for zero-mean
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving-average terms
  int r,        // in      dimension of each x(t)
  int n,        // in      length of each generated series
  int M,        // in      number of replicates to generate
  double X0[],  // in      r×nX0 with nX0 ≥ max(p,q); optional starter or NULL
  int nX0,      // in      length of x0
  RandRng *rng, // in/out  random number generator
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series (or NULL to skip)
  bool *ok      // out     1 if stationary (or X0 non-NULL), else 0
  );

#ifdef __cplusplus
}
#endif
#endif /* VARMASIM_H */
