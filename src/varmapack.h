/// Main header file for varmapack. 

#ifndef VARMAPACK_H
#define VARMAPACK_H
#include "varmapack_config.h"
#ifdef __cplusplus
extern "C" {
#endif

#include "randompack.h"
#include <stdbool.h>
#include <stdio.h>

typedef enum {
  VARMAPACK_OK = 0, VARMAPACK_INVALID_ARGUMENT, VARMAPACK_DIMENSION, VARMAPACK_ALLOCATION,
  VARMAPACK_UNKNOWN_TESTCASE, VARMAPACK_NONSTATIONARY,
  VARMAPACK_NONSTATIONARY_MA, VARMAPACK_SINGULAR,
  VARMAPACK_NOT_POSITIVE_SEMIDEFINITE, VARMAPACK_INTERNAL
} varmapack_error;

const char *varmapack_strerror( varmapack_error error);

varmapack_error varmapack_sim ( // Simulate VARMA time series
  double A[],   // in      r×r×p, autoregressive parameter matrices
  double B[],   // in      r×r×q, moving average parameter matrices
  double Sig[], // in      r×r, covariance of the shock terms ε(t)
  double mu[],  // in      r×nmu, time series means; last supplied mean repeats
  int nmu,      // in      number of supplied mean vectors, 0 means zero mean
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving-average terms
  int r,        // in      dimension of each x(t)
  int n,        // in      length of each generated series
  int M,        // in      number of replicates to generate
  double X0[],  // in      r×nX0 with nX0 ≥ max(p,q); optional starter or NULL
  int nX0,      // in      number of starting values
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series (or NULL to skip)
  randompack_rng *rng // in/out  random number generator
  );

varmapack_error varmapack_simx ( // Simulate VARMAX time series
  double A[],   // in      r×r×p autoregressive parameter matrices
  double B[],   // in      r×r×q moving-average parameter matrices
  double C[],   // in      r×s exogenous coefficient vectors
  double Sig[], // in      r×r covariance of shock terms eps(t)
  double z[],   // in      length-n exogenous scalar sequence
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving-average terms
  int s,        // in      number of exogenous terms
  int r,        // in      dimension of each x(t)
  int n,        // in      length of each generated series
  int M,        // in      number of replicates to generate
  double X0[],  // in      r×h fixed startup values x(0),...,x(h-1)
  int h,        // in      number of fixed startup values
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series (or NULL to skip)
  randompack_rng *rng // in/out  random number generator
  );

double varmapack_specrad( // Spectral radius of model companion matrix
  double A[],   // in   r×r×p, autoregressive parameter matrices
  int r,        // in   dimension of each x(t), row count of A
  int p);       // in   number of autoregressive terms

double varmapack_ma_specrad( // Spectral radius of moving-average companion matrix
  double B[],   // in   r×r×q, moving-average parameter matrices
  int r,        // in   dimension of each shock, row count of B
  int q);       // in   number of moving-average terms

varmapack_error varmapack_acvf( // Theoretical autocovariance function of VARMA model
  double A[],    // in   r×r×p autoregressive matrices
  double B[],    // in   r×r×q moving average matrices
  double Sig[],  // in   r×r shock's covariance
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving-average terms
  int r,         // in   dimension of each x(t)
  double Gamma[],// out  r×r×(maxlag+1) covariance sequence, Γk = Cov(xt, x{t-k})
  int maxlag);   // in   largest lag to compute

varmapack_error varmapack_psi( // VARMA impulse-response coefficients
  double A[],   // in   r×r×p autoregressive parameter matrices
  double B[],   // in   r×r×q moving average parameter matrices
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving-average terms
  int r,        // in   dimension of each x(t)
  int maxlag,   // in   largest coefficient index to compute
  double Psi[]);// out  r×r×(maxlag+1), coefficient sequence, Psi_0,...

varmapack_error varmapack_irf( // Orthogonalized impulse responses
  double A[],     // in   r×r×p autoregressive parameter matrices
  double B[],     // in   r×r×q moving average parameter matrices
  double Sig[],   // in   r×r shock covariance matrix
  int p,          // in   number of autoregressive terms
  int q,          // in   number of moving-average terms
  int r,          // in   dimension of each x(t)
  int maxlag,     // in   largest coefficient index to compute
  double Theta[]);// out  r×r×(maxlag+1), Theta_j = Psi_j*chol(Sig)

varmapack_error varmapack_autocov( // Sample autocovariance of observed data
  const char *transp, // in   "N": X r×n with x(t) in column t; "T": n×r, x(t) in row t
  const char *norm,   // in   "ML" for 1/n scaling, "C" for 1/(n−k) correction
  int r,              // in   dimension of each observation
  int n,              // in   number of observations
  double X[],         // in   data matrix in storage indicated by transp
  int maxlag,         // in   number of lags to compute (≤ n−1)
  double C[]);        // out  r×r×(maxlag+1) array of lag-k covariance matrices

varmapack_error varmapack_testcase( // Create testcase for VARMA calculation
  double A[],    // out     r×r×p, autoregressive parameter matrices (or null)
  double B[],    // out     r×r×q, moving average parameter matrices (or null)
  double Sig[],  // out     r×r, covariance of the shock terms eps(t) (or null)
  char *name,    // out     string w. length >= 12, name of testcase, or ""
  int *pp,       // in/out  Number of autoregressive terms
  int *qp,       // in/out  Number of moving avg. terms
  int *rp,       // in/out  dimension of each x(t)
  int *icase,    // in/out  index of named testcase to create or 0 or -1 to use p,q,r
  double rho,    // in      target spectral radius when name is "rho"
  randompack_rng *rng); // in      random number generator

#ifdef __cplusplus
}
#endif
#endif
