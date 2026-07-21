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

#define VARMAPACK_VERSION "0.1.0"
#define VARMAPACK_TESTCASE_NAME_LEN 32

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
  double X0[],  // in      r×nX0×MX0 starting values, or 0
  int nX0,      // in      number of starting values
  int MX0,      // in      number of X0 starting paths, must be 1 or M
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series, or 0 to skip
  randompack_rng *rng // in/out  random number generator
  );

varmapack_error varmapack_simx ( // Simulate VARMAX time series
  double A[],   // in      r×r×p autoregressive parameter matrices
  double B[],   // in      r×r×q moving-average parameter matrices
  double C[],   // in      r×s exogenous coefficient vectors
  double Sig[], // in      r×r covariance of shock terms eps(t)
  double z[],   // in      n×Mz exogenous input sequences, Mz is 1 or M
  int Mz,       // in      number of z input sequences, must be 1 or M
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving-average terms
  int s,        // in      number of exogenous terms
  int r,        // in      dimension of each x(t)
  int n,        // in      length of each generated series
  int M,        // in      number of replicates to generate
  double X0[],  // in      r×h×MX0 fixed startup values x(0),...,x(h-1)
  int h,        // in      number of fixed startup values, h > max(p,s-1)
  int MX0,      // in      number of X0 startup paths, must be 1 or M
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series, or 0 to skip
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
  char *name,    // in/out  testcase name or output buffer
  int *index,    // in/out  named testcase index, or 0/-1 to use p,q,r
  double rho,    // in      target spectral radius when name is "rho"
  int *pp,       // in/out  number of autoregressive terms
  int *qp,       // in/out  number of moving-average terms
  int *rp,       // in/out  dimension of each x(t)
  double A[],    // out     r×r×p, autoregressive parameter matrices, or 0
  double B[],    // out     r×r×q, moving average parameter matrices, or 0
  double Sig[],  // out     r×r, covariance of the shock terms eps(t), or 0
  randompack_rng *rng); // in      random number generator

#ifdef __cplusplus
}
#endif
#endif
