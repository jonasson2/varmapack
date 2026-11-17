/// Main header file for varmapack.

#ifndef VARMAPACK_H
#define VARMAPACK_H
#include "varmapack_config.h"
#ifdef __cplusplus
extern "C" {
#endif

#include "RandompackGateway.h"
#include <stdbool.h>
#include <stdio.h>

#define VARMAPACK_VERSION "0.1.2"
#define VARMAPACK_TESTCASE_NAME_LEN 32

char *varmapack_last_error(void); // Get last error string for this thread, or 0

typedef struct varmapack_ou_model varmapack_ou_model;

varmapack_ou_model *varmapack_ou_create( // Stationary diagonal-drift OU model
  double tau[], // in   r finite positive persistence times
  double cov[], // in   r*r stationary covariance; lower triangle is used
  double mu[],  // in   r constant means, or 0 for zero means
  int r);       // in   number of series
  // Copies its inputs. cov and D_ij = cov_ij*(1/tau_i + 1/tau_j) must be PSD.
  // Singular covariances are allowed. Returns 0 on failure; see last_error().

void varmapack_ou_free(varmapack_ou_model *model); // Accepts 0

bool varmapack_ou_sim( // Exact OU simulation, with no burn-in or time discretization
  varmapack_ou_model *model, // in   model returned by varmapack_ou_create
  double t[], // in   n finite strictly increasing observation times
  int n,      // in   observations per path, at least 1
  int M,      // in   number of independent paths, at least 1
  double X[], // out  r*n*M values, column-major; must not overlap t
  randompack_rng *rng); // in/out Randompack stream
  // X(t[0]) is drawn from N(mu,cov). Each gap d uses a_i = exp(-d/tau_i)
  // and innovation covariance Q_ij = cov_ij*(1-a_i*a_j).
  // Times and tau have the same units; translating t does not change the law.
  // Returns false on failure. X and the RNG can be partly advanced on failure.

bool varmapack_gamma_times( // Observation times with independent Gamma gaps
  int n,           // in   number of times, at least 1
  double shape,    // in   finite positive shape, or +infinity for equal gaps
  double mean_gap, // in   finite positive expected spacing
  double t[],      // out  n strictly increasing times, starting at zero
  randompack_rng *rng); // in/out stream; may be 0 for equal gaps or n=1
  // Gap scale = mean_gap/shape, CV = 1/sqrt(shape). No sample normalization.
  // Returns false if a positive gap/time cannot be represented as a double.
  // The output and RNG can be partly advanced on failure.

bool varmapack_sim ( // Simulate VARMA time series
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
  // Stationary simulation with X0 requires positive-definite startup covariance.
  // Nonstationary simulation requires X0; with MA terms Sig must be positive definite.
  // See varmapack_sim.c for details.

bool varmapack_simx ( // Simulate VARMAX time series
  double A[],   // in      r×r×p autoregressive parameter matrices
  double B[],   // in      r×r×q moving-average parameter matrices
  double C[],   // in      r×d×s exogenous coefficient matrices C1,...,Cs
  double Sig[], // in      r×r covariance of shock terms eps(t)
  double z[],   // in      d×n×Mz exogenous input sequences, or 0 when s is 0
  int Mz,       // in      number of z input sequences, must be 1 or M
  int p,        // in      number of autoregressive terms
  int q,        // in      number of moving-average terms
  int s,        // in      number of exogenous terms
  int d,        // in      dimension of each exogenous input vector
  int r,        // in      dimension of each x(t)
  int n,        // in      length of each generated series
  int M,        // in      number of replicates to generate
  double X0[],  // in      r×h×MX0 fixed startup values, or 0 when h is 0
  int h,        // in      number of fixed startup values, h >= max(p,q,s-1)
  int MX0,      // in      number of X0 startup paths, must be 1 or M
  double X[],   // out     r×n×M generated series
  double E[],   // out     r×n×M shock series, or 0 to skip
  randompack_rng *rng // in/out  random number generator
  );
  // Fixed-start VARMAX simulation requires positive-definite shock covariance.

double varmapack_specrad( // Spectral radius of model companion matrix
  double A[],   // in   r×r×p, autoregressive parameter matrices
  int r,        // in   dimension of each x(t), row count of A
  int p);       // in   number of autoregressive terms

double varmapack_ma_specrad( // Spectral radius of moving-average companion matrix
  double B[],   // in   r×r×q, moving-average parameter matrices
  int r,        // in   dimension of each shock, row count of B
  int q);       // in   number of moving-average terms

bool varmapack_psi( // VARMA impulse-response coefficients
  double A[],   // in   r×r×p autoregressive parameter matrices
  double B[],   // in   r×r×q moving average parameter matrices
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving-average terms
  int r,        // in   dimension of each x(t)
  int maxlag,   // in   largest coefficient index to compute
  double Psi[]);// out  r×r×(maxlag+1), coefficient sequence, Psi_0,...

bool varmapack_irf( // Orthogonalized impulse responses
  double A[],     // in   r×r×p autoregressive parameter matrices
  double B[],     // in   r×r×q moving average parameter matrices
  double Sig[],   // in   r×r shock covariance matrix; only lower triangle is used
  int p,          // in   number of autoregressive terms
  int q,          // in   number of moving-average terms
  int r,          // in   dimension of each x(t)
  int maxlag,     // in   largest coefficient index to compute
  double Theta[]);// out  r×r×(maxlag+1), Theta_j = Psi_j*L, L*L^T = Sig

bool varmapack_acvf( // Theoretical autocovariance function of VARMA model
  double A[],    // in   r×r×p autoregressive matrices
  double B[],    // in   r×r×q moving average matrices
  double Sig[],  // in   r×r shock's covariance
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving-average terms
  int r,         // in   dimension of each x(t)
  double Gamma[],// out  r×r×(maxlag+1) covariance sequence, Γk = Cov(xt, x{t-k})
  int maxlag);   // in   largest lag to compute
  // See varmapack_acvf.c for more details

bool varmapack_cov2sig( // Recover innovation covariance from model covariance
  double A[],   // in   r*r*p autoregressive matrices; may be 0 when p=0
  double B[],   // in   r*r*q moving-average matrices; may be 0 when q=0
  double cov[], // in   r*r target Cov(x_t,x_t); only the lower triangle is used
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving-average terms
  int r,        // in   dimension of each x(t)
  double Sig[]);// out  r*r symmetric innovation covariance; may equal cov
  // Arrays are column-major. A must define a stationary model. cov and the
  // recovered Sig must be positive semidefinite within floating-point tolerance;
  // singular covariances are allowed.
  // Uses the direct inverse VYW system, without switching to a Lyapunov solver.
  // Returns false for invalid input, an incompatible target, a singular or
  // ill-conditioned inverse system, or a failed forward covariance check.
  // Consult varmapack_last_error() on failure; Sig is then unchanged.
  // Except for Sig == cov, the output must not overlap any input.

bool varmapack_autocov( // Sample autocovariance of observed data
  const char *transp, // in   string beginning "N" for X r×n, or "T" for X n×r
  const char *norm,   // in   string beginning "M" for 1/n, or "C" for 1/(n−k)
  int r,              // in   dimension of each observation
  int n,              // in   number of observations
  double X[],         // in   data matrix in storage indicated by transp
  int maxlag,         // in   number of lags to compute (≤ n−1)
  double C[]);        // out  r×r×(maxlag+1), Ck = Cov(x(t), x(t−k))
  // For "N", x(t) is in column t; for "T", x(t) is in row t. Initial letters are
  // case-insensitive, BLAS-style.

bool varmapack_cov2corr( // Convert a covariance sequence to correlations
  double cov[], // in   r×r×(maxlag+1), with variances on the diagonal of cov(:,:,0)
  int r,        // in   dimension of each observation
  int maxlag,   // in   largest lag in cov
  double corr[]);// out  same shape as cov; may equal cov for in-place conversion
  // Every lag is scaled using the variances in cov(:,:,0). If corr != cov, the
  // arrays must not overlap. Positive-lag values are not clipped to [-1,1]; in
  // particular, corrected sample autocovariances can produce values outside it.

bool varmapack_testcase( // Create testcase for VARMA calculation
  char *name,    // in/out  testcase name, "rho" to use rho, or "max" for max-inquiry
  int *index,    // in/out  named testcase index, or 0/-1 to use p,q,r
  double rho,    // in      target spectral radius when name is "rho"
  int *pp,       // in/out  number of autoregressive terms
  int *qp,       // in/out  number of moving-average terms
  int *rp,       // in/out  dimension of each x(t)
  double A[],    // out     r×r×p, autoregressive parameter matrices, or 0
  double B[],    // out     r×r×q, moving average parameter matrices, or 0
  double Sig[],  // out     r×r, covariance of the shock terms eps(t), or 0
  randompack_rng *rng); // in      random number generator
  // For max-inquiry set name = "max" and receive named-case count in index.
  // For a dimension inquiry set A = B = Sig = 0; "rho" is not supported in this mode.
  // When name returns a testcase name, it must be a writable buffer of at least
  // VARMAPACK_TESTCASE_NAME_LEN characters.
  // For details see varmapack_testcase.c

bool varmapack_testcasex( // Create deterministic exogenous testcase data
  int s,        // in   number of exogenous terms
  int d,        // in   dimension of each exogenous input vector
  int r,        // in   dimension of each x(t)
  int n,        // in   number of exogenous values to create
  double C[],   // out  r×d×s exogenous coefficient matrices
  double z[]);  // out  d×n exogenous input sequence

#ifdef __cplusplus
}
#endif
#endif
