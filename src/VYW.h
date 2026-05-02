// VYW.h:  Declaration of Vector-Yule-Walker functions
#ifndef VYW_H
#define VYW_H

void VYWFactorize ( // Set up and factorize matrix of Yule-Walker equations
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double LU[],  // out  LU factors of the N×N matrix F where N = r*r*p-r*(r-1)/2
  int piv[],    // out  N-vector with pivoting information
  int p,        // in   number of autoregressive terms
  int r,        // in   dimension of each x(t)
  int *info);   // out  nonzero if factorization fails due to singularity

void VYWSolve ( // Solve modified Yule-Walker equations
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double LU[],  // in   LU factors of the N×N matrix F where N = r*r*p-r*(r-1)/2
  double S[],   // out  r×r×nrhs×(p+1), solution to system
  double Y[],   // in   r×r×nrhs×nY, rhs of system to solve
  int nrhs,     // in   number of right hand sides (must be 1 in this version)
  int nY,       // in   number of Y arrays (either q+1 or p+1)
  int piv[],    // in   N-vector with pivoting info from VYWFactorize
  int p,        // in   number of autoregressive terms
  int r);       // in   dimension of each x(t)

#endif
