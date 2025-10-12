#ifndef VYW_H
#define VYW_H
// VYW.h:  Declaration of Vector-Yule-Walker functions

#include "mds.h"
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
  int nrhs,     // in   number of right hand sides (1 or, for derivatives, nPar)
  int nY,       // in   number of Y arrays (either q+1 or p+1)
  int piv[],    // in   N-vector with pivoting info from VYWFactorize
  int p,        // in   number of autoregressive terms
  int r);       // in   dimension of each x(t)

void VYWDeriv ( // Derivative of solution to modified Yule-Walker equations
  double A[],       // in   r×r×p, autoregressive parameter matrices
  struct mds RHS[], // in   RHS in mds-format (G0...Gq or Sig along with derivs)
  double S[],       // in   r×r×(p+1), Si-matrices; Si = cov(x(t),x(t-i))
  double Sd[],      // out  r×r×nPar×(p+1), derivatives of the Si
  double LU[],      // in   N×N, LU-factors of VYW-matrix, N = r·r·p - r·(r-1)/2
  int    piv[],     // in   pivoting vector for the VYW-system
  int    p,         // in   number of autoregressive terms
  int    q,         // in   number of moving average terms
  int    r,         // in   dimension of each x(t)
  int    nPar);     // in   number of model parameters to differentiate w.r.t.

void VYWDerivRHS ( // RHS of equations for derivat. of vector-Yule-Walker soln
  double A[],       // in   r×r×p, autoregressive parameter matrices
  struct mds RHS[], // in   RHS in mds-format (G0...Gq or Sig along with derivs)
  double S[],       // in   r×r×(p+1), Si-matrices; Si = cov(x(t),x(t-i))
  int    p,         // in   number of autoregressive terms
  int    q,         // in   number of moving average terms
  int    r,         // in   dimension of each x(t)
  int    nPar,      // in   number of model parameters to differentiate w.r.t.
  double Sol[]);    // out  r×r×nPar×(p+1), the calculated Sol

#endif
