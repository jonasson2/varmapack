// IsStationary  Check if VARMA model is stationary
#include "allocate.h"
#include "BlasGateway.h"
#include "VarmaMisc.h"
#include "Omega.h"
#include "VYW.h"

int IsStationary ( // Return 1 if model is stationary
  double A[],  // in  r×r×p, autoregressive parameter matrices
  double Sig[],// in  r×r, covariance matrix of shocks
  double LU[], // in  N×N with N = r*r*p-r*(r-1)/2, factors from vyw_factorize
  int piv[],   // in  N-vector, pivots from vyw_factorize
  int p,       // in  number of autoregressive terms
  int r)       // in  dimension of each x(t)
// The model is x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t) where y(t) given by a 
// moving average process, y(t) = B1·y(t-1) + ... + Bq·y(t-q) + eps(t), and 
// eps(t) is N(0,Sig). The stationarity of x(t) does not depend on B so that
// need not be a parameter.
//
// The test is based on checking positive definiteness of Sig and solving the 
// modified Yule-Walker equations for the AR process x(t) = A1·x(t-1) + ... + 
// Ap·x(t-p)+eps(t) with eps(t) ~ N(0,I), but the solution to these equations
// gives a postive definite S0 = Var(x(t)) if and only if all the roots of the 
// polynomial fi(b) = I - A1·b - A2·b^2 - ... - Ap·B^p are outside the unit 
// circle, which characterizes a stationary VARMA model.
  
{
  int i, q=0, info;
  double *I, *S, *LSig, *Su, *Olow=0;
  if (p==0) return 1;
  if (piv[0] == 0) return 0;
  allocate(LSig, r*r);
  copy(r*r, Sig, 1, LSig, 1);
  potrf("Low", r, LSig, r, &info);
  if (info) { freem(LSig); return 0; } // Sig not positive definite
  allocate(I, r*r);
  for (i=0; i<r*r; i+=(r+1)) I[i] = 1.0;
  allocate(S, r*r*(p+1));
  allocate(Su, r*p*r*p);
  VYWSolve(A, LU, S, I, 1, q+1, piv, p, r);
  OmegaBuild(Su, Olow, S, I, I, p, q, r, p);
  potrf("Low", r*p, Su, r*p, &info);
  freem(LSig); freem(I); freem(S); freem(Su);
  return (info==0);
}
