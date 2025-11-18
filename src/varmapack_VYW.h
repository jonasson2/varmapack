// varmapack_VYW.h — shared solver for vector Yule–Walker setup
#ifndef VARMAPACK_VYWSOLVE_H
#define VARMAPACK_VYWSOLVE_H

#include <stdbool.h>

bool varmapack_VYWFactorizeSolve(
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices (ignored when q=0)
  double Sig[], // in   r×r, innovation covariance
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving-average terms
  int r,        // in   dimension of each x(t)
  double S[],   // out  r×r×(p+1), covariance sequence S0…Sp
  double C[],   // out  optional r×r×(q+1) buffer for Ci (can be null)
  double G[]);  // out  optional r×r×(q+1) buffer for Gi (can be null)

#endif
