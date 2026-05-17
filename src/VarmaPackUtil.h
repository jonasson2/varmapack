// VarmaPackUtil.h
// Utility functions for VARMA simulation and covariance computation
// (Part of varmapack / varmasim)

#include <stdio.h>
#include <stdbool.h>
#include "varmapack.h"
#include "varmapack_config.h"

// -----------------------------------------------------------------------------
// Compute the Ci and Gi matrices for a VARMA(p, q) model
// -----------------------------------------------------------------------------
HIDDEN void FindCG( double A[],   // in   r×r×p  autoregressive parameter matrices
  double B[],   // in   r×r×q  moving average parameter matrices
  double Sig[], // in   r×r    covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[],   // out  r×r×(q+1)  Ci = cov(x(t), eps(t−i))
  double G[]    // out  r×r×(q+1)  Gi = cov(y(t), x(t−i))
);

HIDDEN bool FindW(
  double B[],   // in   r×r×q  moving average parameter matrices
  double Sig[], // in   r×r    covariance of the shock terms eps(t)
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double W[]    // out  r×r×(q+1)  Wi = cov(y(t), y(t−i))
);

HIDDEN void WBuild(double Wlag[], int q, int r, int nW, double W[]);

HIDDEN void postmultiplySigmaPrime(double Y[], int ldY, int m, int h, double Sig[],
                                   int r, double Work[]);

HIDDEN void FindPsi(double *A, double *B, double *Psi, int p, int q, int r);
// TODO: doc

HIDDEN void FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h);
// TODO: doc

HIDDEN varmapack_error psdFactor(double Sig[], int r, double L[], bool *triangular);
// TODO: doc
