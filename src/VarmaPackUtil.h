// VarmaPackUtil.h
// Utility functions for VARMA simulation and covariance computation
// (Part of varmapack / varmasim)

#include <stdio.h>
#include <stdbool.h>

// -----------------------------------------------------------------------------
// Compute the Ci and Gi matrices for a VARMA(p, q) model
// -----------------------------------------------------------------------------
void FindCG(
  double A[],   // in   r×r×p  autoregressive parameter matrices
  double B[],   // in   r×r×q  moving average parameter matrices
  double Sig[], // in   r×r    covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[],   // out  r×r×(q+1)  Ci = cov(x(t), eps(t−i))
  double G[]    // out  r×r×(q+1)  Gi = cov(y(t), x(t−i))
);

void flipmat(double *A, double *Aflp, int r, int k);
// [A1 A2...Ak] ––> [Ak...A2 A1]

void FindPsi(double *A, double *B, double *Psi, int p, int q, int r);
// TODO: doc

void FindPsiHat(double *Psi, double *Psi_hat, double *Sig, int r, int h);
// TODO: doc
