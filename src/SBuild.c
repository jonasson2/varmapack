#include "BlasGateway.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"
#include "allocate.h"
void SBuild( // Build covariance matrix of all the values of a VARMA time series
  char *uplo,  // in   Create all of SS (when "A") or lower part only (when "L")
  double S[],  // in   r × r·(p+1), = [S0...Sp], Si = Cov(x(t),x(t-i))
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double G[],  // in   r × r·(q+1), = [G0...Gq], Gi = cov(y(t),x(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n can be any nonnegative value, even 0)
  double SS[]) // out  r·n × r·n, covariance of [x1'...xn']'
{
  //  DESCRIPTION: The time series is given by
  //
  //                  x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)
  //  where
  //                  y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q),
  //
  //  and x(t), y(t) and eps(t) are r-dimensional with eps(t) N(0,Sig). S and G
  //  can (for example) have been obtained with FindCGW. The SS matrix is:
  //
  //                   S0  S1' S2'...Sn-1'
  //                   S1  S0  S1'...Sn-2'
  //                   S2               :
  //                   :                :
  //                   Sn-1 ...... S1  S0
  //
  //  and the Sj are found with the recurrence relation:
  //
  //      Sj = A1*S(j-1) + A2*S(j-2) + ........ + Ap*S(j-p) + Gj
  //
  //  with Gj = 0 for j > q.
  //
  // NOTE: This function may (of course) be used to determine the covariance 
  // matrix of a segment of a timeseries by specifying n smaller then the total
  // length of the series (for example an initial segment, but since the series
  // is stationary all segments of the same length have the same covariance).
  double *Scol, *SSj, *SSi;
  int j, m;
  if (n==0) return;
  m = max(p+1,n);
  allocate(Scol, (r*m)*r);
  SExtend(A, G, S, Scol, p, q, r, m, 0, 0, 0, 0);
  for (j=0; j<n; j++) {
    SSj = SS + j*r*n*r + j*r;
    if (uplo[0] == 'A') {
      SSi = SSj + r*n*r;
      lacpy("All", r*(n-j), r, Scol, r*m, SSj, r*n);
      if (j < n-1) copytranspose(r*(n-j-1), r, Scol+r, r*m, SSi, r*n);
    }
    else {
      lacpy("Low", r*(n-j), r, Scol, r*m, SSj, r*n);
    }
  }
  freem(Scol);
}
