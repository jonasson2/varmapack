#include "BlasGateway.h"
#include "VarmaMisc.h"
#include "VarmaUtilities.h"
#include "xAssert.h"

void CCBuild( // Build covariance between terms and shocks of VARMA time series
  double A[],  // in   r × r × p, A=[A1...Ap], autoregressive parameter matrices
  double C[],  // in   r × r·(q+1), = [C0...Cq], Ci = cov(x(t),eps(t-i))
  int p,       // in   number of autoregressive terms
  int q,       // in   number of moving average terms
  int r,       // in   dimension of xt
  int n,       // in   length of series (n must be <= max(p+q)+1)
  double CC[]) // out  r·n × r·n, cov(x1'...xn',eps1'...epsn')
{
  //  DESCRIPTION: The time series is as shown in SBuild (e.g.). The CC-matrix
  //  is:
  //                   C0  0 ...... 0 
  //                   C1 C0  0 ... 0
  //                   :  C1  .     :
  //                   :         .  0
  //                   C(n-1).. C1 C0
  //
  //  For q < j <= p, Cj is found with the recurrence relation:
  //
  //      Cj = A1*C(j-1) + A2*C(j-2) + ........ + Aj*C0   (*)
  //
  int i,j;
  double *CCj;
  xAssert(n <= max(p,q) + 1);
  // Fill CC with zeros:
  setzero(r*n*r*n, CC);
  // Copy C0..Cq to first q+1 block-rows in first block-column:
  for (j=0; j<=q && j<n; j++) lacpy("All", r, r, C+j*r*r, r, CC + j*r, r*n);
  // Calculate rest of first block-column using (*):
  for (j=q+1; j<n; j++) {
    CCj = CC + j*r;
    for (i=1; i<=j && i<=p; i++)  // (using that CC was set to 0)
      gemm("N", "N", r, r, r, 1.0, A+(i-1)*r*r, r, CC+(j-i)*r, r*n,1.0,CCj,r*n);
  }
  // Copy head of first block-column to other columns
  for (j=1; j<n; j++) {
    CCj = CC + j*r*n*r + j*r;
    lacpy("All", r*(n-j), r, CC, r*n, CCj, r*n);
  }
}
