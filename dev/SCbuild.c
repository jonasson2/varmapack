#include "allocate.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "VarmaMisc.h"
#include "VYW.h"
#include "Tests.h"
void SCbuild(double A[], double B[], double Sig[], int p, int q, int r, int n,
                    double CC[], double SS[]) {
  // CC = cov(x,eps) and SS = cov(x,x) (s = t-i). Used by testing functions.
  // See also CCBuild.c
  int *piv, info, nVYW = r*r*p - r*(r-1)/2, k, j;
  double *C, *G, *W, *S, *vywFactors;
  allocate(C, r*r*(q+1));
  allocate(G, r*r*(q+1));
  allocate(W, r*r*(q+1));
  allocate(S, r*r*(p+1));
  FindCGW(A, B, Sig, p, q, r, C, G, W);
  nVYW = max(0, nVYW);
  allocate(vywFactors, nVYW * nVYW);
  allocate(piv, nVYW);
  VYWFactorize(A, vywFactors, piv, p, r, &info);
  xAssert(info == 0);
  VYWSolve(A, vywFactors, S, G, 1, q+1, piv, p, r);
  SBuild("All", S, A, G, p, q, r, n, SS);
  setzero(r*n * r*n, CC);
  for (k=0; k<n; k++) {
    for (j=0; j<=q && k+j<n; j++) {
      // CC{k+j, k} = C{j};
      lacpy("All", r, r, &C[j*r*r], r, &CC[k*r*r*n + (k+j)*r], r*n); 
    }
  }
  freem(piv); freem(vywFactors); freem(S); freem(W); freem(G); freem(C);
}
