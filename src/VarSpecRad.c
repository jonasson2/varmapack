// VarSpecRad  Spectral radius of companion matrix of a VAR process
#include <math.h>
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "allocate.h"
#include "xAssert.h"

double VarSpecRad(double *A, int r, int p) {
  int i, k, lwork, info;
  int n = r*p;
  double *Ac, *wr, *wi, *work, alwork, hyp, rho, dummy[1];
  if (n == 0) return 0;
  allocate(Ac, n*n);  
  allocate(wr, n);
  allocate(wi, n);
  lacpy("All", r, n, A, r, Ac, n); 
  for (i = k = r; i < n; i++, k += n + 1) Ac[k] = 1;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, &alwork, -1, &info);
  lwork = (int)alwork;
  allocate(work, lwork);
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, work, lwork, &info);
  xAssert(info == 0);
  rho = 0.0;
  for (k=0; k<n; k++) {
    hyp = hypot(wr[k], wi[k]);
    rho = max(rho, hyp);
  }
  freem(work);
  freem(wi);
  freem(wr);
  freem(Ac);
  return rho;
}
