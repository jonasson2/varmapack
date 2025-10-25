#include <math.h>
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "allocate.h"
#include "xAssert.h"

double VarmaSpecRad(double *A, int r, int p) {
  int i, k, lwork;
  int n = r*p;
  double *wr, *wi, *work, info, alwork, hyp;
  allocate(Ac, n*n);  
  allocate(wr, n);
  allocate(wi, n);
  lacpy("All", r, n, A, r, Ac, n); 
  for (i = k = r; i < n; i++, k += n + 1) Ac[k] = 1;
  geev("N", "N", n, Ac, wr, wi, 0, 0, 0, 0, &alwork, -1, &info);
  lwork = alwork;
  alloc(work, lwork);
  geev("N", "N", n, Ac, wr, wi, 0, 0, 0, 0, work, lwork, &info);
  xAssert(info == 0);
  for (k=0; k<n; k++) {
    hyp = hypot(wr[k], wi[k]);
    rho = max(rho, hyp);
  }
  free(work);
  free(wi);
  free(wr);
  free(Ac);
  return rho;
}
