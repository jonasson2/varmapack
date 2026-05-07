// varmapack_specrad  Spectral radius of companion matrix of a VAR process
#include <math.h>
#include <stdbool.h>
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "error.h"

double varmapack_specrad(double *A, int r, int p) {
  int i, k, lwork, info;
  int n = r*p;
  double *Ac, *wr, *wi, *work, alwork, hyp, rho, dummy[1];
  if (p < 0 || r <= 0 || (p > 0 && A == 0)) return NAN;
  if (p == 0) return 0;
  Ac = wr = wi = work = 0;
  if (!ALLOC(Ac, n*n)) goto fail;
  if (!ALLOC(wr, n)) goto fail;
  if (!ALLOC(wi, n)) goto fail;
  lacpy("All", r, n, A, r, Ac, n); 
  for (i = k = r; i < n; i++, k += n + 1) Ac[k] = 1;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, &alwork, -1, &info);
  lwork = (int)alwork;
  if (!ALLOC(work, lwork)) goto fail;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, work, lwork, &info);
  xAssert(info == 0);
  rho = 0.0;
  for (k=0; k<n; k++) {
    hyp = hypot(wr[k], wi[k]);
    rho = fmax(rho, hyp);
  }
  FREE(work);
  FREE(wi);
  FREE(wr);
  FREE(Ac);
  return rho;
fail:
  FREE(work);
  FREE(wi);
  FREE(wr);
  FREE(Ac);
  return NAN;
}
