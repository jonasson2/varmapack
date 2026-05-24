// varmapack_specrad  Spectral radius of companion matrix of a VAR process
#include <math.h>
#include "BlasGateway.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"
#include "error.h"

double varmapack_specrad(double *A, int r, int p) {
  return CompanionSpecrad(A, r, p, 1);
}

double CompanionSpecrad(double P[], int r, int nblocks, double sign) {
  int i, j, lwork, info;
  int n = r*nblocks;
  double *Ac, *wr, *wi, *work, alwork, hyp, rho, dummy[1];
  if (nblocks < 0 || r <= 0 || (nblocks > 0 && P == 0)) return NAN;
  if (nblocks == 0) return 0;
  Ac = wr = wi = work = 0;
  if (!ALLOC(Ac, n*n)) goto fail;
  if (!ALLOC(wr, n)) goto fail;
  if (!ALLOC(wi, n)) goto fail;
  lacpy("All", r, n, P, r, Ac, n);
  if (sign != 1) {
    for (j=0; j<n; j++) scal(r, sign, Ac + j*n, 1);
  }
  for (i = j = r; i < n; i++, j += n + 1) Ac[j] = 1;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, &alwork, -1, &info);
  lwork = (int)alwork;
  if (!ALLOC(work, lwork)) goto fail;
  geev("N", "N", n, Ac, n, wr, wi, dummy, 1, dummy, 1, work, lwork, &info);
  xAssert(info == 0);
  rho = 0;
  for (j=0; j<n; j++) {
    hyp = hypot(wr[j], wi[j]);
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
