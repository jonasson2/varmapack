// varmapack_autocov — sample autocovariance of observed data

#include <stdbool.h>
#include "varmapack.h"
#include "BlasGateway.h"
#include "VarmaUtilities.h"

void varmapack_autocov(const char *transp, const char *norm, int r, int n,
                       double X[], int maxlag, double C[]) {
  xAssert(transp != 0 && norm != 0);
  bool ML = (norm[0] == 'M' || norm[0] == 'm');
  bool CORRECTED = (norm[0] == 'C' || norm[0] == 'c');
  xAssert(ML || CORRECTED);
  bool TRANSP = !(transp[0] == 'N' || transp[0] == 'n');
  xAssert((!TRANSP && (transp[0] == 'N' || transp[0] == 'n')) ||
          (TRANSP && (transp[0] == 'T' || transp[0] == 't')));
  xAssert(maxlag <= n-1);
  double *Xc;
  double *mu;
  xAssert(ALLOC(Xc, r*n));
  xAssert(ALLOC(mu, r));
  copy(r*n, X, 1, Xc, 1);
  if (TRANSP) {
    meanmat("N", n, r, Xc, n, mu);
    for (int i=0; i<n; i++) {
      double *row = Xc + i;
      axpy(r, -1.0, mu, 1, row, n);
    }
  }
  else {
    meanmat("T", r, n, Xc, r, mu);
    for (int j=0; j<n; j++) {
      double *col = Xc + j*r;
      axpy(r, -1.0, mu, 1, col, 1);
    }
  }
  double *Y = Xc;
  for (int k=0; k<=maxlag; k++) {
    double fctr = ML ? 1.0/n : 1.0/(n - k);
    double *Ck = C + r*r*k;
    if (TRANSP) {
      double *Z = Xc + k;
      gemm("T", "N", r, r, n-k, fctr, Y, n, Z, n, 0.0, Ck, r);
    }
    else {
      double *Z = Xc + r*k;
      gemm("N", "T", r, r, n-k, fctr, Y, r, Z, r, 0.0, Ck, r);
    }
  }
  FREE(mu);
  FREE(Xc);
}
