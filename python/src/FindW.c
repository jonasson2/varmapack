#include "BlasGateway.h"
#include "error.h"
#include "varmapack_config.h"
#include "VarmaPackUtil.h"

HIDDEN bool FindW(
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double W[])   // out  r×r×(q+1) with Wi = cov(y(t), y(t-i))
{
  int i, k, rr = r*r;
  double *BSig = 0;
  double *Wk, *Bi, *BipkSig;
  if (!ALLOC(BSig, rr*(q+1))) return false;
  copy(rr, Sig, 1, BSig, 1);
  for (i=1; i<=q; i++) {
    Bi = B + (i-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1, Bi, r, Sig, r, 0, BSig + i*rr, r);
  }
  for (k=0; k<=q; k++) {
    Wk = W + k*rr;
    copy(rr, BSig + k*rr, 1, Wk, 1);
    for (i=1; i<=q-k; i++) {
      Bi = B + (i-1)*rr;
      BipkSig = BSig + (i+k)*rr;
      gemm("NoT", "T", r, r, r, 1, BipkSig, r, Bi, r, 1, Wk, r);
    }
  }
  FREE(BSig);
  return true;
}
