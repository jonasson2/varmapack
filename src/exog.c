#include "BlasGateway.h"
#include "varmapack_config.h"
#include "VarmaPackUtil.h"
#include "VarmaUtilities.h"

HIDDEN void WBuild(double Wlag[], int q, int r, int nW, double W[]) {
  // W is block Toeplitz. With lower-lag blocks Wk = cov(y(t), y(t-k)),
  // the (i,j) lower block is W{i-j}; upper blocks are filled by transpose.
  int i, j, k;
  int rW = r*nW;
  int rr = r*r;
  double *Wij, *Wk;
  laset("All", rW, rW, 0, 0, W, rW);
  for (j=0; j<nW; j++) {
    for (i=j; i<nW; i++) {
      k = i - j;
      if (k > q) continue;
      Wij = W + j*r*rW + i*r;
      Wk = Wlag + k*rr;
      lacpy("All", r, r, Wk, r, Wij, rW);
      if (i > j) {
        Wij = W + i*r*rW + j*r;
        copytranspose(r, r, Wk, r, Wij, rW);
      }
    }
  }
}

HIDDEN void postmultiplySigmaPrime(double Y[], int ldY, int m, int h, double Sig[],
                                   int r, double Work[]) {
  // Replace Y by Y*Sigma', where Sigma' = I_h kron Sig.
  for (int j=0; j<h; j++) {
    double *Yj = Y + j*r*ldY;
    symm("Right", "Low", m, r, 1, Sig, r, Yj, ldY, 0, Work, m);
    lacpy("All", m, r, Work, m, Yj, ldY);
  }
}
