#include "BlasGateway.h"
#include "allocate.h"
#include "VarmaMisc.h"
void FindCGW ( // Calculate the Gi and Wi matrices for VARMALLC and VARMALLM
  double A[],   // in   r×r×p, autoregressive parameter matrices
  double B[],   // in   r×r×q, moving average parameter matrices
  double Sig[], // in   r×r, covariance of the shock terms eps(t)
  int p,        // in   number of autoregressive terms
  int q,        // in   number of moving average terms
  int r,        // in   dimension of each x(t)
  double C[],   // out  r×r×(q+1) with C0, C1...Cq where Ci = cov(x(t),eps(t-i))
  double G[],   // out  r×r×(q+1) with G0, G1...Gq where Gi = cov(y(t),x(t-i))
  double W[])   // out  r×r×(q+1) with W0, W1...Wq where Wi = cov(y(t),y(t-i))
{
  int i, j, rr = r*r;
  double *BSig, *Bj, *Cj, *Gj, *Wj, *Ai, *Bi, *Bimj, *Cimj, *Cjmi,*BjSig,*BiSig;
  //
  // CALCULATE B0·Sig,... Bq·Sig AND C0,... Cq IN BSig AND C:
  allocate(BSig, r*r*(q+1));
  copy(rr, Sig, 1, BSig, 1); // B0Sig = Sig
  copy(rr, Sig, 1, C, 1);    // C0 = Sig
  for (j=1; j<=q; j++) {
    Cj = C + j*rr;
    BjSig = BSig + j*rr; // find Bj·Sig:
    Bj = B + (j-1)*rr;
    gemm("NoT", "NoT", r, r, r, 1.0, Bj, r, Sig, r, 0.0, BjSig, r);
    copy(rr, BjSig, 1, Cj, 1);     // initialize Cj as Bj·Sig
    for (i=1; i<=p && i<=j; i++) { // finish calculation of Cj:
      Cjmi = C + (j-i)*rr;
      Ai = A + (i-1)*rr;
      gemm("NoT", "NoT", r, r, r, 1.0, Ai, r, Cjmi, r, 1.0, Cj, r);
    }
  }
  // NOW TURN ATTENTION TO THE Gj- AND Wj-MATRICES:
  for (j=0; j<=q; j++) {
    Gj = G + j*rr; // initialize Gj as Bj·Sig
    Wj = W + j*rr; // ... and also Wj
    copy(rr, BSig + j*rr, 1, Gj, 1);
    copy(rr, BSig + j*rr, 1, Wj, 1);
    for (i=j+1; i<=q; i++) {
      Bi = B + (i-1)*rr;
      BiSig = BSig + i*rr;
      Cimj = C + (i-j)*rr;
      Bimj = B + (i-j-1)*rr;
      gemm("NoT", "T", r, r, r, 1.0, Bi, r, Cimj, r, 1.0, Gj, r);
      gemm("NoT", "T", r, r, r, 1.0, BiSig, r, Bimj, r, 1.0, Wj, r);
    }
  }
  freem(BSig);
}
