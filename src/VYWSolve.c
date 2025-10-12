#include "xAssert.h"
#include "allocate.h"
#include "VarmaUtilities.h"
#include "BlasGateway.h"
#include "VYW.h"
void VYWSolve ( // Solve modified Yule-Walker equations
  double A[],   // in   r × r × p, autoregressive parameter matrices
  double LU[],  // in   LU factors of the N×N matrix F where N = r*r*p-r*(r-1)/2
  double S[],   // out  r × r × nrhs × (p+1), solution to system
  double Y[],   // in   r × r × nrhs × nY, rhs of system to solve
  int nrhs,     // in   number of right hand sides (1 or, for derivatives, nPar)
  int nY,       // in   number of Y arrays (either q+1 or p+1)
  int piv[],    // in   N-vector with pivoting info from VYWFactorize
  int p,        // in   number of autoregressive terms
  int r)        // in   dimension of each x(t)
{ 
  // If nrhs=1 and Y = G then S = [S0 S1...Sp], with Si = cov(x(t),x(t-i))
  // If nY=q+1 and q<p then Y(q+1)...Y(p) are set to 0.
  int i, j, k, info, rr = r*r;
  double *g, *g1, *Yp, *g0k, *Ykp, *Ap, *gptr, *gii0k, *gj, *Sk, *Sii0k, *Sp;
  double *Skp, *Skpmi, *Ai;
  int mg = r*r*(p+1);         // First 3 dimensions of g
  int N = r*r*p - r*(r-1)/2;  // Order of VYW-system
  // INITIALIZE g0
  allocate(g, mg*nrhs);
  g1 = g + r*(r+1)/2;
  Yp = Y + rr*nrhs*p;
  for (k=0; k<nrhs; k++) {
    g0k = g + mg*k;
    copytranspose(r, r, Y+rr*k, r, g0k, r);
    if (nY>p && p>0) { // g0 = Y0' + Yp·Ap'
      Ykp = Yp + rr*k;
      Ap = A + rr*(p-1);
      gemm("NoT", "T", r, r, r, 1.0, Ykp, r, Ap, r, 1.0, g0k, r);
    }
  }
  // REMOVE ROWS/COLUMNS CORRESPONDING TO UPPER-TRIANGLE FROM g0
  for (k=0; k<nrhs; k++) {
    gptr = g + mg*k;
    for (i=0; i<r; i++) {
      gii0k =  g+(1+r)*i+mg*k;
      for (j=0; j<r-i; j++) *gptr++ = *gii0k++;
    }
    for (i=0; i<r*(r-1)/2; i++) *gptr++=0;
  }
  // LET k-TH COLUMN OF g = Y(:,:,k,:) FOR k>=1
  for (j=1; j<nY && j<p; j++) {
    gj = g1 + (j-1)*rr;
    lacpy("All", rr, nrhs, Y + rr*nrhs*j, rr, gj, mg);
  }
  // SOLVE
  if (p>0) {
    getrs("NoT", N, nrhs, LU, N, piv, g, mg, &info);
    xAssert(info == 0); // only nonzero if an argument is wrong
  }
  // NOW FORM S FROM THE SOLUTION
  if (p>0) {
    for (k=0; k<nrhs; k++) { // S0 = G0 + G0' where G0 = reshape(g0,r,r)
      Sk = S + rr*k;
      gptr = g + mg*k;
      for (i=0; i<r; i++) {
        Sii0k = Sk + (1+r)*i;
        copy(r-i, gptr, 1, Sii0k, 1);       // copy to column
        copy(r-i-1, gptr+1, 1, Sii0k+r, r); // copy to row
        *Sii0k *= 2;                        // double diagonal element
        gptr += r-i;
      }
    }
  }
  for (j=1; j<p; j++) { // Copy g1...gp to S1...Sp
    gj = g1 + (j-1)*rr;
    lacpy("All", rr, nrhs, gj, mg, S + rr*nrhs*j, rr);
  }
  // Finally form Sp:
  Sp = S + rr*nrhs*p;
  if (nY>p) lacpy("All", rr, nrhs, Yp, rr, Sp, rr);
  for (k=0; k<nrhs; k++) {
    Skp = Sp + rr*k;
    for (i=1; i<=p; i++) {
      Skpmi = S + rr*k + rr*nrhs*(p-i);
      Ai = A + rr*(i-1);
      gemm("NoT", "NoT", r, r, r, 1.0, Ai, r, Skpmi, r, 1.0, Skp, r);
    }
  }
  freem(g);
}
