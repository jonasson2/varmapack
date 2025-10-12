#include "BlasGateway.h"
#include "VarmaMisc.h"

static void AGdiff(double[],int, double[], double[],int, double[], int,int,int);

void SExtend ( // Extend Sj matrices to include S(p+1)...S(n-1)
  double A[],    // in   r × r × p, autoregressive parameter matrices
  double G[],    // in   r × r × (q+1), G0, G1,...,Gq
  double S[],    // in   r × r × (p+1), S0, S1,...,Sp
  double Scol[], // out  r·n × r, column of Si-matrices S0, S1,..., S(n-1)
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of x(t)
  int n,         // in   length of time series
  double Gd[],   // in   r × r × nPar × (q+1), derivat. of G (or null if nPar=0)
  double Sd[],   // in   r × r × nPar × (p+1), derivat. of S (or null if nPar=0)
  double Scold[],// out  r·n × r × nPar, derivative of Scol (or null)
  int nPar)      // in   n of params, set to zero to suppress derivative calc
  //
  // It is assumed that Scol and Scold are zero on entry.
{
  int DIFF, iScol = n*r, i, j;
  double *Scolj, *Ai, *Scoli, *Scoldj=0, *Scoldi;
  DIFF = (nPar > 0);
  for (j=0; j<p+1; j++) {
    lacpy("All", r, r, S + j*r*r, r, Scol + j*r, iScol);
    if (DIFF) lacpy("All", r, r*nPar, Sd + j*r*r*nPar, r, Scold + j*r, iScol);    
  }
  for (j=p+1; j<n; j++) {
    Scolj = Scol + j*r;
    if (DIFF) Scoldj = Scold + j*r;
    if (j<=q) {
      lacpy("All", r, r, G + j*r*r, r, Scolj, iScol);
      if (DIFF) lacpy("All", r, r*nPar, Gd + j*r*r*nPar, r, Scoldj, iScol);
    }
    for (i=0; i<p; i++) {
      Ai = A+i*r*r;
      Scoli = Scolj - (i+1)*r;
      gemm("N", "N", r, r, r, 1.0, Ai, r, Scoli, iScol, 1.0, Scolj, iScol);
      if (DIFF) {
        Scoldi = Scoldj - (i+1)*r;
        AGdiff(Ai, i, Scoli, Scoldi, iScol, Scoldj, iScol, r, nPar);
      }
    }
  }
}

static void AGdiff(double Ai[], int i, double C[], double Cd[], int iC, 
                   double Fd[], int iF, int r, int nPar)
// Find F + Ai·C and its derivative.
{
  int k, c, l;
  gemm("N", "N", r, r*nPar, r, 1.0, Ai, r, Cd, iC, 1.0, Fd, iF);
  k = r*r*i;
  for (c=0; c<r; c++) {
    for (l=0; l<r; l++) {
      axpy(r, 1.0, C + c, iC, Fd + l + k*iF*r, iF);
      k++;
    }
  }
}
