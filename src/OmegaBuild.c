#include "BlasGateway.h"
#include "VarmaUtilities.h"
#include "Omega.h"
void OmegaBuild ( // Build Omega = cov(w) with w' = (w1'...wn')
  double Su[],   // out  h·r × h·r, upper left part of Omega, where h = max(p,q)
  double Olow[], // out  (n-h)r × (q+1)r, block-diagonals of lower part of Omega
  double S[],    // in   r × r·(p+1), = [S0 S1...Sp] (lagged) covariances of xt
  double G[],    // in   r × r·(q+1), = [G0...Gq], Gi = cov(y(t),y(t-i))
  double W[],    // in   r × r·(q+1), = [W0...Wq], Wi = cov(y(t),x(t-i))
  int p,         // in   number of autoregressive terms
  int q,         // in   number of moving average terms
  int r,         // in   dimension of xt
  int n)         // in   length of series
{
  // Builds the covariance matrix Omega of w = (w1'...wn')'. Omega is returned
  // in Su and Olow. All matrices are stored in Fortran fashion. See matrix
  // structures in comments at end of routine.
  //
  int j, m, k, t, h, rh;
  double *s, *olow;
  h = max(p,q);
  rh = r*h;    // leading dimension of Su
  m = r*(n-h); // leading dimension of Olow
  //
  // BUILD Su:
  k = 0;
  for (t=0; t<=h-1; t++) {
    // copy to t-th block-subdiagonal:
    s = Su + t*r;    
    for (j=0; j<h-t; j++) {
      if      (j<p-t) lacpy("All", r, r, S+k, r, s, rh);  // Sj
      else if (j<p  ) lacpy("All", r, r, G+k , r, s, rh);  // Gj
      else            lacpy("All", r, r, W+k , r, s, rh);  // Wj
      s += r*rh + r;
    }
    k += r*r;
  }
  // BUILD Olow:
  for (t=h; t<n; t++) {
    olow = Olow + r*(t-h);
    k = r*r*q;
    for (j=t-q; j<=t; j++) { // set Omega(t,j)
      if (j< p) lacpy("All", r, r, G+k, r, olow, m);
      if (j>=p) lacpy("All", r, r, W+k, r, olow, m);
      k-= r*r;
      olow += r*m;
    }
  }
}
// The following shows the structure of Omega for q < p:
//
//  S0 S1'...      Sp-1'|              transpose
//  S1 S0             : |               of lower
//  S2 S1 S0            |                   part    r·p
//   :        ...    S1'|
//  Sp-1             S0 |
//  --------------------+-----------------------
//        Gq ...  G2 G1 | W0           transpose
//           Gq      G2 | W1  W0        of lower
//              ..   G3 | W2  W1  W0        part
//                ..  : |  :   :    ..
//                   Gq | Wq-1        ..            r·(n-p)
//                      | Wq            ..
//                      |     Wq          ..
//                      |        ...        ..
//                      |            Wq ..... W0
//
// For q >= p the lower left, lower right and upper right partitions are as
// above but the upper left partition is:
//
//  S0 S1'...   Sp-1' Gp'.....Gq-1'
//  S1 S0         :    :         : 
//  S2     ..     :    :         :  
//   :        .. S1'  G2'        :  
//  Sp-1      S1 S0   G1'...  Gq-p'   
//  Gp        G2 G1   W0... Wq-p-1'
//  :             :    :  ..     :
//  :             :    :     ..  :
//  Gq-1...     Gq-p Wq-p-1.....W0
//
// The Olow matrix is r·(n-p)×r(q+1) and its structure is:
//
//   Gq ... G2 G1 W0
//   Gq     G2 W1 W0
//   Gq     W2 W1 W0
//   :             :     r·(n-p)
//   Gq Wq-1 ...  W0
//   Wq           W0
//   :             :
//   Wq ...    W1 W0
//         rq
