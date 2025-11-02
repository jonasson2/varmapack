// src/VYW_gateway.c
#include <R.h>
#include <Rinternals.h>
#include "VYW.h"             // declares VYWFactorize, VYWSolve
#include "VarmaUtilities.h"  // if your headers need it
#include "BlasGateway.h"     // if your code uses BLAS via this gateway

/* C cores (already provided in your src):
   void VYWFactorize(double A[], double LU[], int piv[],
                     int p, int r, int *info);

   void VYWSolve(double A[], double LU[], double S[], double Y[],
                 int nrhs, int nY, int piv[], int p, int r);
*/

/* Gateway for factorization
   A_ : r x (r*p) matrix   OR   r x r x p array (either layout is fine)
   p_ : integer
   r_ : integer

   Returns a list(LU = N x N double matrix,
                  piv = integer(N),
                  info = integer(1)),
   where N = r*r*p - r*(r-1)/2.
*/
SEXP VYWFactorize_gateway(SEXP A_, SEXP p_, SEXP r_) {
    int p = asInteger(p_);
    int r = asInteger(r_);
    int N = r * r * p - r * (r - 1) / 2;

    SEXP LU_, piv_, info_;
    PROTECT(LU_   = allocMatrix(REALSXP, N, N));
    PROTECT(piv_  = allocVector(INTSXP, N));
    PROTECT(info_ = allocVector(INTSXP, 1));

    int info = 0;
    VYWFactorize(REAL(A_), REAL(LU_), INTEGER(piv_), p, r, &info);
    INTEGER(info_)[0] = info;

    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, LU_);
    SET_VECTOR_ELT(out, 1, piv_);
    SET_VECTOR_ELT(out, 2, info_);

    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, mkChar("LU"));
    SET_STRING_ELT(nms, 1, mkChar("piv"));
    SET_STRING_ELT(nms, 2, mkChar("info"));
    setAttrib(out, R_NamesSymbol, nms);

    UNPROTECT(5);
    return out;
}

/* Gateway for solve
   A_   : r x (r*p) matrix   OR   r x r x p array
   LU_  : N x N matrix from VYWFactorize
   piv_ : integer(N) from VYWFactorize
   Y_   : r x r x nY array   (nrhs is fixed to 1 here)
   p_, r_ : integers

   Returns S_ as an r x r x (p+1) array.
*/
SEXP VYWSolve_gateway(SEXP A_, SEXP LU_, SEXP piv_, SEXP Y_, SEXP p_, SEXP r_) {
    int p = asInteger(p_);
    int r = asInteger(r_);

    // Infer nY from Y_ (expects a 3D array r x r x nY)
    SEXP Ydim = getAttrib(Y_, R_DimSymbol);
    int nY = 1;
    if (!isNull(Ydim) && LENGTH(Ydim) == 3) {
        nY = INTEGER(Ydim)[2];
    }

    // Allocate S as r x r x (p+1) with a proper dim SEXP first
    SEXP dimS = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dimS)[0] = r;
    INTEGER(dimS)[1] = r;
    INTEGER(dimS)[2] = p + 1;
    SEXP S_ = PROTECT(allocArray(REALSXP, dimS));

    int nrhs = 1; // single RHS path (fits current varmapack needs)

    VYWSolve(REAL(A_), REAL(LU_), REAL(S_), REAL(Y_),
             nrhs, nY, INTEGER(piv_), p, r);

    UNPROTECT(2);
    return S_;
}
