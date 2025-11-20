// src/vpack_FindCG_gateway.c
#include <R.h>
#include <Rinternals.h>
#include "BlasGateway.h"
#include "allocate.h"
#include "VarmaMisc.h"

// void vpack_FindCG(double A[], double B[], double Sig[], int p, int q, int r,
//              double C[], double G[], double W[]);

SEXP vpack_FindCG_gateway(SEXP A_, SEXP B_, SEXP Sig_) {
    // Require: A = r x r x p, B = r x r x q, Sig = r x r
    SEXP Adim = getAttrib(A_, R_DimSymbol);
    SEXP Bdim = getAttrib(B_, R_DimSymbol);
    SEXP Sdim = getAttrib(Sig_, R_DimSymbol);
    //Rprintf("LENGTH(Adim)=%d\n", LENGTH(Adim));

    if (TYPEOF(Adim) != INTSXP || LENGTH(Adim) != 3)
        Rf_error("vpack_FindCG_gateway: A must be a 3D array with dims (r, r, p).");
    if (TYPEOF(Bdim) != INTSXP || LENGTH(Bdim) != 3)
        Rf_error("vpack_FindCG_gateway: B must be a 3D array with dims (r, r, q).");
    if (TYPEOF(Sdim) != INTSXP || LENGTH(Sdim) != 2)
        Rf_error("vpack_FindCG_gateway: Sig must be a 2D matrix with dims (r, r).");

    int rA0 = INTEGER(Adim)[0], rA1 = INTEGER(Adim)[1], p = INTEGER(Adim)[2];
    int rB0 = INTEGER(Bdim)[0], rB1 = INTEGER(Bdim)[1], q = INTEGER(Bdim)[2];
    int rS0 = INTEGER(Sdim)[0], rS1 = INTEGER(Sdim)[1];

    if (rA0 != rA1) Rf_error("vpack_FindCG_gateway: A is not square in first two dims.");
    if (rB0 != rB1) Rf_error("vpack_FindCG_gateway: B is not square in first two dims.");
    if (rS0 != rS1) Rf_error("vpack_FindCG_gateway: Sig is not square.");
    int r = rA0;
    if (rB0 != r || rS0 != r)
        Rf_error("vpack_FindCG_gateway: A, B, Sig must share the same 'r'.");

    // Allocate outputs: r x r x (q+1)
    SEXP dim3 = PROTECT(allocVector(INTSXP, 3));
    INTEGER(dim3)[0] = r;
    INTEGER(dim3)[1] = r;
    INTEGER(dim3)[2] = q + 1;

    SEXP C_ = PROTECT(allocArray(REALSXP, dim3));
    SEXP G_ = PROTECT(allocArray(REALSXP, dim3));

    // Call the core
    vpack_FindCG(REAL(A_), REAL(B_), REAL(Sig_), p, q, r,
            REAL(C_), REAL(G_));

    // Return list(C, G)
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, C_);
    SET_VECTOR_ELT(out, 1, G_);

    SEXP nms = PROTECT(allocVector(STRSXP, 3));
    SET_STRING_ELT(nms, 0, mkChar("C"));
    SET_STRING_ELT(nms, 1, mkChar("G"));
    setAttrib(out, R_NamesSymbol, nms);

    UNPROTECT(6); // dim3, C_, G_, out, nms
    return out;
}
