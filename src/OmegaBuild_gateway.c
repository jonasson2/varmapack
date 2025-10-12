#include <R.h>
#include <Rinternals.h>
#include "Omega.h"          // declares OmegaBuild(...)
#include "VarmaUtilities.h" // whatever your headers need

// S_: r x r*(p+1), G_: r x r*(q+1), W_: r x r*(q+1)
// returns list(Su = (h*r) x (h*r), Olow = ((n-h)*r) x ((q+1)*r))
SEXP OmegaBuild_gateway(SEXP S_, SEXP G_, SEXP W_, SEXP p_, SEXP q_, SEXP r_, SEXP n_) {
    int p = asInteger(p_);
    int q = asInteger(q_);
    int r = asInteger(r_);
    int n = asInteger(n_);
    int h = (p > q) ? p : q;

    SEXP Su_, Olow_;
    PROTECT(Su_   = allocMatrix(REALSXP, h * r, h * r));
    PROTECT(Olow_ = allocMatrix(REALSXP, (n - h) * r, (q + 1) * r));

    // Call the underlying C routine. R matrices are column-major (Fortran order),
    // so REAL() gives the right contiguous memory layout for your routine.
    OmegaBuild(REAL(Su_), REAL(Olow_),
               REAL(S_), REAL(G_), REAL(W_),
               p, q, r, n);

    SEXP out = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, Su_);
    SET_VECTOR_ELT(out, 1, Olow_);

    SEXP names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(names, 0, mkChar("Su"));
    SET_STRING_ELT(names, 1, mkChar("Olow"));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(4);
    return out;
}
