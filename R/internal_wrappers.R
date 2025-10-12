# Internal C gateways
# These functions wrap compiled routines and are used only for testing
# and debugging. They are not exported, documented, or part of the public API.
#' @keywords internal
#' @noRd

omega_build <- function(S, G, W, p, q, r, n) {
  .Call("OmegaBuild_gateway", S, G, W,
        as.integer(p), as.integer(q), as.integer(r), as.integer(n),
        PACKAGE = "varmasim")
}

vyw_factorize <- function(A, p, r) {
  out <- .Call("VYWFactorize_gateway", A, as.integer(p), as.integer(r),
               PACKAGE = "varmasim")
  if (out$info != 0L) stop("VYWFactorize: singular system (info=", out$info, ")")
  out
}

vyw_solve <- function(A, LU, piv, Y, p, r) {
  .Call("VYWSolve_gateway", A, LU, as.integer(piv), Y,
        as.integer(p), as.integer(r), PACKAGE = "varmasim")
}

find_CGW <- function(A, B, Sig) {
  .Call("FindCGW_gateway", A, B, Sig, PACKAGE = "varmasim")
}