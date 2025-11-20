# R interface scaffolding for randompack RNG helpers.
# These functions intentionally mirror the low-level C gateway names
# (rp_create, rp_free, rp_u01) but expose a friendlier R-facing surface.

#' Create a new randompack RNG
#' @param type Optional string selecting the RNG ("Park-Miller", "Xorshift", etc.)
#' @param seed Integer seed; defaults to 0 meaning "randomize"
#' @return An external pointer representing the RNG state.
create <- function(type = NULL, seed = 0L) {
  if (!is.null(type)) {
    stopifnot(is.character(type), length(type) == 1L, !is.na(type))
  }
  .Call("rp_create",
        if (is.null(type)) NULL else as.character(type),
        as.integer(seed),
        PACKAGE = "varmapack")
}

#' Explicitly free an RNG
#' @param rng External pointer previously returned by `create()`
free <- function(rng) {
  .Call("rp_free", rng, PACKAGE = "varmapack")
  invisible(rng)
}

#' Generate U(0,1) deviates from an RNG
#' @param rng External pointer from `create()`
#' @param n Number of variates to draw (non-negative integer)
#' @return Numeric vector of length `n`
u01 <- function(rng, n) {
  stopifnot(length(n) == 1L, is.finite(n))
  .Call("rp_u01", rng, as.integer(n), PACKAGE = "varmapack")
}
