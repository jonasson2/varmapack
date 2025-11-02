#' Generate VARMA test cases (named or random) via C backend
#'
#' @param name Character or NULL. Predefined case name
#'   (one of: "tinyAR","tinyMA","tinyARMA","smallAR","smallMA","smallARMA1",
#'   "smallARMA2","mediumAR","mediumARMA1","mediumARMA2","mediumMA","largeAR").
#' @param p,q,r Integers or NULL. Dimensions for an unnamed random case.
#' @param all Logical. If TRUE, return a named list of all predefined cases.
#' @param seed Optional integer RNG seed (unnamed cases).
#' @return One case as list(A,B,Sig,p,q,r,name), or (if all=TRUE) a named list of such lists.
#' @export
testcase <- function(name = NULL,
                     p = NULL, q = NULL, r = NULL,
                     all = FALSE,
                     seed = NULL) {

  # Canonical set of named testcases (must stay in sync with C's mapping)
  named_cases <- c(
    "tinyAR","tinyMA","tinyARMA",
    "smallAR","smallMA","smallARMA1","smallARMA2",
    "mediumAR","mediumARMA1","mediumARMA2","mediumMA",
    "largeAR"
  )

  # Exclusivity checks
  if (isTRUE(all) && (!is.null(name) || !is.null(p) || !is.null(q) || !is.null(r))) {
    stop("`all = TRUE` cannot be combined with `name` or `p/q/r`.", call. = FALSE)
  }
  if (!is.null(name) && (!is.null(p) || !is.null(q) || !is.null(r))) {
    stop("Supply either `name` OR `p, q, r` (not both).", call. = FALSE)
  }

  seed <- if (is.null(seed)) 42L else as.integer(seed)

  # Return all named cases
  if (isTRUE(all)) {
    out <- lapply(named_cases, function(nm) {
      .Call("Testcase_gateway",
            nm,                # name
            NULL, NULL, NULL,  # p,q,r unused for named
            FALSE,             # park_miller
            seed,
            PACKAGE = "varmapack")
    })
    names(out) <- named_cases
    return(out)
  }

  # Single named case
  if (!is.null(name)) {
    name <- match.arg(name, choices = named_cases)
    return(.Call("Testcase_gateway",
                 name,
                 NULL, NULL, NULL,
                 FALSE,
                 seed,
                 PACKAGE = "varmapack"))
  }

  # Random case: validate p,q,r; delegate to C (which handles stationarity)
  stopifnot(!is.null(p), !is.null(q), !is.null(r))
  stopifnot(.is_integer(p), p >= 0)
  stopifnot(.is_integer(q), q >= 0)
  stopifnot(.is_integer(r), r >= 1)

  tc <- .Call("Testcase_gateway",
              NULL,                          # name = NULL => random
              as.integer(p), as.integer(q), as.integer(r),
              FALSE,                         # park_miller
              seed,
              PACKAGE = "varmapack")

  # Normalize name
  tc$name <- sprintf("random(p=%d,q=%d,r=%d)", tc$p, tc$q, tc$r)
  tc
}

.is_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && floor(x) == x
}
