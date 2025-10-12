#' Simulate VARMA Time Series
#'
#' This function simulates Vector Autoregressive Moving Average (VARMA) time series 
#' using a compiled C function (`VarmaSim`) to ensure high performance and avoid spin-up bias.
#'
#' @details
#' The function interfaces with the C function `VarmaSim`, which generates simulated
#' AR, VAR, ARMA, or VARMA time series with the following parameters:
#' \describe{
#'   \item{A}{A 3D array of dimensions r by r by p for autoregressive parameter matrices.}
#'   \item{B}{A 3D array of dimensions r by r by q for moving average parameter matrices.}
#'   \item{Sig}{An r by r covariance matrix of the shock terms \eqn{\epsilon(t)}.}
#'   \item{mu}{An r-vector specifying the mean of \eqn{x(t)} (or `NULL` for zero-mean series).}
#'   \item{p}{The number of autoregressive terms.}
#'   \item{q}{The number of moving average terms.}
#'   \item{r}{The dimension of each \eqn{x(t)}.}
#'   \item{n}{The length of each generated series.}
#'   \item{M}{The number of replicate series to generate.}
#'   \item{X0}{An r by h matrix, where \eqn{h = \max(p, q)}, providing initial values for the series
#'             (or `NULL` to start from default initial values).}
#'   \item{X}{An r by n by M array containing the generated series (output).}
#'   \item{eps}{An r by n by M array containing the shock series (output).}
#'   \item{ok}{A flag indicating success: `1` if the series is stationary or `X0` is provided; otherwise, `0`.}
#' }
#'
#' The R function performs validation of inputs before passing them to the C function.
#' It also ensures compatibility by converting 2D matrices to 3D arrays where required.
#'
#' @param A A 3D array of dimensions r by r by p for autoregressive parameter
#'   matrices. If a 2D matrix is provided, it will be converted to a 3D array.
#'   Can be specified as NULL for a pure MA model.
#' @param B A 3D array of dimensions r by r by q for moving average parameter
#'   matrices. If a 2D matrix is provided, it will be converted to a 3D array.
#'   Can be specified as NULL for a pure AR model.
#' @param Sig An r by r covariance matrix of the shock terms \eqn{\epsilon(t)}.
#' @param mu An r-vector specifying the mean of \eqn{x(t)}, or `NULL` for
#'   zero-mean series.
#' @param n A positive integer specifying the length of each generated series.
#' @param M A positive integer specifying the number of replicate series to
#'   generate.
#' @param X0 An optional r by h matrix (where \eqn{h = \max(p, q)}) providing
#'   initial values for the series, or `NULL` to start from default initial
#'   values.
#' @param seed An optional positive integer used to seed random number
#'   generation (default NULL to omit seeding).
#' @param park_miller A flag indicating whether the Park Miller random number
#'   generator (RNG) should be used. If FALSE (default), the RNG built into R is
#'   used, but can be set to TRUE to allow comparison with results obtained with
#'   the sister Matlab package vauto, published as Algorithm 878 in ACM TOMS.
#' 
#' @return A list with the following components:
#' \describe{
#'   \item{X}{An r by n by M array containing the generated VARMA time series.}
#'   \item{eps}{An r by n by M array containing the generated shock series.}
#'   When M = 0 X and eps will be r by n matrices.
#' }
#'
#' @examples
#' # Example usage of varma_sim
#' A <- array(c(0.5, -0.2, 0.1, 0.3), dim = c(2, 2, 1))  # AR parameters
#' B <- array(0, dim = c(2, 2, 1))                       # MA parameters
#' Sig <- diag(2)                                        # Covariance matrix
#' mu <- c(0, 0)                                         # Mean
#' n <- 100                                              # Series length
#' M <- 1                                                # Number of replicates
#' result <- varma_sim(A, B, Sig, mu, n, M)
#' str(result$X)  # View the generated series
#'
#' @useDynLib varmasim
#' @importFrom Rcpp sourceCpp
#' @export

varma_sim <- function(A, B, Sig, n, mu = NULL, M = 1, X0 = NULL,
                      seed = NULL, park_miller = FALSE) {
  library(assertthat)
  assert_that(is.matrix(Sig), nrow(Sig) == ncol(Sig))
  r <- nrow(Sig)
  assert_that(is.count(r))
  
  A <- make_array(A, r)
  B <- make_array(B, r)

  assert_that(is.array(A), identical(dim(A)[1:2], c(r, r)))
  assert_that(is.array(B), identical(dim(B)[1:2], c(r, r)))
  
  p = dim(A)[3]
  q = dim(B)[3]

  assert_that(is.count(n), is.count(M))

  if (!is.null(mu))
    assert_that(is.numeric(mu), is.vector(mu), length(mu)==r)

  if (!is.null(X0))
    assert_that(is.matrix(X0), ncol(X0) == max(p,q), nrow(X0) == r)

  if (!is.null(seed))
    assert_that(is.count(seed))
  
  assert_that(is.logical(x), length(x) == 1, !is.na(x))
  
  # Preallocate results
  X <- array(0, dim = if (M==1) c(r, n) else c(r, n, M))
  print('dim(X),X:')
  print(dim(X))
  print(X)
  eps <- array(0, dim = if (M==1) c(r, n) else c(r, n, M))
  ok <- integer(1)  # status flag (e.g., success = 1)

  .Call("VarmaSim_gateway", A, B, Sig, mu, p, q, r, n, M, X0, X, eps, ok, 
        seed, park_miller)
  
  list(X = X, eps = eps, ok = ok)
}

make_array <- function(M, r) {
  if (is.null(M)) M <- array(M, c(r, r, 0))
  if (is.matrix(M)) M <- array(M, c(r, r, 1))
  M
}
