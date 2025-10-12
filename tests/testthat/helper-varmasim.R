# tests/testthat/helper-varmasim.R

# block-matrix helpers
tril   <- function(M) M * lower.tri(M, diag = TRUE)
blk    <- function(...) do.call(cbind, list(...))   # [A B]
rowblk <- function(...) do.call(rbind, list(...))   # [A; B]
split_blocks <- function(Mcat, r) {
  # Mcat is r x (k*r), return list of k blocks (each r x r)
  k <- ncol(Mcat) %/% r
  lapply(seq_len(k), function(j) Mcat[, ((j-1)*r + 1):(j*r), drop = FALSE])
}

# Residual checker like MATLAB find_max_residual
# A_list: {A1..Ap}, S_list: {S0..Sp}, G_list: {G0..Gq}

normF <- function(x) {
  v <- as.numeric(x)          # works for scalar/vector/matrix/array
  sqrt(sum(v * v))
}

max_vyw_residual <- function(A, S, G) {
  p <- if (length(dim(A)) == 3L) dim(A)[3L] else 0L
  q <- if (length(dim(G)) == 3L) dim(G)[3L] - 1L else 0L
  S0 <- S[,,1L]
  maxres <- 0
  for (i in 0:p) {
    sum <- if (i == 0L) {
      S0 - G[,,1L]
    } else if (i <= q) {
      S[,,i + 1L] - G[,,i + 1L]
    } else {
      S[,,i + 1L]
    }
    for (j in seq_len(p)) {
      Aj <- A[,,j]
      if (j <  i) sum <- sum - Aj %*% S[,,(i - j) + 1L]       # A_j * S_{i-j}
      if (j == i) sum <- sum - Aj %*% S0                      # A_i * S0
      if (j >  i) sum <- sum - Aj %*% t(S[,,(j - i) + 1L])    # A_j * S_{j-i}^T
    }
    cur <- normF(sum)
    maxres <- max(maxres, cur)
  }
  maxres
}
# CHECK VYW FOR A RANGE OF DIMENSIONS
check_vyw_grid <- function() {
  combos <- expand.grid(p = 0:3, q = 0:3, r = 1:3)
  for (k in seq_len(nrow(combos))) {
    p <- combos$p[k]; q <- combos$q[k]; r <- combos$r[k]
    cat(sprintf("\np, q, r = %d, %d, %d", p, q, r))
    check_vyw_random(p, q, r)
  }
}  

# CHECK VYW FOR SPECIFIC DIMENSIONS
check_vyw_random <- function(p, q, r, tol = 1e-10,
                             seed = 1 + p + 10*q + 100*r) {
  set.seed(seed)
  A <- if (p == 0) NULL else array(replicate(p, runif(r*r, 0, 1/(4*p*r)), r*r),
                                   dim = c(r, r, p))
  B <- if (q == 0) NULL else array(replicate(q, runif(r*r, 0, 1/(4*q*r)), r*r),
                                   dim = c(r, r, q))
  Sig <- crossprod(matrix(runif(r*r, -0.2, 0.2), r, r)) + diag(r)
  cat('A, Sig = ', A, Sig)
  check_vyw(A, B, Sig, tol = tol)
}

# CHECK VYW FOR A SPECIFIC CASE (ALLOWS COMPARISON WITH MATLAB)
check_vyw <- function(A, B, Sig, tol = 1e-10, tc = NULL) {
  if (!is.null(tc)) {
    if (!is.list(tc) || is.null(tc$A) || is.null(tc$B) || is.null(tc$Sig)) {
      stop("tc must be a list with components A, B, and Sig.")
    }
    A   <- tc$A
    B   <- tc$B
    Sig <- tc$Sig
  }
  r <- nrow(Sig)
  stopifnot(is.matrix(Sig), ncol(Sig) == r)
  
  A <- make_array(A, r)
  B <- make_array(B, r)
  p <- dim(A)[3]
  q <- dim(B)[3]
  
  CGW <- find_CGW(A, B, Sig)
  G   <- CGW$G
  
  if (p == 0L) {
    S <- array(G[,,1, drop = FALSE], dim = c(r, r, 1L))
    maxres <- max_vyw_residual(A, S, G)
    expect_true(maxres < tol,
                sprintf("p=%d q=%d r=%d maxres=%.3e", p, q, r, maxres))
    return(invisible(TRUE))
  }
  
  PLU <- vyw_factorize(A, p, r)
  cat('\nPLU = ', PLU$LU)
  expect_identical(length(PLU), 3L)
  expect_true(is.matrix(PLU$LU))
  expect_true(is.integer(PLU$piv))
  expect_identical(as.integer(PLU$info), 0L)
  
  S <- vyw_solve(A, PLU$LU, PLU$piv, G, p = p, r = r)
  cat('\nS = ', S)
  expect_true(is.array(S))
  expect_identical(dim(S),
                   c(as.integer(r), as.integer(r), as.integer(p) + 1L))
  
  maxres <- max_vyw_residual(A, S, G)
  expect_true(maxres < tol,
              sprintf("p=%d q=%d r=%d maxres=%.3e", p, q, r, maxres))
  
  invisible(TRUE)
}