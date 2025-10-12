# check_vyw_solution(A, B, Sig, verbose = TRUE)
# - Accepts A as r x r, r x (r*p), or r x r x p
# - Accepts B as r x r, r x (r*q), or r x r x q
# - Sig must be r x r
# - Uses varmasim:::find_CGW(), vyw_factorize(), vyw_solve()
# - Returns S plus residual norms; prints a summary when verbose=TRUE

vyw_prufa <- function(A, B, Sig, verbose = TRUE) {
  # --- basic checks ---
  stopifnot(is.matrix(Sig), nrow(Sig) == ncol(Sig))
  r <- nrow(Sig)
  
  to_array3 <- function(M, r) {
    d <- dim(M)
    if (length(d) == 3L) return(M)                          # already r x r x k
    stopifnot(is.matrix(M), nrow(M) == r)
    if (ncol(M) == r) {                                     # r x r  -> k = 1
      return(array(M, dim = c(r, r, 1L)))
    } else {                                                # r x (r*k) blocks
      k <- ncol(M) %/% r
      stopifnot(ncol(M) == r * k, k >= 1L)
      arr <- array(0, dim = c(r, r, k))
      for (j in seq_len(k)) arr[, , j] <- M[, ((j - 1L) * r + 1L):(j * r), drop = FALSE]
      arr
    }
  }
  
  A3 <- to_array3(A, r)
  B3 <- to_array3(B, r)
  p  <- dim(A3)[3L]
  q  <- dim(B3)[3L]
  
  # --- CGW via strict gateway (wrapped by your internal R function) ---
  # (If you didn't make the wrapper, replace with .Call("FindCGW_gateway", A3, B3, Sig, PACKAGE="varmasim"))
  CGW <- varmasim:::find_CGW(A3, B3, Sig)   # returns list(C, G, W)
  Y   <- CGW$G                               # r x r x (q+1)
  
  # --- Build [A1 ... Ap] and solve VYW ---
  A_block <- do.call(cbind, lapply(seq_len(p), function(k) A3[,,k, drop = FALSE][,,1]))
  PLU <- varmasim:::vyw_factorize(A_block, p = p, r = r)
  S   <- varmasim:::vyw_solve(A_block, PLU$LU, PLU$piv, Y, p = p, r = r)  # r x r x (p+1)
  
  # --- residuals for i = 0..p ---
  S_i <- function(i) S[,, i + 1L, drop = FALSE][,,1L]
  G_i <- function(i) Y[,, i + 1L, drop = FALSE][,,1L]
  
  res_norms <- numeric(p + 1L)
  for (i in 0:p) {
    # base term
    sum <- if (i <= q) S_i(i) - G_i(i) else S_i(i)
    
    # A * S terms (modified Yule–Walker structure)
    if (i == 0L) {
      # i = 0: subtract A_j %*% t(S_j) for j=1..p
      for (j in 1:p) sum <- sum - A3[,,j] %*% t(S_i(j))
    } else {
      for (j in 1:p) {
        Aj <- A3[,,j]
        if (j <  i) sum <- sum - Aj %*% S_i(i - j)
        if (j == i) sum <- sum - Aj %*% S_i(0L)
        if (j >  i) sum <- sum - Aj %*% t(S_i(j - i))
      }
    }
    res_norms[i + 1L] <- norm(sum, "F")
  }
  
  max_res <- max(res_norms)
  
  if (isTRUE(verbose)) {
    cat(sprintf("r=%d, p=%d, q=%d\n", r, p, q))
    cat("Residual ||·||_F per i=0..p: ", paste(sprintf("%.3e", res_norms), collapse = " "), "\n", sep = "")
    cat(sprintf("Max residual: %.3e\n", max_res))
  }
  
  invisible(list(
    r = r, p = p, q = q,
    A3 = A3, B3 = B3, Sig = Sig,
    CGW = CGW, G = Y, S = S,
    residuals = res_norms, max_residual = max_res
  ))
}
