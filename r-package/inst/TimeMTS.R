set.seed(1)

n <- 20
burn <- 100

Sigma <- matrix(c(1.0, 0.3,
                  0.3, 2.0), 2, 2)

A1 <- matrix(c(0.5, 0.1,
               0.0, 0.4), 2, 2, byrow=TRUE)

B1 <- matrix(c(0.3, 0.0,
               0.1, 0.2), 2, 2, byrow=TRUE)

cat("\n--- MTS::VARMAsim ---\n")
if (requireNamespace("MTS", quietly=TRUE)) {
  Phi <- array(A1, dim=c(2, 2, 1))
  Theta <- array(B1, dim=c(2, 2, 1))
  z <- MTS::VARMAsim(n, arlags=1, malags=1, phi=Phi, theta=Theta,
                     sigma=Sigma, skip=burn)
  print(head(z$series))
  print(head(z$noises))
}
else {
  cat("Package MTS not installed\n")
}

cat("\n--- portes::varima.sim ---\n")
if (requireNamespace("portes", quietly=TRUE)) {
  phi <- array(A1, dim=c(2, 2, 1))
  theta <- array(B1, dim=c(2, 2, 1))
  z <- portes::varima.sim(n=n, k=2, ar=phi, ma=theta, sigma=Sigma)
  print(head(z))
}
else {
  cat("Package portes not installed\n")
}

cat("\n--- multiwave::varma ---\n")
if (requireNamespace("multiwave", quietly=TRUE)) {
  VAR <- array(A1, dim=c(2, 2, 1))
  VMA <- array(B1, dim=c(2, 2, 1))
  z <- multiwave::varma(N=n, k=2, VAR=VAR, VMA=VMA, cov_matrix=Sigma)
  print(head(z))
}
else {
  cat("Package multiwave not installed\n")
}

cat("\n--- beyondWhittle::sim_varma ---\n")
if (requireNamespace("beyondWhittle", quietly=TRUE)) {
  model <- list(ar=array(A1, dim=c(2, 2, 1)),
                ma=array(B1, dim=c(2, 2, 1)),
                sigma=Sigma)
  z <- try(beyondWhittle::sim_varma(model, n=n, d=2, burnin=burn),
           silent=TRUE)
  print(z)
}
else {
  cat("Package beyondWhittle not installed\n")
}
