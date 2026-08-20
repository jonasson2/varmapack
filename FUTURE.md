# Possible future work

- Support fixed-start VARMA and VARMAX simulation when the relevant covariance
  is singular positive semidefinite. This requires a rank-revealing
  factorization, conditioning within the covariance support, and rejection of
  supplied values outside that support.
- For Varmapack 2, add a bit-exact conditional-MVN mode with deterministic
  rank selection, factorization, matrix kernels, and floating-point settings
  across supported platforms and compilers.
- Compute the log-likelihood of observed series under a supplied model.
- Fit VARMA and VARMAX models by maximum likelihood.
- Use an iterative ARPACK spectral-radius computation for sufficiently large
  companion matrices, retaining LAPACK DGEEV for small models and as a fallback.
