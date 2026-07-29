# Possible future work

- Support fixed-start VARMA and VARMAX simulation when the relevant covariance
  is singular positive semidefinite. This requires a rank-revealing
  factorization, conditioning within the covariance support, and rejection of
  supplied values outside that support.
- Compute the log-likelihood of observed series under a supplied model.
- Fit VARMA and VARMAX models by maximum likelihood.
