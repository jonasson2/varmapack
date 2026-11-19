# varmapack 0.2.0

* Added inverse covariance recovery with `varmapack_cov2sig()`.
* Added stationary OU and Slepian-based models that generalize VAR(1) and
  VMA(1) to unevenly spaced observation times, plus a Gamma-spaced time
  generator.
* Added `acvf(lags)` at fractional lags for both models and `specrad(gap)`
  for the spectral radius of an OU model sampled at a fixed gap.

# varmapack 0.1.2

- Added an empirical example using European stock-index returns and expanded
  repository statistical tests, benchmark analysis, and documentation following
  SoftwareX review.
- Improved portability for builds using link-time optimization (LTO), as
  requested by CRAN.

# varmapack 0.1.1

- Minor maintenance changes for CRAN submission checks. No API or numerical
  changes.
