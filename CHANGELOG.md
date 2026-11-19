# Varmapack Changelog

## 0.2.0 (unreleased)

- Added stationary Ornstein-Uhlenbeck (OU) and Slepian-based models that
  generalize VAR(1) and VMA(1), respectively, to unevenly spaced observation times.
  Also added a Gamma-spaced time generator.
- Added theoretical autocovariance at fractional lags for both models and the
  sampled VAR(1) spectral radius at a specified gap for OU.
- Added `cov2sig` to recover innovation covariance from lag-zero process
  covariance and VARMA coefficients.

## 0.1.2

- Added a European stock-index returns example and expanded statistical
  tests, benchmarks, and documentation following SoftwareX review.
- Fixed build portability issues with link-time optimization reported by CRAN.

## 0.1.1

- Updated the bundled Fortran and documentation for CRAN resubmission.
  No API or numerical changes.

## 0.1.0

- Initial release: burn-in-free simulation of Gaussian VARMA models,
  VARMAX simulation, and tools for autocovariances, spectral radii,
  and impulse responses. C core with Python, R, and MATLAB interfaces.
