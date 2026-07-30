varmapack
================

`varmapack` simulates and analyzes Gaussian VAR, VMA, VARMA, and VARMAX
time-series models. Simulated stationary series have the correct
distribution from the first returned term, without discarding a burn-in
segment. The package also provides model testcases, theoretical and
sample autocovariances, spectral radii, and impulse response functions.

## Documentation

The [Getting Started with Varmapack](vignettes/getting-started.Rmd)
vignette shows model construction, simulation, testcases, and analysis.
The [Mathematical Description of
Varmapack](vignettes/mathematical-description.Rmd) vignette defines the
supported models and describes the simulation method.

After installing the CRAN release, open the vignettes with
`vignette("getting-started", package = "varmapack")` and
`vignette("mathematical-description", package = "varmapack")`. The
installed reference documentation is available through
`help(package = "varmapack")`. For information about the underlying C
library, see the [C
README](https://github.com/jonasson2/varmapack/blob/main/README.md).

## Installation

The package will be installable from CRAN when released:

``` r
install.packages("varmapack")
```

For the development version, install the companion `randompack` package
first:

``` r
remotes::install_github("jonasson2/randompack", subdir = "r-package")
remotes::install_github("jonasson2/varmapack", subdir = "r-package")
```
