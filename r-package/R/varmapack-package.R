#' Varmapack: Exact-Start Simulation of Gaussian VARMA Models
#'
#' An R interface to the Varmapack C library for simulation and analysis of
#' Gaussian VAR, VMA, VARMA, and VARMAX models.
#'
#' @seealso
#' - `vignette("getting-started", package = "varmapack")` for a tutorial.
#' - `vignette("mathematical-description", package = "varmapack")` for the
#'   mathematical description.
#' - [Varmapack on GitHub](https://github.com/jonasson2/varmapack).
#' - [Report bugs](https://github.com/jonasson2/varmapack/issues).
#'
#' @importFrom R6 R6Class
#' @import randompack
#' @useDynLib varmapack, .registration = TRUE
#'
#' @keywords internal
"_PACKAGE"
