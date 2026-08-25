#!/usr/bin/env Rscript

# Compare Varmapack with MTS and ts.extend simulation on the shared testcases.
# Run from r-package/ with: Rscript examples/TimeSimulate.R

ALL_MODELS <- c(
  "tinyAR", "tinyMA", "tinyARMA", "smallAR1", "smallAR2", "smallMA1",
  "smallMA2", "smallARMA1", "smallARMA2", "mediumAR", "mediumMA1",
  "mediumARMA1", "mediumARMA2", "mediumMA2", "largeAR", "largeARMA"
)
PAPER_MODELS <- c(
  "tinyAR", "tinyARMA", "smallAR1", "smallARMA2", "mediumAR",
  "mediumARMA2", "largeAR", "largeARMA"
)
PLATFORM_MODELS <- c(
  "tinyARMA", "smallAR2", "mediumAR", "mediumARMA1", "largeAR"
)

print_help <- function() {
  cat("TimeSimulate.R -- compare Varmapack, MTS::VARMAsim, and ts.extend\n\n")
  cat("Usage: Rscript examples/TimeSimulate.R [options]\n\n")
  cat("MTS uses its default skip = 200 for every simulated path.\n\n")
  cat("Options:\n")
  cat("  -h          show this help\n")
  cat("  -t seconds  timing target per method and testcase (default 0.1)\n")
  cat("  -n length   length of each returned series (default 100)\n")
  cat("  -M count    replicates per Varmapack/ts.extend call (default 1000)\n")
  cat("  -5          time the five selected platform models only\n")
  cat("  -8          time the eight selected paper models only\n")
}

parse_options <- function(args) {
  opts <- list(t = 0.1, n = 100L, M = 1000L, selected = FALSE, platform = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "-h") {
      print_help()
      return(NULL)
    }
    if (arg == "-8") {
      opts$selected <- TRUE
      i <- i + 1L
      next
    }
    if (arg == "-5") {
      opts$platform <- TRUE
      i <- i + 1L
      next
    }
    if (startsWith(arg, "-t") && nchar(arg) > 2L) {
      opts$t <- as.numeric(substr(arg, 3L, nchar(arg)))
      i <- i + 1L
      next
    }
    if (!(arg %in% c("-t", "-n", "-M")) || i == length(args))
      stop("invalid options; use -h for help", call. = FALSE)
    value <- args[[i + 1L]]
    if (arg == "-t") opts$t <- as.numeric(value)
    if (arg == "-n") opts$n <- as.integer(value)
    if (arg == "-M") opts$M <- as.integer(value)
    i <- i + 2L
  }
  if (!is.finite(opts$t) || opts$t <= 0 || is.na(opts$n) || opts$n < 1L ||
      is.na(opts$M) || opts$M < 1L)
    stop("-t, -n, and -M must be positive", call. = FALSE)
  if (opts$selected && opts$platform)
    stop("-5 and -8 cannot be combined", call. = FALSE)
  opts
}

mts_terms <- function(terms, sign) {
  if (is.null(terms)) return(NULL)
  do.call(cbind, lapply(seq_len(dim(terms)[3]), function(i) sign*terms[, , i]))
}

mts_specification <- function(model) {
  # MTS subtracts its MA terms, whereas Varmapack adds them.
  list(
    arlags = if (model$p == 0L) NULL else seq_len(model$p),
    malags = if (model$q == 0L) NULL else seq_len(model$q),
    phi = mts_terms(model$A, 1),
    theta = mts_terms(model$B, -1),
    sigma = model$Sig
  )
}

mts_simulate <- function(specification, n) {
  MTS::VARMAsim(
    nobs = n,
    arlags = specification$arlags,
    malags = specification$malags,
    phi = specification$phi,
    theta = specification$theta,
    sigma = specification$sigma
  )$series
}

time_simulation <- function(draw, values, target) {
  result <- draw()
  if (!all(is.finite(result))) stop("simulation returned non-finite values")
  calls <- 0L
  start <- proc.time()[["elapsed"]]
  elapsed <- 0
  while (elapsed < target) {
    draw()
    calls <- calls + 1L
    elapsed <- proc.time()[["elapsed"]] - start
  }
  1e9*elapsed/(calls*values)
}

time_tsextend <- function(model, opts) {
  ar <- if (is.null(model$A)) numeric() else as.numeric(model$A)
  ma <- if (is.null(model$B)) numeric() else as.numeric(model$B)
  values <- opts$n*opts$M
  set.seed(123L)
  time_simulation(
    function() ts.extend::rGARMA(
      n = opts$M, m = opts$n, errorvar = as.numeric(model$Sig), ar = ar, ma = ma
    ),
    values,
    opts$t
  )
}

time_case <- function(name, model, opts) {
  specification <- mts_specification(model)
  values <- opts$n*opts$M*model$r
  rng <- randompack::randompack_rng()
  rng$seed(123L)
  varmapack_ns <- time_simulation(
    function() model$sim(opts$n, nrep = opts$M, rng = rng), values, opts$t
  )
  set.seed(123L)
  mts_ns <- time_simulation(
    function() mts_simulate(specification, opts$n), opts$n*model$r, opts$t
  )
  tsextend_ns <- if (model$r == 1L) time_tsextend(model, opts) else NA_real_
  rho <- if (model$p == 0L) 0 else model$specrad()
  tsextend_text <- if (is.na(tsextend_ns)) "-" else sprintf("%.0f", tsextend_ns)
  cat(sprintf("%-12s %2d %2d %2d  %5.3f %10.1f %7.0f %9s %6.0f\n", name,
              model$p, model$q, model$r, rho, varmapack_ns, mts_ns, tsextend_text,
              mts_ns/varmapack_ns))
}

select_models <- function(opts) {
  if (opts$selected) return(PAPER_MODELS)
  if (opts$platform) return(PLATFORM_MODELS)
  ALL_MODELS
}

print_setting <- function(label, value) {
  cat(sprintf("%-45s%s\n", label, as.character(value)))
}

main <- function(args) {
  opts <- parse_options(args)
  if (is.null(opts)) return(invisible())
  if (!requireNamespace("MTS", quietly = TRUE))
    stop("install the MTS package before running this benchmark", call. = FALSE)
  if (!requireNamespace("ts.extend", quietly = TRUE))
    stop("install the ts.extend package before running this benchmark", call. = FALSE)
  cat("Varmapack R benchmark and comparison with MTS::VARMAsim and ts.extend\n")
  print_setting("Benchmark unit:", "ns/retained value")
  print_setting("Length per series:", opts$n)
  print_setting("Replicates per Varmapack/ts.extend call:", opts$M)
  print_setting("Replicates per MTS call:", 1)
  print_setting("Target time per method and testcase:", sprintf("%.3g s", opts$t))
  if (opts$selected) print_setting("Models:", "selected paper set")
  if (opts$platform) print_setting("Models:", "selected platform set")
  print_setting("MTS discarded values per path:", 200)
  cat("\n")
  cat(sprintf("%-12s %2s %2s %2s  %5s %10s %7s %9s %6s\n", "Testcase", "p", "q", "r",
              "rho", "Varmapack", "MTS", "ts.extend", "MTS/VP"))
  for (name in select_models(opts)) {
    time_case(name, varmapack::varmapack_testcase(name), opts)
  }
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
