#!/usr/bin/env Rscript

# Controlled comparison of the R interfaces to Varmapack and MTS.

MODELS <- c(
  "tinyAR", "tinyARMA", "smallAR1", "smallARMA2", "mediumAR",
  "mediumARMA2", "largeAR", "largeARMA"
)
THREAD_ENV <- c(
  "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1",
  "BLIS_NUM_THREADS=1", "VECLIB_MAXIMUM_THREADS=1"
)

print_help <- function() {
  cat("AnalyseSpeedup.R -- controlled Varmapack/MTS comparison\n\n")
  cat("Usage: Rscript AnalyseSpeedup.R [options]\n\n")
  cat("Options:\n")
  cat("  -h          show this help\n")
  cat("  -R count    independent timing runs (default 5)\n")
  cat("  -t seconds  target per individual timing (default 0.1)\n")
  cat("  -n length   retained length of each series (default 100)\n")
  cat("  -M count    number of series in a batch (default 1000)\n")
}

script_path <- function() {
  arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(arg) == 0L)
    stop("could not determine the script path", call. = FALSE)
  normalizePath(sub("^--file=", "", arg[[1L]]))
}

ensure_thread_environment <- function(args) {
  if (Sys.getenv("OMP_NUM_THREADS") == "1") return(invisible())
  message("Restarting with one BLAS/OpenMP thread")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c(shQuote(script_path()), vapply(args, shQuote, character(1L))),
    env = THREAD_ENV
  )
  if (!identical(status, 0L))
    stop("restarted benchmark failed with status ", status, call. = FALSE)
  quit(save = "no", status = 0L, runLast = FALSE)
}

parse_options <- function(args) {
  opts <- list(repeats = 5L, target = 0.1, n = 100L, M = 1000L)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "-h") {
      print_help()
      return(NULL)
    }
    if (!(arg %in% c("-R", "-t", "-n", "-M")) ||
        i == length(args))
      stop("invalid options; use -h for help", call. = FALSE)
    value <- args[[i + 1L]]
    if (arg == "-R") opts$repeats <- as.integer(value)
    if (arg == "-t") opts$target <- as.numeric(value)
    if (arg == "-n") opts$n <- as.integer(value)
    if (arg == "-M") opts$M <- as.integer(value)
    i <- i + 2L
  }
  if (is.na(opts$repeats) || opts$repeats < 1L ||
      !is.finite(opts$target) || opts$target <= 0 ||
      is.na(opts$n) || opts$n < 1L || is.na(opts$M) || opts$M < 2L)
    stop("-R, -t, and -n must be positive, and -M must be at least 2",
         call. = FALSE)
  opts
}

mts_terms <- function(terms, sign) {
  if (is.null(terms)) return(NULL)
  do.call(cbind, lapply(seq_len(dim(terms)[3L]), function(i) sign*terms[, , i]))
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

mts_draw <- function(specification, n, skip) {
  function() MTS::VARMAsim(
    nobs = n,
    arlags = specification$arlags,
    malags = specification$malags,
    phi = specification$phi,
    theta = specification$theta,
    sigma = specification$sigma,
    skip = skip
  )
}

vp_draw <- function(model, n, nrep) {
  rng <- randompack::randompack_rng()
  rng$seed(123L)
  function() model$sim(n, nrep = nrep, rng = rng)
}

generation_lengths <- function(n) {
  c(10L*n, 20L*n, 40L*n, 80L*n)
}

finite_result <- function(result) {
  values <- if (is.list(result)) unlist(result, use.names = FALSE) else result
  all(is.finite(values))
}

time_draw <- function(draw, target) {
  if (!finite_result(draw())) stop("simulation returned non-finite values")
  invisible(gc(FALSE))
  calls <- 0L
  start <- proc.time()[["elapsed"]]
  elapsed <- 0
  while (elapsed < target) {
    draw()
    calls <- calls + 1L
    elapsed <- proc.time()[["elapsed"]] - start
  }
  1e9*elapsed/calls
}

time_draw_pair <- function(first, second, target) {
  if (!finite_result(first()) || !finite_result(second()))
    stop("simulation returned non-finite values")
  invisible(gc(FALSE))
  calls <- 0L
  elapsed <- c(first = 0, second = 0)
  samples <- list(first = numeric(), second = numeric())
  while (any(elapsed < target)) {
    order <- if (calls%%2L == 0L) c("first", "second") else c("second", "first")
    for (name in order) {
      start <- proc.time()[["elapsed"]]
      if (name == "first") first() else second()
      duration <- proc.time()[["elapsed"]] - start
      elapsed[[name]] <- elapsed[[name]] + duration
      samples[[name]] <- c(samples[[name]], duration)
    }
    calls <- calls + 1L
  }
  1e9*vapply(samples, mean, numeric(1L), trim = 0.1)
}

time_vp_generation <- function(model, n, target) {
  lengths <- generation_lengths(n)
  times <- vapply(
    lengths,
    function(length) time_draw(vp_draw(model, length, 1L), target),
    numeric(1L)
  )
  generation <- unname(lm.fit(cbind(1, lengths), times)$coefficients[[2L]])
  if (!is.finite(generation) || generation <= 0)
    stop("invalid fitted Varmapack timing; increase -t", call. = FALSE)
  generation
}

time_r_case <- function(name, opts) {
  model <- varmapack::varmapack_testcase(name)
  specification <- mts_specification(model)
  h <- max(model$p, model$q)
  if (opts$n < h)
    stop("series length is shorter than startup for ", name, call. = FALSE)
  set.seed(123L)
  mts <- time_draw_pair(
    mts_draw(specification, opts$n, 0L),
    mts_draw(specification, opts$n, 200L),
    opts$target
  )
  values <- c(
    vp_generation1 = time_vp_generation(model, opts$n, opts$target),
    vp_total1 = time_draw(vp_draw(model, opts$n, 1L), opts$target),
    vp_setup_batch = time_draw(vp_draw(model, h, opts$M), opts$target),
    vp_total_batch = time_draw(
      vp_draw(model, opts$n, opts$M), opts$target
    ),
    mts_skip0 = mts[["first"]],
    mts_skip200 = mts[["second"]]
  )
  data.frame(
    model = name, p = model$p, q = model$q, r = model$r, h = h,
    as.list(values), check.names = FALSE
  )
}

time_r_interfaces <- function(opts) {
  rows <- vector("list", opts$repeats*length(MODELS))
  k <- 1L
  for (trial in seq_len(opts$repeats)) {
    message("R interface run ", trial, " of ", opts$repeats)
    for (model in MODELS) {
      row <- time_r_case(model, opts)
      row$trial <- trial
      rows[[k]] <- row
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

median_r_timings <- function(runs) {
  dimensions <- c("p", "q", "r", "h")
  timings <- setdiff(names(runs), c("trial", "model", dimensions))
  rows <- lapply(MODELS, function(model) {
    selected <- runs[runs$model == model, , drop = FALSE]
    values <- vapply(selected[, timings, drop = FALSE], median, numeric(1L))
    data.frame(
      model = model, selected[1L, dimensions, drop = FALSE], as.list(values),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_timings <- function(r_timings, opts) {
  table <- r_timings
  tail_length <- opts$n - table$h
  intermediate <- table$vp_setup_batch/opts$M +
    tail_length*table$vp_generation1
  table$vp_ns_value <- table$vp_total_batch/(opts$n*table$r*opts$M)
  table$mts_overall <- opts$M*table$mts_skip200/table$vp_total_batch
  table$burnin <- table$mts_skip200/table$mts_skip0
  table$setup_amortization <- table$vp_total1/intermediate
  table$bulk_generation <- intermediate/(table$vp_total_batch/opts$M)
  table$implementation <- table$mts_skip0/table$vp_total1
  table
}

align_markdown <- function(labels, rows) {
  cells <- rbind(labels, do.call(rbind, rows))
  widths <- apply(cells, 2L, function(column) max(nchar(column)))
  format_row <- function(row) {
    values <- vapply(seq_along(row), function(i) {
      flag <- if (i == 1L) "-" else ""
      sprintf(paste0("%", flag, widths[[i]], "s"), row[[i]])
    }, character(1L))
    paste0("| ", paste(values, collapse = " | "), " |")
  }
  separator <- vapply(seq_along(widths), function(i) {
    if (i == 1L) strrep("-", widths[[i]])
    else paste0(strrep("-", widths[[i]] - 1L), ":")
  }, character(1L))
  c(format_row(labels), format_row(separator),
    vapply(rows, format_row, character(1L)))
}

markdown_table <- function(table, opts) {
  labels <- c(
    "Model", "Overall", "Burn-in", "Setup amort.", "Bulk generation",
    "Implementation"
  )
  rows <- lapply(seq_len(nrow(table)), function(i) {
    row <- table[i, ]
    c(
      row$model, sprintf("%.0f", row$mts_overall), sprintf("%.1f", row$burnin),
      sprintf("%.1f", row$setup_amortization),
      sprintf("%.1f", row$bulk_generation), sprintf("%.1f", row$implementation)
    )
  })
  c(
    paste0(
      "Runs: ", opts$repeats, "; target: ", opts$target,
      " s; length: ", opts$n, "; batch: ", opts$M, "."
    ),
    "",
    "Speed-up factor due to:",
    "",
    align_markdown(labels, rows),
    ""
  )
}

main <- function(args) {
  ensure_thread_environment(args)
  opts <- parse_options(args)
  if (is.null(opts)) return(invisible())
  for (assignment in THREAD_ENV) {
    pieces <- strsplit(assignment, "=", fixed = TRUE)[[1L]]
    do.call(Sys.setenv, setNames(list(pieces[[2L]]), pieces[[1L]]))
  }
  for (package in c("varmapack", "randompack", "MTS")) {
    if (!requireNamespace(package, quietly = TRUE))
      stop("install the ", package, " package before running this benchmark",
           call. = FALSE)
  }
  r_runs <- time_r_interfaces(opts)
  table <- summarize_timings(median_r_timings(r_runs), opts)
  cat(markdown_table(table, opts), sep = "\n")
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
