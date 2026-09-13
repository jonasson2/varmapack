#!/usr/bin/env Rscript

# Fit a VAR(1) to European stock-index log returns and simulate from the fit.
returns <- diff(log(EuStockMarkets))
fit <- ar(returns, aic = FALSE, order.max = 1L, method = "yule-walker",
          demean = TRUE)
A <- fit$ar[1L, , ]
Sig <- (fit$var.pred + t(fit$var.pred))/2
model <- varmapack::varmapack_model(A = A, Sig = Sig, mu = fit$x.mean)
rng <- randompack::randompack_rng()
rng$seed(2026L)
X <- model$sim(260L, nrep = 10L, rng = rng)
last <- as.numeric(tail(EuStockMarkets, 1L))
final <- last*exp(rowSums(X[, , 1L]))

cat(sprintf("AR spectral radius: %.6f\n", model$specrad()))
cat("First simulated path, final synthetic index levels:\n")
for (i in seq_along(final))
  cat(sprintf("  %-4s %9.2f\n", colnames(EuStockMarkets)[i], final[i]))
