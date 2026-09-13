# Empirical example: European stock-index returns

This example uses the `EuStockMarkets` dataset distributed with R. It contains
1,860 daily closing prices for the DAX, SMI, CAC, and FTSE stock indices from
1991 to 1998. The example transforms the prices to log returns and fits a
four-dimensional VAR(1) model with the base-R `ar` function using the
multivariate Yule--Walker method.

The R program performs the fit directly. The fitted autoregressive matrix,
innovation covariance, and mean are recorded in the C, Python, and Matlab
programs so that those examples require no separate estimation package or
toolbox. Each program then uses Varmapack to generate ten independent exact
stationary simulations of 260 trading days and converts the first simulated
return path to synthetic index levels starting from the final observed prices.

The four implementations are:

- [C](EuStockMarkets.c)
- [R](../r-package/inst/examples/EuStockMarkets.R)
- [Python](../python/examples/eu_stock_markets.py)
- [Matlab](../matlab/examples/EuStockMarkets.m)

With a configured Meson build, build and run the C example using:

```sh
ninja -C build examples/EuStockMarkets
build/examples/EuStockMarkets
```

The other versions require their corresponding Varmapack and Randompack
interfaces to be installed or available on the language search path.
