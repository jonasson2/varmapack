# Varmapack

Varmapack is a Python interface to the Varmapack C library for simulation and
analysis of Gaussian VAR, VMA, VARMA, and VARMAX time-series models.

It supports burn-in-free simulation, model testcases, theoretical and sample
autocovariances, covariance conversions, innovation-covariance recovery,
spectral radii, and impulse response functions.
VARMAX models support scalar or vector-valued exogenous inputs.
`OUModel` adds exact stationary simulation at uneven observation times;
`gamma_times` generates independent Gamma spacings when needed.

## Quick Example

```python
import varmapack
model = varmapack.Model(
    A=[[0.6, 0.1], [0.0, 0.4]],
    Sig=[[2, 0], [0, 1]],
)
X = model.sim(200)
print(X.shape)
```

For repeated simulations, `Model.sim` accepts a writable, C-contiguous
`float64` NumPy array through `out=`. Reusing an output array avoids an
allocation on each call; its required shape is `(nrep, length, r)`.

The repository also provides a Python version of the
[European stock-index example](https://github.com/jonasson2/varmapack/blob/main/examples/EuStockMarkets.md),
which simulates from a VAR(1) fitted to empirical multivariate returns.

## Uneven Observation Times

```python
import randompack
rng = randompack.Rng()
rng.seed(123)
t = varmapack.gamma_times(100, shape=16, mean_gap=1, rng=rng)
model = varmapack.OUModel(tau=[1, 2], cov=[[1, .5], [.5, 1]])
X = model.sim(t, nrep=10, rng=rng)  # (10, 100, 2)
```

`cov` specifies stationary covariance, not innovation covariance. With
`tau`, the drift is diagonal. For a full matrix use
`OUModel(drift=K, cov=cov)`, with drift `-K @ (X - mu)`. The autoregressive
coefficient matrix over a gap `d` is the matrix exponential `exp(-K*d)`.
Gamma gaps are not normalized to their realized mean;
use `shape=float("inf")` for equal spacing. See the usage guide for covariance
compatibility requirements.

`SlepianModel(theta=[[.5]], cov=[[1.25]])` uses the same `sim` interface
for a Slepian-based moving-average model, which reduces to MA(1) at unit
spacing.

## Timing

`examples/TimeSimulate.py` compares Varmapack with Statsmodels `VARMAX` on all
named testcases and with Statsmodels `VAR` on pure-VAR testcases. The latter
uses a 200-observation burn-in by default. The program uses length 100, 1000
replicates, and a 0.1-second timing target per method and testcase. Its output
can be compared with the C and MATLAB simulation benchmarks. It is a
source-tree benchmark and requires the optional `statsmodels` package:

```sh
    python -m pip install statsmodels
    python examples/TimeSimulate.py
```

Use `-h` to list its timing and workload options.

## Documentation

The full documentation is at:

> https://varmapack.readthedocs.io/

For information about the underlying C library, see the
[C README](https://github.com/jonasson2/varmapack/blob/main/README.md).

For development builds from source, see:

> https://github.com/jonasson2/varmapack/blob/main/DEVELOPMENT.md
