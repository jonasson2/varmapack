varmapack
================

`varmapack` is an R package for simulation and analysis of Gaussian VAR,
VMA, VARMA, and VARMAX time-series models. Simulated stationary series
have the correct distribution from the first returned term, without
discarding a burn-in segment. The package also provides named model
testcases, theoretical and sample autocovariances, spectral radii, and
impulse response functions.

## Installation

The package will be installable from CRAN when released:

``` r
install.packages("varmapack")
```

For development versions, install from the repository root with standard
R package tooling. Varmapack requires the companion `randompack` R
package.

``` r
remotes::install_github("jonasson2/randompack", subdir = "r-package")
remotes::install_github("jonasson2/varmapack", subdir = "r-package")
```

## Quick start

``` r
library(varmapack)
library(randompack)
```

Create a model from coefficient matrices and its innovation covariance.
A single lag may be supplied as a matrix; multiple lags use an `r` by
`r` by lag array.

``` r
A <- matrix(c(0.5, 0.1,
              0.0, 0.3), 2, 2)
B <- matrix(c(0.2, 0.0,
              0.1, 0.1), 2, 2)
model <- varmapack_model(A = A, B = B, Sig = diag(2))

rng <- randompack_rng()
rng$seed(123)
X <- model$sim(100, nrep = 3, rng = rng)
dim(X)
```

The first dimension is the series dimension, the second is time, and the
third selects the replicate. Requesting shocks returns a named list.

``` r
out <- model$sim(20, nrep = 2, rng = rng, return_shocks = TRUE)
names(out)
dim(out$E)
```

Built-in testcases return model objects.

``` r
varmapack_testcases()
test_model <- varmapack_testcase("smallARMA1")
test_model$specrad()
```

Model methods provide theoretical autocovariances, impulse responses,
and the AR and MA spectral radii.

``` r
Gamma <- model$acvf(10)
Psi <- model$psi(10)
Theta <- model$irf(10)
model$specrad()
model$ma_specrad()
```

The `varmapack_autocov()` function computes sample autocovariances for
an observed time-series matrix with observations in rows.

``` r
varmapack_autocov(t(X[, , 1]), maxlag = 5)
```

VARMAX models use exogenous coefficient matrices `C` and input values
`z`. They require fixed starting values `X0`.

``` r
C <- array(c(0.3, -0.2), c(2, 1, 1))
varmax <- varmapack_model(A = A, B = B, C = C, Sig = diag(2))
X0 <- matrix(0, 2, 2)
z <- matrix(sin(seq_len(100)/10), 1, 100)
Xmax <- varmax$sim(100, X0 = X0, z = z, rng = rng)
```

## Model objects

Create a model with:

``` r
varmapack_model(A = NULL, B = NULL, C = NULL, Sig, mu = NULL)
```

`Sig` is required; the other arguments may be omitted. Supplying `C`
creates a VARMAX model, for which `mu` is not supported. A model stores
its coefficient matrices, innovation covariance, optional means, and
dimensions `p`, `q`, `r`, `s`, and `d`. Its methods are:

| Method          | Purpose                                           |
|:----------------|:--------------------------------------------------|
| `$sim()`        | Simulate one or several series.                   |
| `$specrad()`    | Compute the AR spectral radius.                   |
| `$ma_specrad()` | Compute the MA spectral radius.                   |
| `$acvf()`       | Compute theoretical VARMA autocovariances.        |
| `$psi()`        | Compute impulse-response matrices.                |
| `$irf()`        | Compute orthogonalized impulse-response matrices. |

The top-level functions `varmapack_testcase()`, `varmapack_testcases()`,
and `varmapack_autocov()` construct testcases, list named testcases, and
compute sample autocovariances, respectively.

## Array shapes

R arrays use column-major order, which matches the Varmapack C
interface.

- `A` has shape `r` by `r` by `p`, or `r` by `r` for one AR term.
- `B` has shape `r` by `r` by `q`, or `r` by `r` for one MA term.
- `C` has shape `r` by `d` by `s`, or `r` by `d` for one exogenous term.
- `Sig` has shape `r` by `r`.
- `mu` has length `r` or shape `r` by `nmu`, where `nmu` is at most the
  simulated series length; its last value repeats to the series end.
- `X0` has shape `r` by `nX0`, or `r` by `nX0` by `nrep` for different
  starts.
- `z` has shape `d` by length, or `d` by length by `nrep` for different
  paths.
- Simulated series, and returned shocks, have shape `r` by length by
  `nrep`.
- `varmapack_autocov()` accepts observed data with shape length by `r`.

## Mathematical description

### VARMA and VARMAX models

The models considered are either VARMA $(p,q)$:

$$\tag{1}
x_t = \sum_{i=1}^{p} A_i x_{t-i} + \varepsilon_t
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j},
      \quad \varepsilon_t \sim N(0,\Sigma).$$

or VARMAX $(p,q,s)$:

$$\tag{2}
x_t = \varepsilon_t + \sum_{i=1}^{p} A_i x_{t-i}
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j}
      + \sum_{k=1}^{s} C_k z_{t-k+1},
      \quad \varepsilon_t \sim N(0,\Sigma),$$

where in both cases $x_t$ is $r$-dimensional, $\varepsilon_t$ are shocks
or innovations, $A_i$ and $B_j$ are $r$ by $r$ autoregressive and
moving-average coefficient matrices, $C_k$ is an $r$ by $d$ exogenous
coefficient matrix, $z_t$ is a $d$-dimensional exogenous vector, and
$\Sigma$ is the innovation covariance matrix. Varmapack supports a
time-dependent mean path $\mu_t$ for VARMA simulation. Explicit mean
paths are unsupported for VARMAX. For VARMA, the recursion is applied to
the centered series $x_t - \mu_t$, which replaces $x_t$ in equation (1).
The mean may also be fixed, i.e. independent of the time step.

### VAR, VMA, VARX, and VMAX models

A model with $q=s=0$ is a pure vector-autoregressive (VAR) model, and
when $p=s=0$ the model is a pure moving-average (VMA) model. Similarly,
when $q$ or $p$ are zero and $s > 0$ the models are designated VARX or
VMAX.

### Simulation start

There are two possibilities to start VARMA simulation with Varmapack:
(a) by drawing both shocks $\varepsilon_t$ and series values $x_t$ from
the exact joint distribution of $(x,\varepsilon)$ for the initial
segment $t=0,\ldots,h-1$ where $h=\max(p,q)$, and (b) by specifying $h$
initial values of the series and drawing the first $h$ shocks from the
conditional distribution of $(\varepsilon|x)$, where $h\geq\max(p,q)$.
In both cases, for a stationary model, the simulation has the correct
distribution from the first term, so there is no need for discarding a
burn-in start segment. In case (b), the simulated values have the exact
conditional distribution given the supplied starting values. A
nonstationary pure VAR model ($q=0$) can also be simulated from supplied
starting values, but a model with MA terms must be stationary. For a
stationary model with supplied starting values, the covariance of the
startup segment must be positive definite. When that covariance is
singular positive semidefinite, simulation is supported only without
supplied starting values.

For VARMAX simulation, drawing initial values of the series randomly is
not possible. It is necessary to supply the $h$ initial values $x_t$ for
$t=0,\ldots,h-1$, where $h > \max(p,s-1)$, and in addition the whole
sequence of exogenous values $z_t, t=0,\ldots,n-1$ must be specified,
where $n$ is the number of $x_t$ terms to be generated. The initial
shocks are again drawn from the exact conditional distribution of
$(\varepsilon|x,z)$, alleviating the need for burn-in discarding. The
shock covariance $\Sigma$ must be positive definite.

### Autocovariances

Varmapack can compute the theoretical autocovariance function of model
(1),

$$\Gamma_k = \operatorname{Cov}(x_t, x_{t-k}), \qquad k=0,1,2,\ldots,$$

up to a user-specified maximum lag. It can also compute sample
autocovariance matrices from observed data. The default
maximum-likelihood normalization is

$$\widehat{\Gamma}_k =
\frac{1}{n}\sum_{t=k}^{n-1}(x_t-\bar{x})(x_{t-k}-\bar{x})^T,$$

and the corrected normalization replaces the denominator $n$ by $n-k$.

### Spectral radii

The spectral radius $\rho$ of the model (1) is defined as the spectral
radius (maximum absolute eigenvalue) of the autoregressive companion
matrix, which is a square block matrix with $A_1,\ldots,A_p$ on the
first block row and identity matrices on the block-subdiagonal.
Stationary models have $\rho < 1$. The spectral radius of the moving
average companion matrix, which has $-B_1,\ldots,-B_q$ on the top block
row, plays a different role: it provides an invertibility diagnostic.
Models with MA companion spectral radius $<1$ are invertible, in the
sense that the moving-average polynomial of the backshift operator $L$,
$B(L)$ has a convergent inverse,
$B(L)^{-1} = I + D_1L + D_2 L^2 + \ldots$. It follows that the model can
be written as an equivalent infinite pure VAR model

$$x_t = E_1 x_{t-1} + E_2 x_{t-2} + \ldots + \varepsilon_t,$$

where $I-E_1L-E_2L^2-\ldots = (I + D_1L + D_2 L^2 + \ldots) A(L)$.

### Impulse response functions

An important concept in time-series theory is that of impulse response
functions. The impulse response matrix $\Psi_j$ maps a change in the
shock at time $t$ to the resulting change in the process at time
$t + j$. If the time series is stationary, it can be expressed as an
infinite VMA series:

$$x_t = \sum_{j=0}^{\infty} \Psi_j \varepsilon_{t-j},$$

which has the impulse response matrices as coefficients. With
$\Psi_0=I$, the remaining $\Psi$ matrices can be computed with the
recursion

$$\Psi_j = B_j + \sum_{i=1}^{\min(p,j)} A_i\Psi_{j-i},$$

for $j\geq 1$, with $B_j=0$ for $j>q$. Varmapack provides functions to
compute these impulse responses as well as *orthogonalized impulse
responses*, defined as

$$\Theta_j=\Psi_jL,$$

where $L L^T=\Sigma$; for positive-definite $\Sigma$, $L$ is its lower
Cholesky factor. Only the lower triangle of $\Sigma$ is used.
