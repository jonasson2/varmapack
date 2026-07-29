# Mathematical Description

For more details of the mathematics, see the report [1]. (TODO: put in the
right reference).

## VARMA and VARMAX models

The models considered are either VARMA $(p,q)$:

$$
\tag{1}
x_t = \sum_{i=1}^{p} A_i x_{t-i} + \varepsilon_t
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j},
      \quad \varepsilon_t \sim N(0,\Sigma).
$$

or VARMAX $(p,q,s)$:

$$
\tag{2}
x_t = \varepsilon_t + \sum_{i=1}^{p} A_i x_{t-i}
      + \sum_{j=1}^{q} B_j \varepsilon_{t-j}
      + \sum_{k=1}^{s} C_k z_{t-k+1},
      \quad \varepsilon_t \sim N(0,\Sigma),
$$

where in both cases $x_t$ is $r$-dimensional, $\varepsilon_t$ are shocks or
innovations,
$A_i$ and $B_j$ are $r$ by $r$ autoregressive and moving-average coefficient
matrices, $C_k$ is an $r$ by $d$ exogenous coefficient matrix, $z_t$ is a
$d$-dimensional exogenous vector, and $\Sigma$ is the innovation covariance
matrix. Varmapack supports a time-dependent mean path $\mu_t$ for VARMA
simulation. Explicit mean paths are unsupported for VARMAX. For VARMA, the
recursion is applied to the centered series
$x_t - \mu_t$, which replaces $x_t$ in equation (1). The mean may also be
fixed, i.e. independent of the time step.

## VAR, VMA, VARX, and VMAX models

A model with $q=s=0$ is a pure vector-autoregressive (VAR) model, and when
$p=s=0$ the model is a pure moving-average (VMA) model. Similarly, when $q$ or
$p$ are zero and $s > 0$ the models are designated VARX or VMAX.

## Simulation start

There are two possibilities to start VARMA simulation with Varmapack: (a) by
drawing both shocks $\varepsilon_t$ and series values $x_t$ from the exact joint
distribution of $(x,\varepsilon)$ for the initial segment $t=0,\ldots,h-1$ where
$h=\max(p,q)$, and (b) by specifying $h$ initial values of the series and
drawing the first $h$ shocks from the conditional distribution of
$(\varepsilon|x)$,
where $h\geq\max(p,q)$. In both cases, for a stationary model, the simulation
has the correct distribution from the first term, so there is no need for
discarding a burn-in start segment. In case (b), the simulated values have the
exact conditional distribution given the supplied starting values. A
nonstationary pure VAR model ($q=0$) can also be simulated from supplied
starting values, but a model with MA terms must be stationary.
For a stationary model with supplied starting values, the covariance of the
startup segment must be positive definite. When that covariance is singular
positive semidefinite, simulation is supported only without supplied starting
values.

For VARMAX simulation, drawing initial values of the series randomly is not
possible. It is necessary to supply the $h$ initial values $x_t$ for
$t=0,\ldots,h-1$, where $h > \max(p,s-1)$, and in addition the whole sequence
of exogenous values
$z_t, t=0,\ldots,n-1$ must be specified, where $n$ is the number of $x_t$ terms
to be generated. The initial shocks are again drawn from the exact conditional
distribution of $(\varepsilon|x,z)$, alleviating the need for burn-in discarding. The
shock covariance $\Sigma$ must be positive definite.

## Autocovariances

Varmapack can compute the theoretical autocovariance function of model (1),

$$
\Gamma_k = \Cov(x_t, x_{t-k}), \qquad k=0,1,2,\ldots,
$$

up to a user-specified maximum lag. It can also compute sample autocovariance
matrices from observed data. The default maximum-likelihood normalization is

$$
\widehat{\Gamma}_k =
\frac{1}{n}\sum_{t=k}^{n-1}(x_t-\bar{x})(x_{t-k}-\bar{x})^T,
$$

and the corrected normalization replaces the denominator $n$ by $n-k$.

## Spectral radii

The spectral radius $\rho$ of the model (1) is defined as the spectral radius
(maximum absolute eigenvalue) of the autoregressive companion matrix, which is a
square block matrix with $A_1,\ldots,A_p$ on the first block row and identity
matrices on the block-subdiagonal. Stationary models have $\rho < 1$. The
spectral radius of the moving average companion matrix, which has
$-B_1,\ldots,-B_q$ on the top block row, plays a different role: it provides
an invertibility diagnostic. Models with MA companion spectral radius $<1$ are
invertible, in the sense that the moving-average polynomial of the backshift
operator $L$, $B(L)$ has a convergent inverse,
$B(L)^{-1} = I + D_1L + D_2 L^2 + \ldots$. It follows that the model can be
written as an equivalent infinite pure VAR model

$$
x_t = E_1 x_{t-1} + E_2 x_{t-2} + \ldots + \varepsilon_t,
$$

where $I-E_1L-E_2L^2-\ldots = (I + D_1L + D_2 L^2 + \ldots) A(L)$.

## Impulse response functions

An important concept in time-series theory is that of impulse response
functions. The impulse response matrix $\Psi_j$ maps a change in the shock at
time $t$ to the resulting change in the process at time $t + j$. If the time
series is stationary, it can be expressed as an infinite VMA series:

$$
x_t = \sum_{j=0}^{\infty} \Psi_j \varepsilon_{t-j},
$$

which has the impulse response matrices as coefficients. With $\Psi_0=I$, the
remaining $\Psi$ matrices can be computed with the recursion

$$
\Psi_j = B_j + \sum_{i=1}^{\min(p,j)} A_i\Psi_{j-i},
$$

for $j\geq 1$, with $B_j=0$ for $j>q$. Varmapack provides functions to compute
these impulse responses as well as *orthogonalized impulse responses*, defined
as

$$
\Theta_j=\Psi_jL,
$$

where $L L^T=\Sigma$; for positive-definite $\Sigma$, $L$ is its lower
Cholesky factor. Only the lower triangle of $\Sigma$ is used.
