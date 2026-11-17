# Mathematical Description

For more details of the mathematics, see Kristján Jónasson (2026),
["Burn-in-Free Simulation of VARMA Time Series"]
(https://arxiv.org/abs/2608.29199), arXiv:2608.29199.

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
drawing both shocks $\varepsilon_t$ and states $x_t$ from the exact joint
distribution of $(x,\varepsilon)$ for the initial segment $t=0,\ldots,h-1$,
where $h=\max(p,q)$, and (b) by specifying $h \geq\max(p,q)$ initial values of
the series and drawing the first $h$ shocks from the conditional distribution of
$(\varepsilon|x)$. In case (b), the covariance matrix of the startup states must
be positive definite.

For a stationary model, both (a) and (b) are possible, and the simulated series
will have the correct distribution from the first generated term, so that no
burn-in segment needs to be discarded.

For a nonstationary model, starting values must be supplied. For a pure VAR
model ($q=0$), the recurrence simply runs forward from those values. If MA terms
are present, startup innovations are drawn from their theoretical distribution,
conditional on constraints imposed by the model and supplied states. $\Sigma$
must be positive semidefinite and is allowed to be singular except for
nonstationary models with MA terms.

With VARMAX simulation, startup states must be provided: $x_t$ for
$t=0,\ldots,h-1$, where $h\geq\max(p,q,s-1)$, as well as the whole sequence of
exogenous values $z_t$, $t=0,\ldots,n-1$. As for nonstationary VARMA, $\Sigma$
must be positive definite, and the startup innovations are then drawn from their
theoretical distribution, conditional on constraints imposed by the model and
the supplied states and exogenous sequence. When more than $\max(p,q,s-1)$
values of $x_t$ and corresponding $z_t$ are available, they should be supplied,
to improve the information used to draw the innovations.

## Continuous-time OU models

`OUModel` specifies a stationary Gaussian Ornstein--Uhlenbeck process with
diagonal drift. Write $\Gamma=\operatorname{Cov}(X(t),X(t))$ and let
$\tau_i>0$ be component persistence times. Over a gap $d>0$ the exact transition is

$$
X(t+d)-\mu=A_d(X(t)-\mu)+\epsilon_d,\qquad
A_d=\operatorname{diag}(e^{-d/\tau_i}),
$$

where independent transition innovations have covariance

$$
Q_d=\Gamma-A_d\Gamma A_d^T,\qquad
(Q_d)_{ij}=\Gamma_{ij}(1-e^{-d(1/\tau_i+1/\tau_j)}).
$$

Initialization uses $X(t_0)\sim N(\mu,\Gamma)$, with no burn-in. Both
$\Gamma$ and the driving covariance
$D_{ij}=\Gamma_{ij}(1/\tau_i+1/\tau_j)$ must be positive semidefinite.
This permits correlated noise but not cross-variable drift. The lagged
covariance is $\operatorname{Cov}(X(t+d),X(t))=A_d\Gamma$.
See Särkkä and Solin (2019), *Applied Stochastic Differential Equations*,
sections 6.2 and 6.5 ([author-hosted text](https://users.aalto.fi/~asolin/sde-book/sde-book.pdf)).

The separate Gamma time generator draws independent gaps with shape $s$ and
scale $\bar d/s$, hence expected gap $\bar d$ and CV $1/\sqrt{s}$. It returns
$t_0=0$ and cumulative gap sums without normalizing the sample mean. Infinite
shape denotes constant gaps. The simulator also accepts actual observation times.

## Autocovariances

Varmapack can compute the theoretical autocovariance function of model (1),

$$
\Gamma_k = \Cov(x_t, x_{t-k}), \qquad k=0,1,2,\ldots,
$$

up to a user-specified maximum lag.

The inverse problem, recovering $\Sigma$ from a specified $\Gamma_0$, is
described in the separate note *Recovering Innovation Covariance*
(`output/pdf/cov2sig.pdf` in the repository; source: `output/pdf/cov2sig.tex`).

Varmapack can also compute sample autocovariance matrices from observed data.
The default maximum-likelihood normalization is

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
