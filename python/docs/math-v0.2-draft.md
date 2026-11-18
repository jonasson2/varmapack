# Varmapack 0.2: Mathematics of New Features (Draft)

**Draft: not yet proof-read.** This note contains LLM-generated material
describing the mathematics of new features for the Varmapack 0.2 series
(0.2.xxx). It supplements the existing mathematical summary rather than
replacing it.

## Continuous-time OU models

With persistence times, `OUModel` specifies a stationary Gaussian
Ornstein--Uhlenbeck process with diagonal drift. Write
$\Gamma=\operatorname{Cov}(X(t),X(t))$ and let
$\tau_i>0$ be component persistence times. Over a gap $d>0$ the exact transition is

$$
X(t+d)-\mu=A_d(X(t)-\mu)+\epsilon_d,\qquad
A_d=\operatorname{diag}(e^{-d/\tau_i}),
$$

where $A_d$ is the autoregressive coefficient matrix for a gap of length $d$.
The independent transition innovations have covariance

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

## OU models with a general drift matrix

The diagonal model can be extended to allow cross-variable dynamics. Replace
$\operatorname{diag}(1/\tau_i)$ by a real $r\times r$ matrix $K$ whose
eigenvalues have strictly positive real parts. Our sign convention puts
$-K$ in the drift of the stochastic differential equation

$$
dX(t)=-K(X(t)-\mu)\,dt+C\,dW(t),\qquad D=CC^T,
$$

where $W$ is standard multivariate Brownian motion. Unlike diagonal drift,
off-diagonal entries of $K$ allow one component to influence the subsequent
evolution of another. $K$ need not be symmetric or diagonalizable; complex
eigenvalues can produce damped oscillations. There is generally no single
persistence time for each observed component.

### Covariance and transitions

With the stationary covariance $\Gamma$ supplied, the driving covariance is

$$
D=K\Gamma+\Gamma K^T.
$$

Both $\Gamma$ and $D$ must be symmetric positive semidefinite. Stability of
$-K$ and positive semidefiniteness of $\Gamma$ alone do not ensure a valid
model: the implied $D$ must also pass this check. Conversely, if $D$ is
supplied instead of $\Gamma$, the stationary covariance solves this
continuous-time Lyapunov equation. See Särkkä and Solin (2019), section 6.5
([author-hosted text](https://users.aalto.fi/~asolin/sde-book/sde-book.pdf)).

For any gap $d>0$, the exact transition retains the same form as above:

$$
\begin{aligned}
A_d&=\exp(-Kd),\\
X(t+d)-\mu&=A_d(X(t)-\mu)+\epsilon_d,\\
\epsilon_d&\sim N(0,Q_d),\qquad Q_d=\Gamma-A_d\Gamma A_d^T.
\end{aligned}
$$

Here $A_d$ is the autoregressive coefficient matrix over a gap of length $d$,
not the continuous-time drift parameter $K$. The exponential is a matrix
exponential, not entrywise exponentiation.
The innovations for successive transitions are independent of one another
and of the past. Equivalently,

$$
Q_d=\int_0^d \exp(-Ku)D\exp(-K^Tu)\,du,
$$

which also shows why $Q_d$ is positive semidefinite. The lagged covariance is
$\operatorname{Cov}(X(t+d),X(t))=A_d\Gamma$. Initialization remains
$X(t_0)\sim N(\mu,\Gamma)$, so simulation needs no burn-in.

### Relation to discrete-time autoregression

Taking $K=\operatorname{diag}(1/\tau_i)$ recovers the diagonal model.
With a constant observation gap $d$, the sampled process is a VAR(1)
with autoregressive coefficient matrix $A_d$ and innovation covariance
$Q_d$. With unequal gaps, both matrices vary from step to step.

We take $K$, not the unit-step coefficient $A_1$, as input because $A_1$
does not always determine what happens between observations. Recovering
$K$ would require a real matrix logarithm of $A_1$, with a minus sign.
Such a logarithm need not exist: a scalar AR(1) coefficient of $-0.5$
is valid, but no real $K$ gives $e^{-K}=-0.5$. Even when a real logarithm
exists, it need not be unique. For example, two damped rotations whose
angular speeds differ by one full turn per time unit give the same
$A_1$ but different half-step transitions. This is
[system aliasing](https://arxiv.org/abs/1605.06973). Choosing a logarithm
is therefore a modeling assumption, not just a change of notation.

The inputs are $K$, the stationary covariance $\Gamma$, and the mean
$\mu$, with the same observation times and Gamma time generator as for
the diagonal model. The simulator uses a Van Loan block exponential to
evaluate the integral for $Q_d$ without subtracting nearly equal
covariances at short gaps. It computes a short transition using SLICOT
MB05OD and, if necessary, doubles the gap repeatedly using
$A_{2d}=A_d^2$ and $Q_{2d}=Q_d+A_d Q_d A_d^T$. No Lyapunov solve is needed
at each time step. A diagonal $K$ uses the scalar transition formulas.
For general $K$, the simulator checks $A_d\Gamma A_d^T+Q_d=\Gamma$ in
marginal standard-deviation units and reports an error if the computed
transition fails that check.

## Continuous-time Slepian-based MA models

`SlepianModel` is a stationary Gaussian model that can be sampled
at arbitrary observation times and reduces to a vector MA(1) model at unit
time steps. The construction uses Brownian increments over windows of fixed
length one. Thus one unit of time has the same meaning throughout the model;
observation times must not be rescaled without changing that interpretation.

### Windowed noise

Let $W(t)$ be an $r$-dimensional, two-sided Brownian motion whose increments
over disjoint intervals are independent and satisfy

$$
W(t+d)-W(t)\sim N(0,d\Sigma),\qquad d>0,
$$

where $\Sigma$ is positive semidefinite. Define the windowed noise

$$
E(t)=W(t)-W(t-1).
$$

Two such increments share only the noise in the overlap of their windows.
Consequently, writing $k(h)=\max(0,1-|h|)$,

$$
\operatorname{Cov}(E(t+h),E(t))=k(h)\Sigma.
$$

In particular, $E(t+h)$ and $E(t)$ are independent when $|h|\geq1$.
For shorter gaps they are generally correlated. In the scalar case with
$\Sigma=1$, this windowed-noise process is the standard Slepian process;
see Bischoff and Gegg, *Boundary crossing probabilities for
$(q,d)$-Slepian-processes*
([preprint](https://arxiv.org/abs/1607.07260)).

### Moving-average combination

For an arbitrary real $r\times r$ coefficient matrix $\Theta$, set

$$
X(t)=\mu+E(t)+\Theta E(t-1).
$$

This combines two adjacent noise windows, $(t-1,t]$ and $(t-2,t-1]$.
At integer times the $\varepsilon_j=E(j)$ are independent $N(0,\Sigma)$
vectors, so

$$
X(j)=\mu+\varepsilon_j+\Theta\varepsilon_{j-1},
$$

which is an ordinary vector MA(1) model. At noninteger observation times,
the window overlaps define the dependence without rounding the times.
Strictly speaking, $E$ is the Slepian process and $X$ is a filtered,
Slepian-based process; `SlepianModel` names this model family.

With $\Gamma(h)=\operatorname{Cov}(X(t+h),X(t))$, the covariance is

$$
\Gamma(h)=k(h)(\Sigma+\Theta\Sigma\Theta^T)
          +k(h+1)\Sigma\Theta^T+k(h-1)\Theta\Sigma.
$$

Thus $\Gamma(-h)=\Gamma(h)^T$, $\Gamma(1)=\Theta\Sigma$, and observations
are independent at separations of at least two time units. This finite
dependence range contrasts with the exponentially decaying dependence of OU.

### Parameters and exact simulation

Unlike `OUModel`, the Slepian-based model can take the unit-step lag-1
moving-average coefficient matrix $\Theta$ directly: window overlap
determines the dependence at fractional lags, with no matrix logarithm.
The other inputs are the contemporaneous covariance $V=\Gamma(0)$
and the mean $\mu$. The innovation covariance must satisfy

$$
V=\Sigma+\Theta\Sigma\Theta^T.
$$

This is exactly the VMA(1) inverse-covariance problem described below, so
the existing `cov2sig` computation recovers $\Sigma$. As there, an
incompatible target covariance or a singular or poorly conditioned inverse
problem must be reported rather than silently accepted.

For observation times $t_i$, the simulator merges the three sorted streams of
endpoints $t_i$, $t_i-1$, and $t_i-2$. For each consecutive pair of distinct
endpoints $u<v$, it draws an independent increment $N(0,(v-u)\Sigma)$ if
the interval is used by at least one window. Summing the shared
increments over the required windows gives $E(t_i)$ and $E(t_i-1)$, hence
$X(t_i)$. The arbitrary Brownian level at the earliest endpoint cancels in
every difference. Including the windows before the first observation gives
the stationary distribution immediately, with no burn-in. The same Gamma
time generator used for OU can supply the observation times.

## Recovering innovation covariance

Consider a stationary, causal VARMA model

$$
x_t=\sum_{i=1}^{p}A_i x_{t-i}
    +\sum_{j=0}^{q}B_j\varepsilon_{t-j},
\qquad B_0=I,
$$

where $x_t$ is $r$-dimensional and the innovations are temporally uncorrelated
with covariance $\Sigma$.

Suppose the coefficient matrices and the desired contemporaneous covariance

$$
V=\operatorname{Cov}(x_t,x_t)
$$

are given. We wish to determine $\Sigma$.

### The existing VYW system

For $p\geq1$, we use the reduced vector Yule--Walker (VYW) system derived in
Appendix B of Jónasson and Ferrando (2006), cited below. Their lagged covariance
matrices $S_k$ are denoted here by

$$
\Gamma_k=\operatorname{Cov}(x_t,x_{t-k}),
\qquad \Gamma_0=V,\qquad \Gamma_{-k}=\Gamma_k^T.
$$

The appendix eliminates $\Gamma_p$ and exploits the symmetry of $\Gamma_0$.
Let $\operatorname{vec}$ stack matrix columns and $\operatorname{vech}$ stack the
lower triangle column by column. Using the appendix's packing convention, define

$$
\begin{aligned}
v&=\operatorname{vech}(\widehat{\Gamma}_0),&
\sigma&=\operatorname{vech}(\Sigma),\\
h&=\operatorname{col}\bigl(
  \operatorname{vec}\Gamma_1,\ldots,\operatorname{vec}\Gamma_{p-1}\bigr).
\end{aligned}
$$

Here $\widehat{\Gamma}_0$ is lower triangular, with the same strictly lower
triangle as $\Gamma_0$ and half its diagonal, so that
$\Gamma_0=\widehat{\Gamma}_0+\widehat{\Gamma}_0^T$. The diagonal of $\Sigma$ is
*not* halved in $\sigma$. The operator $\operatorname{col}$ stacks vectors;
$h$ is empty when $p=1$.

The reduced system can be written as

$$
F_0v+F_1h=g(\Sigma).
$$

The matrix $F=[\,F_0\;\;F_1\,]$ is the existing VYW coefficient matrix, constructed
as in Appendix B, equations (B.3)--(B.5). It depends only on the AR coefficients.
Writing

$$
m=\frac{r(r+1)}{2},\qquad K=m+(p-1)r^2,
$$

$F_0$ has size $K\times m$ and $F_1$ has size $K\times(K-m)$.

### Linearity of the right-hand side

The existing covariance recurrences are

$$
C_0=\Sigma,\qquad
C_j=B_j\Sigma+\sum_{i=1}^{\min(p,j)}A_iC_{j-i},
\quad j=1,\ldots,q,
$$

and

$$
G_k=\sum_{j=k}^{q}B_jC_{j-k}^T,\quad k=0,\ldots,q,
\qquad G_k=0\quad(k>q).
$$

Thus every $C_j$ and $G_k$ is linear in $\Sigma$. The reduced right-hand side
from Appendix B is

$$
g(\Sigma)=\operatorname{col}\bigl(
  \operatorname{vech}(G_0^T+G_pA_p^T),
  \operatorname{vec}G_1,\ldots,\operatorname{vec}G_{p-1}\bigr)
=R\sigma,
$$

where $R$ is a $K\times m$ matrix determined by the AR and MA coefficients.
The $\operatorname{vech}$ operation in $g$ simply selects lower-triangle entries;
it does not halve diagonal entries.

### Forward and inverse computations

In the usual forward computation, $\Sigma$ is known, and the VYW solve is

$$
\boxed{
\begin{bmatrix}F_0&F_1\end{bmatrix}
\begin{bmatrix}v\\h\end{bmatrix}=R\sigma.
}
$$

For the inverse computation, $\Gamma_0=V$ is known instead. Move its contribution
to the right-hand side and the contribution of $\Sigma$ to the left:

$$
\boxed{
\begin{bmatrix}R&-F_1\end{bmatrix}
\begin{bmatrix}\sigma\\h\end{bmatrix}=F_0v.
}
$$

The two systems have the same dimension $K$. Both $\Gamma_0$ and $\Sigma$ have
$m$ independent entries, so the inverse formulation exchanges which covariance
is prescribed while retaining the auxiliary lagged covariances as unknowns.
Only $\Sigma$ needs to be returned.

No new derivation of the Kronecker-product matrices is required: $F_0$ and $F_1$
come from the existing VYW construction. To assemble $R$, use the symmetric
basis matrices $E_a$ satisfying $\operatorname{vech}(E_a)=e_a$, where $e_a$ is
the $a$th coordinate vector in $\mathbb{R}^m$. Then

$$
R_{\cdot a}=g(E_a),\qquad a=1,\ldots,m.
$$

A diagonal basis matrix has one unit diagonal entry; an off-diagonal basis
matrix has a pair of symmetric unit entries. Each column requires only the
existing right-hand-side construction, not a full covariance solve.

This is an *inverse VYW formulation*, not an adjoint or transpose problem.
The coefficient matrix changes, so a new LU factorization is needed.

### Special cases

For a pure VAR(1) model, $\Gamma_1=A_1V$, giving the explicit formula

$$
\boxed{\Sigma=V-A_1VA_1^T.}
$$

For a pure VMA($q$) model, the equation is

$$
\boxed{V=\Sigma+\sum_{j=1}^{q}B_j\Sigma B_j^T.}
$$

### Existence and validity

A positive-semidefinite target $V$ does not necessarily yield a positive-semidefinite
$\Sigma$. Thus not every target covariance is compatible with fixed VARMA
coefficients.

The inverse calculation must check for singularity or poor conditioning and verify
that the recovered $\Sigma$ is positive semidefinite and reproduces $V$ through
the forward covariance calculation. Both the solve and the checks use rescaled
series: each nonzero marginal standard deviation is brought near one by a
power-of-two change of units. This prevents a large variance in one series
from masking an incompatible covariance in another.

### Reference

Jónasson, K. and Ferrando, S. E. (2006).
*Efficient likelihood evaluation for VARMA processes with missing values*.
Report VHI-01-2006, Faculty of Engineering, University of Iceland.
Appendix B: *Solution of the vector Yule--Walker equations*.
