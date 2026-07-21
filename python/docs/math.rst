Mathematical Description
========================

For more details of the mathematics, see the report [1]. (TODO: put in the
right reference).

VARMA and VARMAX models
-----------------------
The models considered are either VARMA\ :math:`(p,q)`:

.. math::

   \tag 1
   x_t = \sum_{i=1}^{p} A_i x_{t-i} + \eps_t + \sum_{j=1}^{q} B_j \eps_{t-j},
         \quad \eps_t \sim N(0,\Sigma).

or VARMAX\ :math:`(p,q,s)`:

.. math::

   \tag 2
   x_t = \eps_t + \sum_{i=1}^{p} A_i x_{t-i} + \sum_{j=1}^{q} B_j \eps_{t-j}
         + \sum_{k=1}^{s} C_k z_{t-k+1}, \quad \eps_t \sim N(0,\Sigma),

where in both cases :math:`x_t` is :math:`r`-dimensional, :math:`\eps` are
shocks or innovations, :math:`A_i`, :math:`B_j`, and :math:`C_k` are
autoregressive, moving-average, and exogenous coefficient matrices,
respectively, :math:`\Sigma` is the innovation covariance matrix, and
:math:`z_t` are exogenous terms. Varmapack supports a time-dependent mean path
:math:`\mu_t` for VARMA simulation (but not for VARMAX), in which case the
recursion is applied to the centered series :math:`x_t - \mu_t` (which then
replaces :math:`x_t` in (:math:`1`) or (:math:`2`)). The mean may also be fixed,
i.e. independent of the time step.

VAR, VMA, VARX, and VMAX models
-------------------------------
A model with :math:`q=s=0` is a pure vector-autoregressive (VAR) model, and when
:math:`p=s=0` the model is a pure moving-average (VMA) model. Similarly, when
:math:`q` or :math:`p` are zero and :math:`s > 0` the models are designated VARX
or VMAX.

Simulation start
----------------
There are two possibilities to start VARMA simulation with Varmapack: (a) by
drawing both shocks :math:`\eps_t` and series values :math:`x_t` from the exact
joint distribution of :math:`(x,\eps)` for the initial segment
:math:`t=0,\ldots,h-1` where :math:`h=\max(p,q)`, and (b) by specifying
:math:`h` initial values of the series and drawing the first :math:`h` shocks
from the conditional distribution of :math:`(\eps|x)`, where
:math:`h\geq\max(p,q)`. In both cases the simulated series has the correct
distribution from the first term, so there is no need for discarding a burn-in
start segment. In case (a) the model must be stationary, but for case (b) it may
be non-stationary.

For VARMAX simulation, drawing initial values of the series randomly is not
possible. It is necessary to supply initial :math:`x_t` for :math:`t=0,\ldots,h`
where :math:`h > max(p,s-1)`, and in addition the whole sequence of exogenous
values :math:`z_t, t=0,\ldots,n-1` must be specified, where :math:`n` is the
number of :math:`x_t` terms to be generated. The initial shocks are again drawn
from the exact conditional distribution of :math:`(\eps|x,z)`, alleviating the
need for burn-in discarding.

Autocovariances
---------------
Varmapack can compute the theoretical autocovariance function of models (1) and
(2),

.. math::
   \Gamma_k = \Cov(x_t, x_{t-k}), \qquad k=0,1,2,\ldots,

up to a user-specified maximum lag. It can also compute sample autocovariance
matrices from observed data. The default maximum-likelihood normalization is

.. math::         
   \widehat{\Gamma}_k =
   \frac{1}{n}\sum_{t=k}^{n-1}(x_t-\bar{x})(x_{t-k}-\bar{x})^T,

and the corrected normalization replaces the denominator :math:`n` by
:math:`n-k`.

Spectral radii
--------------
The spectral radius :math:`\rho` of the model :math:`(1)` is defined as the
spectral radius (maximum absolute eigenvalue) of the autoregressive companion
matrix, which is a square block matrix with :math:`A_1,\ldots,A_p` on the first
block row and identity matrices on the block-subdiagonal. Stationary models have
:math:`\rho < 1`. The spectral radius of the moving average companion matrix,
which has :math:`B_1,\ldots,B_q` on the top block row, plays a different role:
it provides an invertibility diagnostic. Models with MA companion spectral
radius :math:`<1` are invertible, in the sense that the moving-average
polynomial of the backshift operator :math:`L`, :math:`B(L)` has a convergent
inverse, :math:`B(L)^{-1} = I + D_1L + D_2 L^2 + \ldots`. It follows that the
model can be written as an equivalent infinite pure VAR model

.. math::
   x_t = E_1 x_{t-1} + E_2 x_{t-2} + \ldots + \eps_t,

where :math:`I-E_1L-E_2L^2-\ldots = (I + D_1L + D_2 L^2 + \ldots) A(L)`.

Impulse response functions
--------------------------
An important concept in time-series theory is that of impulse response
functions. The impulse response matrix :math:`\Psi_j` maps a change in the shock
at time :math:`t` to the resulting change in the process at time :math:`t + j`.
If the time series is stationary, it can be expressed as an infinite VMA series:

.. math::
   x_t = \sum_{j=0}^{\infty} \Psi_j \eps_{t-j},

which has the impulse response matrices as coefficients. With :math:`\Psi_0=I`,
the remaining :math:`\Psi` matrices can be computed with the recursion

.. math::
   \Psi_j = B_j + \sum_{i=1}^{\min(p,j)} A_i\Psi_{j-i},

for :math:`j\geq 1`, with :math:`B_j=0` for :math:`j>q`. Varmapack provides
functions to compute these impulse responses as well as *orthogonalized impulse
responses*, defined as

.. math::
   \Theta_j=\Psi_jL,

where :math:`L` is the lower Cholesky factor of the innovation covariance matrix
:math:`\Sigma`.
