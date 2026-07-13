Varmapack
=========
Varmapack is a vector-autoregressive-moving-average (VARMA) time-series
simulation library. Apart from simulation, it provides a variety of model
testcases and supports the computation of theoretical and sample
autocovariances, model spectral radii, and impulse responses. The simulation is
burn-in free: startup shocks are drawn from the correct conditional distribution
prescribed by the model, eliminating the need to discard an initial part of the
simulated series. Exogenous model terms (VARMAX) are also supported, with exact
startup.

The Python interface is built on the C library Varmapack, which is also the
basis for Matlab, R, and Fortran interfaces. Random number generation is
provided by Randompack, a companion package for reproducible random number
streams. For ordinary use, importing ``varmapack`` is sufficient. Users who want
explicit control over seeding or parallel streams can additionally create and
pass a ``randompack.Rng`` object, as shown below.

Short mathematical description
------------------------------
For more details of the mathematics, see the report [1]. (TODO: put in the right
reference).

VARMA and VARMAX models
~~~~~~~~~~~~~~~~~~~~~~~
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

where in both cases :math:`x_t` is :math:`r`-dimensional, :math:`A_i`,
:math:`B_j`, and :math:`C_k` are autoregressive, moving-average, and exogenous
coefficient matrices, respectively, :math:`\Sigma` is the innovation covariance
matrix, and :math:`z_t` are exogenous terms. Varmapack supports a time-dependent
mean path :math:`\mu_t` for VARMA simulation (but not for VARMAX), in which case
the recursion is applied to the centered series :math:`x_t - \mu_t`. The mean
may also be fixed and independent of the time step.

VAR, VMA, VARX, and VMAX models
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A model with :math:`q=s=0` is a pure vector-autoregressive (VAR) model, and when
:math:`p=s=0` the model is a pure moving-average (VMA) model. Similarly, when
:math:`q` or :math:`p` are zero and :math:`s > 0` the models are designated VARX
or VMAX.

Simulation start
~~~~~~~~~~~~~~~~
There are two possibilities to start VARMA simulation with Varmapack (a) by
drawing both shocks :math:`\eps_t` and series values :math:`x_t` from the exact
joint distribution of :math:`(x,\eps)` for the initial segment
:math:`t=0,\ldots,h-1` where :math:`h=\max(p,q)`, and (b) by specifying the
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
~~~~~~~~~~~~~~~
Varmapack can compute both the theoretical autocovariance structure of the
models (1) and (2), as well as sample autocovariance matrices of observed data,
in both cases up to a user-specified maximum lag: :math:`\Cov(x_t, x_{t-k})`
for :math:`k=0,1,2,\ldots`.

Spectral radii
~~~~~~~~~~~~~~
The spectral radius :math:`\rho` of the model :math:`(1)` is defined as the
spectral radius (maximum absolute eigenvalue) of the autoregressive companion
matrix, which is a block matrix with :math:`A_1,\ldots,A_p` on the first block
row and identity matrices on the block-subdiagonal. Stationary models have
:math:`\rho < 1`. The spectral radius of the moving average companion matrix,
which has :math:`B_1,\ldots,B_q` on the top block row plays a different role: it
is an invertibility diagnostic. Models with MA-spectral radius :math:`<1` are
invertible, in the sense that...[TODO]

Impulse response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~
An important concept in normal time-series theory are impulse response
functions. The impulse response :math:`\Psi_j` measures the change in the process at time
:math:`t + j` caused by a shock change at time :math:`t`

.. math::

   x_t = \sum_{j=0}^{\infty} \Psi_j \eps_{t-j}.

Then :math:`\Psi_0=I` and

.. math::

   \Psi_j = B_j + \sum_{i=1}^{\min(p,j)} A_i\Psi_{j-i},

for :math:`j\geq 1`, with :math:`B_j=0` for :math:`j>q`. The function
``Model.psi(maxlag)`` returns these matrices for lags
:math:`0,\ldots,\texttt{maxlag}`. The orthogonalized impulse responses returned
by ``Model.irf(maxlag)`` are :math:`\Theta_j=\Psi_jL`, where :math:`LL'=\Sigma`.

Usage
-----

Quick start
~~~~~~~~~~~
.. code-block:: python

    import varmapack
    Sig = [[2, 0], [0, 1]]
    VAR_model = varmapack.Model(A = [[0.6, 0.1], [0.0, 0.4]], Sig = Sig)
    X = VAR_model.sim(200)
    print("X=", X)

.. code-block:: python

    ... = testcase...
    ...sim...

.. code-block:: python

    import randompack, numpy as np
    A1 = np.array([[0.6, 0.1], [0.0, 0.4]])
    B1 = np.array([[0.4, 0.3], [0.1, 0.2]])
    A2 = 0.3*np.ones((2,2))
    B2 = 0.5*np.eye(2)
    VARMA_model = varmapack.Model(A=[A1,A2], B=B1, Sig=Sig)
    VMA_model = varmapack.Model(B=[B1,B2], Sig=Sig)
    rng = randompack.Rng()
    rng.seed(123)
    X1 = VAR_model.sim(200, nrep=100)
    X2 = VMA_model.sim(200, nrep=100, rng=rng)
    X3,E = VARMA_model.sim(200, nrep=100, rng=rng, return_shocks=True)

.. code-block:: python

   maxlag = 3
   rho = VARMA_model.specrad()
   rhoMA = VARMA_model.ma_specrad()
   Psi = VARMA_model.psi(maxlag)
   Theta = VARMA_model.irf(maxlag)
   Gamma = VARMA_model.acvf(maxlag)

Model objects
~~~~~~~~~~~~~
The central object is :class:`varmapack.Model`. A model stores the AR
coefficient matrices A, MA coefficient matrices B, innovation covariance matrix
Sig, optional time series means mu, and optional exogenous coefficients C. The
model dimensions,

The signature of the constructor is:

.. code-block:: python

    varmapack.Model(*, A=None, B=None, C=None, Sig=None, mu=None)

All arguments are keyword-only. ``Sig`` is required; the other arguments may be
omitted. Supplying ``C`` creates a VARMAX model, in which case ``mu`` is not
accepted. A model has the following methods:

   ``sim``        Create one or several simulated series
   ``specrad``    Compute the model's spectral radius (<1 for stationary models)
   ``ma_specrad`` Compute inverse model's spec.rad. (<1 for invertible models)
   ``acvf``       Return theoretical autocovariance function
   ``irf``        Impulse response function
   ``psi``        Non-orthogonalized impulse response function

Models also expose ``A``, ``B``, ``C``, ``Sig``, ``mu``, ``p``, ``q``, ``r``,
and ``s`` as read-only properties.

Functions
~~~~~~~~~
The Python Varmapack has two functions, ``testcase`` and ``autocov`` The
*testcase* function provides 15 named testcases, describing models of various
complexity ranging from scalar AR(1) and MA(1) models to VARMA(3,3) with
:math:`r=3` and VARMA(4,0) with :math:`r=5`. The function can also give unnamed
testcases with arbitrary user-supplied $p,q,s,r$. Use testcase(xxxx) to display
an overview of the named testcases. [check testcases with specified rho]

The autocov function



Array shapes
~~~~~~~~~~~~
- ``A`` has shape ``(p, r, r)``.
- ``B`` has shape ``(q, r, r)``.
- ``C`` has shape ``(s, r)`` for VARMAX models.
- Simulated series have shape ``(nrep, length, r)``.
- The mean ``mu`` has shape ``(r)`` or ``(nmu,r)`` where ``nmu`` ≤ ``length``;
  the last supplied row repeats to the end.

API Reference
-------------

The complete API reference is generated from the Python docstrings.

.. toctree::
   :maxdepth: 2

   reference/index
