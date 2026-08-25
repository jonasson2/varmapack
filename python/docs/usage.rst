Usage
=====

Quick start
-----------
A runnable version of this sequence is available as ``quickstart.py`` in the
``python/examples`` directory.

.. code-block:: python

    import varmapack
    Sig = [[2, 0], [0, 1]]
    VAR_model = varmapack.Model(A=[[0.6, 0.1], [0.0, 0.4]], Sig=Sig)
    X = VAR_model.sim(200)
    print(X[0, :5])

Returned arrays have shape ``(nrep, length, r)``: the first dimension selects
the replicate, the second is time, and the third is the series dimension.
With ``return_shocks=True``, ``sim`` returns a pair ``(X, E)`` containing the
series and shocks.

Repeated simulations can reuse a preallocated array with ``out``. It must be a
writable, C-contiguous ``float64`` array with shape ``(nrep, length, r)``.

.. code-block:: python

    X = np.empty((100, 200, VAR_model.r))
    VAR_model.sim(200, nrep=100, out=X)

.. code-block:: python

    print(varmapack.testcases())
    model = varmapack.testcase("smallARMA1")
    X = model.sim(200, nrep=10)

.. code-block:: python

    import randompack, numpy as np
    A1 = np.array([[0.6, 0.1], [0.0, 0.4]])
    B1 = np.array([[0.4, 0.3], [0.1, 0.2]])
    A2 = 0.2*np.ones((2, 2))
    B2 = 0.5*np.eye(2)
    z = np.column_stack((np.cos(np.arange(200)/10), np.sin(np.arange(200)/10)))
    VARMA_model = varmapack.Model(A=[A1, A2], B=B1, Sig=Sig)
    VMA_model = varmapack.Model(B=[B1, B2], Sig=Sig)
    VARMAX_model = varmapack.Model(A=A1, B=B1,
                                   C=[[[0.8, 0.2], [-0.3, 0.4]]], Sig=Sig)
    rng = randompack.Rng()
    rng.seed(123)
    X1 = VAR_model.sim(200, nrep=100)
    X2 = VMA_model.sim(200, nrep=100, rng=rng)
    X3, E = VARMA_model.sim(200, nrep=100, rng=rng, return_shocks=True)
    X4 = VARMAX_model.sim(200, nrep=100, X0=np.zeros((2, 2)), z=z, rng=rng)

.. code-block:: python

    maxlag = 3
    rho = VARMA_model.specrad()
    rhoMA = VARMA_model.ma_specrad()
    Psi = VARMA_model.psi(maxlag)
    Theta = VARMA_model.irf(maxlag)
    Gamma = VARMA_model.acvf(maxlag)
    Corr = varmapack.cov2corr(Gamma)
    C = varmapack.autocov(X3[0], maxlag)

Model objects
-------------
The central object is :class:`varmapack.Model`. A model stores the AR
coefficient matrices A, MA coefficient matrices B, innovation covariance matrix
Sig, optional time series means mu, and optional exogenous coefficients C. The
model dimensions are inferred from these arrays.

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
``s``, and ``d`` as read-only properties.

Top-level functions
-------------------
The Python Varmapack has four top-level functions: ``testcase``,
``testcases``, ``autocov``, and ``cov2corr``. The *testcase* function provides
15 named testcases, describing models of various complexity ranging from
scalar AR(1) and MA(1) models to VARMA(3,3) with :math:`r=3` and VARMA(4,0)
with :math:`r=5`. The function can also give unnamed testcases with arbitrary
user-supplied :math:`p`, :math:`q`, and :math:`r`. Use ``testcases()`` to
return a printable overview of the named testcases. [check testcases with
specified rho]

The autocov function computes sample autocovariance matrices, up to a
specified maximum lag, from an observed time series ``X``. The matrices can be
normalized by :math:`1/n` (maximum likelihood) or :math:`1/(n-k)`
(lag-corrected).

The ``cov2corr`` function converts theoretical or sample autocovariances to
correlations by dividing each entry by the product of the corresponding
lag-zero marginal standard deviations. The returned array has the same shape
as the input. With ``out`` omitted, ``cov2corr`` allocates and returns a new
array, leaving ``cov`` unchanged. Use ``out=`` to specify a preallocated output
array and ``out=cov`` for in-place conversion. Lag-zero diagonal entries are
exactly one. Other entries are not clipped to ``[-1,1]``, so correlations
obtained from lag-corrected sample autocovariances may lie outside that
interval.

Array shapes
------------
- ``A`` has shape ``(p, r, r)`` or ``(r, r)`` for a single AR matrix.
- ``B`` has shape ``(q, r, r)`` or ``(r, r)`` for a single MA matrix.
- ``C`` has shape ``(s, d, r)`` for VARMAX models: ``C[k, j, i]`` is element
  ``(i, j)`` of the coefficient matrix at lag ``k+1``. Scalar exogenous input
  also accepts ``(s, r)`` or ``(r,)`` as ``d=1`` conveniences.
- ``Sig`` has shape ``(r, r)``.
- The mean ``mu`` has shape ``(r,)`` or ``(nmu, r)`` where ``nmu`` ≤
  ``length``; the last supplied row repeats to the end.
- Startup values ``X0`` have shape ``(nX0, r)`` or ``(nrep, nX0, r)``.
- Nonstationary VARMA simulation requires ``X0``. With MA terms, ``Sig`` must
  be positive definite and startup shocks are conditioned on the implied
  residual equations.
- For scalar exogenous input, ``z`` has shape ``(length,)`` or
  ``(nrep, length)``. For d-dimensional input, it has shape ``(length, d)``
  or ``(nrep, length, d)``.
- Simulated series, and returned shocks, have shape ``(nrep, length, r)``.
- The ``autocov`` data argument ``X`` has shape ``(n, r)``.
