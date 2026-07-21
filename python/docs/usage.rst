Usage
=====

Quick start
-----------
.. code-block:: python

    import varmapack
    Sig = [[2, 0], [0, 1]]
    VAR_model = varmapack.Model(A=[[0.6, 0.1], [0.0, 0.4]], Sig=Sig)
    X = VAR_model.sim(200)
    print(X[0, :5])

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
    z = np.cos(np.arange(200)/10)
    VARMA_model = varmapack.Model(A=[A1, A2], B=B1, Sig=Sig)
    VMA_model = varmapack.Model(B=[B1, B2], Sig=Sig)
    VARMAX_model = varmapack.Model(A=A1, B=B1, C=[0.8, 0.2], Sig=Sig)
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
and ``s`` as read-only properties.

Top-level functions
-------------------
The Python Varmapack has three top-level functions: ``testcase``,
``testcases``, and ``autocov``. The *testcase* function provides 15 named
testcases, describing models of various complexity ranging from scalar AR(1)
and MA(1) models to VARMA(3,3) with :math:`r=3` and VARMA(4,0) with
:math:`r=5`. The function can also give unnamed testcases with arbitrary
user-supplied :math:`p`, :math:`q`, and :math:`r`. Use ``testcases()`` to
return a printable overview of the named testcases. [check testcases with
specified rho]

The autocov function computes sample autocovariance matrices, up to a
specified maximum lag, from an observed time series ``X``. The matrices can be
normalized by :math:`1/n` (maximum likelihood) or :math:`1/(n-k)`
(lag-corrected).

Array shapes
------------
- ``A`` has shape ``(p, r, r)`` or ``(r, r)`` for a single AR matrix.
- ``B`` has shape ``(q, r, r)`` or ``(r, r)`` for a single MA matrix.
- ``C`` has shape ``(s, r)`` or ``(r,)`` for VARMAX models.
- ``Sig`` has shape ``(r, r)``.
- The mean ``mu`` has shape ``(r,)`` or ``(nmu, r)`` where ``nmu`` ≤
  ``length``; the last supplied row repeats to the end.
- Startup values ``X0`` have shape ``(nX0, r)`` or ``(nrep, nX0, r)``.
- Exogenous inputs ``z`` have shape ``(length,)`` or ``(nrep, length)``.
- Simulated series, and returned shocks, have shape ``(nrep, length, r)``.
- The ``autocov`` data argument ``X`` has shape ``(n, r)``.
