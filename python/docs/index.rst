Varmapack
=========

Varmapack is a VARMA time-series simulation library. It provides exact
stationary startup for Gaussian VARMA models, simulation with deterministic
exogenous inputs, theoretical autocovariances, impulse responses, and sample
autocovariances.

The Python interface is built on the same C core as the C interface. Random
number generation is supplied by Randompack.

Quick Example
-------------

.. code-block:: python

    import numpy as np
    import randompack
    import varmapack

    A = np.array([[[0.6, 0.1], [0.0, 0.4]]])
    B = np.array([[[0.2, 0.0], [0.1, 0.1]]])
    Sig = np.eye(2)

    rng = randompack.Rng()
    rng.seed(123)

    model = varmapack.Model(A=A, B=B, Sig=Sig)
    X = model.sim(200, nrep=100, rng=rng)

Model Objects
-------------

The central object is :class:`varmapack.Model`. A model stores the AR
coefficient matrices, MA coefficient matrices, innovation covariance, optional
time series means, and optional exogenous coefficients.

Python arrays use lag-first shapes:

- ``A`` has shape ``(p, r, r)``.
- ``B`` has shape ``(q, r, r)``.
- ``Sig`` has shape ``(r, r)``.
- ``C`` has shape ``(s, r)`` for VARMAX models.
- Simulated series have shape ``(nrep, length, r)``.

Simulation
----------

For stationary VARMA models without fixed startup values, Varmapack draws the
startup block from its exact stationary distribution. This avoids burn-in for
bulk simulation. Fixed startup values may be supplied with ``X0``.

.. code-block:: python

    X, E = model.sim(100, nrep=10, return_shocks=True)

VARMAX Simulation
-----------------

Supplying ``C`` creates a VARMAX model. VARMAX simulation also requires an
exogenous input sequence ``z`` and fixed startup values ``X0``.

.. code-block:: python

    C = np.array([[0.8, -0.2]])
    z = (-1.0)**np.arange(50)
    X0 = np.zeros((2, 2))

    model = varmapack.Model(A=A, B=B, C=C, Sig=Sig)
    X = model.sim(50, X0=X0, z=z)

Autocovariances and Impulse Responses
-------------------------------------

Theoretical autocovariances and impulse responses are model methods.

.. code-block:: python

    Gamma = model.acvf(10)
    Psi = model.psi(10)
    Theta = model.irf(10)

Sample autocovariances of observed data are computed by
:func:`varmapack.autocov`.

API Reference
-------------

The complete API reference is generated from the Python docstrings.

.. toctree::
   :maxdepth: 2

   reference/index
