"""
Python interface to Varmapack.

Varmapack simulates Gaussian VARMA and VARMAX time-series models and
computes related quantities such as theoretical autocovariances, impulse
responses, and sample autocovariances.

The main entry point is :class:`Model`.

Examples
--------
Create a bivariate VARMA model and simulate 100 independent paths of
length 200.

>>> import numpy as np
>>> import varmapack
>>> A = np.array([[[0.6, 0.1], [0.0, 0.4]]])
>>> B = np.array([[[0.2, 0.0], [0.1, 0.1]]])
>>> Sig = np.eye(2)
>>> model = varmapack.Model(A=A, B=B, Sig=Sig)
>>> X = model.sim(200, nrep=100)
>>> X.shape
(100, 200, 2)

Random number generation is delegated to Randompack. Pass a
``randompack.Rng`` object to control the stream.

>>> import randompack
>>> rng = randompack.Rng()
>>> rng.seed(123)
>>> X = model.sim(20, rng=rng)

Shape conventions
-----------------
Python arrays use lag-first and path-first shapes:

- ``A`` has shape ``(p, r, r)``.
- ``B`` has shape ``(q, r, r)``.
- ``Sig`` has shape ``(r, r)``.
- Simulated series have shape ``(nrep, length, r)``.
"""

try:
    from importlib.metadata import version as _version
except ImportError:  # pragma: no cover
    from importlib_metadata import version as _version

from ._core import Model, VarmapackError, autocov, testcase

try:
    __version__ = _version("varmapack")
except Exception:  # pragma: no cover
    __version__ = "unknown"

__all__ = ["Model", "VarmapackError", "__version__", "autocov", "testcase"]
