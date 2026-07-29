Varmapack
=========

Varmapack is a Python package for simulation and analysis of Gaussian VAR, VMA,
VARMA, and VARMAX time-series models. Its primary purpose is to generate
simulated series from supplied model parameters, but it also provides model
testcases, theoretical and sample autocovariances, spectral radii, and impulse
response functions. In contrast to some other packages, the simulated series
have the correct distribution from the start; they are burn-in (or spin-up)
free. The numerical work is performed by the Varmapack C library, with random
number generation provided by the companion library Randompack.

.. toctree::
   :maxdepth: 2

   install
   math
   usage
   reference
