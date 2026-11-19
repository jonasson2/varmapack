API Reference
=============

.. currentmodule:: varmapack

The public Python objects exposed by ``varmapack`` are listed below.

Classes
-------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   Model
   OUModel
   SlepianModel
   VarmapackError

Top-level functions
-------------------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   autocov
   cov2corr
   cov2sig
   gamma_times
   testcase
   testcases

Model methods
-------------

.. autosummary::
   :toctree: reference/generated
   :template: method.rst
   :nosignatures:

   Model.sim
   Model.acvf
   Model.psi
   Model.irf
   Model.specrad
   Model.ma_specrad

Model attributes
----------------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   Model.p
   Model.q
   Model.r
   Model.s
   Model.A
   Model.B
   Model.C
   Model.Sig
   Model.mu

OU model methods and attributes
-------------------------------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   OUModel.sim
   OUModel.acvf
   OUModel.specrad
   OUModel.tau
   OUModel.drift
   OUModel.cov
   OUModel.mu
   OUModel.r

Slepian model methods and attributes
------------------------------------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   SlepianModel.sim
   SlepianModel.acvf
   SlepianModel.theta
   SlepianModel.cov
   SlepianModel.mu
   SlepianModel.r
