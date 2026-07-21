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
   VarmapackError

Top-level functions
-------------------

.. autosummary::
   :toctree: reference/generated
   :nosignatures:

   autocov
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
