Bivariate VARMA example
=======================

This example considers a zero-mean bivariate VARMA(2,1) model whose components
are denoted by :math:`x_t` and :math:`y_t`, with coefficient matrices

.. math::

   A_1 = \begin{pmatrix}0.75&0.05\\0&0.50\end{pmatrix} \!,\; A_2 = \begin{pmatrix}0.13&0\\0&0.05\end{pmatrix} \!,\; B_1 = \begin{pmatrix}0.40&0.15\\0.05&0.20\end{pmatrix} \!,\; \Sigma = \begin{pmatrix}1&0.99\\0.99&1\end{pmatrix}.

The program simulates five paths with a supplied initial segment and plots the
theoretical cross-correlation.

Program
-------

.. literalinclude:: ../examples/bivariate_varma.py
   :language: python
   :linenos:

Cross-correlation
-----------------

.. figure:: figures/bivariate_varma_cross_correlation.*
   :alt: Theoretical cross-correlation of the two components
   :align: center
   :class: short-caption
   :width: 67%

   Figure 1. The theoretical cross-correlation
   :math:`\rho_{xy}(k)=\operatorname{Corr}(x_t,y_{t-k})`. At positive lags,
   :math:`y_t` leads :math:`x_t`.

Simulated paths
---------------

.. figure:: figures/bivariate_varma_paths.*
   :alt: Five simulated paths from a bivariate VARMA model
   :class: long-caption

   Figure 2. Five simulated replicates of the bivariate VARMA(2,1) model
   specified above, with the fixed initialization path :math:`x_t=y_t=t/2` for
   :math:`t=0,\ldots,10`. Note that :math:`x_t` is much more persistent than
   :math:`y_t`. Looking closely at the paths also reveals that :math:`x_t` and
   :math:`y_t` are noticeably correlated, with :math:`y_t` leading
   :math:`x_t`, consistent with Figure 1.
