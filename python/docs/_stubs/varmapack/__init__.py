"""
Python interface to Varmapack.

Varmapack simulates Gaussian VARMA and VARMAX time-series models and
computes related quantities such as theoretical autocovariances, impulse
responses, and sample autocovariances.
"""

__version__ = "0.1.0"


class VarmapackError(RuntimeError):
  """Exception raised when the Varmapack C library reports an error."""


class Model:
  """
  VARMA or VARMAX model.

  Parameters
  ----------
  A, B, C : array_like, optional
      Autoregressive, moving-average, and exogenous coefficient matrices.
  Sig : array_like
      Shock covariance matrix.
  mu : array_like, optional
      Time-series means.
  """
  def __init__(self, *, A=None, B=None, C=None, Sig=None, mu=None):
    pass

  @property
  def p(self):
    """Number of autoregressive terms."""

  @property
  def q(self):
    """Number of moving-average terms."""

  @property
  def r(self):
    """Dimension of the time series."""

  @property
  def s(self):
    """Number of exogenous terms."""

  @property
  def A(self):
    """Autoregressive coefficient matrices."""

  @property
  def B(self):
    """Moving-average coefficient matrices."""

  @property
  def C(self):
    """Exogenous coefficient matrices."""

  @property
  def Sig(self):
    """Shock covariance matrix."""

  @property
  def mu(self):
    """Time-series means."""

  def sim(self, length, *, nrep=1, X0=None, z=None, return_shocks=False,
          rng=None):
    """
    Simulate independent paths from the model.

    Parameters
    ----------
    length : int
        Length of each generated series.
    nrep : int, optional
        Number of independent simulated paths.
    X0 : array_like, optional
        Starting values.
    z : array_like, optional
        Exogenous input path.
    return_shocks : bool, optional
        Whether to return shock series together with simulated series.
    rng : randompack.Rng, optional
        Random number generator.
    """

  def acvf(self, maxlag):
    """Compute theoretical autocovariances up to ``maxlag``."""

  def psi(self, maxlag):
    """Compute impulse-response coefficient matrices up to ``maxlag``."""

  def irf(self, maxlag):
    """Compute orthogonalized impulse responses up to ``maxlag``."""

  def specrad(self):
    """Compute the autoregressive spectral radius."""

  def ma_specrad(self):
    """Compute the moving-average spectral radius."""


def autocov(X, maxlag, norm="ML"):
  """
  Compute sample autocovariances of one observed series.

  Parameters
  ----------
  X : array_like, shape (n, r)
      Observed time series.
  maxlag : int
      Largest lag to compute.
  norm : {"ML", "C"}, optional
      Normalization.
  """


def testcase(which="random", *, p=None, q=None, r=None, rho=0, rng=None):
  """
  Create a testcase model.

  Parameters
  ----------
  which : str or int, optional
      Testcase selector.
  p, q, r : int, optional
      AR order, MA order, and series dimension.
  rho : float, optional
      Target autoregressive spectral radius when ``which="rho"``.
  rng : randompack.Rng, optional
      Random number generator.
  """


def testcases():
  """Return a printable overview of named VARMA testcases."""


__all__ = ["Model", "VarmapackError", "__version__", "autocov", "testcase",
           "testcases"]
