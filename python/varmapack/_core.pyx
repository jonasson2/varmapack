import operator
import numpy as np
cimport numpy as np

from libc.stddef cimport size_t
from cpython.pycapsule cimport PyCapsule_GetPointer

DTYPE_F64 = np.dtype(np.float64)
RNG_CAPSULE_NAME = b"randompack.randompack_rng"
cdef const char *RNG_CAPSULE_CNAME = "randompack.randompack_rng"

cdef extern from "RandompackPythonGateway.h":
    ctypedef struct randompack_rng
    int randompack_PY_initialize() except -1
    randompack_rng *randompack_PY_create(const char *engine)
    void randompack_PY_free(randompack_rng *rng)

cdef extern from "varmapack.h":
    ctypedef struct varmapack_ou_model
    varmapack_ou_model *varmapack_ou_create(double *tau, double *cov,
                                           double *mu, int r)
    varmapack_ou_model *varmapack_ou_create_drift(double *K, double *cov,
                                                double *mu, int r)
    void varmapack_ou_free(varmapack_ou_model *model)
    bint varmapack_ou_sim(varmapack_ou_model *model, double *t, int n, int M,
                          double *X, randompack_rng *rng)
    ctypedef struct varmapack_slepian_model
    varmapack_slepian_model *varmapack_slepian_create(double *theta, double *cov,
                                                     double *mu, int r)
    void varmapack_slepian_free(varmapack_slepian_model *model)
    bint varmapack_slepian_sim(varmapack_slepian_model *model, double *t, int n,
                               int M, double *X, randompack_rng *rng)
    bint varmapack_gamma_times(int n, double shape, double mean_gap, double *t,
                               randompack_rng *rng)
    enum:
        VARMAPACK_TESTCASE_NAME_LEN
    char *varmapack_last_error()
    bint varmapack_sim(double *A, double *B, double *Sig,
                       double *mu, int nmu, int p, int q, int r, int n,
                       int M, double *X0, int nX0, int MX0, double *X,
                       double *E, randompack_rng *rng)
    bint varmapack_simx(double *A, double *B, double *C, double *Sig,
                        double *z, int Mz, int p, int q, int s, int d, int r, int n,
                        int M, double *X0, int h, int MX0, double *X,
                        double *E, randompack_rng *rng)
    bint varmapack_testcase(char *name, int *index, double rho, int *p, int *q,
                            int *r, double *A, double *B, double *Sig,
                            randompack_rng *rng)
    bint varmapack_acvf(double *A, double *B, double *Sig, int p, int q, int r,
                        double *Gamma, int maxlag)
    bint varmapack_psi(double *A, double *B, int p, int q, int r, int maxlag,
                       double *Psi)
    bint varmapack_irf(double *A, double *B, double *Sig, int p, int q, int r,
                       int maxlag, double *Theta)
    bint varmapack_autocov(const char *transp, const char *norm, int r, int n,
                           double *X, int maxlag, double *C)
    bint varmapack_cov2corr(double *cov, int r, int maxlag, double *corr)
    bint varmapack_cov2sig(double *A, double *B, double *cov,
                          int p, int q, int r, double *Sig)
    double varmapack_specrad(double *A, int r, int p)
    double varmapack_ma_specrad(double *B, int r, int q)

randompack_PY_initialize()
np.import_array()


class VarmapackError(RuntimeError):
    """Exception raised when the Varmapack C library reports an error."""
    pass


cdef int _positive_count(object value, str name) except -1:
    value = operator.index(value)
    if value < 1 or value > np.iinfo(np.int32).max:
        raise ValueError(f"{name} must be a positive 32-bit integer")
    return value


def gamma_times(object n, double shape=16, double mean_gap=1, object rng=None):
    """Generate observation times starting at zero with independent Gamma gaps.

    Parameters
    ----------
    n : int
        Number of times, at least one.
    shape : float, optional
        Positive Gamma shape. ``np.inf`` gives constant spacing.
    mean_gap : float, optional
        Positive expected gap, in the same units as OU persistence times.
    rng : randompack.Rng, optional
        Generator controlling the stream; omitted for a temporary default RNG.

    Returns
    -------
    t : ndarray, shape (n,)
        Strictly increasing times with ``t[0] == 0``. Gaps have scale
        ``mean_gap/shape`` and CV ``1/sqrt(shape)``. Their realized mean is
        not normalized. Unrepresentable gaps or times raise VarmapackError.
    """
    cdef np.ndarray t
    cdef bint owned = False
    cdef randompack_rng *ptr
    cdef int count = _positive_count(n, "n")
    t = np.empty(count, dtype=DTYPE_F64)
    ptr = _rng_from_object(rng, &owned)
    try:
        if not varmapack_gamma_times(count, shape, mean_gap, <double *>t.data, ptr):
            raise VarmapackError(_error_message())
    finally:
        if owned:
            randompack_PY_free(ptr)
    return t


cdef class OUModel:
    """Stationary Gaussian Ornstein--Uhlenbeck model.

    Parameters
    ----------
    tau : array_like, shape (r,), optional
        Finite positive persistence times, in the units used for observation times.
        Specify exactly one of ``tau`` and ``drift``.
    cov : array_like, shape (r, r)
        Stationary process covariance, not innovation covariance. Only its lower
        triangle is used. Both this matrix and the driving covariance
        ``K @ cov + cov @ K.T`` must be positive semidefinite.
        With ``tau``, K is ``diag(1/tau)``. Singular covariances are allowed.
    mu : array_like, shape (r,), optional
        Constant mean vector; defaults to zero.
    drift : array_like, shape (r, r), optional, keyword-only
        Real matrix K in the drift term ``-K @ (X - mu)``. All eigenvalues of K
        must have strictly positive real parts. Its units are inverse time.

    Notes
    -----
    Simulates exact transitions on arbitrary strictly increasing times without
    burn-in. Over a gap of length d, the autoregressive coefficient matrix is
    the matrix exponential ``exp(-K*d)``. At a constant gap, the sampled
    process is a VAR(1). With diagonal drift, component i has autocorrelation
    ``exp(-d/tau[i])``.
    Constructor inputs are copied. Parameter properties return copies.
    """
    cdef varmapack_ou_model *handle
    cdef np.ndarray tauarr, driftarr, covarr, meanarr
    cdef int rval

    def __cinit__(self):
        self.handle = NULL

    def __init__(self, object tau=None, object cov=None, object mu=None, *,
                 object drift=None):
        cdef np.ndarray ta = None, K = None, co, me
        cdef varmapack_ou_model *new_handle
        cdef int r
        if (tau is None) == (drift is None):
            raise ValueError("specify exactly one of tau and drift")
        co = _matrix_to_c(cov, "cov")
        if co is None or co.shape[0] < 1 or co.shape[0] != co.shape[1]:
            raise ValueError("cov must be a non-empty square matrix")
        r = _positive_count(co.shape[0], "model dimension")
        me = np.zeros(r) if mu is None else np.array(mu, dtype=DTYPE_F64)
        if me.ndim != 1 or me.shape[0] != r:
            raise ValueError("mu must be a vector with length r")
        if drift is None:
            ta = np.array(tau, dtype=DTYPE_F64, copy=True)
            if ta.ndim != 1 or ta.shape[0] != r:
                raise ValueError("tau must be a vector with length r")
            new_handle = varmapack_ou_create(<double *>ta.data, <double *>co.data,
                                            <double *>me.data, r)
        else:
            if np.iscomplexobj(drift):
                raise ValueError("drift must be a real matrix")
            K = _matrix_to_c(drift, "drift")
            if K.shape[0] != r or K.shape[1] != r:
                raise ValueError("drift must have shape (r, r)")
            new_handle = varmapack_ou_create_drift(<double *>K.data,
                                                  <double *>co.data,
                                                  <double *>me.data, r)
        if new_handle == NULL:
            raise VarmapackError(_error_message())
        varmapack_ou_free(self.handle)
        self.handle = new_handle
        self.tauarr, self.driftarr, self.covarr, self.meanarr = ta, K, co, me
        self.rval = r

    def __dealloc__(self):
        varmapack_ou_free(self.handle)

    @property
    def r(self):
        """Number of series."""
        return self.rval

    @property
    def tau(self):
        """Persistence times, or None for matrix initialization."""
        return None if self.tauarr is None else self.tauarr.copy()

    @property
    def drift(self):
        """Matrix K in the drift -K @ (X - mu), or None if initialized with tau.

        The autoregressive coefficient matrix over a gap d is exp(-K*d),
        the matrix exponential.
        """
        return None if self.driftarr is None else self.driftarr.T.copy()

    @property
    def cov(self):
        """Stationary process covariance."""
        c = self.covarr.T.copy()
        i, j = np.triu_indices(self.rval, 1)
        c[i, j] = c[j, i]
        return c

    @property
    def mu(self):
        """Constant mean vector."""
        return self.meanarr.copy()

    def sim(self, object t, *, object nrep=1, object rng=None):
        """Simulate independent stationary paths at the supplied times.

        Parameters
        ----------
        t : array_like, shape (n,)
            Finite, strictly increasing observation times.
        nrep : int, optional
            Positive number of independent paths.
        rng : randompack.Rng, optional
            Stream to use instead of a temporary default Randompack generator.

        Returns
        -------
        X : ndarray, shape (nrep, n, r)
            Gaussian paths; the first observation is drawn from N(mu, cov).
        """
        cdef np.ndarray times = np.ascontiguousarray(t, dtype=DTYPE_F64)
        cdef np.ndarray X
        cdef bint owned = False
        cdef randompack_rng *ptr
        cdef int count = _positive_count(nrep, "nrep")
        cdef int n
        if self.handle == NULL:
            raise ValueError("OUModel has not been initialized")
        if times.ndim != 1 or times.shape[0] < 1:
            raise ValueError("t must be a nonempty vector")
        n = _positive_count(times.shape[0], "length of t")
        X = np.empty((count, n, self.rval), dtype=DTYPE_F64)
        ptr = _rng_from_object(rng, &owned)
        try:
            if not varmapack_ou_sim(self.handle, <double *>times.data,
                                    n, count, <double *>X.data, ptr):
                raise VarmapackError(_error_message())
        finally:
            if owned:
                randompack_PY_free(ptr)
        return X


cdef class SlepianModel:
    """Stationary Gaussian model that reduces to a vector MA(1) at unit steps.

    Parameters
    ----------
    theta : array_like, shape (r, r)
        Real unit-step MA(1) coefficient matrix.
    cov : array_like, shape (r, r)
        Stationary process covariance; only the lower triangle is used.
        The driving covariance Sig must be positive semidefinite and
        determined by ``cov = Sig + theta @ Sig @ theta.T``.
    mu : array_like, shape (r,), optional
        Constant mean vector; defaults to zero.

    Notes
    -----
    Adjacent unit-width Brownian noise windows define the model at fractional
    gaps. Paths start in stationarity without burn-in. Parameter arrays
    are copied on input and output.
    """
    cdef varmapack_slepian_model *handle
    cdef np.ndarray thetaarr, covarr, meanarr
    cdef int rval

    def __cinit__(self):
        self.handle = NULL

    def __init__(self, object theta, object cov, object mu=None):
        cdef np.ndarray th, co, me
        cdef varmapack_slepian_model *new_handle
        cdef int r
        if np.iscomplexobj(theta):
            raise ValueError("theta must be a real matrix")
        th = _matrix_to_c(theta, "theta")
        co = _matrix_to_c(cov, "cov")
        if th is None or th.shape[0] < 1 or th.shape[0] != th.shape[1]:
            raise ValueError("theta must be a non-empty square matrix")
        r = _positive_count(th.shape[0], "model dimension")
        if co is None or co.shape[0] != r or co.shape[1] != r:
            raise ValueError("cov must have shape (r, r)")
        me = np.zeros(r) if mu is None else np.array(mu, dtype=DTYPE_F64)
        if me.ndim != 1 or me.shape[0] != r:
            raise ValueError("mu must be a vector with length r")
        new_handle = varmapack_slepian_create(<double *>th.data, <double *>co.data,
                                              <double *>me.data, r)
        if new_handle == NULL:
            raise VarmapackError(_error_message())
        varmapack_slepian_free(self.handle)
        self.handle = new_handle
        self.thetaarr, self.covarr, self.meanarr = th, co, me
        self.rval = r

    def __dealloc__(self):
        varmapack_slepian_free(self.handle)

    @property
    def r(self):
        """Number of series."""
        return self.rval

    @property
    def theta(self):
        """Unit-step MA(1) coefficient matrix."""
        return self.thetaarr.T.copy()

    @property
    def cov(self):
        """Stationary process covariance."""
        c = self.covarr.T.copy()
        i, j = np.triu_indices(self.rval, 1)
        c[i, j] = c[j, i]
        return c

    @property
    def mu(self):
        """Constant mean vector."""
        return self.meanarr.copy()

    def sim(self, object t, *, object nrep=1, object rng=None):
        """Simulate independent stationary paths at the supplied times.

        Parameters
        ----------
        t : array_like, shape (n,)
            Finite, strictly increasing observation times.
        nrep : int, optional
            Positive number of independent paths.
        rng : randompack.Rng, optional
            Stream to use instead of a temporary default Randompack generator.

        Returns
        -------
        X : ndarray, shape (nrep, n, r)
            Stationary Gaussian paths.
        """
        cdef np.ndarray times = np.ascontiguousarray(t, dtype=DTYPE_F64)
        cdef np.ndarray X
        cdef bint owned = False
        cdef randompack_rng *ptr
        cdef int count = _positive_count(nrep, "nrep")
        cdef int n
        if self.handle == NULL:
            raise ValueError("SlepianModel has not been initialized")
        if times.ndim != 1 or times.shape[0] < 1:
            raise ValueError("t must be a nonempty vector")
        n = _positive_count(times.shape[0], "length of t")
        X = np.empty((count, n, self.rval), dtype=DTYPE_F64)
        ptr = _rng_from_object(rng, &owned)
        try:
            if not varmapack_slepian_sim(self.handle, <double *>times.data,
                                         n, count, <double *>X.data, ptr):
                raise VarmapackError(_error_message())
        finally:
            if owned:
                randompack_PY_free(ptr)
        return X


class _TestcaseInfo(dict):
    def __repr__(self):
        if not self:
            return "{}"
        width = max(max(len(name) for name in self), len("name"))
        lines = [f"{'index':>5}  {'name':<{width}}  {'p':>2}  {'q':>2}  {'r':>2}"]
        for name, row in self.items():
            lines.append(
                f"{row['index']:5d}  {name:<{width}}  "
                f"{row['p']:2d}  {row['q']:2d}  {row['r']:2d}")
        return "\n".join(lines)


cdef str _error_message():
    cdef char *msg = varmapack_last_error()
    if msg == NULL:
        return "varmapack error"
    return (<bytes>msg).decode("utf-8", "replace")


cdef np.ndarray _matrix_to_c(object x, str name):
    cdef np.ndarray arr
    if x is None:
        return None
    arr = np.asarray(x, dtype=DTYPE_F64)
    if arr.ndim != 2:
        raise ValueError(f"{name} must have 2 dimensions")
    return arr.T.copy()


cdef np.ndarray _cube_to_c(object x, int r, str name):
    cdef np.ndarray arr
    if x is None:
        return None
    arr = np.asarray(x, dtype=DTYPE_F64)
    if arr.ndim == 2:
        if arr.shape[0] != r or arr.shape[1] != r:
            raise ValueError(f"{name} must have shape (r, r)")
        arr = arr.reshape((1, r, r))
    if arr.ndim != 3:
        raise ValueError(f"{name} must have shape (order, r, r)")
    if arr.shape[1] != r or arr.shape[2] != r:
        raise ValueError(f"{name} must have shape (order, r, r)")
    return arr.transpose(0, 2, 1).copy()


cdef np.ndarray _exog_to_c(object x, int r, str name):
    cdef np.ndarray arr
    if x is None:
        return None
    arr = np.asarray(x, dtype=DTYPE_F64)
    if arr.ndim == 1:
        if arr.shape[0] != r:
            raise ValueError(f"{name} must have length r")
        return arr.reshape((1, 1, r)).copy()
    if arr.ndim == 2:
        if arr.shape[1] != r:
            raise ValueError(f"{name} must have shape (s, r)")
        return arr.reshape((arr.shape[0], 1, r)).copy()
    if arr.ndim != 3 or arr.shape[2] != r:
        raise ValueError(f"{name} must have shape (s, d, r)")
    return arr.copy()


cdef randompack_rng *_rng_from_object(object rng, bint *owned):
    cdef object capsule
    cdef randompack_rng *ptr
    owned[0] = False
    if rng is None:
        ptr = randompack_PY_create(NULL)
        if ptr == NULL:
            raise MemoryError("randompack_create failed")
        owned[0] = True
        return ptr
    if not hasattr(rng, "__randompack_capsule__"):
        raise TypeError("rng must be a randompack.Rng object or None")
    capsule = rng.__randompack_capsule__()
    ptr = <randompack_rng *>PyCapsule_GetPointer(capsule, RNG_CAPSULE_CNAME)
    if ptr == NULL:
        raise RuntimeError("invalid randompack RNG capsule")
    return ptr


def testcase(object which="random", *, object p=None, object q=None,
             object r=None, double rho=0, object rng=None):
    """
    Create a testcase model.

    Parameters
    ----------
    which : str or int, optional
        Testcase selector. Strings may be named testcases, ``"random"``,
        ``"deterministic"``, or ``"rho"``. Integers select named testcases
        by one-based index.
    p, q, r : int, optional
        AR order, MA order, and series dimension. These are required for
        random, deterministic, and rho-controlled testcases.
    rho : float, optional
        Target autoregressive spectral radius when ``which="rho"``. The
        constructed spectral radius is within about ``1e-6`` of this value.
    rng : randompack.Rng, optional
        Randompack generator used for random testcases.

    Returns
    -------
    Model
        Model containing the generated VARMA testcase.
    """
    return _testcase(which, p, q, r, rho, rng)


def testcases():
    """
    List named VARMA testcases.

    Returns
    -------
    dict
        A dictionary-like object mapping testcase names to dictionaries with
        fields ``index``, ``p``, ``q``, and ``r``. When printed, it is displayed
        as a compact table.
    """
    return _testcases()


def autocov(object X, int maxlag, str norm="ML"):
    """
    Compute sample autocovariances of one observed series.

    Parameters
    ----------
    X : array_like, shape (n, r)
        Observed time series.
    maxlag : int
        Largest lag to compute.
    norm : {"ML", "C"}, optional
        Normalization. ``"ML"`` divides by ``n`` and ``"C"`` divides the
        lag-k covariance by ``n-k``.

    Returns
    -------
    C : ndarray, shape (maxlag + 1, r, r)
        Sample autocovariance matrices. ``C[k]`` estimates
        ``Cov(X[t], X[t-k])``.
    """
    return _autocov(X, maxlag, norm)


def cov2corr(object cov, object out=None):
    """
    Convert covariance matrices to correlations.

    Parameters
    ----------
    cov : array_like, shape (r, r) or (maxlag + 1, r, r)
        Covariance matrix or sequence. Lag-zero diagonal entries must be
        positive.
    out : ndarray, optional
        Preallocated output array with the same shape as ``cov``. Use
        ``out=cov`` for in-place conversion.

    Returns
    -------
    corr : ndarray
        Correlations with the same shape as the input. Lag-zero diagonal
        entries are exactly one. Other entries are not clipped to ``[-1,1]``.
    """
    return _cov2corr(cov, out)


def cov2sig(object cov, *, object A=None, object B=None, object out=None):
    """
    Recover innovation covariance from a lag-zero model covariance.

    Parameters
    ----------
    cov : array_like, shape (r, r)
        Target model covariance. Only its lower triangle is read.
    A, B : array_like, shape (order, r, r) or (r, r), optional
        AR and MA coefficient matrices. Omit absent terms. The AR model
        must be stationary.
    out : ndarray, optional
        Writable, aligned, C-contiguous float64 output of shape (r, r).
        Use ``out=cov`` for in-place conversion. Otherwise it must not
        overlap ``cov``; it must never overlap ``A`` or ``B``.

    Returns
    -------
    Sig : ndarray, shape (r, r)
        Symmetric innovation covariance. Both covariances must be positive
        semidefinite; singular covariances are allowed.

    Raises
    ------
    VarmapackError
        If the target is incompatible or the inverse problem cannot be
        solved reliably. The output is unchanged on failure.

    Examples
    --------
    >>> cov2sig([[4/3]], A=[[0.5]])
    array([[1.]])
    """
    return _cov2sig(cov, A, B, out)


cdef object _cov2sig(object cov0, object A0, object B0, object out0):
    cdef np.ndarray cov = _matrix_to_c(cov0, "cov")
    cdef np.ndarray A
    cdef np.ndarray B
    cdef np.ndarray Sig
    cdef np.ndarray out
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef int p = 0, q = 0, r
    if cov is None or cov.shape[0] <= 0 or cov.shape[1] != cov.shape[0]:
        raise ValueError("cov must have shape (r, r)")
    r = cov.shape[0]
    A = _cube_to_c(A0, r, "A")
    B = _cube_to_c(B0, r, "B")
    if A is not None:
        p = A.shape[0]
        Aptr = <double *>np.PyArray_DATA(A)
    if B is not None:
        q = B.shape[0]
        Bptr = <double *>np.PyArray_DATA(B)
    if out0 is not None:
        if not isinstance(out0, np.ndarray):
            raise TypeError("out must be a NumPy array")
        out = out0
        if out.dtype != DTYPE_F64:
            raise TypeError("out must have dtype float64")
        if out.ndim != 2 or out.shape[0] != r or out.shape[1] != r:
            raise ValueError("out must have the same shape as cov")
        if not out.flags.c_contiguous or not out.flags.aligned:
            raise ValueError("out must be C-contiguous and aligned")
        if not out.flags.writeable:
            raise ValueError("out must be writable")
        if (out0 is not cov0 and np.shares_memory(cov0, out)):
            raise ValueError("cov and out must not overlap")
        if ((A0 is not None and np.shares_memory(A0, out)) or
                (B0 is not None and np.shares_memory(B0, out))):
            raise ValueError("A or B and out must not overlap")
    Sig = np.empty((r, r), dtype=DTYPE_F64)
    if not varmapack_cov2sig(Aptr, Bptr, <double *>np.PyArray_DATA(cov),
                             p, q, r, <double *>np.PyArray_DATA(Sig)):
        raise VarmapackError(_error_message())
    if out0 is None:
        return Sig.T.copy()
    out[...] = Sig.T
    return out


cdef object _testcases():
    cdef int p = 0
    cdef int q = 0
    cdef int r = 0
    cdef int icase = 0
    cdef int ncase
    cdef int k
    cdef char namebuf[VARMAPACK_TESTCASE_NAME_LEN]
    cdef bint ok
    cdef bytes raw
    cdef str name
    namebuf[0] = 0
    namebuf[0] = ord("m")
    namebuf[1] = ord("a")
    namebuf[2] = ord("x")
    namebuf[3] = 0
    ok = varmapack_testcase(namebuf, &icase, 0, &p, &q, &r,
                            NULL, NULL, NULL, NULL)
    if not ok:
        raise VarmapackError(_error_message())
    ncase = icase
    out = _TestcaseInfo()
    for k in range(1, ncase + 1):
        p = 0
        q = 0
        r = 0
        icase = k
        namebuf[0] = 0
        ok = varmapack_testcase(namebuf, &icase, 0, &p, &q, &r,
                                NULL, NULL, NULL, NULL)
        if not ok:
            raise VarmapackError(_error_message())
        raw = bytes(namebuf).split(b"\0", 1)[0]
        name = raw.decode("utf-8", "replace")
        out[name] = {"index": k, "p": p, "q": q, "r": r}
    return out


cdef object _autocov(object X0, int maxlag, str norm):
    cdef np.ndarray X
    cdef np.ndarray C
    cdef bytes normbytes = norm.encode("utf-8")
    cdef bint ok
    cdef int n
    cdef int r
    X = np.ascontiguousarray(X0, dtype=DTYPE_F64)
    if X.ndim != 2:
        raise ValueError("X must have shape (n, r)")
    n = X.shape[0]
    r = X.shape[1]
    C = np.empty((maxlag + 1, r, r), dtype=DTYPE_F64)
    ok = varmapack_autocov(
        "N", normbytes, r, n, <double *>np.PyArray_DATA(X), maxlag,
        <double *>np.PyArray_DATA(C))
    if not ok:
        raise VarmapackError(_error_message())
    return C.transpose(0, 2, 1)


cdef object _cov2corr(object cov0, object out0):
    cdef np.ndarray cov
    cdef np.ndarray out
    cdef int axis
    cdef int maxlag
    cdef int r
    cdef bint ok
    cov = np.require(cov0, dtype=DTYPE_F64, requirements=["C", "A"])
    if cov.ndim == 2:
        if cov.shape[0] <= 0 or cov.shape[1] != cov.shape[0]:
            raise ValueError("cov must have shape (r, r)")
        r = cov.shape[0]
        maxlag = 0
    elif cov.ndim == 3:
        if (cov.shape[0] <= 0 or cov.shape[1] <= 0 or
                cov.shape[2] != cov.shape[1]):
            raise ValueError("cov must have shape (maxlag + 1, r, r)")
        maxlag = cov.shape[0] - 1
        r = cov.shape[1]
    else:
        raise ValueError("cov must have shape (r, r) or (maxlag + 1, r, r)")
    if out0 is None:
        out = np.empty_like(cov)
    else:
        if not isinstance(out0, np.ndarray):
            raise TypeError("out must be a NumPy array")
        out = out0
        if out.dtype != DTYPE_F64:
            raise TypeError("out must have dtype float64")
        if out.ndim != cov.ndim:
            raise ValueError("out must have the same shape as cov")
        for axis in range(cov.ndim):
            if out.shape[axis] != cov.shape[axis]:
                raise ValueError("out must have the same shape as cov")
        if not out.flags.c_contiguous or not out.flags.aligned:
            raise ValueError("out must be C-contiguous and aligned")
        if not out.flags.writeable:
            raise ValueError("out must be writable")
        if (np.PyArray_DATA(cov) != np.PyArray_DATA(out) and
                np.shares_memory(cov, out)):
            raise ValueError("cov and out must not overlap")
    ok = varmapack_cov2corr(
        <double *>np.PyArray_DATA(cov), r, maxlag,
        <double *>np.PyArray_DATA(out))
    if not ok:
        raise VarmapackError(_error_message())
    return out


cdef object _testcase(object which, object p0, object q0, object r0,
                      double rho, object rng):
    cdef int p = 0
    cdef int q = 0
    cdef int r = 0
    cdef int icase = 0
    cdef char namebuf[64]
    cdef bytes namebytes = b""
    cdef char *nameptr = namebuf
    cdef bint need_dims = False
    cdef randompack_rng *rngptr
    cdef bint owned
    cdef np.ndarray A
    cdef np.ndarray B
    cdef np.ndarray Sig
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef double dummy = 0
    cdef bint ok
    namebuf[0] = 0
    if isinstance(which, str):
        if which == "random":
            icase = 0
            need_dims = True
        elif which == "deterministic":
            icase = -1
            need_dims = True
        elif which == "rho":
            icase = 0
            need_dims = True
            namebytes = b"rho"
        else:
            namebytes = (<str>which).encode("utf-8")
            if len(namebytes) >= sizeof(namebuf):
                raise ValueError("testcase name is too long")
        if namebytes:
            namebuf[:len(namebytes)] = namebytes
            namebuf[len(namebytes)] = 0
    elif isinstance(which, int):
        icase = <int>which
    else:
        raise TypeError("which must be a testcase name, index, or mode")
    if p0 is not None:
        p = <int>p0
    if q0 is not None:
        q = <int>q0
    if r0 is not None:
        r = <int>r0
    if need_dims:
        if p0 is None or q0 is None or r0 is None:
            raise ValueError("p, q, and r are required for this testcase")
    else:
        ok = varmapack_testcase(nameptr, &icase, rho, &p, &q, &r,
                                NULL, NULL, NULL, NULL)
        if not ok:
            raise VarmapackError(_error_message())
    A = None
    B = None
    if p > 0:
        A = np.empty((p, r, r), dtype=DTYPE_F64)
        Aptr = <double *>np.PyArray_DATA(A)
    else:
        Aptr = &dummy
    if q > 0:
        B = np.empty((q, r, r), dtype=DTYPE_F64)
        Bptr = <double *>np.PyArray_DATA(B)
    else:
        Bptr = &dummy
    Sig = np.empty((r, r), dtype=DTYPE_F64)
    rngptr = _rng_from_object(rng, &owned)
    try:
        ok = varmapack_testcase(nameptr, &icase, rho, &p, &q, &r,
                                Aptr, Bptr,
                                <double *>np.PyArray_DATA(Sig), rngptr)
        if not ok:
            raise VarmapackError(_error_message())
    finally:
        if owned:
            randompack_PY_free(rngptr)
    return _make_model_internal(A, B, Sig, None, p, q, r)


cdef Model _make_model_internal(np.ndarray A, np.ndarray B, np.ndarray Sig,
                                np.ndarray mu, int p, int q, int r):
    cdef Model model = Model.__new__(Model)
    model.Aarr = A
    model.Barr = B
    model.Carr = None
    model.Sigarr = Sig
    model.muarr = mu
    model.pval = p
    model.qval = q
    model.sval = 0
    model.dval = 0
    model.rval = r
    model.nmu = 0
    return model


cdef class Model:
    """
    Gaussian VARMA or VARMAX model.

    Parameters
    ----------
    A : array_like, shape (p, r, r) or (r, r), optional
        Autoregressive coefficient matrices.
    B : array_like, shape (q, r, r) or (r, r), optional
        Moving-average coefficient matrices.
    C : array_like, shape (s, d, r), optional
        Exogenous coefficient matrices. ``C[k, j, i]`` is element ``(i, j)``
        of the r-by-d matrix multiplying the d-vector ``z[t-k]``. Supplying
        ``C`` creates a VARMAX model and requires ``z`` and ``X0`` when
        simulating.
    Sig : array_like, shape (r, r)
        Innovation covariance matrix.
    mu : array_like, shape (r,) or (nmu, r), optional
        Time series means for VARMA simulation. If more than one mean vector
        is supplied, the last supplied mean repeats to the end of the series.
        ``mu`` is not supported for VARMAX models.

    Notes
    -----
    Without exogenous terms, the model is

    ``x_t - mu_t = eps_t + sum_i A_i (x_{t-i} - mu_{t-i})
    + sum_j B_j eps_{t-j}``.

    With exogenous terms, Varmapack simulates

    ``x_t = eps_t + sum_i A_i x_{t-i} + sum_j B_j eps_{t-j}
    + sum_k C_k z_{t-k+1}``.

    Python input arrays use lag-first shapes. Simulated series are returned
    with shape ``(nrep, length, r)``.
    """
    cdef np.ndarray Aarr
    cdef np.ndarray Barr
    cdef np.ndarray Carr
    cdef np.ndarray Sigarr
    cdef np.ndarray muarr
    cdef int pval
    cdef int qval
    cdef int sval
    cdef int dval
    cdef int rval
    cdef int nmu

    def __init__(self, *, object A=None, object B=None, object C=None,
                 object Sig=None, object mu=None):
        cdef np.ndarray muarr
        self.Sigarr = _matrix_to_c(Sig, "Sig")
        self.muarr = None
        self.nmu = 0
        self.pval = 0
        self.qval = 0
        self.sval = 0
        self.dval = 0
        if self.Sigarr is None:
            raise ValueError("Sig must be supplied")
        if self.Sigarr.shape[0] != self.Sigarr.shape[1]:
            raise ValueError("Sig must be square")
        self.rval = self.Sigarr.shape[0]
        self.Aarr = _cube_to_c(A, self.rval, "A")
        self.Barr = _cube_to_c(B, self.rval, "B")
        self.Carr = _exog_to_c(C, self.rval, "C")
        if self.Aarr is not None:
            self.pval = self.Aarr.shape[0]
        if self.Barr is not None:
            self.qval = self.Barr.shape[0]
        if self.Carr is not None:
            self.sval = self.Carr.shape[0]
            self.dval = self.Carr.shape[1]
            if mu is not None:
                raise ValueError("mu is not supported for VARMAX models")
        if mu is not None:
            muarr = np.asarray(mu, dtype=DTYPE_F64)
            if muarr.ndim == 1:
                if muarr.shape[0] != self.rval:
                    raise ValueError("mu must have length r")
                self.muarr = muarr.copy()
                self.nmu = 1
            elif muarr.ndim == 2:
                if muarr.shape[1] != self.rval:
                    raise ValueError("mu must have shape (nmu, r)")
                self.muarr = muarr.copy()
                self.nmu = muarr.shape[0]
            else:
                raise ValueError("mu must have 1 or 2 dimensions")

    @property
    def p(self):
        """Number of autoregressive coefficient matrices."""
        return self.pval

    @property
    def q(self):
        """Number of moving-average coefficient matrices."""
        return self.qval

    @property
    def r(self):
        """Dimension of each observation vector."""
        return self.rval

    @property
    def s(self):
        """Number of exogenous coefficient matrices."""
        return self.sval

    @property
    def d(self):
        """Dimension of each exogenous input vector."""
        return self.dval

    @property
    def A(self):
        """Autoregressive coefficient matrices, or ``None``."""
        if self.Aarr is None:
            return None
        return self.Aarr.transpose(0, 2, 1).copy()

    @property
    def B(self):
        """Moving-average coefficient matrices, or ``None``."""
        if self.Barr is None:
            return None
        return self.Barr.transpose(0, 2, 1).copy()

    @property
    def C(self):
        """Exogenous coefficient matrices, or ``None``."""
        if self.Carr is None:
            return None
        if self.dval == 1:
            return self.Carr[:, 0, :].copy()
        return self.Carr.copy()

    @property
    def Sig(self):
        """Innovation covariance matrix."""
        return self.Sigarr.T.copy()

    @property
    def mu(self):
        """Time series means supplied at construction, or ``None``."""
        if self.muarr is None:
            return None
        if self.nmu == 1:
            return self.muarr.copy()
        return self.muarr.copy()

    def sim(self, int length, *, int nrep=1, object X0=None, object z=None,
            object rng=None, bint return_shocks=False, object out=None):
        """
        Simulate time series from the model.

        Parameters
        ----------
        length : int
            Number of observations to return in each simulated path.
        nrep : int, optional
            Number of independent simulated paths.
        X0 : array_like, optional
            Fixed startup values. For VARMA, shapes ``(nX0, r)`` and
            ``(nrep, nX0, r)`` are accepted. For VARMAX, ``X0`` is required
            when ``max(p, q, s - 1)`` is positive and must contain at least
            that many fixed startup values. Nonstationary VARMA models also
            require ``X0``; with MA terms, ``Sig`` must be positive definite.
        z : array_like, optional
            Exogenous input for VARMAX models. Scalar inputs use shapes
            ``(length,)`` and ``(nrep, length)``. Vector inputs use
            ``(length, d)`` and ``(nrep, length, d)``.
        rng : randompack.Rng, optional
            Randompack generator. If omitted, a temporary default generator is
            created for the call.
        return_shocks : bool, optional
            If true, return both simulated series and innovations.
        out : ndarray, optional
            Writable, C-contiguous float64 array with shape
            ``(nrep, length, r)`` in which to store the simulated series.

        Returns
        -------
        X : ndarray, shape (nrep, length, r)
            Simulated series.
        E : ndarray, shape (nrep, length, r)
            Simulated innovations. Returned only when ``return_shocks=True``.
        """
        return _sim_model(self, length, nrep, X0, z, rng, return_shocks, out)

    def acvf(self, int maxlag):
        """
        Compute theoretical autocovariances of the VARMA part.

        Parameters
        ----------
        maxlag : int
            Largest lag to compute.

        Returns
        -------
        Gamma : ndarray, shape (maxlag + 1, r, r)
            Autocovariance sequence. ``Gamma[k]`` is
            ``Cov(x_t, x_{t-k})`` for the stationary VARMA model.
        """
        return _acvf_model(self, maxlag)

    def psi(self, int maxlag):
        """
        Compute VARMA impulse-response coefficient matrices.

        Parameters
        ----------
        maxlag : int
            Largest impulse-response lag to compute.

        Returns
        -------
        Psi : ndarray, shape (maxlag + 1, r, r)
            Coefficients satisfying
            ``x_t = sum_j Psi[j] eps_{t-j}`` for the zero-mean VARMA model.
        """
        return _psi_model(self, maxlag)

    def irf(self, int maxlag):
        """
        Compute orthogonalized impulse-response matrices.

        Parameters
        ----------
        maxlag : int
            Largest impulse-response lag to compute.

        Returns
        -------
        Theta : ndarray, shape (maxlag + 1, r, r)
            Orthogonalized responses ``Theta[j] = Psi[j] L``, where
            ``Sig = L L.T``. Positive semidefinite ``Sig`` is allowed.
        """
        return _irf_model(self, maxlag)

    def specrad(self):
        """
        Compute the spectral radius of the autoregressive companion matrix.

        Returns
        -------
        rho : float
            Spectral radius of the VAR companion matrix. Pure VMA and white
            noise models have spectral radius zero.
        """
        return _specrad_model(self)

    def ma_specrad(self):
        """
        Compute the spectral radius of the moving-average companion matrix.

        Returns
        -------
        rho : float
            Spectral radius of the moving-average companion matrix. Pure VAR
            and white noise models have moving-average spectral radius zero.
        """
        return _ma_specrad_model(self)


cdef object _acvf_model(Model model, int maxlag):
    cdef np.ndarray Gamma
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef double dummy = 0
    cdef bint ok
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    else:
        Aptr = &dummy
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    else:
        Bptr = &dummy
    Gamma = np.empty((maxlag + 1, model.rval, model.rval), dtype=DTYPE_F64)
    ok = varmapack_acvf(Aptr, Bptr, <double *>np.PyArray_DATA(model.Sigarr),
                         model.pval, model.qval, model.rval,
                         <double *>np.PyArray_DATA(Gamma), maxlag)
    if not ok:
        raise VarmapackError(_error_message())
    return Gamma.transpose(0, 2, 1).copy()


cdef object _psi_model(Model model, int maxlag):
    cdef np.ndarray Psi
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef bint ok
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    Psi = np.empty((maxlag + 1, model.rval, model.rval), dtype=DTYPE_F64)
    ok = varmapack_psi(Aptr, Bptr, model.pval, model.qval, model.rval, maxlag,
                       <double *>np.PyArray_DATA(Psi))
    if not ok:
        raise VarmapackError(_error_message())
    return Psi.transpose(0, 2, 1).copy()


cdef object _irf_model(Model model, int maxlag):
    cdef np.ndarray Theta
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef bint ok
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    Theta = np.empty((maxlag + 1, model.rval, model.rval), dtype=DTYPE_F64)
    ok = varmapack_irf(Aptr, Bptr, <double *>np.PyArray_DATA(model.Sigarr),
                        model.pval, model.qval, model.rval, maxlag,
                        <double *>np.PyArray_DATA(Theta))
    if not ok:
        raise VarmapackError(_error_message())
    return Theta.transpose(0, 2, 1).copy()


cdef double _specrad_model(Model model):
    cdef double *Aptr = NULL
    cdef double rho
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    rho = varmapack_specrad(Aptr, model.rval, model.pval)
    if rho != rho:
        raise VarmapackError(_error_message())
    return rho


cdef double _ma_specrad_model(Model model):
    cdef double *Bptr = NULL
    cdef double rho
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    rho = varmapack_ma_specrad(Bptr, model.rval, model.qval)
    if rho != rho:
        raise VarmapackError(_error_message())
    return rho


cdef object _sim_model(Model model, int length, int nrep, object X0, object z,
                       object rng, bint return_shocks, object out0):
    cdef np.ndarray X0arr
    cdef np.ndarray zarr
    cdef np.ndarray X
    cdef np.ndarray E
    cdef np.ndarray X0public
    cdef np.ndarray zpublic
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef double *Cptr = NULL
    cdef double *muptr = NULL
    cdef double *X0ptr = NULL
    cdef double *zptr = NULL
    cdef double *Eptr = NULL
    cdef randompack_rng *rngptr
    cdef bint owned
    cdef bint ok
    cdef int nX0 = 0
    cdef int MX0 = 1
    cdef int Mz = 1
    cdef int nmu = 0
    cdef int h
    cdef int min_h
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    if model.Carr is not None:
        Cptr = <double *>np.PyArray_DATA(model.Carr)
    if model.muarr is not None:
        muptr = <double *>np.PyArray_DATA(model.muarr)
        nmu = model.nmu
    if X0 is not None:
        X0public = np.asarray(X0, dtype=DTYPE_F64)
        if X0public.ndim == 2:
            if X0public.shape[1] != model.rval:
                raise ValueError("X0 must have shape (nX0, r)")
            nX0 = X0public.shape[0]
            X0arr = X0public.copy()
        elif X0public.ndim == 3:
            if X0public.shape[0] != nrep or X0public.shape[2] != model.rval:
                raise ValueError("X0 must have shape (nrep, nX0, r)")
            nX0 = X0public.shape[1]
            MX0 = nrep
            X0arr = X0public.copy()
        else:
            raise ValueError("X0 must have shape (nX0, r) or (nrep, nX0, r)")
        X0ptr = <double *>np.PyArray_DATA(X0arr)
    if z is not None:
        zpublic = np.asarray(z, dtype=DTYPE_F64)
        if model.dval == 1 and zpublic.ndim == 1:
            if zpublic.shape[0] < length:
                raise ValueError("z must have at least length entries")
            zarr = zpublic[:length].copy()
        elif model.dval == 1 and zpublic.ndim == 2:
            if zpublic.shape[0] != nrep or zpublic.shape[1] < length:
                raise ValueError("z must have shape (nrep, length)")
            Mz = nrep
            zarr = zpublic[:, :length].copy()
        elif zpublic.ndim == 2:
            if zpublic.shape[0] < length or zpublic.shape[1] != model.dval:
                raise ValueError("z must have shape (length, d)")
            zarr = zpublic[:length, :].copy()
        elif zpublic.ndim == 3:
            if (zpublic.shape[0] != nrep or zpublic.shape[1] < length or
                    zpublic.shape[2] != model.dval):
                raise ValueError("z must have shape (nrep, length, d)")
            Mz = nrep
            zarr = zpublic[:, :length, :].copy()
        else:
            raise ValueError("z has an invalid shape")
        zptr = <double *>np.PyArray_DATA(zarr)
    if model.Carr is not None:
        if zptr == NULL:
            raise ValueError("z must be supplied for VARMAX simulation")
        min_h = max(model.pval, model.qval, model.sval - 1)
        if X0ptr == NULL and min_h > 0:
            raise ValueError("X0 must be supplied for VARMAX simulation")
        if nX0 < min_h:
            raise ValueError("X0 has fewer than max(p, q, s - 1) values")
    elif zptr != NULL:
        raise ValueError("z can only be supplied when C is present")
    if out0 is None:
        X = np.empty((nrep, length, model.rval), dtype=DTYPE_F64, order="C")
    else:
        if not isinstance(out0, np.ndarray):
            raise TypeError("out must be a NumPy array")
        X = out0
        if X.dtype != DTYPE_F64:
            raise TypeError("out must have dtype float64")
        if (X.ndim != 3 or X.shape[0] != nrep or X.shape[1] != length or
                X.shape[2] != model.rval):
            raise ValueError("out must have shape (nrep, length, r)")
        if not X.flags.c_contiguous or not X.flags.aligned:
            raise ValueError("out must be C-contiguous and aligned")
        if not X.flags.writeable:
            raise ValueError("out must be writable")
    if return_shocks:
        E = np.empty((nrep, length, model.rval), dtype=DTYPE_F64, order="C")
        Eptr = <double *>np.PyArray_DATA(E)
    else:
        E = None
    rngptr = _rng_from_object(rng, &owned)
    try:
        if model.Carr is not None:
            h = nX0
            ok = varmapack_simx(Aptr, Bptr, Cptr,
                                 <double *>np.PyArray_DATA(model.Sigarr),
                                 zptr, Mz, model.pval, model.qval, model.sval, model.dval,
                                 model.rval, length, nrep, X0ptr, h, MX0,
                                 <double *>np.PyArray_DATA(X), Eptr, rngptr)
        else:
            ok = varmapack_sim(Aptr, Bptr,
                                <double *>np.PyArray_DATA(model.Sigarr), muptr,
                                nmu, model.pval, model.qval, model.rval, length,
                                nrep, X0ptr, nX0, MX0,
                                <double *>np.PyArray_DATA(X), Eptr, rngptr)
        if not ok:
            raise VarmapackError(_error_message())
    finally:
        if owned:
            randompack_PY_free(rngptr)
    if return_shocks:
        return X, E
    return X
