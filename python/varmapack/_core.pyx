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
    double varmapack_specrad(double *A, int r, int p)
    double varmapack_ma_specrad(double *B, int r, int q)

randompack_PY_initialize()
np.import_array()


class VarmapackError(RuntimeError):
    """Exception raised when the Varmapack C library reports an error."""
    pass


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
