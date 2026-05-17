import numpy as np
cimport numpy as np

from libc.stddef cimport size_t
from cpython.pycapsule cimport PyCapsule_GetPointer

DTYPE_F64 = np.dtype(np.float64)
RNG_CAPSULE_NAME = b"randompack.randompack_rng"
cdef const char *RNG_CAPSULE_CNAME = "randompack.randompack_rng"

cdef extern from "randompack.h":
    ctypedef struct randompack_rng
    randompack_rng *randompack_create(const char *engine)
    void randompack_free(randompack_rng *rng)

cdef extern from "varmapack.h":
    ctypedef enum varmapack_error:
        VARMAPACK_OK
        VARMAPACK_INVALID_ARGUMENT
        VARMAPACK_DIMENSION
        VARMAPACK_ALLOCATION
        VARMAPACK_UNKNOWN_TESTCASE
        VARMAPACK_NONSTATIONARY
        VARMAPACK_NONSTATIONARY_MA
        VARMAPACK_SINGULAR
        VARMAPACK_NOT_POSITIVE_SEMIDEFINITE
        VARMAPACK_INTERNAL
    const char *varmapack_strerror(varmapack_error error)
    varmapack_error varmapack_sim(double *A, double *B, double *Sig,
                                  double *mu, int nmu, int p, int q, int r, int n,
                                  int M, double *X0, int nX0,
                                  double *X, double *E, randompack_rng *rng)
    varmapack_error varmapack_testcase(double *A, double *B, double *Sig,
                                       char *name, int *pp, int *qp, int *rp,
                                       int *icase, double rho, randompack_rng *rng)
    varmapack_error varmapack_acvf(double *A, double *B, double *Sig,
                                   int p, int q, int r, double *Gamma, int maxlag)
    varmapack_error varmapack_psi(double *A, double *B, int p, int q,
                                  int r, int h, double *Psi)
    varmapack_error varmapack_irf(double *A, double *B, double *Sig,
                                  int p, int q, int r, int h, double *Theta)
    varmapack_error varmapack_autocov(const char *transp, const char *norm,
                                      int r, int n, double *X, int maxlag, double *C)

np.import_array()


class VarmapackError(RuntimeError):
    pass


cdef str _error_message(varmapack_error error):
    cdef const char *msg = varmapack_strerror(error)
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
    if arr.ndim != 3:
        raise ValueError(f"{name} must have 3 dimensions")
    if arr.shape[1] != r or arr.shape[2] != r:
        raise ValueError(f"{name} must have shape (order, r, r)")
    return arr.transpose(0, 2, 1).copy()


cdef randompack_rng *_rng_from_object(object rng, bint *owned):
    cdef object capsule
    cdef randompack_rng *ptr
    owned[0] = False
    if rng is None:
        ptr = randompack_create(NULL)
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
    Create a Varmapack testcase and return it as a Model.

    `which` may be a named testcase, a one-based integer testcase index,
    "random", "deterministic", or "rho". The dimension arguments `p`, `q`,
    and `r` are needed for random, deterministic, and rho testcases.
    """
    return _testcase(which, p, q, r, rho, rng)


def autocov(object X, int maxlag, str norm="ML"):
    """
    Compute sample autocovariances of one series.

    `X` has shape `(n, r)`. The result has shape `(maxlag+1, r, r)`.
    """
    return _autocov(X, maxlag, norm)


cdef object _autocov(object X0, int maxlag, str norm):
    cdef np.ndarray X
    cdef np.ndarray Xc
    cdef np.ndarray C
    cdef bytes normbytes = norm.encode("utf-8")
    cdef varmapack_error error
    cdef int n
    cdef int r
    X = np.asarray(X0, dtype=DTYPE_F64)
    if X.ndim != 2:
        raise ValueError("X must have shape (n, r)")
    n = X.shape[0]
    r = X.shape[1]
    Xc = X.T.copy()
    C = np.empty((maxlag + 1, r, r), dtype=DTYPE_F64)
    error = varmapack_autocov(
        "T", normbytes, r, n, <double *>np.PyArray_DATA(Xc), maxlag,
        <double *>np.PyArray_DATA(C))
    if error != VARMAPACK_OK:
        raise VarmapackError(_error_message(error))
    return C.transpose(0, 2, 1).copy()


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
    cdef varmapack_error error
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
        error = varmapack_testcase(NULL, NULL, NULL, nameptr, &p, &q, &r,
                                   &icase, rho, NULL)
        if error != VARMAPACK_OK:
            raise VarmapackError(_error_message(error))
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
        error = varmapack_testcase(Aptr, Bptr, <double *>np.PyArray_DATA(Sig),
                                   nameptr, &p, &q, &r, &icase, rho, rngptr)
        if error != VARMAPACK_OK:
            raise VarmapackError(_error_message(error))
    finally:
        if owned:
            randompack_free(rngptr)
    return _make_model_internal(A, B, Sig, None, p, q, r)


cdef Model _make_model_internal(np.ndarray A, np.ndarray B, np.ndarray Sig,
                                np.ndarray mu, int p, int q, int r):
    cdef Model model = Model.__new__(Model)
    model.Aarr = A
    model.Barr = B
    model.Sigarr = Sig
    model.muarr = mu
    model.p = p
    model.q = q
    model.r = r
    model.nmu = 0
    return model


cdef class Model:
    """
    VARMA model.

    Parameters use Python's lag-first convention. `A` has shape `(p, r, r)`,
    `B` has shape `(q, r, r)`, and `Sig` has shape `(r, r)`.
    """
    cdef np.ndarray Aarr
    cdef np.ndarray Barr
    cdef np.ndarray Sigarr
    cdef np.ndarray muarr
    cdef int p
    cdef int q
    cdef int r
    cdef int nmu

    def __init__(self, object A=None, object B=None, object Sig=None, object mu=None):
        cdef np.ndarray muarr
        self.Sigarr = _matrix_to_c(Sig, "Sig")
        self.muarr = None
        self.nmu = 0
        self.p = 0
        self.q = 0
        if self.Sigarr is None:
            raise ValueError("Sig must be supplied")
        if self.Sigarr.shape[0] != self.Sigarr.shape[1]:
            raise ValueError("Sig must be square")
        self.r = self.Sigarr.shape[0]
        self.Aarr = _cube_to_c(A, self.r, "A")
        self.Barr = _cube_to_c(B, self.r, "B")
        if self.Aarr is not None:
            self.p = self.Aarr.shape[0]
        if self.Barr is not None:
            self.q = self.Barr.shape[0]
        if mu is not None:
            muarr = np.asarray(mu, dtype=DTYPE_F64)
            if muarr.ndim == 1:
                if muarr.shape[0] != self.r:
                    raise ValueError("mu must have length r")
                self.muarr = muarr.copy()
                self.nmu = 1
            elif muarr.ndim == 2:
                if muarr.shape[1] != self.r:
                    raise ValueError("mu must have shape (nmu, r)")
                self.muarr = muarr.copy()
                self.nmu = muarr.shape[0]
            else:
                raise ValueError("mu must have 1 or 2 dimensions")

    @property
    def order(self):
        return self.p, self.q

    @property
    def dimension(self):
        return self.r

    @property
    def A(self):
        if self.Aarr is None:
            return None
        return self.Aarr.transpose(0, 2, 1).copy()

    @property
    def B(self):
        if self.Barr is None:
            return None
        return self.Barr.transpose(0, 2, 1).copy()

    @property
    def Sig(self):
        return self.Sigarr.T.copy()

    @property
    def mu(self):
        if self.muarr is None:
            return None
        if self.nmu == 1:
            return self.muarr.copy()
        return self.muarr.copy()

    def sim(self, int length, *, int nrep=1, object X0=None, object rng=None,
            bint return_shocks=False):
        """
        Simulate from this VARMA model.

        Optional `X0` has shape `(nX0, r)`. The returned `X` has shape
        `(nrep, length, r)`. If `return_shocks` is true, return `(X, E)`.
        """
        return _sim_model(self, length, nrep, X0, rng, return_shocks)

    def acvf(self, int maxlag):
        """
        Compute theoretical autocovariances.

        The result has shape `(maxlag+1, r, r)`.
        """
        return _acvf_model(self, maxlag)

    def psi(self, int maxlag):
        """
        Compute impulse-response coefficient matrices.

        The result has shape `(maxlag+1, r, r)`.
        """
        return _psi_model(self, maxlag)

    def irf(self, int maxlag):
        """
        Compute orthogonalized impulse-response matrices.

        The result has shape `(maxlag+1, r, r)`.
        """
        return _irf_model(self, maxlag)


cdef object _acvf_model(Model model, int maxlag):
    cdef np.ndarray Gamma
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef double dummy = 0
    cdef varmapack_error error
    cdef int call_maxlag = maxlag
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if call_maxlag < model.p:
        call_maxlag = model.p
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    else:
        Aptr = &dummy
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    else:
        Bptr = &dummy
    Gamma = np.empty((call_maxlag + 1, model.r, model.r), dtype=DTYPE_F64)
    error = varmapack_acvf( Aptr, Bptr, <double *>np.PyArray_DATA(model.Sigarr), model.p,
        model.q, model.r, <double *>np.PyArray_DATA(Gamma), call_maxlag)
    if error != VARMAPACK_OK:
        raise VarmapackError(_error_message(error))
    return Gamma[:maxlag + 1].transpose(0, 2, 1).copy()


cdef object _psi_model(Model model, int maxlag):
    cdef np.ndarray Psi
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef varmapack_error error
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    Psi = np.empty((maxlag + 1, model.r, model.r), dtype=DTYPE_F64)
    error = varmapack_psi(Aptr, Bptr, model.p, model.q, model.r, maxlag,
                          <double *>np.PyArray_DATA(Psi))
    if error != VARMAPACK_OK:
        raise VarmapackError(_error_message(error))
    return Psi.transpose(0, 2, 1).copy()


cdef object _irf_model(Model model, int maxlag):
    cdef np.ndarray Theta
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef varmapack_error error
    if maxlag < 0:
        raise ValueError("maxlag must be nonnegative")
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    Theta = np.empty((maxlag + 1, model.r, model.r), dtype=DTYPE_F64)
    error = varmapack_irf(Aptr, Bptr, <double *>np.PyArray_DATA(model.Sigarr),
                          model.p, model.q, model.r, maxlag,
                          <double *>np.PyArray_DATA(Theta))
    if error != VARMAPACK_OK:
        raise VarmapackError(_error_message(error))
    return Theta.transpose(0, 2, 1).copy()


cdef object _sim_model(Model model, int length, int nrep, object X0, object rng,
                       bint return_shocks):
    cdef np.ndarray X0arr
    cdef np.ndarray X
    cdef np.ndarray E
    cdef double *Aptr = NULL
    cdef double *Bptr = NULL
    cdef double *muptr = NULL
    cdef double *X0ptr = NULL
    cdef double *Eptr = NULL
    cdef randompack_rng *rngptr
    cdef bint owned
    cdef varmapack_error error
    cdef int nX0 = 0
    cdef int nmu = 0
    cdef np.ndarray X0public
    if model.Aarr is not None:
        Aptr = <double *>np.PyArray_DATA(model.Aarr)
    if model.Barr is not None:
        Bptr = <double *>np.PyArray_DATA(model.Barr)
    if model.muarr is not None:
        muptr = <double *>np.PyArray_DATA(model.muarr)
        nmu = model.nmu
    if X0 is not None:
        X0public = np.asarray(X0, dtype=DTYPE_F64)
        if X0public.ndim != 2 or X0public.shape[1] != model.r:
            raise ValueError("X0 must have shape (nX0, r)")
        nX0 = X0public.shape[0]
        X0arr = X0public.T.copy()
        X0ptr = <double *>np.PyArray_DATA(X0arr)
    X = np.empty((nrep, length, model.r), dtype=DTYPE_F64, order="C")
    if return_shocks:
        E = np.empty((nrep, length, model.r), dtype=DTYPE_F64, order="C")
        Eptr = <double *>np.PyArray_DATA(E)
    else:
        E = None
    rngptr = _rng_from_object(rng, &owned)
    try:
        error = varmapack_sim( Aptr, Bptr, <double *>np.PyArray_DATA(model.Sigarr), muptr,
            nmu, model.p, model.q, model.r, length, nrep, X0ptr, nX0,
            <double *>np.PyArray_DATA(X), Eptr, rngptr)
        if error != VARMAPACK_OK:
            raise VarmapackError(_error_message(error))
    finally:
        if owned:
            randompack_free(rngptr)
    if return_shocks:
        return X, E
    return X
