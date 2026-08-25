import numpy as np
import pytest
import randompack
import varmapack


def test_model_argument_errors():
    Sig = np.eye(2)
    A = np.zeros((1, 2, 2))
    B = np.zeros((1, 2, 2))
    C = np.zeros((1, 2))
    with pytest.raises(TypeError):
        varmapack.Model(A, B, Sig)
    with pytest.raises(ValueError):
        varmapack.Model()
    with pytest.raises(ValueError):
        varmapack.Model(Sig=np.ones((2, 3)))
    with pytest.raises(ValueError):
        varmapack.Model(A=np.zeros((3, 3)), Sig=Sig)
    with pytest.raises(ValueError):
        varmapack.Model(A=np.zeros((1, 3, 3)), Sig=Sig)
    with pytest.raises(ValueError):
        varmapack.Model(Sig=Sig, mu=np.ones(3))
    with pytest.raises(ValueError):
        varmapack.Model(Sig=Sig, mu=np.ones((2, 3)))
    with pytest.raises(ValueError):
        varmapack.Model(C=C, Sig=Sig, mu=np.zeros(2))


def test_simulation_argument_errors():
    model = varmapack.Model(
        A=np.zeros((1, 2, 2)), B=np.zeros((1, 2, 2)), Sig=np.eye(2))
    rng = randompack.Rng()
    with pytest.raises(ValueError):
        model.sim(4, X0=np.zeros((1, 3)), rng=rng)
    with pytest.raises(ValueError):
        model.sim(4, nrep=2, X0=np.zeros((3, 1, 2)), rng=rng)
    with pytest.raises(ValueError):
        model.sim(4, z=np.zeros(4), rng=rng)
    with pytest.raises(TypeError):
        model.sim(4, rng=object())
    with pytest.raises((ValueError, varmapack.VarmapackError)):
        model.sim(0, rng=rng)
    with pytest.raises((ValueError, varmapack.VarmapackError)):
        model.sim(4, nrep=0, rng=rng)
    with pytest.raises(TypeError):
        model.sim(4, nrep=2, rng=rng, out=np.empty((2, 4, 2), dtype=np.float32))
    with pytest.raises(TypeError):
        model.sim(4, nrep=2, rng=rng, out=[[[0.0, 0.0]]*4]*2)
    with pytest.raises(ValueError):
        model.sim(4, nrep=2, rng=rng, out=np.empty((1, 4, 2)))
    with pytest.raises(ValueError):
        model.sim(4, nrep=2, rng=rng, out=np.empty((2, 4, 4))[:, :, ::2])
    readonly = np.empty((2, 4, 2))
    readonly.flags.writeable = False
    with pytest.raises(ValueError):
        model.sim(4, nrep=2, rng=rng, out=readonly)


def test_derived_quantity_errors():
    model = varmapack.Model(Sig=np.eye(2))
    with pytest.raises(ValueError):
        model.acvf(-1)
    with pytest.raises(ValueError):
        model.psi(-1)
    with pytest.raises(ValueError):
        model.irf(-1)
    unstable = varmapack.Model(A=np.array([[[1.1]]]), Sig=np.ones((1, 1)))
    with pytest.raises(varmapack.VarmapackError):
        unstable.sim(4)
    with pytest.raises(varmapack.VarmapackError):
        unstable.acvf(0)
    bad_ar = varmapack.Model(A=np.array([[[np.nan]]]), Sig=np.ones((1, 1)))
    bad_ma = varmapack.Model(B=np.array([[[np.inf]]]), Sig=np.ones((1, 1)))
    with pytest.raises(varmapack.VarmapackError):
        bad_ar.specrad()
    with pytest.raises(varmapack.VarmapackError):
        bad_ma.ma_specrad()
    non_psd = varmapack.Model(Sig=np.array([[1.0, 2.0], [2.0, 1.0]]))
    with pytest.raises(varmapack.VarmapackError):
        non_psd.irf(0)


def test_autocovariance_errors():
    with pytest.raises(ValueError):
        varmapack.autocov(np.zeros((2, 2, 1)), 1)
    with pytest.raises(varmapack.VarmapackError):
        varmapack.autocov(np.zeros((3, 2)), 1, "bad")
    with pytest.raises(varmapack.VarmapackError):
        varmapack.autocov(np.zeros((3, 2)), -1)
    with pytest.raises(varmapack.VarmapackError):
        varmapack.autocov(np.zeros((3, 2)), 3)


def test_covariance_to_correlation_errors():
    cov = np.eye(2)
    with pytest.raises(TypeError):
        varmapack.cov2corr(cov, out=[[0.0, 0.0], [0.0, 0.0]])
    with pytest.raises(TypeError):
        varmapack.cov2corr(cov, out=np.empty((2, 2), dtype=np.float32))
    with pytest.raises(ValueError):
        varmapack.cov2corr(cov, out=np.empty((2, 3)))
    with pytest.raises(ValueError):
        varmapack.cov2corr(cov, out=np.empty((2, 2)).T)
    readonly = np.empty((2, 2))
    readonly.flags.writeable = False
    with pytest.raises(ValueError):
        varmapack.cov2corr(cov, out=readonly)
    storage = np.empty(5)
    overlapping_cov = storage[:4].reshape(2, 2)
    overlapping_out = storage[1:].reshape(2, 2)
    with pytest.raises(ValueError):
        varmapack.cov2corr(overlapping_cov, out=overlapping_out)
    with pytest.raises(varmapack.VarmapackError):
        varmapack.cov2corr(np.array([[0.0, 0.0], [0.0, 1.0]]))


def test_varmax_argument_errors():
    Sig = np.eye(2)
    A = np.zeros((1, 2, 2))
    B = np.zeros((1, 2, 2))
    C = np.zeros((1, 2))
    model = varmapack.Model(A=A, B=B, C=C, Sig=Sig)
    rng = randompack.Rng()
    with pytest.raises(ValueError):
        model.sim(4, X0=np.zeros((2, 2)), rng=rng)
    with pytest.raises(ValueError):
        model.sim(4, z=np.zeros(4), rng=rng)
    with pytest.raises(ValueError):
        model.sim(4, X0=np.zeros((2, 2)), z=np.zeros(3), rng=rng)
    with pytest.raises(ValueError):
        model.sim(4, nrep=2, X0=np.zeros((2, 2)), z=np.zeros((3, 4)), rng=rng)


def test_testcase_errors():
    with pytest.raises(varmapack.VarmapackError):
        varmapack.testcase("does-not-exist")
    with pytest.raises(TypeError):
        varmapack.testcase(object())
    with pytest.raises(ValueError):
        varmapack.testcase("rho")
    with pytest.raises(ValueError):
        varmapack.testcase("x"*64)
    with pytest.raises(TypeError):
        varmapack.testcase("random", p=1, q=1, r=1, rng=object())
