import numpy as np
import randompack
import varmapack


def _rng(seed):
    rng = randompack.Rng()
    rng.seed(seed)
    return rng


def test_white_noise_and_mean_path():
    Sig = np.array([[1.0, 0.25], [0.25, 2.0]])
    white = varmapack.Model(Sig=Sig)
    X, E = white.sim(4, return_shocks=True)
    assert X.shape == (1, 4, 2)
    assert E.shape == (1, 4, 2)
    assert np.allclose(X, E)
    mu = np.array([[10.0, 20.0], [11.0, 21.0], [12.0, 22.0]])
    white_mu = varmapack.Model(Sig=Sig, mu=mu)
    X, E = white_mu.sim(5, return_shocks=True)
    expected_mu = np.vstack([mu, mu[-1], mu[-1]])
    assert np.allclose(X[0] - E[0], expected_mu)
    assert np.array_equal(white_mu.mu, mu)


def test_var_model_properties_and_simulation():
    A = np.array([[[0.4]]])
    Sig = np.array([[1.0]])
    model = varmapack.Model(A=A, Sig=Sig)
    assert model.p == 1
    assert model.q == 0
    assert model.r == 1
    assert model.s == 0
    assert model.d == 0
    assert model.B is None
    assert model.C is None
    assert model.mu is None
    assert np.isclose(model.specrad(), 0.4)
    assert np.isclose(model.ma_specrad(), 0)
    rng = _rng(123)
    X1 = model.sim(5, rng=rng)
    assert X1.shape == (1, 5, 1)
    assert np.isfinite(X1).all()
    rng.seed(123)
    X2 = model.sim(5, rng=rng)
    assert np.allclose(X1, X2)
    X3 = model.sim(length=5, nrep=3, rng=rng)
    assert X3.shape == (3, 5, 1)


def test_shared_and_replicate_starting_values():
    A = np.array([[[0.1, 0.2], [0.3, 0.4]]])
    B = np.array([[[0.5, 0.6], [0.7, 0.8]]])
    Sig = np.array([[1.0, 0.25], [0.25, 2.0]])
    model = varmapack.Model(A=A, B=B, Sig=Sig)
    rng = _rng(123)
    X0shared = np.array([[1.0, 2.0], [3.0, 4.0]])
    Xshared = model.sim(4, nrep=2, X0=X0shared, rng=rng)
    assert np.allclose(Xshared[:, :2, :], X0shared)
    X0multi = np.array([
        [[2.0, 3.0]], [[4.0, 5.0]], [[6.0, 7.0]],
    ])
    Xmulti = model.sim(5, nrep=3, X0=X0multi, rng=rng)
    assert Xmulti.shape == (3, 5, 2)
    assert np.allclose(Xmulti[:, 0, :], X0multi[:, 0, :])


def test_model_accepts_single_matrices_and_returns_copies():
    A = np.array([[0.1, 0.2], [0.3, 0.4]])
    B = np.array([[0.5, 0.6], [0.7, 0.8]])
    Sig = np.array([[1.0, 0.25], [0.25, 2.0]])
    model = varmapack.Model(A=A, B=B, Sig=Sig)
    assert model.p == 1
    assert model.q == 1
    assert np.allclose(model.A[0], A)
    assert np.allclose(model.B[0], B)
    assert np.allclose(model.Sig, Sig)
    A[0, 0] = 99
    B[0, 0] = 99
    Sig[0, 0] = 99
    assert np.isclose(model.A[0, 0, 0], 0.1)
    assert np.isclose(model.B[0, 0, 0], 0.5)
    assert np.isclose(model.Sig[0, 0], 1)
    returned_A = model.A
    returned_A[0, 0, 0] = 88
    assert np.isclose(model.A[0, 0, 0], 0.1)


def test_nonstationary_varma_with_starting_values():
    model = varmapack.Model(
        A=np.array([[[1.25]]]), B=np.array([[[0.5]]]), Sig=np.array([[1.0]]))
    X0 = np.array([[2.0], [3.0], [4.0]])
    X, E = model.sim(6, X0=X0, rng=_rng(46), return_shocks=True)
    assert np.allclose(X[0, :3, 0], X0[:, 0])
    for t in range(1, 6):
        expected = 1.25*X[0, t - 1, 0] + E[0, t, 0] + 0.5*E[0, t - 1, 0]
        assert np.allclose(X[0, t, 0], expected)
