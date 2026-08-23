import numpy as np
import randompack
import varmapack


def test_white_noise_moments():
    Sig = np.array([[2.0]])
    model = varmapack.Model(Sig=Sig)
    Gamma = model.acvf(3)
    assert np.allclose(Gamma[0], Sig)
    assert np.allclose(Gamma[1:], 0)
    assert np.isclose(model.specrad(), 0)
    assert np.isclose(model.ma_specrad(), 0)


def test_scalar_ar_moments():
    model = varmapack.Model(A=np.array([[[0.5]]]), Sig=np.array([[1.0]]))
    Gamma = model.acvf(4)
    expected = np.array([
        1/(1 - 0.25), 0.5/(1 - 0.25), 0.25/(1 - 0.25),
        0.125/(1 - 0.25), 0.0625/(1 - 0.25),
    ])
    assert np.allclose(Gamma[:, 0, 0], expected)
    assert np.allclose(model.psi(4)[:, 0, 0], [1.0, 0.5, 0.25, 0.125, 0.0625])
    assert np.allclose(model.irf(2)[:, 0, 0], [1.0, 0.5, 0.25])


def test_scalar_ma_moments():
    model = varmapack.Model(B=np.array([[[0.25]]]), Sig=np.array([[4.0]]))
    Gamma = model.acvf(3)
    assert np.allclose(Gamma[:, 0, 0], [4.25, 1.0, 0.0, 0.0])
    assert np.allclose(model.psi(3)[:, 0, 0], [1.0, 0.25, 0.0, 0.0])
    assert np.isclose(model.ma_specrad(), 0.25)


def test_multivariate_impulse_responses():
    A = np.array([[[0.1, 0.2], [0.3, 0.4]]])
    B = np.array([[[0.5, 0.6], [0.7, 0.8]]])
    Sig = np.array([[1.0, 0.25], [0.25, 2.0]])
    model = varmapack.Model(A=A, B=B, Sig=Sig)
    assert model.specrad() > 0
    assert model.ma_specrad() > 0
    Psi = model.psi(2)
    expected = np.array([
        [[1.0, 0.0], [0.0, 1.0]],
        [[0.6, 0.8], [1.0, 1.2]],
        [[0.26, 0.32], [0.58, 0.72]],
    ])
    assert Psi.shape == (3, 2, 2)
    assert np.allclose(Psi, expected)
    Theta = model.irf(2)
    L = np.linalg.cholesky(Sig)
    assert Theta.shape == (3, 2, 2)
    assert np.allclose(Theta, np.array([value @ L for value in expected]))


def test_irf_accepts_singular_covariance():
    A = np.array([[[0.1, 0.2], [0.3, 0.4]]])
    B = np.array([[[0.5, 0.6], [0.7, 0.8]]])
    for Sig in (
        np.array([[1.0, 2.0], [2.0, 4.0]]),
        np.array([[1.0, 1.0], [1.0, 1.0]]),
    ):
        model = varmapack.Model(A=A, B=B, Sig=Sig)
        Theta = model.irf(2)
        assert np.allclose(Theta[0] @ Theta[0].T, Sig)


def test_singular_white_noise_simulation():
    rng = randompack.Rng()
    rng.seed(321)
    model = varmapack.Model(Sig=np.array([[1.0, 1.0], [1.0, 1.0]]))
    X, E = model.sim(6, nrep=3, rng=rng, return_shocks=True)
    assert np.allclose(X, E)
    assert np.allclose(X[:, :, 0], X[:, :, 1])
    assert np.isfinite(X).all()


def test_larger_model_smoke_test():
    r = 20
    p = 5
    A = np.zeros((p, r, r))
    A[0] = 0.1*np.eye(r)
    model = varmapack.Model(A=A, Sig=np.eye(r))
    rng = randompack.Rng()
    rng.seed(321)
    Gamma = model.acvf(1)
    X = model.sim(6, nrep=2, rng=rng)
    assert Gamma.shape == (2, r, r)
    assert X.shape == (2, 6, r)
    assert np.isfinite(Gamma).all()
    assert np.isfinite(X).all()
