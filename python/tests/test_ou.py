import numpy as np
import pytest
import randompack
import varmapack as vp


def rng(seed=423):
    r = randompack.Rng()
    r.seed(seed)
    return r


def test_ou_stationary_and_lagged_covariance():
    tau = np.array([1., 2., 3.])
    cov = np.array([[1., .2, -.1], [.2, 2., .3], [-.1, .3, 1.]])
    mu = np.array([2., -1., .5])
    model = vp.OUModel(tau, cov, mu)
    t = [-3, -2.9, -1, 3]
    X = model.sim(t, nrep=60000, rng=rng())
    assert X.shape == (60000, 4, 3)
    np.testing.assert_allclose(X.mean(axis=0), np.tile(mu, (4, 1)), atol=.025)
    for k in range(4):
        for h in range(k + 1):
            empirical = (X[:, k] - mu).T @ (X[:, h] - mu) / len(X)
            expected = np.exp(-(t[k] - t[h]) / tau)[:, None] * cov
            np.testing.assert_allclose(empirical, expected, atol=.04)
    assert model.r == 3
    np.testing.assert_equal(model.tau, tau)
    np.testing.assert_equal(model.cov, cov)
    np.testing.assert_equal(model.mu, mu)


def test_ou_copy_and_reproducibility():
    tau = np.array([1., 2.])
    cov = np.array([[1., np.nan], [.2, 1.]])
    mu = np.array([1., 2.])
    model = vp.OUModel(tau, cov, mu)
    tau[:] = 99
    cov[:] = 99
    mu[:] = 99
    model.tau[:] = 9
    np.testing.assert_equal(model.tau, [1, 2])
    np.testing.assert_equal(model.cov, [[1, .2], [.2, 1]])
    np.testing.assert_equal(model.mu, [1, 2])
    t = np.array([0., .1, 1., 100.])
    a = model.sim(t, nrep=4, rng=rng())
    b = model.sim(t + 1, nrep=4, rng=rng())
    np.testing.assert_allclose(a, b, atol=1e-14)
    np.testing.assert_equal(a, model.sim(t, nrep=4, rng=rng()))
    assert model.sim(t[:1], rng=rng()).shape == (1, 1, 2)
    assert model.sim(t, rng=rng()).dtype == np.float64
    assert model.sim(t).shape == (1, 4, 2)


def test_ou_semidefinite_and_extreme_gaps():
    t = [0, 1e-15, 1e6]
    X = vp.OUModel([2, 2], np.ones((2, 2))).sim(t, nrep=10, rng=rng())
    np.testing.assert_allclose(X[:, :, 0], X[:, :, 1], atol=1e-12)
    X = vp.OUModel([1, 2], np.zeros((2, 2)), [3, 4]).sim(t, rng=rng())
    np.testing.assert_equal(X, [[[3, 4], [3, 4], [3, 4]]])
    # Equal-step VAR(1) innovation covariance C - A C A.T.
    cov = np.array([[2., .3], [.3, 1.]])
    a = np.diag(np.exp(-.4 / np.array([1., 2.])))
    sig = vp.cov2sig(cov, A=a)
    np.testing.assert_allclose(sig, cov * -np.expm1(-.4 *
        (1 / np.array([1., 2.])[:, None] + 1 / np.array([1., 2.])[None, :])))


def test_ou_errors():
    for tau, cov, mu in [([], [], None), ([1], [[1, 0], [0, 1]], None),
                         ([1], [[1]], [1, 2]), ([0], [[1]], None),
                         ([np.inf], [[1]], None), ([1], [[-1]], None),
                         ([1], [[1]], [np.nan]), ([1, 2], [[1, 1], [1, 1]], None),
                         ([1, 1], [[1, 2], [2, 1]], None)]:
        with pytest.raises((ValueError, vp.VarmapackError)):
            vp.OUModel(tau, cov, mu)
    model = vp.OUModel([1], [[1]])
    for t in [[], [[1, 2]], [0, 0], [1, 0], [0, np.nan], [0, np.inf]]:
        with pytest.raises((ValueError, vp.VarmapackError)):
            model.sim(t)
    with pytest.raises(ValueError):
        model.sim([0], nrep=0)
    with pytest.raises(TypeError):
        model.sim([0], rng=object())
    with pytest.raises(TypeError):
        model.sim([0], nrep=1.5)


def test_ou_matrix_drift():
    K = np.array([[1., -2.], [2., 1.]])
    model = vp.OUModel(cov=np.eye(2), drift=K)
    assert model.tau is None and model.r == 2
    K[:] = 99
    model.drift[:] = 99
    np.testing.assert_equal(model.drift, [[1, -2], [2, 1]])
    h = .5
    X = model.sim([0, h, h + .1], nrep=30000, rng=rng())
    np.testing.assert_allclose(X[:, 1].T @ X[:, 1] / len(X), np.eye(2), atol=.03)
    c, s = np.cos(2*h), np.sin(2*h)
    expected = np.exp(-h) * np.array([[c, s], [-s, c]])
    np.testing.assert_allclose(X[:, 1].T @ X[:, 0] / len(X), expected, atol=.03)
    with pytest.raises(ValueError, match="exactly one"):
        vp.OUModel(tau=[1, 1], cov=np.eye(2), drift=np.eye(2))
    with pytest.raises(ValueError, match="shape"):
        vp.OUModel(cov=np.eye(2), drift=[[1]])
    with pytest.raises(ValueError, match="real"):
        vp.OUModel(cov=np.eye(2), drift=1j * np.eye(2))
    with pytest.raises(vp.VarmapackError, match="incompatible"):
        vp.OUModel(cov=np.eye(2), drift=[[1, 3], [0, 1]])


def test_gamma_times():
    for shape in [.5, 1, 4, 16]:
        t = vp.gamma_times(60001, shape=shape, mean_gap=2, rng=rng())
        d = np.diff(t)
        assert t[0] == 0 and np.all(d > 0)
        assert abs(d.mean() - 2) < 7 * np.sqrt(4 / shape / len(d))
        assert abs(d.var() - 4 / shape) < .08 * 4 / shape
        np.testing.assert_equal(t, vp.gamma_times(len(t), shape, 2, rng()))
    np.testing.assert_equal(vp.gamma_times(5, np.inf, 2), [0, 2, 4, 6, 8])
    np.testing.assert_equal(vp.gamma_times(1), [0])
    with pytest.raises(TypeError):
        vp.gamma_times(1.5)
    r = rng()
    vp.gamma_times(4, np.inf, rng=r)
    np.testing.assert_equal(vp.gamma_times(4, rng=r), vp.gamma_times(4, rng=rng()))
    for n, shape, gap in [(0, 4, 1), (2, 0, 1), (2, np.nan, 1),
                          (2, -np.inf, 1), (2, 4, 0), (2, 4, np.inf),
                          (4, np.inf, 1e308)]:
        with pytest.raises((ValueError, vp.VarmapackError)):
            vp.gamma_times(n, shape, gap)
def test_ou_specrad():
    diagonal = vp.OUModel([1, 2], np.eye(2))
    rotation = vp.OUModel(drift=[[.5, -3], [3, .5]], cov=np.eye(2))
    assert diagonal.specrad(.25) == pytest.approx(np.exp(-.125))
    assert rotation.specrad(2.) == pytest.approx(np.exp(-1))
    for gap in (0, -1, np.inf, np.nan, [1], 1j):
        with pytest.raises(ValueError, match="gap"):
            rotation.specrad(gap)
