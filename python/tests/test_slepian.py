import numpy as np
import pytest
import randompack
import varmapack as vp


def seeded_rng():
    rng = randompack.Rng()
    rng.seed(127)
    return rng


def test_slepian_moments_and_parameters():
    theta = np.array([[0., .5], [.25, 0.]])
    expected_theta = theta.copy()
    model = vp.SlepianModel(theta, [[1.25, np.nan], [0., 1.0625]], [2., -1.])
    theta[:] = 99
    model.theta[:] = 99
    np.testing.assert_equal(model.theta, expected_theta)
    np.testing.assert_equal(model.cov, np.diag([1.25, 1.0625]))
    assert model.r == 2
    t = [0, .5, 1, 2]
    X = model.sim(t, nrep=30000, rng=seeded_rng()) - model.mu
    assert X.shape == (30000, 4, 2)
    expected_lags = [model.cov, .5*(model.cov + expected_theta),
                     expected_theta, np.zeros((2, 2))]
    for j, expected in enumerate(expected_lags):
        lagcov = X[:, j].T @ X[:, 0] / len(X)
        np.testing.assert_allclose(lagcov, expected, atol=.03)
    np.testing.assert_allclose(X.mean(axis=0), 0, atol=.025)


def test_slepian_rng_and_validation():
    model = vp.SlepianModel([[.5]], [[1.25]])
    X = model.sim([.1, .5, 1], nrep=3, rng=seeded_rng())
    np.testing.assert_equal(X, model.sim([.1, .5, 1], nrep=3, rng=seeded_rng()))
    assert model.sim([0]).shape == (1, 1, 1)
    for t in [[], [0, 0], [1, 0], [0, np.inf]]:
        with pytest.raises((ValueError, vp.VarmapackError)):
            model.sim(t)
    with pytest.raises(ValueError, match="square"):
        vp.SlepianModel([[1, 2]], [[1]])
    with pytest.raises(vp.VarmapackError):
        vp.SlepianModel([[0, 1], [0, 0]], [[1, 0], [0, 2]])


def test_slepian_mixed_units():
    for scale in [1., 1e8]:
        with pytest.raises(vp.VarmapackError):
            vp.SlepianModel([[0., 0.], [1 / scale, 0.]], np.diag([scale**2, .25]))
        model = vp.SlepianModel(np.zeros((3, 3)), np.diag([scale**2, .25, 0.]))
        X = model.sim([0.], nrep=20000, rng=seeded_rng())[:, 0]
        np.testing.assert_allclose(X.var(axis=0), [scale**2, .25, 0.], rtol=.05)
