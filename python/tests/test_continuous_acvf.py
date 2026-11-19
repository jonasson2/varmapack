import numpy as np
import pytest
import varmapack as vp


def test_ou_acvf():
    drift = [[1., -1.], [0., 1.]]
    cov = np.eye(2)
    model = vp.OUModel(drift=drift, cov=cov)
    lags = np.array([1., -.5, 0., 2., .5, -1.])
    result = model.acvf(lags)
    unit = np.exp(-1) * np.array([[1., 1.], [0., 1.]])
    discrete = vp.Model(A=unit, Sig=cov - unit @ cov @ unit.T).acvf(2)
    np.testing.assert_allclose(result[[2, 0, 3]], discrete, atol=1e-13)
    half = np.exp(-.5) * np.array([[1., .5], [0., 1.]])
    np.testing.assert_allclose(result[[1, 4, 5]], [half.T, half, unit.T], atol=1e-13)
    tau = np.array([1., 2.])
    cov = np.array([[1., .5], [.5, 2.]])
    model = vp.OUModel(tau=tau, cov=cov)
    np.testing.assert_allclose(model.acvf(.25), np.exp(-.25 / tau)[:, None] * cov)
    np.testing.assert_equal(model.acvf([]).shape, (0, 2, 2))


def test_slepian_acvf():
    theta = np.array([[0., .5], [.25, 0.]])
    cov = np.eye(2) + theta @ theta.T
    model = vp.SlepianModel(theta, cov)
    discrete = vp.Model(B=theta, Sig=np.eye(2)).acvf(2)
    np.testing.assert_allclose(model.acvf([0, 1, 2]), discrete, atol=1e-13)
    np.testing.assert_allclose(model.acvf(.25), .75 * cov + .25 * theta)
    np.testing.assert_allclose(model.acvf(-1.75), .25 * theta.T)
    np.testing.assert_allclose(model.acvf(1e-20), cov + 1e-20 * theta)
    np.testing.assert_equal(model.acvf([-3, 3]), np.zeros((2, 2, 2)))
    np.testing.assert_equal(model.acvf([]).shape, (0, 2, 2))


@pytest.mark.parametrize("lags", [np.inf, [np.nan], [1j], [[0]]])
def test_invalid_lags(lags):
    for model in [vp.OUModel([1], [[1]]), vp.SlepianModel([[0]], [[1]])]:
        with pytest.raises(ValueError, match="lags"):
            model.acvf(lags)
