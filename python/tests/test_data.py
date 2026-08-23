import numpy as np
import varmapack


def test_autocovariance_orientation_and_ml_normalization():
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 8.0]])
    Gamma = varmapack.autocov(X, 1)
    centered = X - X.mean(axis=0)
    assert Gamma.shape == (2, 2, 2)
    assert np.allclose(Gamma[0], centered.T @ centered/len(X))
    assert np.allclose(Gamma[1], centered[1:].T @ centered[:-1]/len(X))


def test_autocovariance_corrected_normalization():
    X = np.array([[1.0], [3.0], [5.0], [7.0]])
    centered = X - X.mean(axis=0)
    ml = varmapack.autocov(X, 2, "ML")
    corrected = varmapack.autocov(X, 2, "C")
    for lag in range(3):
        product = centered[lag:].T @ centered[:len(X) - lag]
        assert np.allclose(ml[lag], product/len(X))
        assert np.allclose(corrected[lag], product/(len(X) - lag))


def test_covariance_to_correlation():
    cov = np.array([
        [[4.0, 3.0], [3.0, 9.0]],
        [[2.0, -3.0], [6.0, 1.5]],
    ])
    original = cov.copy()
    expected = np.array([
        [[1.0, 0.5], [0.5, 1.0]],
        [[0.5, -0.5], [1.0, 1/6]],
    ])
    corr = varmapack.cov2corr(cov)
    assert corr.shape == cov.shape
    assert np.allclose(corr, expected)
    assert np.array_equal(cov, original)
    assert np.allclose(varmapack.cov2corr(cov[0]), expected[0])


def test_covariance_to_correlation_output_array():
    cov = np.array([
        [[4.0, 3.0], [3.0, 9.0]],
        [[2.0, -3.0], [6.0, 1.5]],
    ])
    expected = varmapack.cov2corr(cov)
    out = np.empty_like(cov)
    assert varmapack.cov2corr(cov, out=out) is out
    assert np.allclose(out, expected)
    inplace = cov.copy()
    assert varmapack.cov2corr(inplace, out=inplace) is inplace
    assert np.allclose(inplace, expected)
