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


def test_innovation_covariance_scalar_and_white_noise():
    assert np.allclose(varmapack.cov2sig([[4.0]], A=[[0.5]]), [[3.0]])
    assert np.allclose(varmapack.cov2sig([[2.125]], B=[[0.25]]), [[2.0]])
    assert np.allclose(varmapack.cov2sig([[3.5]], A=[[0.5]], B=[[0.25]]), [[2.0]])
    cov = np.array([[2, 0.3], [0.3, 1]])
    assert np.allclose(varmapack.cov2sig(cov), cov)


def test_innovation_covariance_orientation_and_storage():
    A = np.array([[0.4, 0.2], [-0.1, 0.3]])
    cov = np.array([[2.0, 0.3], [0.3, 1.0]])
    expected = cov - A @ cov @ A.T
    original = cov.copy()
    storage = np.empty((4, 4))
    storage[::2, ::2] = cov
    for target in (cov, np.asfortranarray(cov), storage[::2, ::2]):
        actual_target = np.tril(target) + np.tril(target, -1).T
        Sig = varmapack.cov2sig(target, A=A)
        assert np.allclose(Sig, actual_target - A @ actual_target @ A.T)
        assert Sig.shape == (2, 2)
        assert Sig.flags.c_contiguous
    assert np.array_equal(cov, original)
    lower = cov.copy()
    lower[0, 1] = np.nan
    assert np.allclose(varmapack.cov2sig(lower, A=A), expected)
    assert np.isnan(lower[0, 1])
    out = np.empty_like(cov)
    assert varmapack.cov2sig(cov, A=A, out=out) is out
    assert np.allclose(out, expected)
    assert varmapack.cov2sig(cov, A=A, out=cov) is cov
    assert np.allclose(cov, expected)


def test_innovation_covariance_named_round_trips():
    for name in varmapack.testcases():
        model = varmapack.testcase(name)
        cov = model.acvf(0)[0]
        A, B = model.A, model.B
        original = cov.copy()
        Sig = varmapack.cov2sig(cov, A=A, B=B)
        assert np.allclose(Sig, model.Sig, rtol=1e-9, atol=1e-9), name
        recovered = varmapack.Model(A=A, B=B, Sig=Sig).acvf(0)[0]
        assert np.allclose(recovered, cov, rtol=1e-9, atol=1e-9), name
        assert np.array_equal(cov, original), name
        if A is not None:
            assert np.array_equal(model.A, A), name
        if B is not None:
            assert np.array_equal(model.B, B), name


def test_innovation_covariance_semidefinite():
    for cov in (np.ones((2, 2)), np.zeros((2, 2))):
        assert np.allclose(varmapack.cov2sig(cov, A=0.5*np.eye(2)), 0.75*cov)
