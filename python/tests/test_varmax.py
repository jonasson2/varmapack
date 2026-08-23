import numpy as np
import randompack
import varmapack


def _rng(seed):
    rng = randompack.Rng()
    rng.seed(seed)
    return rng


def _assert_scalar_recurrence(X, E, z, start, A, B, C):
    for path in range(len(X)):
        path_z = z if z.ndim == 1 else z[path]
        for t in range(start, X.shape[1]):
            expected = E[path, t, 0]
            if A:
                expected += A*X[path, t - 1, 0]
            if B:
                expected += B*E[path, t - 1, 0]
            expected += C[0]*path_z[t]
            if len(C) > 1:
                expected += C[1]*path_z[t - 1]
            assert np.allclose(X[path, t, 0], expected)


def _assert_vector_recurrence(model, X, E, z, start):
    for path in range(len(X)):
        path_z = z if z.ndim == 2 else z[path]
        for t in range(start, X.shape[1]):
            expected = E[path, t] + model.A[0] @ X[path, t - 1]
            expected += model.B[0] @ E[path, t - 1]
            expected += model.C[0].T @ path_z[t]
            expected += model.C[1].T @ path_z[t - 1]
            assert np.allclose(X[path, t], expected)


def test_scalar_shared_start_and_exogenous_path():
    model = varmapack.Model(
        A=np.array([[[0.4]]]), B=np.array([[[0.2]]]), C=np.array([[0.8]]),
        Sig=np.array([[1.0]]))
    X0 = np.array([[0.25], [0.5]])
    z = np.array([1.0, -1.0, 1.0, -1.0, 1.0])
    X, E = model.sim(5, nrep=3, X0=X0, z=z, rng=_rng(42), return_shocks=True)
    assert model.p == 1
    assert model.q == 1
    assert model.s == 1
    assert model.d == 1
    assert np.allclose(model.C, np.array([[0.8]]))
    assert X.shape == (3, 5, 1)
    assert E.shape == (3, 5, 1)
    assert np.allclose(X[:, :2, :], X0)
    _assert_scalar_recurrence(X, E, z, 2, 0.4, 0.2, [0.8])


def test_scalar_replicate_starting_values():
    model = varmapack.Model(
        A=np.array([[[0.4]]]), B=np.array([[[0.2]]]), C=np.array([[0.8]]),
        Sig=np.array([[1.0]]))
    X0 = np.array([[[0.0], [0.1]], [[1.0], [1.1]], [[2.0], [2.1]]])
    z = np.array([1.0, -1.0, 1.0, -1.0, 1.0])
    X, E = model.sim(5, nrep=3, X0=X0, z=z, rng=_rng(42), return_shocks=True)
    assert np.allclose(X[:, :2, :], X0)
    assert not np.allclose(X[0], X[1])
    _assert_scalar_recurrence(X, E, z, 2, 0.4, 0.2, [0.8])


def test_scalar_replicate_exogenous_paths():
    model = varmapack.Model(
        A=np.array([[[0.4]]]), B=np.array([[[0.2]]]), C=np.array([[0.8]]),
        Sig=np.array([[1.0]]))
    X0 = np.array([[[0.0], [0.1]], [[1.0], [1.1]], [[2.0], [2.1]]])
    z = np.array([
        [1.0, -1.0, 1.0, -1.0, 1.0],
        [1.2, -0.8, 1.2, -0.8, 1.2],
        [1.4, -0.6, 1.4, -0.6, 1.4],
    ])
    X, E = model.sim(5, nrep=3, X0=X0, z=z, rng=_rng(42), return_shocks=True)
    assert np.allclose(X[:, :2, :], X0)
    assert not np.allclose(X[0], X[1])
    _assert_scalar_recurrence(X, E, z, 2, 0.4, 0.2, [0.8])


def test_minimum_starting_path():
    C = [0.4, -0.1]
    model = varmapack.Model(
        A=np.array([[[0.4]]]), B=np.array([[[0.2]]]),
        C=np.array([[[C[0]]], [[C[1]]]]), Sig=np.array([[1.0]]))
    X0 = np.array([[0.25]])
    z = np.array([1.0, -1.0, 1.0, -1.0, 1.0])
    X, E = model.sim(5, nrep=3, X0=X0, z=z, rng=_rng(43), return_shocks=True)
    assert np.allclose(X[:, 0, :], X0)
    _assert_scalar_recurrence(X, E, z, 1, 0.4, 0.2, C)


def test_no_starting_values_required_when_orders_are_zero():
    z = np.array([1.0, -1.0, 1.0, -1.0, 1.0])
    model = varmapack.Model(C=np.array([[0.5]]), Sig=np.array([[1.0]]))
    X, E = model.sim(5, nrep=3, z=z, rng=_rng(44), return_shocks=True)
    assert np.allclose(X[:, :, 0], E[:, :, 0] + 0.5*z)


def test_scalar_exogenous_series_with_vector_states():
    C = np.array([[0.5, -0.25]])
    model = varmapack.Model(
        A=np.zeros((1, 2, 2)), B=np.zeros((1, 2, 2)), C=C, Sig=np.eye(2))
    X = model.sim(
        4, nrep=2, X0=np.zeros((2, 2)), z=np.arange(4.0), rng=_rng(7))
    assert X.shape == (2, 4, 2)
    assert np.isfinite(X).all()


def test_vector_exogenous_series():
    C = np.array([
        [[0.5, -0.4], [-0.2, 0.2]],
        [[0.1, 0.2], [0.3, 0.1]],
    ])
    model = varmapack.Model(
        A=np.array([[[0.3, 0.1], [-0.1, 0.2]]]),
        B=np.array([[[0.2, 0.0], [0.1, -0.2]]]), C=C, Sig=np.eye(2))
    z = np.array([
        [0.6, -0.2], [-0.2, 0.8], [0.8, 0.0], [0.0, -1.0], [-1.0, 0.4],
        [0.4, 0.2],
    ])
    X0 = np.zeros((2, 2, 2))
    X, E = model.sim(6, nrep=2, X0=X0, z=z, rng=_rng(7), return_shocks=True)
    assert model.s == 2
    assert model.d == 2
    assert np.allclose(model.C, C)
    _assert_vector_recurrence(model, X, E, z, 2)


def test_replicate_vector_exogenous_series():
    C = np.array([
        [[0.5, -0.4], [-0.2, 0.2]],
        [[0.1, 0.2], [0.3, 0.1]],
    ])
    model = varmapack.Model(
        A=np.array([[[0.3, 0.1], [-0.1, 0.2]]]),
        B=np.array([[[0.2, 0.0], [0.1, -0.2]]]), C=C, Sig=np.eye(2))
    z = np.array([
        [0.6, -0.2], [-0.2, 0.8], [0.8, 0.0], [0.0, -1.0], [-1.0, 0.4],
        [0.4, 0.2],
    ])
    z = np.stack([z, z + 0.1])
    X, E = model.sim(
        6, nrep=2, X0=np.zeros((2, 2, 2)), z=z, rng=_rng(7),
        return_shocks=True)
    _assert_vector_recurrence(model, X, E, z, 2)
