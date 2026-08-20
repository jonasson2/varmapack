import numpy as np
import randompack
import varmapack


rng = randompack.Rng()
A = np.array([[[0.4]]])
B = np.array([[[0.2]]])
C = np.array([[0.8]])
Sig = np.array([[1.0]])
model = varmapack.Model(A=A, B=B, C=C, Sig=Sig)
assert model.p == 1
assert model.q == 1
assert model.s == 1
assert np.allclose(model.C, C)

length = 5
nrep = 3
X0shared = np.array([[0.25], [0.5]])
z_shared = np.array([1.0, -1.0, 1.0, -1.0, 1.0])
rng.seed(42)
X1, E1 = model.sim(length, nrep=nrep, X0=X0shared, z=z_shared, rng=rng,
                   return_shocks=True)
assert X1.shape == (nrep, length, 1)
assert E1.shape == (nrep, length, 1)
assert np.allclose(X1[:, :2, :], X0shared)

X0multi = np.array([[[0.0], [0.1]], [[1.0], [1.1]], [[2.0], [2.1]]])
rng.seed(42)
X2 = model.sim(length, nrep=nrep, X0=X0multi, z=z_shared, rng=rng)
assert X2.shape == (nrep, length, 1)
assert np.allclose(X2[:, :2, :], X0multi)
assert not np.allclose(X2[0], X2[1])

z_multi = np.array([
    [1.0, -1.0, 1.0, -1.0, 1.0],
    [1.2, -0.8, 1.2, -0.8, 1.2],
    [1.4, -0.6, 1.4, -0.6, 1.4],
])
rng.seed(42)
X3 = model.sim(length, nrep=nrep, X0=X0shared, z=z_multi, rng=rng)
assert X3.shape == (nrep, length, 1)
assert np.allclose(X3[:, :2, :], X0shared)
assert not np.allclose(X3[0], X3[1])

rng.seed(42)
X4, E4 = model.sim(length, nrep=nrep, X0=X0multi, z=z_multi, rng=rng,
                   return_shocks=True)
assert X4.shape == (nrep, length, 1)
assert E4.shape == (nrep, length, 1)
assert np.allclose(X4[:, :2, :], X0multi)

C_lagged = np.array([[[0.4]], [[-0.1]]])
minimum_model = varmapack.Model(A=A, B=B, C=C_lagged, Sig=Sig)
X0minimum = np.array([[0.25]])
rng.seed(43)
Xminimum, Eminimum = minimum_model.sim(
    length, nrep=nrep, X0=X0minimum, z=z_shared, rng=rng, return_shocks=True)
assert np.allclose(Xminimum[:, 0, :], X0minimum)
for j in range(nrep):
    for t in range(1, length):
        expected = (0.4*Xminimum[j, t - 1, 0] + Eminimum[j, t, 0] +
                    0.2*Eminimum[j, t - 1, 0] + 0.4*z_shared[t] -
                    0.1*z_shared[t - 1])
        assert np.allclose(Xminimum[j, t, 0], expected)

zero_model = varmapack.Model(C=np.array([[0.5]]), Sig=Sig)
rng.seed(44)
Xzero, Ezero = zero_model.sim(length, nrep=nrep, z=z_shared, rng=rng,
                              return_shocks=True)
for j in range(nrep):
    assert np.allclose(Xzero[j, :, 0], Ezero[j, :, 0] + 0.5*z_shared)

C2 = np.array([[0.5, -0.25]])
Sig2 = np.eye(2)
model2 = varmapack.Model(A=np.zeros((1, 2, 2)), B=np.zeros((1, 2, 2)), C=C2,
                         Sig=Sig2)
rng.seed(7)
X5 = model2.sim(4, nrep=2, X0=np.zeros((2, 2)), z=np.arange(4.0), rng=rng)
assert X5.shape == (2, 4, 2)
assert np.all(np.isfinite(X5))

C3 = np.array([
    [[0.5, -0.4], [-0.2, 0.2]],
    [[0.1, 0.2], [0.3, 0.1]],
])
model3 = varmapack.Model(A=np.array([[[0.3, 0.1], [-0.1, 0.2]]]),
                          B=np.array([[[0.2, 0.0], [0.1, -0.2]]]),
                          C=C3, Sig=Sig2)
z_vector = np.array([
    [0.6, -0.2], [-0.2, 0.8], [0.8, 0.0], [0.0, -1.0], [-1.0, 0.4],
    [0.4, 0.2],
])
z_vector_multi = np.stack([z_vector, z_vector + 0.1])
X0vector = np.zeros((2, 2, 2))
rng.seed(7)
Xcommon, Ecommon = model3.sim(6, nrep=2, X0=X0vector, z=z_vector, rng=rng,
                              return_shocks=True)
assert model3.s == 2
assert model3.d == 2
assert np.allclose(model3.C, C3)
for j in range(2):
    for t in range(2, 6):
        expected = (Ecommon[j, t] + model3.A[0] @ Xcommon[j, t - 1] +
                    model3.B[0] @ Ecommon[j, t - 1] + C3[0].T @ z_vector[t] +
                    C3[1].T @ z_vector[t - 1])
        assert np.allclose(Xcommon[j, t], expected)
rng.seed(7)
X6, E6 = model3.sim(6, nrep=2, X0=X0vector, z=z_vector_multi, rng=rng,
                     return_shocks=True)
for j in range(2):
    for t in range(2, 6):
        expected = (E6[j, t] + model3.A[0] @ X6[j, t - 1] +
                    model3.B[0] @ E6[j, t - 1] + C3[0].T @ z_vector_multi[j, t] +
                    C3[1].T @ z_vector_multi[j, t - 1])
        assert np.allclose(X6[j, t], expected)
