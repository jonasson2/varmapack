import numpy as np
import randompack
import varmapack


rng = randompack.Rng()
A = np.array([[[0.4]]])
B = np.array([[[0.2]]])
C = np.array([[0.8]])
Sig = np.array([[1.0]])
model = varmapack.Model(A=A, B=B, C=C, Sig=Sig)
assert model.order == (1, 1)
assert model.exog_order == 1
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

C2 = np.array([[0.5, -0.25]])
Sig2 = np.eye(2)
model2 = varmapack.Model(A=np.zeros((1, 2, 2)), B=np.zeros((1, 2, 2)), C=C2,
                         Sig=Sig2)
rng.seed(7)
X5 = model2.sim(4, nrep=2, X0=np.zeros((2, 2)), z=np.arange(4.0), rng=rng)
assert X5.shape == (2, 4, 2)
assert np.all(np.isfinite(X5))
