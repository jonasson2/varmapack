import numpy as np
import randompack
import varmapack


Sig = np.array([[2.0]])
white = varmapack.Model(Sig=Sig)
G = white.acvf(3)
assert np.allclose(G[0], Sig)
assert np.allclose(G[1:], 0)
assert np.isclose(white.specrad(), 0)
assert np.isclose(white.ma_specrad(), 0)

A = np.array([[[0.5]]])
ar = varmapack.Model(A=A, Sig=np.array([[1.0]]))
G = ar.acvf(4)
expected = np.array([1/(1 - 0.25), 0.5/(1 - 0.25), 0.25/(1 - 0.25),
                     0.125/(1 - 0.25), 0.0625/(1 - 0.25)])
assert np.allclose(G[:, 0, 0], expected)
assert np.allclose(ar.psi(4)[:, 0, 0], [1.0, 0.5, 0.25, 0.125, 0.0625])
assert np.allclose(ar.irf(2)[:, 0, 0], [1.0, 0.5, 0.25])

B = np.array([[[0.25]]])
ma = varmapack.Model(B=B, Sig=np.array([[4.0]]))
G = ma.acvf(3)
assert np.allclose(G[:, 0, 0], [4.25, 1.0, 0.0, 0.0])
assert np.allclose(ma.psi(3)[:, 0, 0], [1.0, 0.25, 0.0, 0.0])
assert np.isclose(ma.ma_specrad(), 0.25)

X = np.array([[1.0], [3.0], [5.0], [7.0]])
Cml = varmapack.autocov(X, 2, "ML")
Cc = varmapack.autocov(X, 2, "C")
Xc = X - X.mean(axis=0)
assert np.allclose(Cml[0], Xc.T @ Xc / 4)
assert np.allclose(Cml[1], Xc[:-1].T @ Xc[1:] / 4)
assert np.allclose(Cc[1], Xc[:-1].T @ Xc[1:] / 3)

rng = randompack.Rng()
rng.seed(321)
model = varmapack.Model(Sig=np.array([[1.0, 1.0], [1.0, 1.0]]))
X, E = model.sim(6, nrep=3, rng=rng, return_shocks=True)
assert np.allclose(X, E)
assert np.allclose(X[:, :, 0], X[:, :, 1])
assert np.all(np.isfinite(X))
