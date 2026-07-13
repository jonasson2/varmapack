import numpy as np
import randompack
import varmapack


Sig = np.array([[1.0, 0.25], [0.25, 2.0]])
white = varmapack.Model(Sig=Sig)
X, E = white.sim(4, return_shocks=True)
assert X.shape == (1, 4, 2)
assert E.shape == (1, 4, 2)
assert np.allclose(X, E)
mu = np.array([[10.0, 20.0], [11.0, 21.0], [12.0, 22.0]])
white_mu = varmapack.Model(Sig=Sig, mu=mu)
X, E = white_mu.sim(5, return_shocks=True)
expected_mu = np.vstack([mu, mu[-1], mu[-1]])
assert np.allclose(X[0] - E[0], expected_mu)

rng = randompack.Rng()
rng.seed(123)
A = np.empty((1, 1, 1))
A[0, 0, 0] = 0.4
Sig1 = np.array([[1.0]])
model = varmapack.Model(A=A, Sig=Sig1)
assert model.p == 1
assert model.q == 0
assert model.r == 1
assert np.isclose(model.specrad(), 0.4)
assert np.isclose(model.ma_specrad(), 0.0)
X1 = model.sim(5, rng=rng)
assert X1.shape == (1, 5, 1)
assert np.isfinite(X1).all()

rng.seed(123)
X2 = model.sim(5, rng=rng)
assert np.allclose(X1, X2)

X3 = model.sim(length=5, nrep=3, rng=rng)
assert X3.shape == (3, 5, 1)
X0multi = np.array([[[2.0]], [[4.0]], [[6.0]]])
X3start = model.sim(length=5, nrep=3, X0=X0multi, rng=rng)
assert X3start.shape == (3, 5, 1)
assert np.allclose(X3start[:, 0, :], X0multi[:, 0, :])

G = white.acvf(2)
assert G.shape == (3, 2, 2)
assert np.allclose(G[0], Sig)
assert np.allclose(G[1:], 0)

Xsmall = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 8.0]])
C = varmapack.autocov(Xsmall, 1)
Xc = Xsmall - Xsmall.mean(axis=0)
assert C.shape == (2, 2, 2)
assert np.allclose(C[0], Xc.T @ Xc / Xsmall.shape[0])
assert np.allclose(C[1], Xc[:-1].T @ Xc[1:] / Xsmall.shape[0])

A2 = np.array([[[0.1, 0.2], [0.3, 0.4]]])
B2 = np.array([[[0.5, 0.6], [0.7, 0.8]]])
model2 = varmapack.Model(A=A2, B=B2, Sig=Sig)
assert np.allclose(model2.A, A2)
assert np.allclose(model2.B, B2)
assert np.allclose(model2.Sig, Sig)
assert model2.specrad() > 0
assert model2.ma_specrad() > 0
model2_single = varmapack.Model(A=A2[0], B=B2[0], Sig=Sig)
assert model2_single.p == 1
assert model2_single.q == 1
assert np.allclose(model2_single.A, A2)
assert np.allclose(model2_single.B, B2)
X0shared = np.array([[1.0, 2.0], [3.0, 4.0]])
Xshared = model2.sim(4, nrep=2, X0=X0shared, rng=rng)
assert np.allclose(Xshared[:, :2, :], X0shared)
Psi2 = model2.psi(2)
Psi2_expected = np.array([
    [[1.0, 0.0], [0.0, 1.0]], [[0.6, 0.8], [1.0, 1.2]], [[0.26, 0.32], [0.58, 0.72]], ])
assert Psi2.shape == (3, 2, 2)
assert np.allclose(Psi2, Psi2_expected)
Theta2 = model2.irf(2)
L = np.linalg.cholesky(Sig)
Theta2_expected = np.array([Psi2_expected[j] @ L for j in range(3)])
assert Theta2.shape == (3, 2, 2)
assert np.allclose(Theta2, Theta2_expected)
model3 = varmapack.Model(A=A2, B=B2, Sig=np.array([[1.0, 2.0], [2.0, 4.0]]))
Theta3 = model3.irf(2)
assert np.allclose(Theta3[0] @ Theta3[0].T, model3.Sig)
model4 = varmapack.Model(A=A2, B=B2, Sig=np.array([[1.0, 1.0], [1.0, 1.0]]))
Theta4 = model4.irf(2)
assert np.allclose(Theta4[0] @ Theta4[0].T, model4.Sig)

C1 = np.array([[0.8]])
varmax = varmapack.Model(A=A, B=np.array([[[0.2]]]), C=C1, Sig=Sig1)
assert varmax.s == 1
assert np.allclose(varmax.C, C1)
z = np.array([[1.0, -1.0, 1.0, -1.0], [1.5, -0.5, 1.5, -0.5]])
X0x = np.array([[[0.25], [0.5]], [[-0.25], [-0.5]]])
Xx, Ex = varmax.sim(4, nrep=2, X0=X0x, z=z, rng=rng, return_shocks=True)
assert Xx.shape == (2, 4, 1)
assert Ex.shape == (2, 4, 1)
assert np.allclose(Xx[:, :2, :], X0x)

tc = varmapack.testcase("tinyAR")
assert tc.p == 1
assert tc.q == 0
assert tc.r == 1
assert tc.A.shape == (1, 1, 1)
assert tc.B is None
assert tc.Sig.shape == (1, 1)
assert np.isclose(tc.A[0, 0, 0], 0.5)
assert np.isclose(tc.specrad(), 0.5)

tc2 = varmapack.testcase(8)
assert tc2.p == 1
assert tc2.q == 1
assert tc2.r == 2
assert tc2.A.shape == (1, 2, 2)
assert tc2.B.shape == (1, 2, 2)
X4 = tc2.sim(5, nrep=2, rng=rng)
assert X4.shape == (2, 5, 2)

tc3 = varmapack.testcase("rho", p=3, q=1, r=2, rho=0.5)
assert tc3.p == 3
assert tc3.q == 1
assert tc3.r == 2
assert tc3.A.shape == (3, 2, 2)
assert tc3.B.shape == (1, 2, 2)
G0 = tc3.acvf(0)
G4 = tc3.acvf(4)
assert G0.shape == (1, 2, 2)
assert G4.shape == (5, 2, 2)
assert np.allclose(G0[0], G4[0])

rng.seed(123)
tc4 = varmapack.testcase("random", p=1, q=1, r=2, rng=rng)
rng.seed(123)
tc5 = varmapack.testcase("random", p=1, q=1, r=2, rng=rng)
assert np.allclose(tc4.A, tc5.A)
assert np.allclose(tc4.B, tc5.B)
assert np.allclose(tc4.Sig, tc5.Sig)
