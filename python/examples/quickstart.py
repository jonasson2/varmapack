import randompack
import numpy as np
import varmapack

# Basic VAR simulation
Sig = [[2, 0], [0, 1]]
VAR_model = varmapack.Model(A=[[0.6, 0.1], [0.0, 0.4]], Sig=Sig)
X = VAR_model.sim(200)
print("VAR simulation, first five values of its first replicate:")
print(X[0, :5])

# Named testcases
print("\nNamed testcases:")
print(varmapack.testcases())
model = varmapack.testcase("smallARMA1")
X = model.sim(200, nrep=10)
print("smallARMA1 simulation shape:", X.shape)

# VAR, VMA, VARMA, and VARMAX simulation
A1 = np.array([[0.6, 0.1], [0.0, 0.4]])
B1 = np.array([[0.4, 0.3], [0.1, 0.2]])
A2 = 0.2*np.ones((2, 2))
B2 = 0.5*np.eye(2)
z = np.column_stack((np.cos(np.arange(200)/10), np.sin(np.arange(200)/10)))
VARMA_model = varmapack.Model(A=[A1, A2], B=B1, Sig=Sig)
VMA_model = varmapack.Model(B=[B1, B2], Sig=Sig)
VARMAX_model = varmapack.Model(A=A1, B=B1,
                               C=[[[0.8, 0.2], [-0.3, 0.4]]], Sig=Sig)
rng = randompack.Rng()
rng.seed(123)
X1 = VAR_model.sim(200, nrep=100)
X2 = VMA_model.sim(200, nrep=100, rng=rng)
X3, E = VARMA_model.sim(200, nrep=100, rng=rng, return_shocks=True)
X4 = VARMAX_model.sim(200, nrep=100, X0=np.zeros((2, 2)), z=z, rng=rng)
print("\nSimulation shapes:")
print("VAR:   ", X1.shape)
print("VMA:   ", X2.shape)
print("VARMA: ", X3.shape, "with shocks", E.shape)
print("VARMAX:", X4.shape)

# Model analysis
maxlag = 3
rho = VARMA_model.specrad()
rhoMA = VARMA_model.ma_specrad()
Psi = VARMA_model.psi(maxlag)
Theta = VARMA_model.irf(maxlag)
Gamma = VARMA_model.acvf(maxlag)
C = varmapack.autocov(X3[0], maxlag)
print("\nVARMA model analysis:")
print(f"spectral radius: {rho:.6f}")
print(f"MA spectral radius: {rhoMA:.6f}")
print("Psi, Theta, Gamma, and sample autocovariance shapes:")
print(Psi.shape, Theta.shape, Gamma.shape, C.shape)
