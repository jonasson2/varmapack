# Varmapack

Varmapack is a Python interface to the Varmapack C library for simulation and
analysis of Gaussian VAR, VMA, VARMA, and VARMAX time-series models.

It supports burn-in-free simulation, model testcases, theoretical and sample
autocovariances, spectral radii, and impulse response functions.

## Quick Example

```python
import varmapack
model = varmapack.Model(
    A=[[0.6, 0.1], [0.0, 0.4]],
    Sig=[[2, 0], [0, 1]],
)
X = model.sim(200)
print(X.shape)
```

## Documentation

The full documentation is at:

> https://varmapack.readthedocs.io/

For development builds from source, see:

> https://github.com/jonasson2/varmapack/blob/main/DEVELOPMENT.md
