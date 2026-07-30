import subprocess
import sys

import numpy as np
import randompack
import varmapack


model = varmapack.Model(A=np.array([[[0.4]]]), Sig=np.ones((1, 1)))
rng = randompack.Rng()
rng.seed(314)
X1 = model.sim(8, nrep=3, rng=rng)
tail1 = rng.unif(5)
rng.seed(314)
X2 = model.sim(8, nrep=3, rng=rng)
tail2 = rng.unif(5)
assert np.array_equal(X1, X2)
assert np.array_equal(tail1, tail2)
assert np.isfinite(rng.normal(3)).all()
assert np.isfinite(model.sim(4)).all()

random_case = varmapack.testcase("random", p=1, q=1, r=2)
assert random_case.A.shape == (1, 2, 2)
assert random_case.B.shape == (1, 2, 2)
case15 = varmapack.testcase(15)
assert case15.p == 4
assert case15.r == 5

missing_api = """
import randompack._core as core
del core._C_API
import varmapack
"""
result = subprocess.run(
    [sys.executable, "-c", missing_api],
    capture_output=True,
    text=True,
)
assert result.returncode != 0
assert "Randompack Python C API version" in result.stderr

wrong_version = """
import ctypes
import randompack._core as core

class ApiPrefix(ctypes.Structure):
    _fields_ = [
        ("abi_version", ctypes.c_uint32),
        ("struct_size", ctypes.c_size_t),
    ]

api = ApiPrefix(99, ctypes.sizeof(ApiPrefix))
name = b"randompack._core._C_API"
new_capsule = ctypes.pythonapi.PyCapsule_New
new_capsule.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_void_p]
new_capsule.restype = ctypes.py_object
core._C_API = new_capsule(ctypes.addressof(api), name, None)
import varmapack
"""
result = subprocess.run(
    [sys.executable, "-c", wrong_version],
    capture_output=True,
    text=True,
)
assert result.returncode != 0
assert "incompatible" in result.stderr
