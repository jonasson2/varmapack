import numpy as np
import randompack
import varmapack


def assert_raises(exc, func):
    try:
        func()
    except exc:
        return
    raise AssertionError(f"expected {exc.__name__}")


Sig = np.eye(2)
A = np.zeros((1, 2, 2))
B = np.zeros((1, 2, 2))
C = np.zeros((1, 2))
model = varmapack.Model(A=A, B=B, Sig=Sig)
rng = randompack.Rng()
rng.seed(123)

assert_raises(ValueError, lambda: varmapack.Model(Sig=np.ones((2, 3))))
assert_raises(ValueError, lambda: varmapack.Model(A=np.zeros((2, 2)), Sig=Sig))
assert_raises(ValueError, lambda: varmapack.Model(A=np.zeros((1, 3, 3)), Sig=Sig))
assert_raises(ValueError, lambda: varmapack.Model(Sig=Sig, mu=np.ones(3)))
assert_raises(ValueError, lambda: varmapack.Model(Sig=Sig, mu=np.ones((2, 3))))
assert_raises(ValueError, lambda: varmapack.Model(C=C, Sig=Sig, mu=np.zeros(2)))

assert_raises(ValueError, lambda: model.sim(4, X0=np.zeros((1, 3)), rng=rng))
assert_raises(ValueError, lambda: model.sim(4, nrep=2, X0=np.zeros((3, 1, 2)), rng=rng))
assert_raises(ValueError, lambda: model.sim(4, z=np.zeros(4), rng=rng))
assert_raises(ValueError, lambda: model.acvf(-1))
assert_raises(ValueError, lambda: model.psi(-1))
assert_raises(ValueError, lambda: model.irf(-1))
assert_raises(ValueError, lambda: varmapack.autocov(np.zeros((2, 2, 1)), 1))

varmax = varmapack.Model(A=A, B=B, C=C, Sig=Sig)
assert_raises(ValueError, lambda: varmax.sim(4, X0=np.zeros((2, 2)), rng=rng))
assert_raises(ValueError, lambda: varmax.sim(4, z=np.zeros(4), rng=rng))
assert_raises(ValueError, lambda: varmax.sim(4, X0=np.zeros((2, 2)), z=np.zeros(3),
                                             rng=rng))
assert_raises(ValueError, lambda: varmax.sim(4, nrep=2, X0=np.zeros((2, 2)),
                                             z=np.zeros((3, 4)), rng=rng))

assert_raises(varmapack.VarmapackError, lambda: varmapack.testcase("does-not-exist"))
