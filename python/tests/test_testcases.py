import numpy as np
import randompack
import varmapack


def test_named_testcase_and_catalog():
    model = varmapack.testcase("tinyAR")
    assert model.p == 1
    assert model.q == 0
    assert model.r == 1
    assert model.A.shape == (1, 1, 1)
    assert model.B is None
    assert model.Sig.shape == (1, 1)
    assert np.isclose(model.A[0, 0, 0], 0.5)
    assert np.isclose(model.specrad(), 0.5)
    cases = varmapack.testcases()
    assert isinstance(cases, dict)
    assert len(cases) == 16
    assert cases["tinyAR"] == {"index": 1, "p": 1, "q": 0, "r": 1}
    assert cases["largeAR"] == {"index": 15, "p": 5, "q": 0, "r": 7}
    assert cases["largeARMA"] == {"index": 16, "p": 3, "q": 3, "r": 7}
    assert "tinyAR" in repr(cases)


def test_indexed_testcase():
    model = varmapack.testcase(8)
    assert model.p == 1
    assert model.q == 1
    assert model.r == 2
    assert model.A.shape == (1, 2, 2)
    assert model.B.shape == (1, 2, 2)
    rng = randompack.Rng()
    rng.seed(123)
    assert model.sim(5, nrep=2, rng=rng).shape == (2, 5, 2)


def test_spectral_radius_testcase():
    model = varmapack.testcase("rho", p=3, q=1, r=2, rho=0.5)
    assert model.p == 3
    assert model.q == 1
    assert model.r == 2
    assert model.A.shape == (3, 2, 2)
    assert model.B.shape == (1, 2, 2)
    Gamma0 = model.acvf(0)
    Gamma4 = model.acvf(4)
    assert Gamma0.shape == (1, 2, 2)
    assert Gamma4.shape == (5, 2, 2)
    assert np.allclose(Gamma0[0], Gamma4[0])


def test_random_testcase_reproducibility():
    rng = randompack.Rng()
    rng.seed(123)
    first = varmapack.testcase("random", p=1, q=1, r=2, rng=rng)
    rng.seed(123)
    second = varmapack.testcase("random", p=1, q=1, r=2, rng=rng)
    assert np.allclose(first.A, second.A)
    assert np.allclose(first.B, second.B)
    assert np.allclose(first.Sig, second.Sig)
