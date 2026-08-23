#!/usr/bin/env python3
# TimeSimulate.py: compare Varmapack and statsmodels simulation (ns/value).

import argparse
import time
import warnings
from collections.abc import Callable

import numpy as np
import randompack
import varmapack
from statsmodels.tools.sm_exceptions import EstimationWarning
from statsmodels.tsa.statespace.varmax import VARMAX

ALL_MODELS = (
    "tinyAR", "tinyMA", "tinyARMA", "smallAR1", "smallAR2", "smallMA1",
    "smallMA2", "smallARMA1", "smallARMA2", "mediumAR", "mediumMA1",
    "mediumARMA1", "mediumARMA2", "mediumMA2", "largeAR", "largeARMA"
)
PAPER_MODELS = {
    "tinyAR", "tinyARMA", "smallAR1", "smallARMA2", "mediumAR",
    "mediumARMA2", "largeAR", "largeARMA"
}
PLATFORM_MODELS = {
    "tinyARMA", "smallAR2", "mediumAR", "mediumARMA1", "largeAR"
}


def varmax_params(model: varmapack.Model) -> np.ndarray:
  parts = []
  if model.p > 0:
    ar = np.concatenate([model.A[lag] for lag in range(model.p)], axis=1)
    parts.append(ar.reshape(-1))
  if model.q > 0:
    ma = np.concatenate([model.B[lag] for lag in range(model.q)], axis=1)
    parts.append(ma.reshape(-1))
  L = np.linalg.cholesky(model.Sig)
  parts.append(L[np.tril_indices(model.r)])
  return np.concatenate(parts)


def make_varmax(model: varmapack.Model) -> VARMAX:
  dummy = np.zeros((1, model.r))
  return VARMAX(dummy, order=(model.p, model.q), trend="n",
                error_cov_type="unstructured", enforce_stationarity=False,
                enforce_invertibility=False)


def time_simulation(fn: Callable[[], np.ndarray], values: int,
                    bench_time: float) -> float:
  fn()
  calls = 0
  t0 = time.perf_counter()
  elapsed = 0
  while elapsed < bench_time:
    fn()
    calls += 1
    elapsed = time.perf_counter() - t0
  return 1e9*elapsed/(calls*values)


def time_model(name: str, model: varmapack.Model, length: int, nrep: int,
               bench_time: float) -> None:
  varmax = make_varmax(model)
  params = varmax_params(model)
  rho = 0 if model.p == 0 else model.specrad()
  values = length*nrep*model.r
  rng = randompack.Rng()
  rng.seed(123)
  statsmodels_rng = np.random.default_rng(123)
  ns_varmapack = time_simulation(
      lambda: model.sim(length, nrep=nrep, rng=rng), values, bench_time)
  ns_varmax = time_simulation(
      lambda: varmax.simulate(params, nsimulations=length, repetitions=nrep,
                              transformed=True, random_state=statsmodels_rng),
      values, bench_time)
  print(f"{name:<12} {model.p:>2} {model.q:>2} {model.r:>2} {rho:>5.3f} "
        f"{ns_varmapack:>10.1f} {ns_varmax:>7.0f} "
        f"{ns_varmax/ns_varmapack:>6.0f}")


def main() -> None:
  parser = argparse.ArgumentParser(
      description="Compare Varmapack and statsmodels VARMAX simulation")
  parser.add_argument("-M", type=int, default=1000, dest="nrep",
                      help="replicates per call (default 1000)")
  parser.add_argument("-n", type=int, default=100, dest="length",
                      help="length of each series (default 100)")
  parser.add_argument("-t", type=float, default=0.1, dest="bench_time",
                      help="benchmark time per implementation and testcase (seconds)")
  selection = parser.add_mutually_exclusive_group()
  selection.add_argument("-8", action="store_true", dest="selected",
                         help="time the eight selected paper models only")
  selection.add_argument("-5", action="store_true", dest="platform",
                         help="time the five selected platform models only")
  args = parser.parse_args()
  if args.nrep < 1 or args.length < 1 or args.bench_time <= 0:
    parser.error("-M, -n, and -t must be positive")
  warnings.filterwarnings("ignore", category=EstimationWarning)
  print("Varmapack Python benchmark and comparison with statsmodels VARMAX simulate")
  print(f"{'benchmark unit:':<25}ns/value")
  print(f"{'length per series:':<25}{args.length}")
  print(f"{'replicates per call:':<25}{args.nrep}")
  print(f"{'benchmark time per case:':<25}{args.bench_time:.1f} s")
  if args.selected:
    print(f"{'models:':<25}selected paper set")
  if args.platform:
    print(f"{'models:':<25}selected platform set")
  print()
  print(f"{'Testcase':<12} {'p':>2} {'q':>2} {'r':>2} {'rho':>5} {'Varmapack':>10} "
        f"{'VARMAX':>7} {'Ratio':>6}")
  names = ALL_MODELS
  if args.selected:
    names = [name for name in ALL_MODELS if name in PAPER_MODELS]
  if args.platform:
    names = [name for name in ALL_MODELS if name in PLATFORM_MODELS]
  for name in names:
    time_model(name, varmapack.testcase(name), args.length, args.nrep,
               args.bench_time)


if __name__ == "__main__":
  main()
