import argparse
import time
import warnings

import numpy as np
import randompack
import varmapack
from statsmodels.tools.sm_exceptions import EstimationWarning
from statsmodels.tsa.statespace.varmax import VARMAX

def varmax_params(model):
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


def make_varmax(model):
  dummy = np.zeros((max(model.p + model.q + 10, 20), model.r))
  return VARMAX(dummy, order=(model.p, model.q), trend="n",
                error_cov_type="unstructured", enforce_stationarity=False,
                enforce_invertibility=False)


def compute_reps(values):
  return max(1, 1_000_000//values)


def time_simulation(fn, values, bench_time):
  calls = 0
  reps = compute_reps(values)
  t0 = time.perf_counter()
  elapsed = 0
  while elapsed < bench_time:
    for _ in range(reps):
      X = fn()
      sink = X.flat[-1]
    calls += reps
    elapsed = time.perf_counter() - t0
  return 1e9*elapsed/(calls*values)


def time_model(name, model, length, nrep, bench_time):
  varmax = make_varmax(model)
  params = varmax_params(model)
  rho = 0 if model.p == 0 else model.specrad()
  values = length*nrep*model.r
  rng = randompack.Rng()
  rng.seed(123)
  ns_varmapack = time_simulation(
      lambda: model.sim(length, nrep=nrep, rng=rng), values, bench_time)
  ns_varmax = time_simulation(
      lambda: varmax.simulate(params, nsimulations=length, repetitions=nrep,
                              transformed=True), values, bench_time)
  print(f"{name:<12} {model.p:>2} {model.q:>2} {model.r:>2} {rho:>5.2f} "
        f"{ns_varmapack:>10.1f} {ns_varmax:>7.0f} "
        f"{ns_varmax/ns_varmapack:>6.0f}")


def main():
  parser = argparse.ArgumentParser(
      description="Compare Varmapack and statsmodels VARMAX simulation")
  parser.add_argument("-M", type=int, default=1000, dest="nrep",
                      help="replicates per call (default 1000)")
  parser.add_argument("-n", type=int, default=100, dest="length",
                      help="length of each series (default 100)")
  parser.add_argument("-t", type=float, default=0.1, dest="bench_time",
                      help="benchmark time per implementation and testcase (seconds)")
  args = parser.parse_args()
  if args.nrep < 1 or args.length < 1 or args.bench_time <= 0:
    parser.error("-M, -n, and -t must be positive")
  warnings.filterwarnings("ignore", category=EstimationWarning)
  print("Varmapack Python benchmark and comparison with statsmodel's VARMAX simulate")
  print(f"{'benchmark unit:':<25}ns/value")
  print(f"{'length per series:':<25}{args.length}")
  print(f"{'replicates per call:':<25}{args.nrep}")
  print(f"{'benchmark time per case:':<25}{args.bench_time:.1f} s")
  print()
  print(f"{'Testcase':<12} {'p':>2} {'q':>2} {'r':>2} {'rho':>5} {'Varmapack':>10} "
        f"{'VARMAX':>7} {'Ratio':>6}")
  for name in [
      "tinyAR", "tinyMA", "tinyARMA", "smallAR1", "smallAR2", "smallMA1",
      "smallMA2", "smallARMA1", "smallARMA2", "mediumAR", "mediumMA1",
      "mediumARMA1", "mediumARMA2", "mediumMA2", "largeAR"
  ]:
    time_model(name, varmapack.testcase(name), args.length, args.nrep,
               args.bench_time)
  model = varmapack.testcase("rho", p=3, q=3, r=10, rho=0.98)
  time_model("Unnamed", model, args.length, args.nrep, args.bench_time)


if __name__ == "__main__":
  main()
