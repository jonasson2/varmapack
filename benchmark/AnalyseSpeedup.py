#!/usr/bin/env python3

"""Controlled comparison of the Python interfaces to Varmapack and VARMAX."""

import argparse
import gc
import os
import sys
import time
import warnings
from collections.abc import Callable

THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
)
for variable in THREAD_VARIABLES:
  os.environ[variable] = "1"

import numpy as np
import randompack
import varmapack
from statsmodels.tools.sm_exceptions import EstimationWarning
from statsmodels.tsa.statespace.varmax import VARMAX

MODELS = (
    "tinyAR", "tinyARMA", "smallAR1", "smallARMA2", "mediumAR",
    "mediumARMA2", "largeAR", "largeARMA",
)


def parse_options() -> argparse.Namespace:
  parser = argparse.ArgumentParser(
      description="Controlled Varmapack/Statsmodels VARMAX comparison"
  )
  parser.add_argument("-R", type=int, default=5, dest="repeats",
                      help="independent timing runs (default 5)")
  parser.add_argument("-t", type=float, default=0.1, dest="target",
                      help="target per individual timing (default 0.1 s)")
  parser.add_argument("-n", type=int, default=100, dest="length",
                      help="length of each series (default 100)")
  parser.add_argument("-M", type=int, default=1000, dest="nrep",
                      help="number of series in a batch (default 1000)")
  options = parser.parse_args()
  if (options.repeats < 1 or options.target <= 0 or options.length < 1 or
      options.nrep < 2):
    parser.error("-R, -t, and -n must be positive, and -M must be at least 2")
  return options


def varmax_params(model: varmapack.Model) -> np.ndarray:
  parts = []
  if model.p > 0:
    ar = np.concatenate([model.A[lag] for lag in range(model.p)], axis=1)
    parts.append(ar.reshape(-1))
  if model.q > 0:
    ma = np.concatenate([model.B[lag] for lag in range(model.q)], axis=1)
    parts.append(ma.reshape(-1))
  factor = np.linalg.cholesky(model.Sig)
  parts.append(factor[np.tril_indices(model.r)])
  return np.concatenate(parts)


def make_varmax(model: varmapack.Model) -> VARMAX:
  dummy = np.zeros((1, model.r))
  return VARMAX(
      dummy, order=(model.p, model.q), trend="n",
      error_cov_type="unstructured", enforce_stationarity=False,
      enforce_invertibility=False,
  )


def vp_draw(model: varmapack.Model, length: int, nrep: int,
            in_place: bool) -> Callable[[], np.ndarray]:
  rng = randompack.Rng()
  rng.seed(123)
  out = None
  if in_place:
    out = np.empty((nrep, length, model.r), dtype=np.float64)
  return lambda: model.sim(length, nrep=nrep, rng=rng, out=out)


def varmax_draw(model: VARMAX, params: np.ndarray, length: int,
                nrep: int) -> Callable[[], np.ndarray]:
  rng = np.random.default_rng(123)
  return lambda: model.simulate(
      params, nsimulations=length, repetitions=nrep, transformed=True,
      random_state=rng,
  )


def finite_result(result: np.ndarray) -> bool:
  return bool(np.all(np.isfinite(result)))


def time_draw(draw: Callable[[], np.ndarray], target: float,
              validate: bool = True) -> float:
  if validate and not finite_result(draw()):
    raise RuntimeError("simulation returned non-finite values")
  if not validate:
    draw()
  gc.collect()
  enabled = gc.isenabled()
  gc.disable()
  calls = 0
  start = time.perf_counter()
  elapsed = 0
  try:
    while elapsed < target:
      draw()
      calls += 1
      elapsed = time.perf_counter() - start
  finally:
    if enabled:
      gc.enable()
  return 1e9*elapsed/calls


def generation_lengths(length: int) -> np.ndarray:
  return length*np.array((10, 20, 40, 80))


def repetition_counts(nrep: int) -> np.ndarray:
  counts = np.rint(nrep*np.array((0.01, 0.02, 0.04, 0.08))).astype(int)
  counts = np.unique(np.maximum(1, counts))
  if len(counts) < 4:
    counts = np.array((1, 2, 4, 8))
  return counts


def fit_slope(x: np.ndarray, y: np.ndarray, description: str) -> float:
  design = np.column_stack((np.ones(len(x)), x))
  slope = float(np.linalg.lstsq(design, y, rcond=None)[0][1])
  if not np.isfinite(slope) or slope <= 0:
    raise RuntimeError(f"invalid fitted {description}; increase -t")
  return slope


def time_vp_generation(model: varmapack.Model, length: int,
                       target: float) -> float:
  lengths = generation_lengths(length)
  times = np.array([
      time_draw(vp_draw(model, int(n), 1, False), target) for n in lengths
  ])
  return fit_slope(lengths, times, "Varmapack generation time")


def time_varmax_path(model: VARMAX, params: np.ndarray, length: int,
                     nrep: int, target: float) -> float:
  counts = repetition_counts(nrep)
  times = np.array([
      time_draw(varmax_draw(model, params, length, int(m)), target, False)
      for m in counts
  ])
  return fit_slope(counts, times, "VARMAX path time")


def time_case(name: str, options: argparse.Namespace) -> dict[str, float] | None:
  model = varmapack.testcase(name)
  varmax = make_varmax(model)
  params = varmax_params(model)
  if not finite_result(varmax_draw(varmax, params, options.length, 1)()):
    print(f"{name}: VARMAX returned non-finite values; skipped", file=sys.stderr)
    return None
  startup = max(model.p, model.q)
  if options.length < startup:
    raise RuntimeError(f"series length is shorter than startup for {name}")
  values = {
      "vp_generation1": time_vp_generation(
          model, options.length, options.target
      ),
      "varmax_path": time_varmax_path(
          varmax, params, options.length, options.nrep, options.target
      ),
      "vp_total1": time_draw(
          vp_draw(model, options.length, 1, False), options.target
      ),
      "vp_setup_batch": time_draw(
          vp_draw(model, startup, options.nrep, False), options.target
      ),
      "vp_return_batch": time_draw(
          vp_draw(model, options.length, options.nrep, False), options.target
      ),
      "vp_in_place_batch": time_draw(
          vp_draw(model, options.length, options.nrep, True), options.target
      ),
      "varmax_batch": time_draw(
          varmax_draw(varmax, params, options.length, options.nrep),
          options.target, False,
      ),
  }
  return {
      "model": name, "p": model.p, "q": model.q, "r": model.r,
      "startup": startup, **values,
  }


def run_benchmarks(options: argparse.Namespace) -> tuple[list[dict[str, float]],
                                                         set[str]]:
  rows = []
  failures = set()
  for trial in range(options.repeats):
    print(
        f"Python interface run {trial + 1} of {options.repeats}",
        file=sys.stderr,
    )
    for name in MODELS:
      row = time_case(name, options)
      if row is None:
        failures.add(name)
        continue
      row["trial"] = trial + 1
      rows.append(row)
  return rows, failures


def median_timings(runs: list[dict[str, float]]) -> list[dict[str, float]]:
  dimensions = ("p", "q", "r", "startup")
  timings = tuple(
      key for key in runs[0] if key not in ("trial", "model", *dimensions)
  )
  rows = []
  names = [name for name in MODELS if any(row["model"] == name for row in runs)]
  for name in names:
    selected = [row for row in runs if row["model"] == name]
    row = {"model": name}
    row.update({key: selected[0][key] for key in dimensions})
    row.update({
        key: float(np.median([sample[key] for sample in selected]))
        for key in timings
    })
    rows.append(row)
  return rows


def summarize(timings: list[dict[str, float]],
              options: argparse.Namespace) -> list[dict[str, float]]:
  rows = []
  for timing in timings:
    tail = options.length - timing["startup"]
    intermediate = (
        timing["vp_setup_batch"]/options.nrep +
        tail*timing["vp_generation1"]
    )
    row = dict(timing)
    row["overall"] = timing["varmax_batch"]/timing["vp_in_place_batch"]
    row["allocation"] = (
        timing["vp_return_batch"]/timing["vp_in_place_batch"]
    )
    row["setup_amortization"] = (
        timing["vp_total1"]/intermediate*
        timing["varmax_batch"]/options.nrep/timing["varmax_path"]
    )
    row["bulk_generation"] = (
        intermediate/(timing["vp_return_batch"]/options.nrep)
    )
    row["implementation"] = timing["varmax_path"]/timing["vp_total1"]
    rows.append(row)
  return rows


def align_markdown(labels: tuple[str, ...], rows: list[tuple[str, ...]]) -> str:
  widths = [
      max(len(row[i]) for row in [labels, *rows]) for i in range(len(labels))
  ]

  def format_row(row: tuple[str, ...]) -> str:
    cells = [row[0].ljust(widths[0])]
    cells.extend(row[i].rjust(widths[i]) for i in range(1, len(row)))
    return "| " + " | ".join(cells) + " |"

  separator = tuple(
      "-"*widths[0] if i == 0 else "-"*(widths[i] - 1) + ":"
      for i in range(len(widths))
  )
  return "\n".join([
      format_row(labels), format_row(separator),
      *(format_row(row) for row in rows),
  ])


def print_table(rows: list[dict[str, float]], options: argparse.Namespace,
                failures: set[str]) -> None:
  labels = (
      "Model", "Overall", "Allocation", "Setup/amort.", "Bulk generation",
      "Implementation",
  )
  results = {row["model"]: row for row in rows}
  formatted = []
  for name in MODELS:
    if name in failures:
      formatted.append((name, "failed", "--", "--", "--", "--"))
      continue
    row = results[name]
    formatted.append(
        (
            name, f'{row["overall"]:.0f}', f'{row["allocation"]:.2f}',
            f'{row["setup_amortization"]:.1f}',
            f'{row["bulk_generation"]:.1f}', f'{row["implementation"]:.1f}',
        )
    )
  print(
      f"Runs: {options.repeats}; target: {options.target} s; "
      f"length: {options.length}; batch: {options.nrep}.\n"
  )
  print("Speed-up factor due to:\n")
  print(align_markdown(labels, formatted))


def main() -> None:
  warnings.filterwarnings("ignore", category=EstimationWarning)
  options = parse_options()
  runs, failures = run_benchmarks(options)
  print_table(summarize(median_timings(runs), options), options, failures)


if __name__ == "__main__":
  main()
