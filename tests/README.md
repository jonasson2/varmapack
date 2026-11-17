# Tests

Varmapack has automated C tests, run after a Meson build with:

```sh
meson test -C build --print-errorlogs
```

`RunTests` is the main test executable. Run it directly with `-v` to show its
test groups, or with `-C` for the longer statistical tests (200,000 rather than
10,000 independent replicates per model). `AgainstMatlab` compares results with
checked-in MATLAB reference fixtures; MATLAB is not needed to run it.

The reported counts are logical checks, not individual array elements or model
cases. Each array comparison counts once while still examining every element;
a failure reports the number of mismatches and the first failing index.
An expected API outcome also counts once, checking both the return value and
the error state. Related elementwise checks are grouped by property: one mean
vector, covariance matrix, simulation recurrence, or support condition. Every
element is still checked with the original tolerance. Distinct models, times,
lags, goodness-of-fit tests, and invalid-input cases remain separate checks.

Allocations, RNG setup, and fixture-file creation are prerequisites rather than
counted checks; failure still stops the test. Empty arrays or property groups
pass as one check. `TestAssertions` verifies these counting rules, including
deliberately failing checks and preservation of the first failure's context.

## Statistical validation

`TestDistribution.c` provides distribution-level tests for 23 VARMA models. It
covers the named test models, four AR persistence levels through $\rho=0.999$,
a $(p,q,r)=(3,3,20)$ model, a near-invertible MA model with spectral radius
0.999, and a model with an ill-conditioned positive-definite innovation
covariance (condition number $10^{10}$). For each model, independently simulated
paths are checked against theoretical quantities at times or lags 0, 1, 2, and 10:

- means, same-time covariance matrices, and lagged cross-covariances are
  compared with their theoretical values using sampling-error bounds;
- normal probability-integral transforms (PITs) are checked by both binned
  chi-square and Kolmogorov--Smirnov goodness-of-fit tests;
- selected pair and triple linear combinations across components and times are
  checked with KS tests of their theoretical normal laws.

Thus the tests cover stationary marginal distributions, autocovariance
structures, and cross-dependencies across a range of orders, dimensions, and
persistence levels. The $\rho=0.999$ cases are explicit near-stationary stress
tests; the MA spectral-radius-0.999 model is a near-invertibility stress test.

## Startup and covariance robustness

The startup covariance $S$ is the covariance of the stacked initial states;
the innovation covariance $\Sigma$ is the covariance of each shock vector.
For unconditional stationary simulation, either may be singular. Conditioning
on supplied stationary initial states requires $S$ to be numerically positive
definite: if its Cholesky factorization fails, an error is returned; otherwise
simulation proceeds without a condition-number warning.

The edge-case and conditional-covariance groups document the supported and
rejected cases directly in executable tests:

- A singular positive-semidefinite innovation covariance is accepted for
  unconditional stationary simulation; the generated shocks and series are
  checked to be finite and to obey the VARMA recursion.
- A singular stationary startup covariance is accepted for unconditional
  simulation, but a supplied stationary initial state is rejected when its
  startup covariance is not positive definite.
- Rank-deficient conditional startup covariances are tested for several models.
  The resulting draws are checked to lie on the correct conditional support.
- Indefinite innovation covariances, invalid dimensions, unavailable required
  startup values, and oversized problems are rejected.
- Nonstationary models with supplied states, including MA models and presample
  shocks, are checked separately. Their required positive-definiteness
  conditions are enforced.

These checks are principally in `TestSimEdgeCases.c` and `TestPsdCondCov.c`.
They complement deterministic comparisons with MATLAB, theoretical
autocovariance tests, and tests of the VYW and Lyapunov covariance solvers.

## Other coverage

`TestCov2sig.c` checks innovation-covariance recovery against scalar formulas,
the explicit VAR(1) formula, independently summed impulse responses, the
Lyapunov solver, and all named VARMA testcases. It also covers lower-triangle
input, in-place conversion, singular positive-semidefinite covariances,
very small and large covariance scales, near-unit persistence, incompatible
targets, singular and ill-conditioned inverse systems, and invalid inputs.
Failed calls must leave the output unchanged.

The remaining `RunTests` groups exercise named testcase construction, spectral
radii, theoretical and sample autocovariances, covariance-to-correlation
conversion, impulse responses, VARMAX simulation, error handling,
reproducibility with seeded random-number generators, and the internal
Lyapunov solver.

## OU and Gamma Times

`TestOU.c` checks stationary means, covariance at startup and later times, and
lagged cross-covariance on irregular times, using 100000 independent paths.
It also checks singular/zero covariance, incompatible persistence/covariance,
and input validation. Gamma-time tests check gap means and variances for shapes
below, at, and above one, plus constant spacing and invalid inputs. Language
tests exercise the same model, shape, and RNG conventions through each wrapper.
