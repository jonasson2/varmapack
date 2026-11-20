# Varmapack Development

This file describes local development of the C library and its MATLAB, Python,
and R interfaces. For ordinary installation, see [README.md](README.md).

## C library

### Configure and build

Varmapack development uses an installed Randompack release. Install it
system-wide or under `~/.local` as described in the root README, including the
`PKG_CONFIG_PATH` setting for a user-local installation. Before configuring
Varmapack, verify that `pkg-config` finds it:

```sh
pkg-config --modversion randompack
pkg-config --cflags --libs randompack
```

Meson obtains the Randompack compiler and linker flags through `pkg-config`; no
Randompack paths need to be passed to `meson setup`. From the repository root,
configure and build an optimized development tree:

```sh
meson setup build --buildtype=release
ninja -C build
```

For an existing build directory whose Meson configuration has changed, use:

```sh
meson setup build --reconfigure
```

### Test

Run all C tests with:

```sh
meson test -C build --print-errorlogs
```

This runs `RunTests`, including the Lyapunov-solver checks, and `AgainstMatlab`.
`AgainstMatlab` normally uses the checked-in MATLAB reference fixtures, so
MATLAB is not needed for the ordinary C test run.

### Compiler checks

Use `scripts/compiler-check.sh` to recreate a release build for an explicit C
compiler, Fortran compiler, and BLAS provider, then run the C test suite:

```sh
scripts/compiler-check.sh -c clang -f gfortran -b openblas
```

The build directory is named from that compiler combination. `-s` instead
creates a debug build with C address and undefined-behavior sanitizers. On Elja
and Spark1, source the ignored machine-local preparation helper first; use
`scripts/elja/prepare.sh -h` or `scripts/spark1/prepare.sh -h` for its compiler
and library options.

## Internal matrix exponential

`src/mb05odGateway.h` declares the internal C call
`mb05odGateway(A, delta, n, E, accuracy)`, computing `E = exp(delta*A)` for
column-major real matrices. It preserves `A` (unless `E == A`) and leaves `E`
unchanged on failure. The optional `mb05od_accuracy` receives the SLICOT
`MDIG`, `IDIG`, and `IWARN` values. A successful return with a nonzero warning
must be assessed by the caller; it does not imply full accuracy.

The gateway uses SLICOT MB05OD with balancing and Padé order 9. MB05OD and
its two SLICOT dependencies, MB04MD and MB05OY, are bundled unchanged in
`src/mb05od-complete.F`, with provenance, license, and SciPy-OpenBLAS symbol
mapping in the preamble. No system SLICOT installation is needed. Matrix-drift
OU uses `mb05odExponential` without balancing on a scaled Van Loan block matrix;
the diagonal persistence-time path does not use the gateway.
With GNU Fortran, `meson test` also compiles this source using
`-std=f2008 -pedantic-errors -Werror`, independently of the ordinary build's
warning suppression for bundled Fortran. This retains the upstream labelled
DO loops, which Fortran 2008 permits.

## MATLAB reference implementation

The independent MATLAB reference functions are in `matlab-reference/`. They use
small test-only MEX gateways for the needed Randompack functions. The gateways
are platform-specific and not tracked by Git.

Before running the MATLAB reference tests on a new machine, build the Randompack
gateways. The Randompack source tree must already have been built in its `build`
directory:

```matlab
mex -setup C                         % only needed once per MATLAB installation
addpath("tests/matlab")
build_randompack_mex("/path/to/randompack")
```

From the Varmapack repository root, run the MATLAB reference tests and refresh
the comparison fixtures with:

```matlab
cd("matlab-reference/tests")
run_reference_tests
cd("../../tests")
matlabcompare
```

Then run the C comparison from the repository root:

```sh
meson test -C build AgainstMatlab --print-errorlogs
```

`matlabcompare` rewrites `tests/matlabcompare.txt` and
`tests/matlabcompare-rolling.txt`. Review and commit deliberate fixture changes.

## Python package

### First-time setup

Use a virtual environment in the Varmapack repository root. These commands use
the current Randompack checkout; replace its path as needed.

```sh
uv venv
source .venv/bin/activate
uv pip install meson-python meson ninja cython numpy pytest
uv pip install -r python/docs/requirements.txt
uv pip install -e /path/to/randompack/python --no-build-isolation
scripts/syncpy.sh
uv pip install -e python --no-build-isolation
python -c "import varmapack; print(varmapack.__file__)"
```

For the Python build, Randompack supplies its C API through the installed
Python package. A separately installed Randompack C library is not required.

### Randompack interoperability

The Randompack extension exposes a private, versioned function table through
the `randompack._core._C_API` capsule. Varmapack imports that table when its
extension is loaded, and therefore does not link against `librandompack`.

When `Model.sim()` or `testcase()` receives a `randompack.Rng`, Varmapack uses
the borrowed `randompack_rng *` supplied by that object. Randompack retains
ownership, so the Python object must remain alive for the duration of the call.
When no generator is supplied, Varmapack creates and frees a temporary RNG
through the same function table.

### After changing C or Cython sources

The Python package compiles the synchronized source snapshot in `python/src/`,
not the top-level `src/` directory. After changing `src/`, synchronize it before
testing or importing a rebuilt extension:

```sh
scripts/syncpy.sh
```

Restart an active Python or Jupyter session after rebuilding the extension.

### Test and document

Run the Python tests with:

```sh
python -m pytest -q python/tests
```

Build and open the Python documentation with:

```sh
cd python/docs
make html
make open
```

## R package

The R package in `r-package/` depends on the current Randompack R package.
Install that package first, then synchronize and install Varmapack from the
repository root:

```sh
scripts/syncR.sh
R CMD INSTALL --preclean r-package
```

`syncR.sh` copies the C source snapshot to `r-package/src/` and regenerates the
mathematical-description vignette include from `python/docs/math.md`. Run it
after every change to `src/` before rebuilding the R package. Do not edit the
synchronized files in `r-package/src/` directly.

Run the R tests and package check with:

```sh
Rscript -e 'testthat::test_local("r-package")'
R CMD build r-package
R CMD check --as-cran varmapack_*.tar.gz
```

R interface documentation is written as roxygen comments in `r-package/R/`.
After changing it, regenerate the manual pages, namespace, and README (requires
the `devtools` package):

```sh
Rscript -e 'devtools::document("r-package")'
Rscript -e 'devtools::build_readme("r-package")'
```

## Source synchronization and versions

The C sources are canonical. `scripts/syncpy.sh` and `scripts/syncR.sh` copy
them to the Python and R package trees respectively. Run the relevant script
before building or testing an interface after changing `src/`.

Set the shared C, Python, and R package version with:

```sh
scripts/set_version.sh 0.1.0
```

The script updates `VERSION`, the C library version, Python package metadata,
and `r-package/DESCRIPTION`.

## Benchmarking

The regular timing programs are described in the root README:

- `TimeSimulate` compares the common named-testcase workload across C, MATLAB,
  and Python.
- `TimeScalability` varies one model characteristic at a time around a reference
  workload.

Two specialized programs investigate the startup covariance solvers. Build them
from the repository root with:

```sh
ninja -C build benchmark/TimeSetup benchmark/TimeBreakEven
```

`TimeSetup` measures complete startup work, reporting the VYW and SLICOT
Lyapunov solver portions separately. Run it with:

```sh
build/benchmark/TimeSetup
```

`TimeBreakEven` compares only those two solvers over a grid of model orders and
dimensions. Summarize its crossover tables with:

```sh
build/benchmark/TimeBreakEven -t 0.05 | \
  awk -f benchmark/SummarizeBreakEven.awk
```

Cross-platform crossover results, including the averaged cutoffs used for
automatic solver selection, are recorded in `vyw-slicot-breakeven.txt`.

All benchmark programs accept `-h` for workload and timing options.

## Releases

### 1. Prepare the candidate

Work from the repository root. Set `version` to the intended release, verify
the required Randompack versions are published on PyPI and CRAN, and update
`CHANGELOG.md`, `r-package/NEWS.md`, and the Python lockfile as needed.

```sh
scripts/set_version.sh "$version"   # only if a bump is needed
scripts/syncpy.sh
scripts/syncR.sh
scripts/generate-readme.sh
Rscript -e 'devtools::document("r-package")'
Rscript -e 'devtools::build_readme("r-package")'
git diff --check
```

Commit the candidate and test from a clean checkout against published dependencies:
C, Python, R, the MATLAB bindings, and the documentation examples affected by
the release. Push the tested commit. Keep release-specific SHAs, run links,
tarball checksums, and outstanding results in the local `codex-plans/PUBLISH.md`.

### 2. Verify and publish on GitHub

Dispatch a trial build on the candidate branch, for example
`gh workflow run python-wheels.yml --ref main -f build='cp313-*'`.
Choose a supported selector from the current [wheel workflow](.github/workflows/python-wheels.yml);
that file defines the platforms, runners, and any per-job overrides.
Use `gh run view "$trial_run" --json headSha,conclusion,jobs` to verify the SHA
and that every expected job passed. If a fix is needed, push it and repeat.

Tag the **tested SHA**, not whatever `main` points to later:

```sh
git tag -a "v$version" "$sha" -m "Varmapack $version"
git push origin "refs/tags/v$version"
gh release create "v$version" --verify-tag --title "Varmapack $version" --notes-file "$notes"
```

Use the version's changelog entry for `$notes`. Publishing the GitHub release,
not pushing the tag alone, starts the full wheel matrix. This does not upload to PyPI.

### 3. Publish to PyPI

Wait for the full **release** run to pass and verify its SHA against the tag.
From a clean checkout of that tag, using a Python environment with `build` and
`twine` installed, collect the artifacts in a fresh directory:

```sh
work=$(mktemp -d)
mkdir "$work/dist"
python -m build --sdist --outdir "$work/dist" python/
gh run download "$release_run" --pattern 'wheels-*' --dir "$work/wheels"
find "$work/wheels" -name '*.whl' -exec cp {} "$work/dist/" \;
python -m twine check "$work"/dist/*
```

Verify the directory contains one sdist and a wheel for every supported target;
`twine check` does not check matrix completeness. Upload to TestPyPI, then
smoke-test in a fresh environment. Install runtime dependencies from regular
PyPI (see [pyproject.toml](python/pyproject.toml)) and only Varmapack from TestPyPI:

```sh
python -m twine upload --repository testpypi "$work"/dist/*
python -m venv "$work/smoke"
"$work/smoke/bin/python" -m pip install numpy randompack
"$work/smoke/bin/python" -m pip install --index-url https://test.pypi.org/simple/ \
    --no-deps --only-binary=:all: "varmapack==$version"
"$work/smoke/bin/python" -m pip check
"$work/smoke/bin/python" -I -c 'import numpy as np, varmapack as vp; np.testing.assert_allclose(vp.cov2corr(np.eye(2)), np.eye(2))'
```

Only after the smoke test passes, upload the **unchanged** files with
`python -m twine upload "$work"/dist/*`. See the [PyPA publishing guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/) for current authentication instructions.

### 4. Submit to CRAN

This can proceed independently of the Python wheels. Check the current
[CRAN submission policy](https://cran.r-project.org/web/packages/policies.html#Submission)
before building. Use current R-release or R-patched for the build, and current
R-devel for the check; set the two variables below to the respective executables.

```sh
"$R_RELEASE" CMD build r-package
"$R_DEVEL" CMD check --as-cran "varmapack_${version}.tar.gz"
```

Keep that tarball and its SHA-256. For native-code changes, check the same
archive on Windows R-devel even for a package update. A [Win-Builder](https://win-builder.r-project.org/)
upload is **pending**, not a pass: retrieve `00check.log` from the link emailed
to the package maintainer. A delayed report holds up CRAN, not GitHub or PyPI;
a defect in shared code warrants reassessing all affected builds.

Record results, any explained notes, and reverse-dependency checks in
`r-package/cran-comments.md`. Resolve errors, warnings, and significant notes.
If the tarball is rebuilt, recheck the new file. Submit the exact checked
archive through the [CRAN form](https://cran.r-project.org/submit.html) and
confirm the maintainer email.
