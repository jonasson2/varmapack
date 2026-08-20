# Varmapack Development

This file describes local development of the C library and its MATLAB, Python,
and R interfaces. For ordinary installation, API examples, testing, and timing,
see the repository [README.md](README.md).

## C library

### Configure and build

Varmapack uses Meson and Ninja. Its C dependency, Randompack, must already be
installed where `pkg-config` can find it. Follow the root README installation
instructions until this succeeds:

```sh
pkg-config --modversion randompack
```

From the repository root, configure and build an optimized development tree:

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

This runs `RunTests`, `AgainstMatlab`, and `TestLyapunov`. `AgainstMatlab`
normally uses the checked-in MATLAB reference fixtures, so MATLAB is not needed
for the ordinary C test run.

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
uv pip install meson-python meson ninja cython numpy
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

The current Python tests are executable scripts:

```sh
for t in python/tests/test_*.py; do python "$t" || exit 1; done
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

The release steps below are repository-specific. Before any release, work from a
clean tree, choose a new version `X.Y.Z`, and regenerate synchronized and
generated files:

```sh
scripts/set_version.sh X.Y.Z
scripts/syncpy.sh
scripts/syncR.sh
scripts/generate-readme.sh
Rscript -e 'devtools::document("r-package")'
Rscript -e 'devtools::build_readme("r-package")'
git diff --check
```

Run the C, Python, and R checks above. Commit the resulting files, push the
commit, and tag that exact commit as `vX.Y.Z` before publishing artifacts.

### PyPI

PyPI releases should contain an sdist and wheels for every supported Python,
operating-system, and architecture combination. Do not publish a source-only
release: building the extension requires a BLAS/LAPACK implementation and a
Fortran compiler.

Publish the required Randompack Python version first. Then update
`python/uv.lock` and build the Varmapack artifacts against that release.

The `Python wheels` GitHub Actions workflow builds the supported wheel matrix.
macOS wheels use Accelerate. Linux and Windows wheels use the SciPy OpenBLAS
runtime and repair the resulting wheel so that its native libraries are
included. Linux ARM wheels for CPython 3.10 through 3.13 use `manylinux2014`;
CPython 3.14 uses `manylinux_2_28` because that is the oldest platform tag
available for its NumPy build dependency.

For a trial run, select one Python version when manually dispatching the
workflow:

```sh
gh workflow run python-wheels.yml -f build='cp313-*'
```

A published GitHub release runs the complete matrix. Start with a clean
`python/dist/`, build the source distribution, then collect the wheels produced
by the same commit:

```sh
python3 -m venv /tmp/varmapack-release-tools
/tmp/varmapack-release-tools/bin/python -m pip install -U build twine
rm -rf python/dist
/tmp/varmapack-release-tools/bin/python -m build --sdist python
# Download or copy the CI-built wheels to python/dist/.
/tmp/varmapack-release-tools/bin/python -m twine check python/dist/*
```

First upload and test the exact artifacts on TestPyPI, using a fresh virtual
environment. Then upload the unchanged artifacts to PyPI:

```sh
/tmp/varmapack-release-tools/bin/python -m twine upload --repository testpypi \
  python/dist/*
# Install and smoke-test varmapack from TestPyPI in a fresh environment.
/tmp/varmapack-release-tools/bin/python -m twine upload python/dist/*
```

TestPyPI and PyPI have separate accounts and cannot replace an already uploaded
version. See the [PyPA packaging guide][pypa-packaging] for credentials and
TestPyPI installation details.

### CRAN

CRAN submission requires the Randompack R package to be available on CRAN,
because Varmapack lists it in `Imports` and `LinkingTo`. Build the exact source
tarball and check it before submission:

```sh
R CMD build r-package
R CMD check --as-cran varmapack_X.Y.Z.tar.gz
```

Resolve all warnings and significant notes. For a new submission, also obtain a
Windows R-devel check through Win-Builder. Submit the unchanged tarball through
the [CRAN submission form][cran-submit] and confirm the email sent by CRAN. The
[CRAN submission checklist][cran-checklist] and repository policy cover the
current requirements.

[pypa-packaging]: https://packaging.python.org/en/latest/tutorials/packaging-projects/
[cran-submit]: https://CRAN.R-project.org/submit.html
[cran-checklist]: https://cran.r-project.org/web/packages/submission_checklist.html
