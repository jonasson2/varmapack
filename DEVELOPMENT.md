# Varmapack Development

## C build

Varmapack uses Meson/Ninja. Randompack must be installed so that `pkg-config`
can find `randompack`.

```sh
meson setup build --buildtype=release
ninja -C build
meson test -C build --print-errorlogs
```

Use the `build` directory for normal optimized development builds and tests.

## Setup benchmarks

`TimeSetup` measures complete startup preparation using both the modified
Yule-Walker (VYW) and SLICOT Lyapunov covariance solvers. It reports each
solver's time separately from the remaining preparation work that turns its
result into a startup draw. Build and run it with:

```sh
ninja -C build benchmark/TimeSetup
build/benchmark/TimeSetup
```

`TimeBreakEven` times only the VYW and Lyapunov solvers over a grid of model
orders and dimensions. Its output identifies where the Lyapunov route becomes
faster; `SummarizeBreakEven.awk` produces first-winning and stable break-even
tables:

```sh
ninja -C build benchmark/TimeBreakEven
build/benchmark/TimeBreakEven -t 0.05 | \
    awk -f benchmark/SummarizeBreakEven.awk
```

## MATLAB reference tests

The MATLAB reference implementation lives in `matlab-reference/`.
`tests/matlabcompare.m` writes fixtures for the `AgainstMatlab` C comparison
test.

In a MATLAB session, change the current folder to the repository root, then run:

```matlab
cd("matlab-reference/tests")
run_reference_tests
cd("../../tests")
matlabcompare
```

From the repository root, run the C comparison with:

```sh
meson test -C build AgainstMatlab --print-errorlogs
```

## Python package

Use a virtual environment in the repository root for Python development.

```sh
uv venv
uvactivate
uv pip install meson-python meson ninja cython pytest numpy
uv pip install -r python/docs/requirements.txt
uv pip install /Users/jonasson/randompack/python
uv pip install -e python --no-build-isolation
python -c "import varmapack; print(varmapack.__file__)"
```

The editable install assumes the C Randompack library is installed so that
`pkg-config` can find `randompack`. The Python `randompack` package should be
installed in the same virtual environment.

Run the Python smoke tests directly:

```sh
for t in python/tests/test_*.py; do python "$t" || exit 1; done
```

Build the Python documentation with:

```sh
cd python/docs
make html
```

## Source synchronization

The Python package keeps a snapshot of C sources in `python/src`.

```sh
scripts/syncpy.sh
```

Version strings are updated with:

```sh
scripts/set_version.sh 0.1.0
```

R packaging is kept separate for now.
