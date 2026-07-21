# Varmapack Development

## C build

Varmapack uses Meson/Ninja. Randompack must be installed so that `pkg-config`
can find `randompack`.

```sh
meson setup release --buildtype=release
ninja -C release
meson test -C release --print-errorlogs
```

Use the `release` build directory for normal development builds and tests.

## MATLAB reference tests

The MATLAB reference implementation lives in `matlab-reference/`. The C tests
compare against files written by `tests/matlabcompare.m`.

```sh
matlab -batch "cd('matlab-reference/tests'); run_reference_tests"
matlab -batch "addpath('tests'); matlabcompare"
meson test -C release RunTests --print-errorlogs
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
