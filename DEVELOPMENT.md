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

The Python interface is built from the `python/` subtree.

```sh
cd python
python -m pip install -e . --no-build-isolation
python -c "import varmapack; print(varmapack.__file__)"
```

The editable install assumes the C Randompack library and the Python
`randompack` package are both available in the active environment.

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
