#!/usr/bin/env bash
set -eu

python -m pip install -r python/docs/requirements.txt
python -m pip install meson meson-python ninja cython numpy scipy-openblas32
python -m pip install "randompack @ git+https://github.com/jonasson2/randompack.git#subdirectory=python"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

git clone --depth 1 https://github.com/jonasson2/randompack.git "$tmpdir/randompack"
meson setup "$tmpdir/randompack-build" "$tmpdir/randompack" \
  --prefix="$READTHEDOCS_VIRTUALENV_PATH" -Dbuildtype=release -Dblas=builtin
meson install -C "$tmpdir/randompack-build"

scipy_pc=$(python - <<'PY'
import pathlib
import scipy_openblas32

print(pathlib.Path(scipy_openblas32.__file__).parent / "lib" / "pkgconfig")
PY
)

export PKG_CONFIG_PATH="$READTHEDOCS_VIRTUALENV_PATH/lib/pkgconfig:$scipy_pc:${PKG_CONFIG_PATH:-}"
python -m pip install ./python
