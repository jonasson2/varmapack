#!/bin/sh
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(dirname "$SCRIPT_DIR")
cd "$REPO_ROOT"

echo "Syncing core sources to r-package/src..."
rsync -av --delete \
  --exclude='*_R.c' \
    --exclude='init.c' \
    --exclude='Makevars' \
    --exclude='Makevars.in' \
    --exclude='RandompackRGateway.h' \
    --exclude='varmapack_R.c' \
  --exclude='.DS_Store' \
  --exclude='meson.build' \
  --exclude='printX.c' \
  --exclude='tests/' \
  --exclude='varmapack.def' \
  src/ \
  r-package/src/

echo "Copying license..."
cp -f LICENSE r-package/inst/THIRD-PARTY-NOTICES

echo "Sync complete."
