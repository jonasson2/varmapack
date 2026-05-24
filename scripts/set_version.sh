#!/bin/sh
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(dirname "$SCRIPT_DIR")
cd "$REPO_ROOT"

if [ -z "${1:-}" ]; then
  echo "Usage: $0 NEW_VERSION"
  exit 1
fi

VER="$1"

echo "$VER" > VERSION
sed -i.bak -E "s/^project\\('varmapack', 'c', 'fortran', version : '[^']*'/project('varmapack', 'c', 'fortran', version : '$VER'/" meson.build
rm -f meson.build.bak
sed -i.bak -E "s/^#define[[:space:]]+VARMAPACK_VERSION[[:space:]]+\".*\"/#define VARMAPACK_VERSION \"$VER\"/" src/varmapack.h
rm -f src/varmapack.h.bak
sed -i.bak -E "s/^  version: '[^']*'/  version: '$VER'/" python/meson.build
rm -f python/meson.build.bak
sed -i.bak -E "s/^version[[:space:]]*=[[:space:]]*\".*\"/version = \"$VER\"/" python/pyproject.toml
rm -f python/pyproject.toml.bak

echo "Version set to $VER"
