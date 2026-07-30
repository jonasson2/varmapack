#!/bin/sh
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(dirname "$SCRIPT_DIR")

exec python3 "$SCRIPT_DIR/process-markdown.py" \
  --root "$REPO_ROOT" \
  --include-heading-offset 1 \
  "$REPO_ROOT/README.in" "$REPO_ROOT/README.md"
