#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
python3 - <<'PY'
import pathlib
import re

root = pathlib.Path.cwd()
include_re = re.compile(r"^\s*<!--\s*include\s+([^>]+?)\s*-->\s*$")
macros = {
    r"\Cov": r"\operatorname{Cov}",
    r"\Var": r"\operatorname{Var}",
    r"\eps": r"\varepsilon",
}


def demote_heading(line):
    match = re.match(r"^(#{1,6})(\s+.*)$", line)
    if not match or len(match.group(1)) == 6:
        return line
    return "#" + line


def read_with_includes(path, stack, demote):
    if path in stack:
        chain = " -> ".join(str(p) for p in stack + [path])
        raise SystemExit(f"recursive README include: {chain}")
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        match = include_re.match(line)
        if match:
            include = root / match.group(1).strip()
            if not include.is_file():
                raise SystemExit(f"README include not found: {include}")
            lines.append(read_with_includes(include, stack + [path], True))
        else:
            if demote:
                line = demote_heading(line)
            lines.append(line)
    return "".join(lines)


def expand_macros(text):
    for name, replacement in macros.items():
        text = text.replace(name, replacement)
    return text


readme = expand_macros(read_with_includes(root / "README.in", [], False))
(root / "README.md").write_text(readme, encoding="utf-8")
PY
