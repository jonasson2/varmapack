#!/usr/bin/env python3
"""Expand Markdown includes, adjust headings, and translate shared macros."""

import argparse
import pathlib
import re
import sys


INCLUDE_RE = re.compile(r"^\s*<!--\s*include\s+([^>]+?)\s*-->\s*$")
HEADING_RE = re.compile(r"^(#{1,6})([ \t]+.*)$")
MACROS = {
    r"\Cov": r"\operatorname{Cov}",
    r"\Var": r"\operatorname{Var}",
    r"\eps": r"\varepsilon",
}


def shift_heading(line, offset):
    match = HEADING_RE.match(line)
    if not match or offset == 0:
        return line
    level = min(max(len(match.group(1)) + offset, 1), 6)
    return "#"*level + match.group(2) + ("\n" if line.endswith("\n") else "")


def read_document(path, root, stack, offset, include_offset, drop_first_heading):
    if path in stack:
        chain = " -> ".join(str(item) for item in stack + [path])
        raise ValueError(f"recursive include: {chain}")
    lines = []
    drop_heading = drop_first_heading
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        match = INCLUDE_RE.match(line)
        if match:
            include = root / match.group(1).strip()
            if not include.is_file():
                raise ValueError(f"include not found: {include}")
            lines.append(read_document(include, root, stack + [path],
                                       offset + include_offset, include_offset, False))
            continue
        if drop_heading and HEADING_RE.match(line):
            drop_heading = False
            continue
        lines.append(shift_heading(line, offset))
    return "".join(lines)


def expand_macros(text):
    for name, replacement in MACROS.items():
        text = text.replace(name, replacement)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path.cwd(),
                        help="repository root used to resolve include paths")
    parser.add_argument("--heading-offset", type=int, default=0,
                        help="heading-level adjustment for the input document")
    parser.add_argument("--include-heading-offset", type=int, default=0,
                        help="additional heading-level adjustment for each include")
    parser.add_argument("--drop-first-heading", action="store_true",
                        help="omit the first ATX heading in the input document")
    parser.add_argument("input", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    args = parser.parse_args()
    root = args.root.resolve()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    try:
        text = read_document(input_path, root, [], args.heading_offset,
                             args.include_heading_offset, args.drop_first_heading)
    except ValueError as error:
        sys.exit(str(error))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(expand_macros(text), encoding="utf-8")


if __name__ == "__main__":
    main()
