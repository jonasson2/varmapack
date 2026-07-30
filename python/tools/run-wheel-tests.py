#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: run-wheel-tests.py TEST-DIRECTORY")
    test_directory = Path(sys.argv[1]).resolve()
    tests = sorted(test_directory.glob("test_*.py"))
    if not tests:
        raise SystemExit(f"no tests found in {test_directory}")
    for test in tests:
        print(f"Running {test.name}", flush=True)
        subprocess.run([sys.executable, str(test)], check=True)


if __name__ == "__main__":
    main()
