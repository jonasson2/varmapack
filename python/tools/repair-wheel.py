#!/usr/bin/env python3

import os
from pathlib import Path
import subprocess
import sys

import scipy_openblas32


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: repair-wheel.py WHEEL DESTINATION")
    wheel = Path(sys.argv[1]).resolve()
    destination = Path(sys.argv[2]).resolve()
    libdir = Path(scipy_openblas32.get_lib_dir()).resolve()
    environment = os.environ.copy()
    if sys.platform.startswith("linux"):
        old_path = environment.get("LD_LIBRARY_PATH", "")
        paths = [str(libdir)]
        if old_path:
            paths.append(old_path)
        environment["LD_LIBRARY_PATH"] = os.pathsep.join(paths)
        command = ["auditwheel", "repair", "-w", str(destination), str(wheel)]
    elif sys.platform == "win32":
        command = [
            sys.executable,
            "-m",
            "delvewheel",
            "repair",
            "--add-path",
            str(libdir),
            "-w",
            str(destination),
            str(wheel),
        ]
    else:
        raise SystemExit(f"unsupported repair platform: {sys.platform}")
    subprocess.run(command, check=True, env=environment)


if __name__ == "__main__":
    main()
