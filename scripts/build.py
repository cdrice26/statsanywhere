#!/usr/bin/env python3
"""Cross-platform builder for the library target (statsanywhere).

Usage:
  python scripts/build.py          # build statsanywhere in Release
  python scripts/build.py -t name  # build a different target
  python scripts/build.py --clean  # remove build/ and re-configure
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
DEFAULT_TARGET = "statsanywhere"


def run(cmd, cwd=None):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def main():
    p = argparse.ArgumentParser(description="Build the library target using CMake")
    p.add_argument("-t", "--target", default=DEFAULT_TARGET, help="CMake target to build")
    p.add_argument("-c", "--config", default="Release", help="Build configuration")
    p.add_argument("--clean", action="store_true", help="Remove build/ before configuring")
    args = p.parse_args()

    if not shutil.which("cmake"):
        print("Error: 'cmake' not found in PATH", file=sys.stderr)
        sys.exit(1)

    if args.clean and BUILD_DIR.exists():
        print("Removing build/ directory")
        shutil.rmtree(BUILD_DIR)

    # Configure
    run(["cmake", "-S", str(ROOT), "-B", str(BUILD_DIR)])

    # Build requested target
    run(["cmake", "--build", str(BUILD_DIR), "--target", args.target, "--config", args.config])


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print("Command failed:", e, file=sys.stderr)
        sys.exit(e.returncode)
