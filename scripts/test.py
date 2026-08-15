#!/usr/bin/env python3
"""Cross-platform test runner for this repository.

Replaces scripts/test.sh with a Python script so other scripts stay Python-based
and execution works on Windows, macOS and Linux.
"""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"


def run(cmd, cwd=None):
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=cwd)


def main():
    if not shutil.which("cmake"):
        print("Error: 'cmake' not found in PATH", file=sys.stderr)
        sys.exit(1)

    # Configure
    run(["cmake", "-S", str(ROOT), "-B", str(BUILD_DIR)])

    # Build the test runner target
    run(["cmake", "--build", str(BUILD_DIR), "--target", "test_runner", "--config", "Release"]) 

    # Run tests with ctest if available; otherwise try RUN_TESTS target as a fallback
    if shutil.which("ctest"):
        run(["ctest", "-V"], cwd=str(BUILD_DIR))
    else:
        print("'ctest' not found in PATH; attempting to run tests via CMake target 'RUN_TESTS'", file=sys.stderr)
        run(["cmake", "--build", str(BUILD_DIR), "--target", "RUN_TESTS", "--config", "Release"]) 


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        print("Command failed:", e, file=sys.stderr)
        sys.exit(e.returncode)
