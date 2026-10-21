#!/usr/bin/env python3
"""Check the fatal configuration diagnostic in its own process."""
import subprocess
import sys

result = subprocess.run(
    [sys.argv[1], "[.MultigridLimitError]"], capture_output=True, text=True, timeout=30
)
expected = "MGLEVEL must be between 0 and 10"
if result.returncode == 0 or expected not in result.stdout + result.stderr:
    raise RuntimeError(
        "Unsupported MGLEVEL did not fail with its diagnostic:\n"
        + result.stdout
        + result.stderr
    )
print("PASS: MGLEVEL=11 is rejected with its supported range")
