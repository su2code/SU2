#!/usr/bin/env python3
"""Run the hidden CConfig tests that must stop with an error (SU2_MPI::Error terminates the process)."""

import argparse
import subprocess

# Each hidden test tag, with the text its error message has to contain.
CASES = (
    ("HeatSourceHeatEquation", "HEAT_SOURCE= YES requires a RADIATION_MODEL"),
    ("HeatSourceIncompressible", "HEAT_SOURCE= YES requires a RADIATION_MODEL"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_driver")
    parser.add_argument(
        "--launcher",
        nargs=argparse.REMAINDER,
        default=[],
        help="Optional MPI launcher and arguments, for example mpiexec -n 2",
    )
    args = parser.parse_args()
    for tag, expected in CASES:
        command = args.launcher + [args.test_driver, "[." + tag + "]"]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
        if result.returncode == 0 or expected not in result.stdout:
            raise RuntimeError(
                "Config check "
                + tag
                + " did not stop with the expected error:\n"
                + result.stdout
            )
        print("PASS: " + tag + " stops with: " + expected, flush=True)


if __name__ == "__main__":
    main()
