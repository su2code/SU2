#!/usr/bin/env python3
"""Check periodic configuration errors in isolated processes (MPI errors terminate)."""
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_driver")
    parser.add_argument("--streamwise", action="store_true")
    parser.add_argument("--launcher", nargs=argparse.REMAINDER, default=[])
    args = parser.parse_args()
    cases = (
        (
            ("Rotation", "must be purely translational"),
            ("Zero", "must have a nonzero translation"),
            (
                "Convection",
                "MARKER_HEATTRANSFER and MARKER_CHT_INTERFACE are unsupported",
            ),
            ("CHT", "MARKER_HEATTRANSFER and MARKER_CHT_INTERFACE are unsupported"),
        )
        if args.streamwise
        else (
            ("Continuous", "Continuous adjoints do not implement MARKER_PERIODIC"),
            ("Radiation", "RADIATION_MODEL does not implement MARKER_PERIODIC"),
            ("Structure", "SOLVER= ELASTICITY does not implement MARKER_PERIODIC"),
        )
    )
    prefix = "StreamwisePeriodicSupport" if args.streamwise else "PeriodicSupport"
    for tag, message in cases:
        result = subprocess.run(
            args.launcher + [args.test_driver, "[." + prefix + tag + "]"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
        if result.returncode == 0 or message not in result.stdout:
            raise RuntimeError(
                tag
                + " did not report the expected configuration error:\n"
                + result.stdout
            )
        print("PASS: " + tag, flush=True)


if __name__ == "__main__":
    main()
