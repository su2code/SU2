#!/usr/bin/env python3
"""Run isolated convergence-migration probes (SU2_MPI::Error terminates the process)."""

import argparse
import subprocess


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
    cases = (
        ("MAX", "MAX_DENSITY"),
        ("BGS", "BGS_DENSITY"),
        ("RelativeMAX", "REL_MAX_DENSITY"),
        ("RelativeBGS", "REL_BGS_DENSITY"),
        ("ZonedMAX", "MAX_DENSITY[0]"),
        ("ZonedBGS", "BGS_DENSITY[0]"),
        ("ZonedRelativeMAX", "REL_MAX_DENSITY[0]"),
        ("ZonedRelativeBGS", "REL_BGS_DENSITY[0]"),
        ("Mixed", "MAX_DENSITY"),
    )
    for tag, field in cases:
        command = args.launcher + [args.test_driver, "[.NEMOObsolete" + tag + "]"]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=30,
            check=False,
        )
        base, separator, zone = field.partition("[")
        replacement = base + "_0" + separator + zone
        expected = (
            "Obsolete NEMO CONV_FIELD '" + field + "'",
            "'" + replacement + "' for species 0",
            "the obsolete field cannot be ignored",
        )
        if result.returncode == 0 or any(
            message not in result.stdout for message in expected
        ):
            raise RuntimeError(
                "Migration probe "
                + tag
                + " did not fail with its migration diagnostic:\n"
                + result.stdout
            )
        print(
            "PASS: " + tag + " rejects " + field + " with explicit-species guidance",
            flush=True,
        )


if __name__ == "__main__":
    main()
