#!/usr/bin/env python3

## \file run_aeroacoustics.py
#  \brief Regression driver for aeroacoustic pressure-history analysis.
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import argparse
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace


su2_run = os.environ.get("SU2_RUN")
if su2_run:
    sys.path.insert(0, su2_run)

from aeroacoustics import analyze  # noqa: E402


SEGMENT_LENGTH = 2048


def main():
    parser = argparse.ArgumentParser(
        description="Run the square-cylinder aeroacoustics regression analysis."
    )
    parser.add_argument("history", help="SU2 pressure-history CSV")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as directory:
        summary = analyze(
            SimpleNamespace(
                input=Path(args.history),
                pressure=["mic1"],
                sample_rate=None,
                time="Cur_Time",
                pressure_scale=1.0,
                reference_pressure=20.0e-6,
                skip_samples=0,
                segment_length=SEGMENT_LENGTH,
                overlap=0.5,
                window="hann",
                detrend="mean",
                output=Path(directory) / "spectrum.csv",
                summary=None,
            )
        )

    metrics = summary["signals"]["mic1"]
    print("\n------------------------------ Begin Solver -----------------------------\n")
    print(
        "| 0 | {:.10f} | {:.10f} |".format(
            metrics["oaspl_db"], metrics["dominant_frequency_hz"]
        )
    )
    print(
        "Analyzed {} pressure samples using {} Welch segments.".format(
            summary["samples"], metrics["segments"]
        )
    )


if __name__ == "__main__":
    main()
