#!/usr/bin/env python3
"""Check that an Euler wall corner blocks velocity through both faces."""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    solver = str(Path(sys.argv[1]).resolve())
    fixture = Path(__file__).resolve().parent
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OMPI_ALLOW_RUN_AS_ROOT="1", OMPI_ALLOW_RUN_AS_ROOT_CONFIRM="1")
    for ranks in (1, 2):
        with tempfile.TemporaryDirectory(prefix="sharp_corner_", dir=Path.cwd()) as work:
            work = Path(work)
            for name in ("config.cfg", "mesh.su2"):
                shutil.copy(fixture / name, work / name)
            command = [solver, "config.cfg"]
            if ranks > 1:
                command = ["mpirun", "-n", str(ranks)] + command
            run = subprocess.run(command, cwd=work, env=env, text=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, timeout=120)
            if run.returncode:
                raise RuntimeError(f"SU2 failed on {ranks} rank(s):\n{run.stdout}")
            with (work / "restart_flow.csv").open(newline="") as f:
                rows = list(csv.DictReader(f, skipinitialspace=True))
            corner = next(row for row in rows if int(row["PointID"]) == 12)
            momenta = [float(corner["Momentum_x"]), float(corner["Momentum_y"])]
            scale = max(math.hypot(float(row["Momentum_x"]), float(row["Momentum_y"])) for row in rows)
            if not math.isfinite(scale) or scale <= 0 or any(not math.isfinite(x) or abs(x) > 1e-10 * scale for x in momenta):
                raise AssertionError(f"{ranks} rank(s): corner momentum {momenta}; reference {scale}")
            print(f"{ranks} rank(s): corner momentum {momenta}; reference {scale}")


if __name__ == "__main__":
    main()
