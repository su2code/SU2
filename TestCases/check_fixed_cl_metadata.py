#!/usr/bin/env python3
"""Regression for #2937: screen frequency must not change fixed-CL metadata.

Run with a serial/OpenMP SU2_CFD and the NACA0012 mesh from the TestCases
repository, for example:
  python3 check_fixed_cl_metadata.py --su2-cfd /path/to/SU2_CFD \
      --mesh fixed_cl/naca0012/mesh_NACA0012_inv.su2

All runs start from scratch in temporary directories. No reference derivatives
are hard-coded: identical solves must produce identical metadata at frequencies
1, 10, 50 and 1000, including the pre-finite-difference restart iteration.
"""

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def set_option(config, name, value):
    pattern = rf"(?m)^\s*{name}\s*=.*$"
    line = f"{name}= {value}"
    return re.sub(pattern, line, config) if re.search(pattern, config) else config + "\n" + line


def read_metadata(path):
    return {name.strip(): float(value) for name, value in
            (line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--su2-cfd", required=True, type=Path)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--timeout", default=180, type=int)
    args = parser.parse_args()
    binary, mesh = args.su2_cfd.resolve(), args.mesh.resolve()
    env = dict(os.environ, OMP_NUM_THREADS="1")
    config = (Path(__file__).parent / "fixed_cl/naca0012/inv_NACA0012.cfg").read_text()
    config = set_option(config, "RESTART_SOL", "NO")
    config = set_option(config, "MESH_FILENAME", "mesh.su2")
    reference = None
    with tempfile.TemporaryDirectory(prefix="su2-fixed-cl-") as tmp:
        for frequency in (1, 10, 50, 1000):
            run = Path(tmp) / str(frequency)
            run.mkdir()
            shutil.copyfile(mesh, run / "mesh.su2")
            (run / "case.cfg").write_text(set_option(config, "SCREEN_WRT_FREQ_INNER", frequency))
            result = subprocess.run([str(binary), "case.cfg"], cwd=run, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    timeout=args.timeout, env=env)
            if result.returncode:
                raise AssertionError(f"SU2 failed at frequency {frequency}:\n{result.stdout[-4000:]}")
            metadata = read_metadata(run / "flow.meta")
            for key in ("ITER", "AOA", "DCD_DCL_VALUE", "DCMZ_DCL_VALUE"):
                assert key in metadata and math.isfinite(metadata[key]), f"Invalid {key} at {frequency}"
            if reference is None:
                reference = metadata
            else:
                assert metadata == reference, f"Metadata changed at frequency {frequency}: {metadata} != {reference}"
            print(f"PASS screen frequency {frequency}: {metadata}", flush=True)


if __name__ == "__main__":
    main()
