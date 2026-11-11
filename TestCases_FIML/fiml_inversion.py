#!/usr/bin/env python3
"""Classic FIML field inversion: L-BFGS-B on beta at every point, gradients from the SU2 discrete adjoint.

beta lives in the Beta_FIML column of the ASCII flow restart. Each evaluation runs the direct solver
(objective = ComboObj in the history file) and the discrete adjoint (Sensitivity_Beta_FIML in the adjoint restart).
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize


def write_config(base, path, overrides):
    """Copy a config file, replacing or appending the given options."""
    done = set()
    with open(base) as f_in, open(path, "w") as f_out:
        for line in f_in:
            key = line.split("=")[0].strip()
            if "=" in line and not line.lstrip().startswith("%") and key in overrides:
                f_out.write(f"{key}= {overrides[key]}\n")
                done.add(key)
            else:
                f_out.write(line)
        for key, value in overrides.items():
            if key not in done:
                f_out.write(f"{key}= {value}\n")


def read_csv(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    header = [h.strip().strip('"') for h in rows[0]]
    return header, rows[1:]


def read_column(path, name):
    """Column of an SU2 ASCII restart, ordered by PointID."""
    header, rows = read_csv(path)
    ids = np.array([int(r[0]) for r in rows])
    values = np.array([float(r[header.index(name)]) for r in rows])
    out = np.empty_like(values)
    out[ids] = values
    return out


def write_beta(path, beta):
    """Replace the Beta_FIML column of an SU2 ASCII restart."""
    with open(path) as f:
        rows = list(csv.reader(f))
    col = [h.strip().strip('"') for h in rows[0]].index("Beta_FIML")
    for r in rows[1:]:
        r[col] = f"{beta[int(r[0])]:.15e}"
    with open(path, "w") as f:
        f.write(",".join(rows[0]) + "\n")
        for r in rows[1:]:
            f.write(",".join(r) + "\n")


def last_history(path):
    header, rows = read_csv(path)
    return {h: float(v) for h, v in zip(header, rows[-1])}


class FieldInversion:
    def __init__(self, args):
        self.args = args
        self.n_eval = 0
        self.direct_cfg = "inv_direct.cfg"
        self.adjoint_cfg = "inv_adjoint.cfg"
        write_config(args.config, self.direct_cfg, {"MATH_PROBLEM": "DIRECT", "RESTART_SOL": "YES"})
        write_config(args.config, self.adjoint_cfg,
                     {"MATH_PROBLEM": "DISCRETE_ADJOINT", "RESTART_SOL": "NO", "ITER": str(args.adjoint_iter),
                      "CONV_FIELD": "RMS_ADJ_PRESSURE", "CONV_RESIDUAL_MINVAL": str(args.adjoint_tol),
                      "OUTPUT_FILES": "RESTART_ASCII"})
        self.log = open("inversion_history.csv", "w")
        self.log.write("eval,J,CL,CD,beta_min,beta_max\n")

    def run(self, exe, cfg, log):
        cmd = ["mpirun", "-n", str(self.args.np), exe, cfg]
        with open(log, "w") as f:
            if subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT) != 0:
                sys.exit(f"{' '.join(cmd)} failed, see {log}")

    def evaluate(self, beta):
        write_beta(self.args.restart, beta)
        self.run(self.args.su2_cfd, self.direct_cfg, "log_direct.txt")
        hist = last_history("history.csv")
        self.run(self.args.su2_cfd_ad, self.adjoint_cfg, "log_adjoint.txt")
        grad = read_column(self.args.adjoint_restart, "Sensitivity_Beta_FIML")

        shutil.copy(self.args.restart, f"restart_eval_{self.n_eval:03d}.csv")
        self.log.write(f"{self.n_eval},{hist['ComboObj']:.10e},{hist['CL']:.8f},{hist['CD']:.8f},"
                       f"{beta.min():.6f},{beta.max():.6f}\n")
        self.log.flush()
        print(f"eval {self.n_eval}: J = {hist['ComboObj']:.6e}, CL = {hist['CL']:.6f}, |grad| = "
              f"{np.linalg.norm(grad):.3e}", flush=True)
        self.n_eval += 1
        return self.args.scale * hist["ComboObj"], self.args.scale * grad


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-f", "--config", default="config_inversion.cfg")
    p.add_argument("-n", "--np", type=int, default=8, help="MPI ranks")
    p.add_argument("--su2_cfd", default="SU2_CFD")
    p.add_argument("--su2_cfd_ad", default="SU2_CFD_AD")
    p.add_argument("--restart", default="restart_flow.csv", help="flow restart holding Beta_FIML")
    p.add_argument("--adjoint_restart", default="restart_adj.csv")
    p.add_argument("--adjoint_iter", type=int, default=5000)
    p.add_argument("--adjoint_tol", type=float, default=-11)
    p.add_argument("--maxiter", type=int, default=20)
    p.add_argument("--scale", type=float, default=1.0, help="objective scale seen by the optimizer")
    p.add_argument("--beta_bounds", type=float, nargs=2, default=(0.0, 2.0))
    args = p.parse_args()

    inv = FieldInversion(args)
    beta0 = read_column(args.restart, "Beta_FIML")
    result = minimize(inv.evaluate, beta0, jac=True, method="L-BFGS-B",
                      bounds=[tuple(args.beta_bounds)] * beta0.size,
                      options={"maxiter": args.maxiter, "ftol": 1e-12, "gtol": 1e-14})
    print(result.message)
    np.savetxt("beta_inverted.dat", result.x)


if __name__ == "__main__":
    main()
