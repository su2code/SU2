#!/usr/bin/env python3
"""FIML inversion with L-BFGS-B and gradients from the SU2 discrete adjoint.

--mode field:   the design variables are beta at every point (KIND_TRAIN_NN= NONE). beta lives in the
                Beta_FIML column of the ASCII flow restart, the gradient is Sensitivity_Beta_FIML.
--mode weights: the design variables are the network weights in nn_weights.dat (KIND_TRAIN_NN= WEIGHTS),
                the gradient is nn_weights_grad.dat.

Each evaluation runs the direct solver from the previous restart, then the discrete adjoint, which restarts
from the previous adjoint solution when there is one. The objective is ComboObj of the last history
line, normally a CUSTOM_OBJFUNC expression.
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


def write_column(path, name, values):
    """Replace a column of an SU2 ASCII restart."""
    with open(path) as f:
        rows = list(csv.reader(f))
    col = [h.strip().strip('"') for h in rows[0]].index(name)
    for r in rows[1:]:
        r[col] = f"{values[int(r[0])]:.15e}"
    with open(path, "w") as f:
        f.write(",".join(rows[0]) + "\n")
        for r in rows[1:]:
            f.write(",".join(r) + "\n")


def last_history(path):
    header, rows = read_csv(path)
    return {h: float(v) for h, v in zip(header, rows[-1])}


def read_weights(path):
    """Header, input names (if present) and scaling lines, and the weights of an nn_weights.dat file."""
    with open(path) as f:
        lines = f.read().splitlines()
    n_head = 3 if lines[1].startswith("Inputs:") else 2
    return lines[:n_head], np.array([float(v) for v in lines[n_head:] if v.strip()])


def write_weights(path, head, w):
    with open(path, "w") as f:
        f.write("\n".join(head) + "\n")
        f.writelines(f"{v:.15e}\n" for v in w)


class Inversion:
    def __init__(self, args):
        self.args = args
        self.n_eval = 0
        self.direct_cfg = "inv_direct.cfg"
        self.adjoint_cfg = "inv_adjoint.cfg"
        # The residual alone can already be converged after a change of beta, so ComboObj must settle too.
        write_config(args.config, self.direct_cfg, {
            "MATH_PROBLEM": "DIRECT", "RESTART_SOL": "YES", "ITER": str(args.direct_iter),
            "CONV_FIELD": "RMS_PRESSURE, COMBO", "CONV_STARTITER": "10",
            "CONV_CAUCHY_ELEMS": "100", "CONV_CAUCHY_EPS": str(args.cauchy_eps)})
        adjoint = {"MATH_PROBLEM": "DISCRETE_ADJOINT", "RESTART_SOL": "NO", "ITER": str(args.adjoint_iter),
                   "CONV_FIELD": "RMS_ADJ_PRESSURE", "CONV_RESIDUAL_MINVAL": str(args.adjoint_tol),
                   "OUTPUT_FILES": "RESTART_ASCII", "CONV_FILENAME": "history_adjoint",
                   "SCREEN_OUTPUT": "INNER_ITER, RMS_ADJ_PRESSURE, RMS_ADJ_NU_TILDE"}
        if args.adjoint_cfl:
            adjoint.update({"CFL_NUMBER": str(args.adjoint_cfl), "CFL_ADAPT": "NO"})
        self.adjoint = adjoint
        if args.mode == "weights":
            self.weight_head, _ = read_weights(args.weights)
        self.log = open(args.log, "w")
        self.log.write("eval,J," + ",".join(args.fields) + ",beta_min,beta_max,grad_norm,adjoint_rms_p\n")

    def run_adjoint(self):
        """Run the adjoint, from the previous adjoint solution with --warm_adjoint, and return the gradient
        and the final rms[A_P]. Stops if the adjoint residual did not fall."""
        overrides = dict(self.adjoint)
        if self.args.warm_adjoint and os.path.exists(self.args.adjoint_restart):
            overrides["RESTART_SOL"] = "YES"
        write_config(self.args.config, self.adjoint_cfg, overrides)
        self.run(self.args.su2_cfd_ad, self.adjoint_cfg, f"log_adjoint_{self.n_eval:03d}.txt")
        header, rows = read_csv("history_adjoint.csv")
        col = header.index("rms[A_P]")
        rms_start, rms_end = float(rows[0][col]), float(rows[-1][col])
        if self.args.mode == "field":
            grad = read_column(self.args.adjoint_restart, "Sensitivity_Beta_FIML")
        else:
            grad = np.loadtxt("nn_weights_grad.dat")
        if not (np.all(np.isfinite(grad)) and rms_end < rms_start):
            sys.exit(f"eval {self.n_eval}: the adjoint diverged (rms[A_P] {rms_start:.2f} -> {rms_end:.2f}), "
                     "lower --adjoint_cfl or the step of the optimizer")
        return grad, rms_end

    def run(self, exe, cfg, log):
        cmd = ["mpirun", "-n", str(self.args.np), exe, cfg]
        with open(log, "w") as f:
            if subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT) != 0:
                sys.exit(f"{' '.join(cmd)} failed, see {log}")

    def evaluate(self, x):
        if self.args.mode == "field":
            write_column(self.args.restart, "Beta_FIML", x)
        else:
            write_weights(self.args.weights, self.weight_head, x)

        self.run(self.args.su2_cfd, self.direct_cfg, f"log_direct_{self.n_eval:03d}.txt")
        hist = last_history("history.csv")
        if hist["Inner_Iter"] + 1 >= self.args.direct_iter:
            print(f"eval {self.n_eval}: the direct solver did not converge in {self.args.direct_iter} iterations",
                  flush=True)
        shutil.copy("history.csv", f"history_eval_{self.n_eval:03d}.csv")
        J = hist["ComboObj"]
        grad, adj_rms = self.run_adjoint()

        beta = read_column(self.args.restart, "Beta_FIML")
        np.save(f"x_eval_{self.n_eval:03d}.npy", x)
        terms = ",".join(f"{hist[k]:.10e}" for k in self.args.fields)
        self.log.write(f"{self.n_eval},{J:.10e},{terms},{beta.min():.6f},{beta.max():.6f},"
                       f"{np.linalg.norm(grad):.6e},{adj_rms:.3f}\n")
        self.log.flush()
        print(f"eval {self.n_eval}: J = {J:.6e}, " +
              ", ".join(f"{k} = {hist[k]:.6g}" for k in self.args.fields), flush=True)
        self.n_eval += 1
        return self.args.scale * J, self.args.scale * grad


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-f", "--config", default="config_inversion.cfg")
    p.add_argument("-n", "--np", type=int, default=8, help="MPI ranks")
    p.add_argument("--mode", choices=("field", "weights"), default="field")
    p.add_argument("--su2_cfd", default="SU2_CFD")
    p.add_argument("--su2_cfd_ad", default="SU2_CFD_AD")
    p.add_argument("--restart", default="restart_flow.csv", help="ASCII flow restart holding Beta_FIML")
    p.add_argument("--adjoint_restart", default="restart_adj_custom.csv")
    p.add_argument("--weights", default="nn_weights.dat")
    p.add_argument("--direct_iter", type=int, default=3000)
    p.add_argument("--cauchy_eps", type=float, default=1e-5,
                   help="mean relative change of ComboObj over 100 iterations at which the direct solver stops")
    p.add_argument("--adjoint_iter", type=int, default=5000)
    p.add_argument("--adjoint_tol", type=float, default=-11)
    p.add_argument("--adjoint_cfl", type=float, help="fixed CFL of the adjoint, which does not adapt it")
    p.add_argument("--warm_adjoint", action="store_true", help="start each adjoint from the previous one")
    p.add_argument("--maxiter", type=int, default=20)
    p.add_argument("--scale", type=float, default=1.0, help="objective scale seen by the optimizer")
    p.add_argument("--x_scale", type=float, default=1.0,
                   help="the optimizer works on x / x_scale, so its first step has length x_scale")
    p.add_argument("--beta_bounds", type=float, nargs=2, default=(0.0, 2.0))
    p.add_argument("--fields", nargs="*", default=["CL", "CD", "Beta_Reg"], help="history columns to log")
    p.add_argument("--log", default="inversion_history.csv")
    p.add_argument("--output", default="beta_inverted.dat")
    args = p.parse_args()

    inv = Inversion(args)
    if args.mode == "field":
        x0 = read_column(args.restart, "Beta_FIML")
        bounds = [tuple(args.beta_bounds)] * x0.size
    else:
        _, x0 = read_weights(args.weights)
        bounds = None
    s = args.x_scale
    if bounds:
        bounds = [(lo / s, hi / s) for lo, hi in bounds]

    def scaled(z):
        J, grad = inv.evaluate(s * z)
        return J, s * grad

    result = minimize(scaled, x0 / s, jac=True, method="L-BFGS-B", bounds=bounds,
                      options={"maxiter": args.maxiter, "ftol": 1e-12, "gtol": 1e-14})
    print(result.message)
    result.x = s * result.x

    # The restart must hold the flow of the optimum, not of the last line-search trial.
    last = np.load(f"x_eval_{inv.n_eval - 1:03d}.npy")
    if not np.array_equal(last, result.x):
        if args.mode == "field":
            write_column(args.restart, "Beta_FIML", result.x)
        else:
            write_weights(args.weights, inv.weight_head, result.x)
        inv.run(args.su2_cfd, inv.direct_cfg, "log_direct_final.txt")
    if args.mode == "field":
        np.savetxt(args.output, result.x)
    else:
        write_weights(args.weights, inv.weight_head, result.x)


if __name__ == "__main__":
    main()
