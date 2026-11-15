#!/usr/bin/env python3
"""Objective of the Morse field inversion: probes at mesh nodes with the measured velocities as targets.

Writes config_inversion.cfg (morse_sa.cfg with SA_FIML, the probes and CUSTOM_OBJFUNC), config_predict.cfg (the same
with beta from the network in nn_weights.dat, inputs --inputs) and morse_targets.csv.
The targets are the centreline velocity for x/D in --centreline and the axial and swirl velocity profiles at
--stations. Each target sits at a mesh node of the first azimuthal plane; the measurement is interpolated to the
radius of that node. The objective is

    J = sum ((U - U_exp) / sigma)^2 + sum ((W - W_exp) / sigma)^2 + lambda * BETA_REG

split into the history fields Misfit_Ucl (centreline), Misfit_U and Misfit_W (profiles).
"""

import argparse
import os
import sys

import numpy as np

import morse_data

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fiml_inversion import write_config  # noqa: E402


TERMS_PER_FUNCTION = 30


def read_mesh_points(path):
    """Coordinates of the points of an SU2 mesh."""
    with open(path) as f:
        for line in f:
            if line.startswith("NPOIN"):
                npoin = int(line.split("=")[1].split()[0])
                return np.loadtxt(f, max_rows=npoin, usecols=(0, 1, 2))
    sys.exit(f"{path} has no NPOIN section")


def nearest(values, target):
    return int(np.argmin(np.abs(values - target)))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data_dir", default=".", help="directory with morse0.dat ... morse6.dat")
    p.add_argument("--config", default="morse_sa.cfg", help="SA configuration the inversion starts from")
    p.add_argument("--mesh", default="morse_wedge.su2")
    p.add_argument("--output", default="config_inversion.cfg")
    p.add_argument("--targets", default="morse_targets.csv")
    p.add_argument("--stations", default="2,3,4,5,6", help="profile stations (morse<n>.dat), comma separated")
    p.add_argument("--centreline", type=float, nargs=2, default=(1.0, 15.0), help="x/D range of U_cl targets")
    p.add_argument("--sigma", type=float, default=0.02 * morse_data.U_EXIT, help="velocity uncertainty, m/s")
    p.add_argument("--beta_reg", type=float, default=0.05, help="lambda, weight of BETA_REG")
    p.add_argument("--no_swirl", action="store_true", help="leave the swirl velocity out of the objective")
    p.add_argument("--predict", default="config_predict.cfg")
    p.add_argument("--inputs", default="CHI,VISC_RATIO,STRAIN_VORT,Q_CRIT",
                   help="FIML_INPUTS of the network, comma separated")
    a = p.parse_args()

    D = morse_data.DIAMETER
    xyz = read_mesh_points(a.mesh)
    theta = np.arctan2(xyz[:, 2], xyz[:, 1])
    plane = np.flatnonzero(np.isclose(theta, theta.min(), atol=1e-6))
    xs = np.unique(np.round(xyz[plane, 0] / D, 9))
    th = theta.min()

    def node(x_over_d, r_over_d):
        """Index of the plane point nearest to (x/D, r/D)."""
        px, pr = xyz[plane, 0] / D, np.hypot(xyz[plane, 1], xyz[plane, 2]) / D
        column = np.isclose(px, xs[nearest(xs, x_over_d)])
        cand = plane[column]
        return cand[nearest(pr[column], r_over_d)]

    targets = []  # (group, name, point, quantity, value)
    exp = morse_data.read_development(a.data_dir)
    in_range = (exp["x"] >= a.centreline[0]) & (exp["x"] <= a.centreline[1])
    for x in exp["x"][in_range]:
        i = node(x, 0.0)
        u = np.interp(xyz[i, 0] / D, exp["x"], exp["U_cl"])
        targets.append(("Ucl", f"ucl_{len(targets):03d}", i, "U", u))

    has_w = np.isfinite(exp["W_max"])
    for st in (int(s) for s in a.stations.split(",")):
        x_st = morse_data.STATIONS[st]
        prof = morse_data.read_profile(a.data_dir, st)
        order = np.argsort(prof["r"])
        r_m, u_m, w_m = prof["r"][order], prof["U"][order], prof["W"][order]
        u_cl = np.interp(x_st, exp["x"], exp["U_cl"])
        w_max = np.interp(x_st, exp["x"][has_w], exp["W_max"][has_w])
        used = set()
        for r in r_m[r_m > 0]:
            i = node(x_st, r)
            if i in used:
                continue
            used.add(i)
            r_node = np.hypot(xyz[i, 1], xyz[i, 2]) / D
            targets.append(("U", f"u{st}_{len(used):02d}", i, "U", u_cl * np.interp(r_node, r_m, u_m)))
            if not a.no_swirl:
                targets.append(("W", f"w{st}_{len(used):02d}", i, "W", w_max * np.interp(r_node, r_m, w_m)))

    # Swirl velocity in the plane theta = th: W = -v sin(th) + w cos(th).
    expr = {"U": "VELOCITY_X", "W": f"{-np.sin(th):.12f} * VELOCITY_Y + {np.cos(th):.12f} * VELOCITY_Z"}
    probes = [f"{name} : Probe{{{expr[q]}}}[{xyz[i, 0]:.12e}, {xyz[i, 1]:.12e}, {xyz[i, 2]:.12e}]"
              for _, name, i, q, _ in targets]
    groups = {}
    for g, name, _, _, value in targets:
        groups.setdefault(g, []).append(f"pow({name} - {value:.6f}, 2)")
    # An SU2 expression holds at most 255 nodes, so each misfit is a sum of partial sums.
    misfits = []
    for g, terms in groups.items():
        parts = [f"part_{g}_{k}" for k in range(0, len(terms), TERMS_PER_FUNCTION)]
        misfits += [f"{part} : Function{{{' + '.join(terms[k:k + TERMS_PER_FUNCTION])}}}"
                    for part, k in zip(parts, range(0, len(terms), TERMS_PER_FUNCTION))]
        misfits.append(f"Misfit_{g} : Function{{({' + '.join(parts)}) / {a.sigma**2:.8f}}}")

    objective = " + ".join(f"Misfit_{g}" for g in groups) + f" + {a.beta_reg:g} * BETA_REG"
    overrides = {
        "KIND_TURB_MODEL": "SA_FIML",
        "KIND_TRAIN_NN": "NONE",
        "TRAIN_NN": "NO",
        "FIML_SHIELD": "NONE",
        "RESTART_SOL": "YES",
        "OBJECTIVE_FUNCTION": "CUSTOM_OBJFUNC",
        "OBJECTIVE_WEIGHT": "1.0",
        "CUSTOM_OUTPUTS": "'" + "; ".join(probes + misfits) + "'",
        "CUSTOM_OBJFUNC": f"'{objective}'",
        "READ_BINARY_RESTART": "NO",
        "SOLUTION_FILENAME": "restart_flow",
        "RESTART_FILENAME": "restart_flow",
        "SOLUTION_ADJ_FILENAME": "restart_adj",
        "RESTART_ADJ_FILENAME": "restart_adj",
        "CONV_FILENAME": "history",
        "SCREEN_OUTPUT": "INNER_ITER, RMS_PRESSURE, RMS_NU_TILDE, LINSOL_ITER, AVG_CFL, "
                         + ", ".join(f"Misfit_{g}" for g in groups) + ", COMBO",
        "HISTORY_OUTPUT": "ITER, RMS_RES, CFL_NUMBER, CUSTOM, FIML, COMBO",
        "VOLUME_OUTPUT": "COORDINATES, SOLUTION, PRIMITIVE, FIML",
        "OUTPUT_FILES": "RESTART_ASCII",
        "WRT_RESTART_COMPACT": "NO",
        "OUTPUT_WRT_FREQ": "100000",
    }
    write_config(a.config, a.output, overrides)
    write_config(a.output, a.predict, {"KIND_TRAIN_NN": "WEIGHTS", "ITER_START_NN": "1",
                                       "FIML_INPUTS": f"({a.inputs.replace(',', ', ')})"})

    with open(a.targets, "w") as f:
        f.write("name,group,quantity,x/D,r/D,target\n")
        for g, name, i, q, value in targets:
            f.write(f"{name},{g},{q},{xyz[i, 0] / D:.6f},{np.hypot(xyz[i, 1], xyz[i, 2]) / D:.6f},{value:.6f}\n")
    counts = ", ".join(f"{len(t)} {g}" for g, t in groups.items())
    print(f"{a.output}: {len(targets)} probes ({counts}), sigma = {a.sigma:.3f} m/s, lambda = {a.beta_reg:g}")
    print(f"{a.predict}: network inputs {a.inputs}")
    print(f"{a.targets}: probe positions and targets")


if __name__ == "__main__":
    main()
