"""Readers and runners shared by the steps of the S809 field inversion."""

import argparse
import csv
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from fiml_inversion import write_config  # noqa: E402

BASE_CONFIG = "config_inversion.cfg"


def read_airfoil_marker(mesh, marker="airfoil"):
    """Global ids, coordinates and dual lengths (half of the adjacent edges) of the wall nodes."""
    with open(mesh) as f:
        lines = f.read().splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("NPOIN"))
    npoin = int(lines[i].split("=")[1].split()[0])
    coords = np.array([l.split()[:2] for l in lines[i + 1:i + 1 + npoin]], dtype=float)
    i = next(k for k, l in enumerate(lines) if l.replace(" ", "") == f"MARKER_TAG={marker}")
    nedge = int(lines[i + 1].split("=")[1])
    edges = np.array([l.split()[1:3] for l in lines[i + 2:i + 2 + nedge]], dtype=int)

    ids = np.unique(edges)
    pos = {p: k for k, p in enumerate(ids)}
    length = np.zeros(len(ids))
    for a, b in edges:
        half = 0.5 * np.linalg.norm(coords[a] - coords[b])
        length[pos[a]] += half
        length[pos[b]] += half
    return ids, coords[ids], length


def upper_surface(xy, airfoil):
    """True for points above the camber line of a Selig coordinate file; xy in chords."""
    d = np.loadtxt(airfoil, skiprows=1)
    i_le = np.argmin(d[:, 0])
    upper, lower = d[:i_le + 1][::-1], d[i_le:]
    camber = 0.5 * (np.interp(xy[:, 0], upper[:, 0], upper[:, 1]) + np.interp(xy[:, 0], lower[:, 0], lower[:, 1]))
    return xy[:, 1] >= camber


def read_table(path):
    """Columns of an SU2 CSV file (restart, surface or history) as a dict of arrays."""
    with open(path) as f:
        rows = list(csv.reader(f))
    header = [h.strip().strip('"') for h in rows[0]]
    data = np.array([[float(v) for v in r] for r in rows[1:] if r], dtype=float)
    return {h: data[:, k] for k, h in enumerate(header)}


def read_surface(path, ids, field="Pressure_Coefficient"):
    """A field of a surface CSV file at the given global ids."""
    t = read_table(path)
    lookup = dict(zip(t["PointID"].astype(int), t[field]))
    return np.array([lookup[i] for i in ids])


def read_restart(path):
    """Columns of an ASCII restart, ordered by PointID."""
    t = read_table(path)
    order = np.argsort(t["PointID"])
    return {k: v[order] for k, v in t.items()}


def last_history(path):
    """Last line of an SU2 history file as a dict."""
    return {k: v[-1] for k, v in read_table(path).items()}


def run_parser(description):
    """Argument parser with the options to run SU2."""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("-n", "--np", type=int, default=8, help="MPI ranks")
    p.add_argument("--su2_cfd", default="SU2_CFD")
    p.add_argument("--su2_cfd_ad", default="SU2_CFD_AD")
    return p


def run_su2(args, name, overrides, exe=None):
    """Write name.cfg from the base configuration with the overrides and run it, logging to log_name.txt."""
    cfg = f"{name}.cfg"
    write_config(BASE_CONFIG, cfg, overrides)
    cmd = ["mpirun", "-n", str(args.np), exe or args.su2_cfd, cfg]
    with open(f"log_{name}.txt", "w") as f:
        if subprocess.call(cmd, stdout=f, stderr=subprocess.STDOUT) != 0:
            sys.exit(f"{' '.join(cmd)} failed, see log_{name}.txt")


def keep(src, dst):
    """Copy an output file to the name a later step reads."""
    if not os.path.exists(src):
        sys.exit(f"{src} was not written")
    shutil.copy(src, dst)


def require(*paths):
    for path in paths:
        if not os.path.exists(path):
            sys.exit(f"{path} is missing, run the previous step first")


def report(label, hist, max_iter=None):
    """Print the lift, drag and objective terms of the last history line."""
    it = int(hist["Inner_Iter"])
    note = "" if max_iter is None else (" (stopped at ITER)" if it + 1 >= max_iter else " (converged)")
    print(f"{label}: {it + 1} iterations{note}, CL = {hist['CL']:.4f}, CD = {hist['CD']:.4f}, "
          f"Cp_Diff = {hist['Cp_Diff']:.5e}, Beta_Reg = {hist['Beta_Reg']:.4e}, J = {hist['ComboObj']:.4f}")
