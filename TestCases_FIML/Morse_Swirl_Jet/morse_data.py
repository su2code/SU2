"""Readers for the Morse swirling-jet measurements and for SU2 solutions on the wedge mesh.

Used by plot_centreline.py, plot_profiles.py and plot_stresses.py.
"""

import csv
import os

import numpy as np

DIAMETER = 0.0246
U_EXIT = 34.46
STATIONS = {1: 0.5, 2: 1.0, 3: 2.0, 4: 4.0, 5: 6.0, 6: 10.0}
PROFILE_COLUMNS = ("r/x", "r/r05", "U", "V", "W", "uu", "vv", "ww", "uv", "uw", "vw")


def read_development(data_dir):
    """Centreline velocity, maximum swirl velocity [m/s] and half-width r_1/2 / D against x/D, from morse0.dat."""
    rows = []
    with open(os.path.join(data_dir, "morse0.dat")) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith("#"):
                vals = [float(v) for v in line.split()]
                rows.append(vals + [np.nan] * (5 - len(vals)))
    m = np.array(rows)
    return {"x": m[:, 0], "U_cl": m[:, 1], "W_max": m[:, 3], "r05": m[:, 4]}


def read_profile(data_dir, station):
    """Profiles at one station: r/D, U/U_cl, V/U_cl, W/W_max and the Reynolds stresses divided by U_cl^2."""
    d = np.loadtxt(os.path.join(data_dir, f"morse{station}.dat"), comments="#")
    out = dict(zip(PROFILE_COLUMNS, d.T))
    out["r"] = out["r/x"] * STATIONS[station]
    return out


def read_solution(path, diameter=DIAMETER):
    """Solution on the (x/D, r/D) grid of the first azimuthal plane of the wedge.

    Returns xs, rs and a dict of 2D arrays: axial (U), radial (V) and swirl (W) velocity and kinematic eddy viscosity.
    """
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [c.strip().strip('"') for c in rows[0]]
    a = np.array([[float(v) for v in r] for r in rows[1:] if r])
    col = {k: a[:, i] for i, k in enumerate(head)}
    y, z = col["y"], col["z"]
    th = np.arctan2(z, y)
    plane = np.isclose(th, th.min(), atol=1e-6)
    x, r, t = col["x"][plane] / diameter, np.hypot(y, z)[plane] / diameter, th[plane]
    u, v, w = col["Velocity_x"][plane], col["Velocity_y"][plane], col["Velocity_z"][plane]
    fields = {"U": u, "V": v * np.cos(t) + w * np.sin(t), "W": -v * np.sin(t) + w * np.cos(t),
              "nu_t": col["Eddy_Viscosity"][plane] / col["Density"][plane]}
    xs, rs = np.unique(np.round(x, 9)), np.unique(np.round(r, 9))
    ix, ir = np.searchsorted(xs, np.round(x, 9)), np.searchsorted(rs, np.round(r, 9))
    grid = {}
    for k, val in fields.items():
        grid[k] = np.full((len(xs), len(rs)), np.nan)
        grid[k][ix, ir] = val
    if np.nansum(grid["W"]) < 0:
        grid["W"] = -grid["W"]
    return xs, rs, grid


def development(xs, rs, g):
    """Centreline velocity, maximum swirl velocity and half-width r_1/2 / D along the jet."""
    u_cl = g["U"][:, 0]
    w_max = np.abs(g["W"]).max(axis=1)
    r05 = np.array([np.interp(0.5 * u[0], u[::-1], rs[::-1]) for u in g["U"]])
    return u_cl, w_max, r05


def shear_stresses(xs, rs, g, diameter=DIAMETER):
    """Boussinesq shear stresses u'v' = -nu_t (dU/dr + dV/dx) and v'w' = -nu_t r d(W/r)/dr [m^2/s^2]."""
    x, r = xs * diameter, rs * diameter
    dU_dr = np.gradient(g["U"], r, axis=1)
    dV_dx = np.gradient(g["V"], x, axis=0)
    dWr_dr = np.gradient(g["W"] / r[None, :], r, axis=1)
    return -g["nu_t"] * (dU_dr + dV_dx), -g["nu_t"] * r[None, :] * dWr_dr


def parse_cases(items):
    """Label and solution of each --case label=restart.csv argument."""
    cases = {}
    for item in items:
        label, path = item.rsplit("=", 1) if "=" in item else (os.path.basename(item), item)
        cases[label] = read_solution(path)
    return cases
