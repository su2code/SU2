#!/usr/bin/env python3
"""Postprocessing of the Morse (1980) swirling jet: SU2 against the measurements.

Samples the volume solution on lines in the plane z = 0, writes the line data to csv and plots the jet
development, the velocity profiles and the shear stresses.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# ---------------------------------------------------------------- user input
CASES = {"SA": "flow.vtu"}          # label: SU2 volume solution, several cases are plotted together
DATA_DIR = "test1"                  # directory with morse0.dat ... morse6.dat
CSV_DIR = "lines"                   # output directory of the line data
OUT_DEVELOPMENT = "morse_centreline.png"
OUT_PROFILES = "morse_profiles.png"
OUT_STRESSES = "morse_stresses.png"

DIAMETER = 0.0246                   # nozzle diameter [m]
U_EXIT = 34.46                      # nozzle exit velocity [m/s]
R_AXIS = 0.02                       # r/D of the slip wall that replaces the axis
R_MAX = 6.0                         # r/D end of the radial lines
N_RADIAL = 600                      # points on a radial line
X_DEV = (0.5, 15.0)                 # x/D range of the jet development
N_DEV = 300                         # radial lines for the jet development
PLOT_STATIONS = (2, 3, 4, 5, 6)     # measurement stations in the profile figures
ERROR_STATIONS = (1, 2, 4, 6, 10, 15)
# ----------------------------------------------------------------------------

STATIONS = {1: 0.5, 2: 1.0, 3: 2.0, 4: 4.0, 5: 6.0, 6: 10.0}
PROFILE_COLUMNS = ("r/x", "r/r05", "U", "V", "W", "uu", "vv", "ww", "uv", "uw", "vw")
LINE_COLUMNS = ("x/D", "r/D", "U", "V", "W", "nu_t", "uv", "vw")


def read_development():
    """Measured U_cl, W_max [m/s] and r_1/2 / D against x/D."""
    rows = []
    with open(os.path.join(DATA_DIR, "morse0.dat")) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith("#"):
                vals = [float(v) for v in line.split()]
                rows.append(vals + [np.nan] * (5 - len(vals)))
    m = np.array(rows)
    return {"x": m[:, 0], "U_cl": m[:, 1], "W_max": m[:, 3], "r05": m[:, 4]}


def read_profile(station):
    """Measured profiles at one station, velocities and stresses normalised as in the data file."""
    d = np.loadtxt(os.path.join(DATA_DIR, f"morse{station}.dat"), comments="#")
    out = dict(zip(PROFILE_COLUMNS, d.T))
    out["r"] = out["r/x"] * STATIONS[station]
    return out


def load_case(path):
    """Volume solution with the velocity gradient and the swirl direction for which W > 0."""
    mesh = pv.read(path).compute_derivative(scalars="Velocity", gradient="grad_vel")
    sign = 1.0 if np.mean(mesh.point_data["Velocity"][:, 2]) >= 0 else -1.0
    return mesh, sign


def sample_points(mesh, sign, x, r):
    """Fields at (x/D, r/D) in the plane z = 0, where y is radial and z is azimuthal.

    Shear stresses follow from the eddy viscosity: u'v' = -nu_t (dU/dr + dV/dx), v'w' = -nu_t (dW/dr - W/r).
    """
    pts = np.column_stack([x * DIAMETER, r * DIAMETER, np.zeros_like(x)])
    s = pv.PolyData(pts).sample(mesh)
    vel, g = s.point_data["Velocity"], s.point_data["grad_vel"]
    nu_t = s.point_data["Eddy_Viscosity"] / s.point_data["Density"]
    w = sign * vel[:, 2]
    out = {"x/D": x, "r/D": r, "U": vel[:, 0], "V": vel[:, 1], "W": w, "nu_t": nu_t,
           "uv": -nu_t * (g[:, 1] + g[:, 3]), "vw": -nu_t * (sign * g[:, 7] - w / (r * DIAMETER))}
    valid = s.point_data["vtkValidPointMask"].astype(bool)
    for k in LINE_COLUMNS[2:]:
        out[k] = np.where(valid, out[k], np.nan)
    return out


def radial_lines(mesh, sign, xs):
    """Radial lines at each x/D in xs, as 2D arrays of shape (len(xs), N_RADIAL)."""
    r = np.linspace(R_AXIS * 1.001, R_MAX, N_RADIAL)
    X, R = np.meshgrid(xs, r, indexing="ij")
    s = sample_points(mesh, sign, X.ravel(), R.ravel())
    return {k: v.reshape(X.shape) for k, v in s.items()}


def development(lines):
    """Centreline velocity, maximum swirl velocity and half-width r_1/2 / D along the jet."""
    u, r = lines["U"], lines["r/D"]
    u_cl = u[:, 0]
    w_max = np.nanmax(np.abs(lines["W"]), axis=1)
    r05 = np.full(len(u_cl), np.nan)
    for i in range(len(u_cl)):
        j = np.argmax(u[i] < 0.5 * u_cl[i])
        if j > 0:
            r05[i] = np.interp(0.5 * u_cl[i], [u[i, j], u[i, j - 1]], [r[i, j], r[i, j - 1]])
    return {"x/D": lines["x/D"][:, 0], "U_cl": u_cl, "W_max": w_max, "r05": r05}


def write_csv(path, data, columns):
    np.savetxt(path, np.column_stack([np.ravel(data[c]) for c in columns]), delimiter=",",
               header=",".join(columns), comments="", fmt="%.8e")


def process(label, path):
    """Sample one case, write its line data to csv and return the development and station profiles."""
    mesh, sign = load_case(path)
    dev = development(radial_lines(mesh, sign, np.linspace(*X_DEV, N_DEV)))
    xs = np.array([STATIONS[s] for s in sorted(STATIONS)])
    st = radial_lines(mesh, sign, xs)
    stations = {s: {k: v[i] for k, v in st.items()} for i, s in enumerate(sorted(STATIONS))}

    write_csv(os.path.join(CSV_DIR, f"{label}_development.csv"), dev, ("x/D", "U_cl", "W_max", "r05"))
    for s, prof in stations.items():
        write_csv(os.path.join(CSV_DIR, f"{label}_station{s}_xD{STATIONS[s]:g}.csv"), prof, LINE_COLUMNS)
    return dev, stations


def print_errors(label, dev, exp):
    """Errors at ERROR_STATIONS: U_cl and r_1/2 in percent, W_max in m/s."""
    has_w, has_r = np.isfinite(exp["W_max"]), np.isfinite(exp["r05"])

    def err(x, xe, ye, ys, relative):
        sim, ref = np.interp(x, dev["x/D"], ys), np.interp(x, xe, ye)
        return (sim / ref - 1) * 100 if relative else sim - ref

    rows = {"U_cl [%]": [err(x, exp["x"], exp["U_cl"], dev["U_cl"], True) for x in ERROR_STATIONS],
            "W_max [m/s]": [err(x, exp["x"][has_w], exp["W_max"][has_w], dev["W_max"], False)
                            if x <= exp["x"][has_w].max() else np.nan for x in ERROR_STATIONS],
            "r_1/2 [%]": [err(x, exp["x"][has_r], exp["r05"][has_r], dev["r05"], True)
                          if x <= exp["x"][has_r].max() else np.nan for x in ERROR_STATIONS]}
    for name, vals in rows.items():
        print(f"{label:>12} {name:>14} " + " ".join(f"{v:+9.1f}" if np.isfinite(v) else f"{'':>9}" for v in vals))


def plot_development(results, exp):
    has_w, has_r = np.isfinite(exp["W_max"]), np.isfinite(exp["r05"])
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(exp["x"], exp["U_cl"] / U_EXIT, "ko", ms=4, mfc="none", label="Morse (1980)")
    ax[1].plot(exp["x"][has_w], exp["W_max"][has_w] / U_EXIT, "ko", ms=5, mfc="none")
    ax[2].plot(exp["x"][has_r], exp["r05"][has_r], "ko", ms=5, mfc="none")
    for label, (dev, _) in results.items():
        ax[0].plot(dev["x/D"], dev["U_cl"] / U_EXIT, label=label)
        ax[1].plot(dev["x/D"], dev["W_max"] / U_EXIT)
        ax[2].plot(dev["x/D"], dev["r05"])
    ax[0].set(xlabel="x/D", ylabel="$U_{cl}/U_{exit}$", xlim=(0, 15), ylim=(0, 1.05))
    ax[1].set(xlabel="x/D", ylabel="$W_{max}/U_{exit}$", xlim=(0, 10), ylim=(0, None))
    ax[2].set(xlabel="x/D", ylabel="$r_{1/2}/D$", xlim=(0, 10), ylim=(0, 1.5 * np.nanmax(exp["r05"])))
    ax[0].legend()
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DEVELOPMENT, dpi=120)
    print(f"wrote {OUT_DEVELOPMENT}")


def plot_stations(results, rows, output):
    """Profiles against r/D at PLOT_STATIONS; rows holds (key, ylabel, normalisation)."""
    n = len(PLOT_STATIONS)
    fig, ax = plt.subplots(len(rows), n, figsize=(3.2 * n, 2.9 * len(rows)), sharey="row", squeeze=False)
    for col, st in enumerate(PLOT_STATIONS):
        exp = read_profile(st)
        for row, (key, ylabel, norm) in enumerate(rows):
            panel = ax[row, col]
            panel.plot(exp["r"], exp[key], "ko", ms=4, mfc="none", label="Morse (1980)")
            for label, (_, stations) in results.items():
                p = stations[st]
                panel.plot(p["r/D"], p[key] / norm(p), label=label)
            panel.set(xlim=(0, 1.25 * exp["r"].max()))
            panel.grid(alpha=0.3)
            if row == 0:
                panel.set_title(f"x/D = {STATIONS[st]:g}")
            if row == len(rows) - 1:
                panel.set_xlabel("r/D")
            if col == 0:
                panel.set_ylabel(ylabel)
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=120)
    print(f"wrote {output}")


def main():
    os.makedirs(CSV_DIR, exist_ok=True)
    exp = read_development()
    results = {label: process(label, path) for label, path in CASES.items()}

    print(f"{'case':>12} {'quantity':>14} " + " ".join(f"{f'x/D={x:g}':>9}" for x in ERROR_STATIONS))
    for label, (dev, _) in results.items():
        print_errors(label, dev, exp)

    u_cl = lambda p: p["U"][0]
    u_cl2 = lambda p: p["U"][0] ** 2
    w_max = lambda p: np.nanmax(np.abs(p["W"]))
    plot_development(results, exp)
    plot_stations(results, (("U", "$U/U_{cl}$", u_cl), ("V", "$V/U_{cl}$", u_cl), ("W", "$W/W_{max}$", w_max)),
                  OUT_PROFILES)
    plot_stations(results, (("uv", "$\\overline{u'v'}/U_{cl}^2$", u_cl2),
                            ("vw", "$\\overline{v'w'}/U_{cl}^2$", u_cl2)), OUT_STRESSES)
    print(f"wrote line data to {CSV_DIR}/")


if __name__ == "__main__":
    main()
