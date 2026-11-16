#!/usr/bin/env python3
"""Inverted beta and the eddy viscosity ratio to the baseline in the (x/D, r/D) plane.

    python3 plot_beta.py --inverted restart_inverted.csv --baseline restart_baseline.csv --targets morse_targets.csv

Writes morse_beta.png; the probe positions of the objective are marked.
"""

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm

import morse_data


def read_plane(path, names):
    """Columns of an ASCII restart on the (x/D, r/D) grid of the first azimuthal plane."""
    with open(path) as f:
        rows = list(csv.reader(f))
    head = [c.strip().strip('"') for c in rows[0]]
    a = np.array([[float(v) for v in r] for r in rows[1:] if r])
    col = {k: a[:, i] for i, k in enumerate(head)}
    th = np.arctan2(col["z"], col["y"])
    plane = np.isclose(th, th.min(), atol=1e-6)
    x = np.round(col["x"][plane] / morse_data.DIAMETER, 9)
    r = np.round(np.hypot(col["y"], col["z"])[plane] / morse_data.DIAMETER, 9)
    xs, rs = np.unique(x), np.unique(r)
    ix, ir = np.searchsorted(xs, x), np.searchsorted(rs, r)
    grid = {}
    for k in names:
        grid[k] = np.full((len(xs), len(rs)), np.nan)
        grid[k][ix, ir] = col[k][plane]
    return xs, rs, grid


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inverted", default="restart_inverted.csv")
    p.add_argument("--baseline", default="restart_baseline.csv")
    p.add_argument("--targets", default="morse_targets.csv")
    p.add_argument("--xmax", type=float, default=16.0, help="x/D range of the plot")
    p.add_argument("--rmax", type=float, default=3.0, help="r/D range of the plot")
    p.add_argument("--output", default="morse_beta.png")
    a = p.parse_args()

    xs, rs, inv = read_plane(a.inverted, ("Beta_FIML", "Eddy_Viscosity"))
    _, _, base = read_plane(a.baseline, ("Eddy_Viscosity",))
    tgt = np.genfromtxt(a.targets, delimiter=",", names=True, dtype=None, encoding=None)

    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    beta = inv["Beta_FIML"].T
    im = ax[0].pcolormesh(xs, rs, beta, cmap="RdBu_r", shading="nearest",
                          norm=TwoSlopeNorm(vcenter=1.0, vmin=min(np.nanmin(beta), 0.99), vmax=max(np.nanmax(beta), 1.01)))
    fig.colorbar(im, ax=ax[0], label=r"$\beta$")
    ratio = (inv["Eddy_Viscosity"] / base["Eddy_Viscosity"]).T
    lim = max(np.nanmax(ratio), 1.0 / np.nanmin(ratio))
    im = ax[1].pcolormesh(xs, rs, ratio, cmap="RdBu_r", shading="nearest", norm=LogNorm(vmin=1 / lim, vmax=lim))
    fig.colorbar(im, ax=ax[1], label=r"$\mu_t / \mu_{t,SA}$")

    for axis in ax:
        for group, marker in (("Ucl", "s"), ("U", "o")):
            sel = tgt["group"] == group
            axis.plot(tgt["xD"][sel], tgt["rD"][sel], marker, ms=2.5, mfc="none", mec="k", mew=0.5, ls="none")
        axis.set(xlim=(0, a.xmax), ylim=(0, a.rmax), ylabel="r/D")
        axis.set_aspect("equal")
    ax[0].set_title(r"inverted $\beta$, probes of the objective")
    ax[1].set_title("eddy viscosity, inverted / SA")
    ax[1].set_xlabel("x/D")
    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
