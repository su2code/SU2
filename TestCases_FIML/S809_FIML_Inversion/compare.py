#!/usr/bin/env python3
"""Table of lift, drag and objective terms for every step that has run, and a plot of the upper-surface Cp."""

import os

import numpy as np

from fiml_case import last_history, read_airfoil_marker, read_surface, upper_surface

LABELS = ("baseline", "inverted", "predicted", "direct")


def main():
    print(f"{'':10s} {'CL':>8s} {'CD':>8s} {'Cp_Diff':>11s} {'Beta_Reg':>10s} {'J':>10s}")
    for label in LABELS:
        if os.path.exists(f"history_{label}.csv"):
            h = last_history(f"history_{label}.csv")
            print(f"{label:10s} {h['CL']:8.4f} {h['CD']:8.4f} {h['Cp_Diff']:11.4e} {h['Beta_Reg']:10.2f} "
                  f"{h['ComboObj']:10.3f}")
    print("experiment     1.0546   (target)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    ids, xy, _ = read_airfoil_marker("../S809_FIML_Forward/s809_cmesh.su2")
    up = upper_surface(xy / 0.6, "../S809_FIML_Forward/s809.dat")
    order = np.argsort(xy[up, 0])
    exp = np.loadtxt("s809_cp_upper_exp.dat")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(exp[:, 0], exp[:, 1], "ko", ms=4, label="experiment (upper)")
    for label in LABELS:
        if os.path.exists(f"surface_{label}.csv"):
            cp = read_surface(f"surface_{label}.csv", ids)
            ax.plot(xy[up, 0][order] / 0.6, cp[up][order], label=label)
    ax.invert_yaxis()
    ax.set_xlabel("x/c")
    ax.set_ylabel("Cp, upper surface")
    ax.legend()
    fig.savefig("cp_upper.png", dpi=120, bbox_inches="tight")
    print("wrote cp_upper.png")


if __name__ == "__main__":
    main()
