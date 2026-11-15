#!/usr/bin/env python3
"""Shear stresses u'v' and v'w' divided by U_cl^2 against r/D at the measurement stations, SU2 and Morse (1980).

    python3 plot_stresses.py --data_dir <dir> --case SA=restart_flow.csv [--case ...]

SU2 stresses follow from the eddy viscosity (Boussinesq). Writes morse_stresses.png.
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import morse_data

ROWS = (("uv", "$\\overline{u'v'}/U_{cl}^2$"), ("vw", "$\\overline{v'w'}/U_{cl}^2$"))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", action="append", required=True, help="label=restart.csv, may be repeated")
    p.add_argument("--data_dir", default=".", help="directory with morse1.dat ... morse6.dat")
    p.add_argument("--stations", type=int, nargs="+", default=[2, 3, 4, 5, 6], help="station numbers 1-6")
    p.add_argument("--output", default="morse_stresses.png")
    a = p.parse_args()

    cases = morse_data.parse_cases(a.case)
    stresses = {label: morse_data.shear_stresses(xs, rs, g) for label, (xs, rs, g) in cases.items()}
    n = len(a.stations)
    fig, ax = plt.subplots(len(ROWS), n, figsize=(3.2 * n, 5.8), sharey="row", squeeze=False)
    for col, st in enumerate(a.stations):
        x_st = morse_data.STATIONS[st]
        exp = morse_data.read_profile(a.data_dir, st)
        r_lim = 1.25 * exp["r"].max()
        for row, (key, ylabel) in enumerate(ROWS):
            panel = ax[row, col]
            panel.plot(exp["r"], exp[key], "ko", ms=4, mfc="none", label="Morse (1980)")
            for label, (xs, rs, g) in cases.items():
                i = int(np.argmin(abs(xs - x_st)))
                uv, vw = stresses[label]
                panel.plot(rs, (uv if key == "uv" else vw)[i] / g["U"][i, 0] ** 2, label=label)
            panel.set(xlim=(0, r_lim))
            panel.grid(alpha=0.3)
            if row == 0:
                panel.set_title(f"x/D = {x_st:g}")
            if row == len(ROWS) - 1:
                panel.set_xlabel("r/D")
            if col == 0:
                panel.set_ylabel(ylabel)
    ax[0, 0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.output, dpi=120)
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
