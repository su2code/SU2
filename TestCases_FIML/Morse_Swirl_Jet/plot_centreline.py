#!/usr/bin/env python3
"""Jet development: centreline velocity, maximum swirl velocity and half-width against x/D, SU2 and Morse (1980).

    python3 plot_centreline.py --data_dir <dir> --case SA=restart_flow.csv [--case ...]

Prints the errors at the measurement stations and writes morse_centreline.png.
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import morse_data


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--case", action="append", required=True, help="label=restart.csv, may be repeated")
    p.add_argument("--data_dir", default=".", help="directory with morse0.dat")
    p.add_argument("--output", default="morse_centreline.png")
    a = p.parse_args()

    exp = morse_data.read_development(a.data_dir)
    has_w, has_r = np.isfinite(exp["W_max"]), np.isfinite(exp["r05"])
    cases = morse_data.parse_cases(a.case)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(exp["x"], exp["U_cl"] / morse_data.U_EXIT, "ko", ms=4, mfc="none", label="Morse (1980)")
    ax[1].plot(exp["x"][has_w], exp["W_max"][has_w] / morse_data.U_EXIT, "ko", ms=5, mfc="none")
    ax[2].plot(exp["x"][has_r], exp["r05"][has_r], "ko", ms=5, mfc="none")

    print(f"{'case':>12} {'quantity':>14} " + " ".join(f"{f'x/D={x:g}':>9}" for x in (1, 2, 4, 6, 10, 15)))
    for label, (xs, rs, g) in cases.items():
        u_cl, w_max, r05 = morse_data.development(xs, rs, g)
        ax[0].plot(xs, u_cl / morse_data.U_EXIT, label=label)
        ax[1].plot(xs, w_max / morse_data.U_EXIT)
        ax[2].plot(xs, r05)

        def err(xe, ye, ys, relative):
            sim = np.interp(xe, xs, ys)
            return (sim / ye - 1) * 100 if relative else sim - ye

        rows = {"U_cl [%]": [err(x, np.interp(x, exp["x"], exp["U_cl"]), u_cl, True) for x in (1, 2, 4, 6, 10, 15)],
                "W_max [m/s]": [err(x, np.interp(x, exp["x"][has_w], exp["W_max"][has_w]), w_max, False)
                                if x <= 6 else np.nan for x in (1, 2, 4, 6, 10, 15)],
                "r_1/2 [%]": [err(x, np.interp(x, exp["x"][has_r], exp["r05"][has_r]), r05, True)
                              if x <= 10 else np.nan for x in (1, 2, 4, 6, 10, 15)]}
        for name, vals in rows.items():
            print(f"{label:>12} {name:>14} " + " ".join(f"{v:+9.1f}" if np.isfinite(v) else f"{'':>9}" for v in vals))

    ax[0].set(xlabel="x/D", ylabel="$U_{cl}/U_{exit}$", xlim=(0, 15), ylim=(0, 1.05))
    ax[1].set(xlabel="x/D", ylabel="$W_{max}/U_{exit}$", xlim=(0, 10), ylim=(0, None))
    ax[2].set(xlabel="x/D", ylabel="$r_{1/2}/D$", xlim=(0, 10), ylim=(0, 1.5 * np.nanmax(exp["r05"])))
    ax[0].legend()
    for axis in ax:
        axis.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(a.output, dpi=120)
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
