#!/usr/bin/env python3
"""Step 1: interpolate the measured upper-surface Cp onto the wall nodes and write TargetCp.dat.

Only upper-surface nodes between the first and the last pressure tap get a target, so the lower surface
does not contribute to INVERSE_DESIGN_PRESSURE. With --check, the script also prints the Cp misfit of a
surface CSV file, sum of (Cp - Cp_target)^2 times the vertex length, the definition used by SU2.
"""

import argparse

import numpy as np

from fiml_case import read_airfoil_marker, read_surface, upper_surface


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mesh", default="../S809_FIML_Forward/s809_cmesh.su2")
    p.add_argument("--airfoil", default="../S809_FIML_Forward/s809.dat")
    p.add_argument("--exp", default="s809_cp_upper_exp.dat")
    p.add_argument("--chord", type=float, default=0.6)
    p.add_argument("--output", default="TargetCp.dat")
    p.add_argument("--check", help="surface CSV whose Cp misfit is printed")
    args = p.parse_args()

    ids, xy, length = read_airfoil_marker(args.mesh)
    exp = np.loadtxt(args.exp)
    xc = xy[:, 0] / args.chord
    on_taps = upper_surface(xy / args.chord, args.airfoil) & (xc >= exp[0, 0]) & (xc <= exp[-1, 0])
    cp_target = np.interp(xc, exp[:, 0], exp[:, 1])

    with open(args.output, "w") as f:
        f.write('"PointID","x","y","Pressure","Pressure_Coefficient"\n')
        for i in np.flatnonzero(on_taps):
            f.write(f"{ids[i]}, {xy[i, 0]:.10e}, {xy[i, 1]:.10e}, 0.0, {cp_target[i]:.8f}\n")
    print(f"{on_taps.sum()} of {len(ids)} wall nodes have a target Cp, "
          f"x/c from {xc[on_taps].min():.4f} to {xc[on_taps].max():.4f}")

    if args.check:
        cp = read_surface(args.check, ids)
        misfit = np.sum(length[on_taps] * (cp[on_taps] - cp_target[on_taps]) ** 2)
        print(f"Cp misfit of {args.check}: {misfit:.6e}")


if __name__ == "__main__":
    main()
