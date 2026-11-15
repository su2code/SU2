#!/usr/bin/env python3
"""Inlet profile of the Morse swirling jet at x/D = 0.5 for the SU2 wedge mesh.

Fills the template that SU2 writes (example_inlet_jet.dat) with the measured axial, radial and swirl velocity of
station 1 (morse1.dat: r/x, r/r05, U/Ucl, V/Ucl, W/Wmax, ...), converted to Cartesian components on the wedge
face. The eddy viscosity follows from the measured shear stress, nu_t = -u'v' / (dU/dr), or is uniform with
--nut_ratio; nu_tilde follows from nu_t = nu_tilde f_v1.
"""

import argparse

import numpy as np


def nu_tilde_from_nu_t(nu_t, nu):
    """Invert nu_t = nu_tilde chi^3 / (chi^3 + 7.1^3) for nu_tilde by Newton iteration on chi."""
    target = nu_t / nu
    chi = np.maximum(target, 1.0)
    for _ in range(50):
        f = chi**4 / (chi**3 + 7.1**3) - target
        df = (chi**6 + 4 * 7.1**3 * chi**3) / (chi**3 + 7.1**3) ** 2
        chi = np.maximum(chi - f / df, 1e-6)
    return chi * nu


def nu_t_from_uv(r, U, uv, strain_fraction):
    """Eddy viscosity -u'v' / (dU/dr) at the points where |dU/dr| exceeds strain_fraction of its maximum.

    Returns those radii and eddy viscosities; dV/dx is neglected.
    """
    dU_dr = np.gradient(U, r)
    sheared = np.abs(dU_dr) >= strain_fraction * np.abs(dU_dr).max()
    return r[sheared], -uv[sheared] / dU_dr[sheared]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", default="morse1.dat")
    p.add_argument("--template", default="example_inlet_jet.dat")
    p.add_argument("--output", default="inlet_jet.dat")
    p.add_argument("--diameter", type=float, default=0.0246, help="D in m")
    p.add_argument("--x_over_d", type=float, default=0.5, help="station of the profile")
    p.add_argument("--u_cl", type=float, default=34.10, help="centreline velocity at the station, m/s")
    p.add_argument("--w_max", type=float, default=16.50, help="maximum swirl velocity at the station, m/s")
    p.add_argument("--nu", type=float, default=1.82e-5 / 1.2, help="kinematic viscosity, m^2/s")
    p.add_argument("--nut_ratio", type=float, default=None,
                   help="uniform nu_t / nu; without it nu_t comes from the measured u'v'")
    p.add_argument("--strain_fraction", type=float, default=0.1,
                   help="u'v' is used where |dU/dr| exceeds this fraction of its maximum")
    p.add_argument("--swirl_sign", type=float, default=1.0, help="+1: swirl in the +theta direction about +x")
    a = p.parse_args()

    d = np.loadtxt(a.data, comments="#")
    r_data = d[:, 0] * a.x_over_d * a.diameter
    U_data, V_data, W_data = d[:, 2] * a.u_cl, d[:, 3] * a.u_cl, d[:, 4] * a.w_max
    uv_data = d[:, 8] * a.u_cl**2

    with open(a.template) as f:
        lines = f.read().splitlines()
    head = [l for l in lines if not l[:1].isdigit() and not l[:1] == "-"]
    rows = np.array([[float(v) for v in l.split()] for l in lines if l[:1].isdigit() or l[:1] == "-"])
    y, z = rows[:, 1], rows[:, 2]
    r, th = np.hypot(y, z), np.arctan2(z, y)

    U = np.interp(r, r_data, U_data)
    V = np.interp(r, r_data, V_data)
    W = a.swirl_sign * np.interp(r, r_data, W_data)
    vel = np.column_stack([U, V * np.cos(th) - W * np.sin(th), V * np.sin(th) + W * np.cos(th)])
    mag = np.linalg.norm(vel, axis=1)
    rows[:, 4] = mag
    rows[:, 5:8] = vel / mag[:, None]
    if a.nut_ratio is None:
        r_fit, nu_t_fit = nu_t_from_uv(r_data, U_data, uv_data, a.strain_fraction)
        # Outside the sheared points nu_t is held at the nearest fitted value.
        nu_t = np.maximum(np.interp(r, r_fit, nu_t_fit), a.nu)
        for rf, nf in zip(r_fit, nu_t_fit):
            print(f"  r/D = {rf / a.diameter:5.3f}  nu_t/nu = {nf / a.nu:6.1f}")
    else:
        nu_t = np.full_like(r, a.nut_ratio * a.nu)
    rows[:, 8] = nu_tilde_from_nu_t(nu_t, a.nu)

    with open(a.output, "w") as f:
        f.write("\n".join(head) + "\n")
        for row in rows:
            f.write("\t".join(f"{v:.15e}" for v in row) + "\n")
    print(f"{a.output}: {len(rows)} points, |U| in [{mag.min():.2f}, {mag.max():.2f}] m/s, "
          f"nu_t/nu in [{nu_t.min() / a.nu:.1f}, {nu_t.max() / a.nu:.1f}]")


if __name__ == "__main__":
    main()
