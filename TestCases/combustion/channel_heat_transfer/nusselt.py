"""Local Nusselt number of the half channel from an SU2 Tecplot ASCII volume file.

Nu = q_w D_h / (k (T_w - T_b)), with D_h = 4 h, the bulk temperature T_b weighted with the mass flux and k
from Cantera at T_b. The Graetz coordinate is x* = x / (D_h Re Pr).
"""
import argparse
import re
from collections import defaultdict

import cantera as ct


def read_tecplot(filename):
    lines = open(filename).read().split("\n")
    names = re.findall(r'"([^"]+)"', lines[1])
    npoin = int(re.search(r"NODES= *(\d+)", lines[2]).group(1))
    rows = [list(map(float, lines[3 + i].split())) for i in range(npoin)]
    return {n: [r[k] for r in rows] for k, n in enumerate(names)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("volume", nargs="?", default="flow.dat")
    p.add_argument("--mechanism", default="h2o2.yaml")
    p.add_argument("--composition", default="O2:0.2329, N2:0.7671")
    p.add_argument("--wall_y", type=float, default=1.25e-4)
    p.add_argument("--wall_T", type=float, default=None, help="isothermal wall temperature, else from the data")
    p.add_argument("--csv", default="nusselt.csv")
    a = p.parse_args()

    data = read_tecplot(a.volume)
    gas = ct.Solution(a.mechanism)
    columns = defaultdict(list)
    for i, x in enumerate(data["x"]):
        columns[round(x, 12)].append(i)

    dh = 4.0 * a.wall_y
    out = []
    for x in sorted(columns):
        idx = sorted(columns[x], key=lambda i: data["y"][i])
        y = [data["y"][i] for i in idx]
        mflux = [data["Density"][i] * data["Velocity_x"][i] for i in idx]
        T = [data["Temperature"][i] for i in idx]
        m = sum(0.5 * (mflux[j] + mflux[j + 1]) * (y[j + 1] - y[j]) for j in range(len(y) - 1))
        mT = sum(0.5 * (mflux[j] * T[j] + mflux[j + 1] * T[j + 1]) * (y[j + 1] - y[j]) for j in range(len(y) - 1))
        if m <= 0.0:
            continue
        Tb = mT / m
        wall = idx[-1]
        Tw = a.wall_T if a.wall_T is not None else data["Temperature"][wall]
        qw = abs(data["Heat_Flux"][wall])
        gas.TPY = Tb, data["Pressure"][wall], a.composition
        k, mu, cp = gas.thermal_conductivity, gas.viscosity, gas.cp_mass
        Re = (m / a.wall_y) * dh / mu
        Pr = cp * mu / k
        if abs(Tw - Tb) < 1e-8:
            continue
        out.append((x, x / (dh * Re * Pr), Tb, Tw, qw, qw * dh / (k * abs(Tw - Tb)), Re))

    with open(a.csv, "w") as f:
        f.write("x,x_star,T_bulk,T_wall,q_wall,Nu,Re\n")
        for r in out:
            f.write(",".join(f"{v:.8g}" for v in r) + "\n")
    for r in out[:: max(1, len(out) // 16)]:
        print("x=%.4f x*=%.4f Tb=%.3f Tw=%.3f q=%.1f Nu=%.3f Re=%.1f" % r)


if __name__ == "__main__":
    main()
