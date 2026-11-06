"""Energy balance of the channel: the enthalpy flux drop of the fluid against the heat through the wall."""
import argparse
import csv

import cantera as ct

from nusselt import read_tecplot


def enthalpy_flux(data, gas, composition, x):
    """Mass flow and total enthalpy flux per unit width through the column at x."""
    idx = sorted((i for i, xi in enumerate(data["x"]) if abs(xi - x) < 1e-12), key=lambda i: data["y"][i])
    m = H = 0.0
    for a, b in zip(idx[:-1], idx[1:]):
        dy = data["y"][b] - data["y"][a]
        for i in (a, b):
            gas.TPY = data["Temperature"][i], data["Pressure"][i], composition
            mflux = data["Density"][i] * data["Velocity_x"][i]
            h0 = gas.enthalpy_mass + 0.5 * data["Velocity_x"][i] ** 2
            m += 0.5 * mflux * dy
            H += 0.5 * mflux * h0 * dy
    return m, H


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("volume", nargs="?", default="flow_0.dat")
    p.add_argument("--history", default="cht_channel.csv")
    p.add_argument("--mechanism", default="h2o2.yaml")
    p.add_argument("--composition", default="O2:0.2329, N2:0.7671")
    a = p.parse_args()

    data = read_tecplot(a.volume)
    gas = ct.Solution(a.mechanism)
    x_in, x_out = min(data["x"]), max(data["x"])
    m_in, H_in = enthalpy_flux(data, gas, a.composition, x_in)
    m_out, H_out = enthalpy_flux(data, gas, a.composition, x_out)
    print(f"mass flow per width  inlet {m_in:.6e}  outlet {m_out:.6e} kg/(m s)")
    print(f"enthalpy flux drop of the fluid  {H_in - H_out:.4f} W/m")

    rows = list(csv.DictReader(open(a.history)))
    last = {k.strip().strip('"'): v for k, v in rows[-1].items()}
    for key in (k for k in last if not k.startswith("bgs") and "Iter" not in k):
        print(f"{key}  {float(last[key]):.4f} W/m")


if __name__ == "__main__":
    main()
