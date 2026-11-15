"""Mass flow and enthalpy flux along the channel, compressible (channel_lowmach.cfg) against incompressible
(channel_lowmach_inc.cfg) solver. The incompressible pressure is relative to the thermodynamic pressure."""
import argparse

import cantera as ct

from energy_balance import enthalpy_flux
from nusselt import read_tecplot


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("compressible")
    p.add_argument("incompressible")
    p.add_argument("--pressure_thermodynamic", type=float, default=101325.0)
    p.add_argument("--mechanism", default="h2o2.yaml")
    p.add_argument("--composition", default="O2:0.2329, N2:0.7671")
    p.add_argument("--stations", type=float, nargs="+", default=[0.00125, 0.0025, 0.005, 0.01, 0.0199])
    a = p.parse_args()

    gas = ct.Solution(a.mechanism)
    comp = read_tecplot(a.compressible)
    inc = read_tecplot(a.incompressible)
    inc["Pressure"] = [q + a.pressure_thermodynamic for q in inc["Pressure"]]

    print("   x [mm]   m_comp      m_inc     diff%   H_comp [W/m]  H_inc [W/m]  diff%")
    for x in a.stations:
        xc = min(set(comp["x"]), key=lambda v: abs(v - x))
        mc, hc = enthalpy_flux(comp, gas, a.composition, xc)
        mi, hi = enthalpy_flux(inc, gas, a.composition, xc)
        print(f"{xc * 1e3:8.3f}  {mc:.5e}  {mi:.5e}  {100 * (mc - mi) / mi:6.2f}  {hc:12.4f}  {hi:12.4f}"
              f"  {100 * (hc - hi) / hi:6.2f}")


if __name__ == "__main__":
    main()
