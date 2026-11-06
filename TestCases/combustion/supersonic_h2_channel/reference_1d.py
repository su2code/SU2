"""Steady inviscid constant-area reacting flow (Rayleigh flow with chemistry), supersonic branch."""
import sys
import numpy as np
import cantera as ct
from scipy.integrate import solve_ivp


def inflow(gas, T, P, mach):
    gas.TP = T, P
    a = np.sqrt(gas.cp_mass / gas.cv_mass * P / gas.density)
    return mach * a


def solve(mech, X, T, P, mach, length, npts=400):
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    u0 = inflow(gas, T, P, mach)
    m = gas.density * u0
    mom = P + m * u0
    H = gas.enthalpy_mass + 0.5 * u0**2
    state = {"u": u0, "T": T}

    def thermo(Y):
        """Newton on (u, T) for p = rho R T and h = H - u^2/2 at fixed m, momentum and H."""
        u, Tk = state["u"], state["T"]
        for _ in range(50):
            gas.TPY = Tk, max(mom - m * u, 1.0), Y
            rho = m / u
            p = mom - m * u
            R = ct.gas_constant / gas.mean_molecular_weight
            f1 = p - rho * R * Tk
            f2 = gas.enthalpy_mass + 0.5 * u**2 - H
            J = np.array([[-m + m / u**2 * R * Tk, -rho * R], [u, gas.cp_mass]])
            du, dT = np.linalg.solve(J, [-f1, -f2])
            u, Tk = u + du, Tk + dT
            if abs(du) < 1e-10 * u and abs(dT) < 1e-8:
                break
        state["u"], state["T"] = u, Tk
        gas.TDY = Tk, m / u, Y
        return u

    def rhs(x, Y):
        u = thermo(Y)
        return gas.net_production_rates * gas.molecular_weights / (gas.density * u)

    x = np.linspace(0.0, length, npts)
    sol = solve_ivp(rhs, (0.0, length), gas.Y.copy(), method="BDF", t_eval=x, rtol=1e-8, atol=1e-12)
    rows = []
    for i in range(sol.y.shape[1]):
        u = thermo(sol.y[:, i])
        a = np.sqrt(gas.cp_mass / gas.cv_mass * gas.P / gas.density)
        rows.append([sol.t[i], gas.T, gas.P, gas.density, u, u / a] + list(gas.Y))
    return gas, np.array(rows)


if __name__ == "__main__":
    mech, X = "h2o2.yaml", "H2:0.4, O2:1, N2:3.76"
    gas, data = solve(mech, X, T=1200.0, P=ct.one_atm, mach=2.5, length=0.2)
    header = "x,T,P,rho,u,Mach," + ",".join(gas.species_names)
    np.savetxt(sys.argv[1] if len(sys.argv) > 1 else "reference_1d.csv", data, delimiter=",", header=header,
               comments="")
    for row in data[:: len(data) // 12]:
        print(f"x={row[0]:.4f} T={row[1]:.1f} P={row[2]:.0f} u={row[4]:.1f} M={row[5]:.3f}")
    print("inlet Y:", dict(zip(gas.species_names, data[0, 6:].round(6))))
