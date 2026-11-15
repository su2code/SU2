"""Centerline state of an SU2 compressible Cantera restart, reconstructed with Cantera."""
import csv, sys
import numpy as np
import cantera as ct

names = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR"]
gas = ct.Solution(sys.argv[2] if len(sys.argv) > 2 else "h2o2.yaml")
rows = list(csv.DictReader(open(sys.argv[1])))
key = {k.strip('" '): k for k in rows[0]}
get = lambda r, n: float(r[key[n]])
ymid = sorted({get(r, "y") for r in rows})[len({get(r, "y") for r in rows}) // 2]
line = sorted((r for r in rows if abs(get(r, "y") - ymid) < 1e-9), key=lambda r: get(r, "x"))
out = []
for r in line:
    rho = get(r, "Density")
    u = get(r, "Momentum_x") / rho
    e = get(r, "Energy") / rho - 0.5 * u**2
    Y = {n: max(get(r, "Species_" + n), 0.0) for n in names}
    Y["N2"] = max(1.0 - sum(Y.values()), 0.0)
    gas.UVY = e, 1.0 / rho, Y
    out.append([get(r, "x"), gas.T, gas.P, rho, u, gas["H2"].Y[0], gas["OH"].Y[0], gas["H2O"].Y[0]])
out = np.array(out)
np.savetxt(sys.argv[3] if len(sys.argv) > 3 else "centerline.csv", out, delimiter=",",
           header="x,T,P,rho,u,H2,OH,H2O", comments="")
for o in out[:: max(1, len(out) // 20)]:
    print("x=%.3f T=%.1f P=%.0f rho=%.4f u=%.1f H2=%.2e OH=%.2e H2O=%.3e" % tuple(o))
