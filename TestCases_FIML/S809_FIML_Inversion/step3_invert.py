#!/usr/bin/env python3
"""Step 3: field inversion of beta at every point, starting from the baseline.

Runs ../fiml_inversion.py in field mode on restart_flow.csv and keeps the result as restart_inverted.csv,
surface_inverted.csv and history_inverted.csv. inversion_history.csv holds J and its terms per evaluation.
"""

import os
import subprocess
import sys

import numpy as np

from fiml_case import keep, last_history, read_restart, report, require, run_parser, write_config, BASE_CONFIG


def main():
    p = run_parser(__doc__)
    p.add_argument("--maxiter", type=int, default=15, help="L-BFGS-B iterations")
    p.add_argument("--adjoint_iter", type=int, default=5000)
    p.add_argument("--adjoint_cfl", type=float, default=100.0)
    p.add_argument("--warm_adjoint", action="store_true", help="start each adjoint from the previous one")
    p.add_argument("--resume", action="store_true", help="continue from the beta in restart_flow.csv")
    args = p.parse_args()

    require("TargetCp.dat", "restart_baseline.csv")
    if not args.resume:
        keep("restart_baseline.csv", "restart_flow.csv")
    if not args.resume and os.path.exists("restart_adj_custom.csv"):
        os.remove("restart_adj_custom.csv")
    write_config(BASE_CONFIG, "invert.cfg", {})
    driver = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fiml_inversion.py")
    cmd = [sys.executable, driver, "-f", "invert.cfg", "-n", str(args.np), "--mode", "field",
           "--su2_cfd", args.su2_cfd, "--su2_cfd_ad", args.su2_cfd_ad, "--maxiter", str(args.maxiter),
           "--adjoint_iter", str(args.adjoint_iter), "--adjoint_cfl", str(args.adjoint_cfl),
           "--fields", "CL", "CD", "Cp_Diff", "Beta_Reg"]
    if args.warm_adjoint:
        cmd.append("--warm_adjoint")
    if subprocess.call(cmd) != 0:
        sys.exit("the inversion failed")

    keep("restart_flow.csv", "restart_inverted.csv")
    keep("surface_flow.csv", "surface_inverted.csv")
    keep("history.csv", "history_inverted.csv")
    report("baseline", last_history("history_baseline.csv"))
    report("inverted", last_history("history_inverted.csv"))
    beta = read_restart("restart_inverted.csv")["Beta_FIML"]
    print(f"beta: min {beta.min():.3f}, max {beta.max():.3f}, |beta - 1| > 0.1 at {(abs(beta - 1) > 0.1).sum()} "
          f"of {beta.size} points")


if __name__ == "__main__":
    main()
