#!/usr/bin/env python3
"""Step 6: optimize the network weights directly with the adjoint (KIND_TRAIN_NN= WEIGHTS).

Starts from the trained weights (nn_weights.dat, kept as nn_weights_trained.dat) and the predicted solution.
The objective is the same as in the field inversion; beta = 1 + tanh(network) stays a function of the features.
Writes nn_weights_direct.dat, restart_direct.csv, surface_direct.csv and history_direct.csv.
"""

import os
import shutil
import subprocess
import sys

from fiml_case import keep, last_history, report, require, run_parser, write_config, BASE_CONFIG


def main():
    p = run_parser(__doc__)
    p.add_argument("--maxiter", type=int, default=10, help="L-BFGS-B iterations")
    p.add_argument("--adjoint_iter", type=int, default=5000)
    p.add_argument("--adjoint_cfl", type=float, default=100.0)
    p.add_argument("--adjoint_iter_warm", type=int, default=2000)
    p.add_argument("--x_scale", type=float, default=0.01, help="length of the first L-BFGS-B step in weight space")
    p.add_argument("--resume", action="store_true",
                   help="continue from nn_weights.dat, restart_flow.csv and the adjoint restart of an earlier run")
    p.add_argument("--log", default="direct_history.csv")
    args = p.parse_args()

    require("nn_weights.dat", "restart_predicted.csv")
    if not os.path.exists("nn_weights_trained.dat"):
        shutil.copy("nn_weights.dat", "nn_weights_trained.dat")
    if not args.resume:
        shutil.copy("nn_weights_trained.dat", "nn_weights.dat")
        keep("restart_predicted.csv", "restart_flow.csv")
    write_config(BASE_CONFIG, "direct.cfg", {"KIND_TRAIN_NN": "WEIGHTS", "ITER_START_NN": "1"})
    driver = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fiml_inversion.py")
    cmd = [sys.executable, driver, "-f", "direct.cfg", "-n", str(args.np), "--mode", "weights",
           "--su2_cfd", args.su2_cfd, "--su2_cfd_ad", args.su2_cfd_ad, "--maxiter", str(args.maxiter),
           "--adjoint_iter", str(args.adjoint_iter), "--adjoint_cfl", str(args.adjoint_cfl), "--adjoint_iter_warm", str(args.adjoint_iter_warm), "--x_scale", str(args.x_scale),
           "--fields", "CL", "CD", "Cp_Diff", "Beta_Reg", "--log", args.log]
    if subprocess.call(cmd) != 0:
        sys.exit("the weight optimization failed")

    keep("nn_weights.dat", "nn_weights_direct.dat")
    keep("restart_flow.csv", "restart_direct.csv")
    keep("surface_flow.csv", "surface_direct.csv")
    keep("history.csv", "history_direct.csv")
    for label in ("baseline", "inverted", "predicted", "direct"):
        report(label, last_history(f"history_{label}.csv"))


if __name__ == "__main__":
    main()
