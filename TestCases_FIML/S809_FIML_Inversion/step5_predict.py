#!/usr/bin/env python3
"""Step 5: SA-FIML with the trained network (KIND_TRAIN_NN= WEIGHTS), from the baseline solution.

Writes restart_predicted.csv, surface_predicted.csv and history_predicted.csv, and compares the predicted
beta with the inverted beta.
"""

import numpy as np

from fiml_case import keep, last_history, read_restart, report, require, run_parser, run_su2

ITER = 15000


def main():
    args = run_parser(__doc__).parse_args()
    require("nn_weights.dat", "restart_baseline.csv")
    keep("restart_baseline.csv", "restart_flow.csv")
    run_su2(args, "predict", {"KIND_TRAIN_NN": "WEIGHTS", "ITER_START_NN": "1", "ITER": str(ITER)})
    keep("restart_flow.csv", "restart_predicted.csv")
    keep("surface_flow.csv", "surface_predicted.csv")
    keep("history.csv", "history_predicted.csv")

    for label in ("baseline", "inverted"):
        try:
            report(label, last_history(f"history_{label}.csv"))
        except FileNotFoundError:
            pass
    report("predicted", last_history("history_predicted.csv"), ITER)
    try:
        b_inv = read_restart("restart_inverted.csv")["Beta_FIML"]
    except FileNotFoundError:
        return
    b_nn = read_restart("restart_predicted.csv")["Beta_FIML"]
    changed = abs(b_inv - 1.0) > 0.05
    print(f"where |beta_inverted - 1| > 0.05 ({changed.sum()} points): "
          f"RMS(beta_predicted - beta_inverted) = {np.sqrt(np.mean((b_nn - b_inv)[changed] ** 2)):.4f}")


if __name__ == "__main__":
    main()
