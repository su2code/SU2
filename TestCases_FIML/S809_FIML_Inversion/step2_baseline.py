#!/usr/bin/env python3
"""Step 2: SA-FIML with beta = 1 (plain SA) from free stream.

Writes restart_baseline.csv, surface_baseline.csv and history_baseline.csv and prints the lift, drag and
objective terms. The run should stop on the Cauchy test, well before ITER.
"""

from fiml_case import keep, last_history, report, require, run_parser, run_su2

ITER = 15000


def main():
    args = run_parser(__doc__).parse_args()
    require("TargetCp.dat")
    run_su2(args, "baseline", {"RESTART_SOL": "NO", "ITER": str(ITER)})
    keep("restart_flow.csv", "restart_baseline.csv")
    keep("surface_flow.csv", "surface_baseline.csv")
    keep("history.csv", "history_baseline.csv")
    report("baseline", last_history("history_baseline.csv"), ITER)


if __name__ == "__main__":
    main()
