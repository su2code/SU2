#!/usr/bin/env python3
"""Compare per-iteration convergence history between two SU2_CFD builds.

Runs a case set that, between them, exercises every solver on the new
scalar edge-flux numerics framework (SA, SST, LM transition, species
transport, flamelet, solid heat conduction, flow-coupled/CHT heat), then
diffs the resulting history.csv files row by row and column by column.

Typical workflow: build SU2_CFD once per branch (e.g. develop and this
branch) into two separate prefixes, then

    ./compare_numerics_migration.py run --exe /path/to/develop/SU2_CFD --leg before
    ./compare_numerics_migration.py run --exe /path/to/branch/SU2_CFD  --leg after
    ./compare_numerics_migration.py compare

`run` can be invoked once per leg (e.g. from two different checkouts, at
different times) or you can pass both --before-exe/--after-exe to
`run-both` to do it in one shot. `compare` only reads back what `run`
already produced under --work-dir; it never launches SU2_CFD itself.

Large mesh/data files are usually not tracked in this repo's TestCases/
tree; if a case's data file isn't found next to its .cfg, it is pulled
from --data-root (default: the sibling TestCases repo checkout).

IMPORTANT: set the same OMP_NUM_THREADS (and mpirun rank count) for both
legs before invoking this script. Parallel reductions sum in an order
that depends on the thread count, so even the *same* binary run twice
with different OMP_NUM_THREADS will show spurious differences of up to
~1e-2 within the first few iterations that have nothing to do with any
code change.

FIXED (previously KNOWN ISSUE): flamelet cases used to crash ("free():
invalid pointer" inside an OpenMP parallel region) whenever OMP_NUM_THREADS
> 1, independent of mpirun rank count. Root cause: CSpeciesFlameletSolver::
Preprocessing accumulated a thread-local counter and then called
SU2_MPI::Reduce unguarded, so every OpenMP thread issued the collective
call concurrently. Fixed in CSpeciesFlameletSolver.{hpp,cpp} (shared,
atomically-accumulated counter; MPI_Reduce now runs on the master thread
only). Needs the same fix applied to whatever "before" binary you build,
or flamelet_ch4/flamelet_h2_prefdiff will still crash there at
OMP_NUM_THREADS > 1.
"""
from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

TESTCASES_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA_ROOT = "/home/pedro/Documents/SoftwareProjects/TestCases"


@dataclass
class Case:
    tag: str
    numerics: str
    cfg_dir: str
    cfg_file: str
    # {cfg filename (within cfg_dir): {KEY: value}}, applied to a copy of that file.
    overrides: dict = field(default_factory=dict)
    # Other cfg files (multizone sub-configs) copied alongside cfg_file, unmodified
    # unless also present in `overrides`.
    extra_cfg_files: list = field(default_factory=list)
    # Extra input files (mesh, inlet profiles, LUTs, ...), paths relative to
    # TESTCASES_ROOT. Resolved locally first, then under --data-root.
    data_files: list = field(default_factory=list)
    multizone: bool = False
    np: int = 2


def _disable_early_convergence(iters):
    return {
        "ITER": str(iters),
        "CONV_RESIDUAL_MINVAL": "-50",
        "CONV_STARTITER": "999999999",
        "TABULAR_FORMAT": "CSV",
    }


def build_cases(iters, outer_iters):
    single_zone_common = lambda hist: {**_disable_early_convergence(iters), "HISTORY_OUTPUT": hist}

    cases = [
        Case(
            tag="sa_flatplate",
            numerics="SA turbulence (CScalarFlux_SA)",
            cfg_dir="rans/flatplate",
            cfg_file="turb_SA_flatplate.cfg",
            overrides={"turb_SA_flatplate.cfg": single_zone_common("(ITER, RMS_RES, LINSOL, AERO_COEFF)")},
            data_files=["rans/flatplate/mesh_flatplate_turb_137x97.su2"],
        ),
        Case(
            tag="sst_rae2822",
            numerics="SST turbulence (CScalarFlux_SST)",
            cfg_dir="rans/rae2822",
            cfg_file="turb_SST_RAE2822.cfg",
            overrides={"turb_SST_RAE2822.cfg": single_zone_common("(ITER, RMS_RES, LINSOL, AERO_COEFF)")},
            data_files=["rans/rae2822/mesh_RAE2822_turb.su2"],
        ),
        Case(
            tag="lm_flatplate",
            numerics="LM transition (CScalarFlux_TransLM)",
            cfg_dir="transition/flatplate_lm",
            cfg_file="turb_LM_flatplate.cfg",
            overrides={"turb_LM_flatplate.cfg": single_zone_common("(ITER, RMS_RES, LINSOL, AERO_COEFF)")},
            data_files=[
                "rans/flatplate/mesh_flatplate_turb_137x97.su2",
                "rans/flatplate/inlet.dat",
            ],
        ),
        Case(
            tag="species_flatplate",
            numerics="species transport (CScalarFlux_Species)",
            cfg_dir="rans/flatplate",
            cfg_file="turb_SA_flatplate_species.cfg",
            overrides={
                "turb_SA_flatplate_species.cfg": single_zone_common(
                    "(ITER, RMS_RES, LINSOL, AERO_COEFF, SPECIES_COEFF, SPECIES_COEFF_SURF)"
                )
            },
            data_files=[
                "rans/flatplate/mesh_flatplate_turb_137x97.su2",
                "rans/flatplate/inlet.dat",
            ],
        ),
        Case(
            tag="flamelet_ch4",
            numerics="flamelet scalars (CScalarFlux_Flamelet)",
            cfg_dir="flamelet/01_laminar_premixed_ch4_flame_cfd",
            cfg_file="lam_prem_ch4_cfd.cfg",
            overrides={
                "lam_prem_ch4_cfd.cfg": single_zone_common(
                    "(ITER, RMS_RES, LINSOL, AERO_COEFF, FLOW_COEFF, FLOW_COEFF_SURF)"
                )
            },
            data_files=[
                "flamelet/01_laminar_premixed_ch4_flame_cfd/mesh_structured.cgns",
                "flamelet/01_laminar_premixed_ch4_flame_cfd/fgm_ch4.drg",
                "flamelet/01_laminar_premixed_ch4_flame_cfd/solution.dat",
            ],
        ),
        Case(
            tag="flamelet_h2_prefdiff",
            numerics="flamelet scalars with preferential diffusion (CScalarFlux_Flamelet, SetPreferentialDiffusionScalars)",
            cfg_dir="flamelet/07_laminar_premixed_h2_flame_cfd",
            cfg_file="laminar_premixed_h2_flame_cfd.cfg",
            overrides={
                "laminar_premixed_h2_flame_cfd.cfg": {
                    **single_zone_common("(ITER, RMS_RES, LINSOL, FLOW_COEFF, FLOW_COEFF_SURF)"),
                    # Original case defines no MARKER_ANALYZE, so FLOW_COEFF (outlet avg temperature,
                    # mass flow, ...) would otherwise be empty.
                    "MARKER_ANALYZE": "(outlet)",
                }
            },
            data_files=[
                "flamelet/07_laminar_premixed_h2_flame_cfd/2Dhex_BL.su2",
                "flamelet/07_laminar_premixed_h2_flame_cfd/MLP_TD.mlp",
                "flamelet/07_laminar_premixed_h2_flame_cfd/MLP_PD.mlp",
                "flamelet/07_laminar_premixed_h2_flame_cfd/MLP_PPV.mlp",
                "flamelet/07_laminar_premixed_h2_flame_cfd/solution.dat",
            ],
        ),
        Case(
            tag="heat_solid",
            numerics="heat, solid conduction (CScalarFlux_Heat, flow_solver=null)",
            cfg_dir="solid_heat_conduction/periodic_pins",
            cfg_file="configSolid.cfg",
            overrides={
                "configSolid.cfg": single_zone_common("(ITER, RMS_RES, HEAT, LINSOL)"),
            },
            data_files=["solid_heat_conduction/periodic_pins/solid.su2"],
        ),
        Case(
            tag="cht_multizone",
            numerics="heat, flow-coupled CHT (CScalarFlux_Heat, BC_ConjugateHeat_Interface)",
            cfg_dir="coupled_cht/incomp_2d",
            cfg_file="cht_2d_3cylinders.cfg",
            multizone=True,
            extra_cfg_files=["flow_cylinder.cfg", "solid_cylinder1.cfg", "solid_cylinder2.cfg", "solid_cylinder3.cfg"],
            overrides={
                "cht_2d_3cylinders.cfg": {
                    "OUTER_ITER": str(outer_iters),
                    "HISTORY_OUTPUT": (
                        "(ITER, BGS_RES[0], BGS_RES[1], BGS_RES[2], BGS_RES[3], "
                        "RMS_RES[0], RMS_RES[1], RMS_RES[2], RMS_RES[3], "
                        "LINSOL[0], LINSOL[1], LINSOL[2], LINSOL[3], "
                        "HEAT[0], HEAT[1], HEAT[2], HEAT[3])"
                    ),
                },
                "flow_cylinder.cfg": {"HISTORY_OUTPUT": "(ITER, RMS_RES, HEAT, LINSOL)", "TABULAR_FORMAT": "CSV"},
                "solid_cylinder1.cfg": {"HISTORY_OUTPUT": "(ITER, RMS_RES, HEAT, LINSOL)", "TABULAR_FORMAT": "CSV"},
                "solid_cylinder2.cfg": {"HISTORY_OUTPUT": "(ITER, RMS_RES, HEAT, LINSOL)", "TABULAR_FORMAT": "CSV"},
                "solid_cylinder3.cfg": {"HISTORY_OUTPUT": "(ITER, RMS_RES, HEAT, LINSOL)", "TABULAR_FORMAT": "CSV"},
            },
            data_files=["coupled_cht/incomp_2d/mesh_cht_3cyl.su2"],
        ),
    ]
    return {c.tag: c for c in cases}


def resolve_data_file(rel_path, data_root):
    local = os.path.join(TESTCASES_ROOT, rel_path)
    if os.path.exists(local):
        return local
    fallback = os.path.join(data_root, rel_path)
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(
        f"Could not find data file '{rel_path}' locally ({local}) or under --data-root ({fallback})"
    )


def apply_overrides(path, overrides):
    if not overrides:
        return
    with open(path) as f:
        lines = f.readlines()

    remaining = dict(overrides)
    out = []
    for line in lines:
        stripped = line.strip()
        key = stripped.split("=", 1)[0].strip() if ("=" in stripped and not stripped.startswith("%")) else None
        if key in remaining:
            out.append(f"{key}= {remaining.pop(key)}\n")
        else:
            out.append(line)
    for key, value in remaining.items():
        out.append(f"{key}= {value}\n")

    with open(path, "w") as f:
        f.writelines(out)


def prepare_case(case: Case, dest_dir, data_root):
    if os.path.exists(dest_dir):
        shutil.rmtree(dest_dir)
    os.makedirs(dest_dir)

    src_dir = os.path.join(TESTCASES_ROOT, case.cfg_dir)
    for cfg in [case.cfg_file] + case.extra_cfg_files:
        shutil.copy2(os.path.join(src_dir, cfg), os.path.join(dest_dir, cfg))

    for rel_path in case.data_files:
        src = resolve_data_file(rel_path, data_root)
        shutil.copy2(src, os.path.join(dest_dir, os.path.basename(rel_path)))

    for cfg, overrides in case.overrides.items():
        apply_overrides(os.path.join(dest_dir, cfg), overrides)


def run_case(case: Case, exe, work_dir, leg, data_root, np):
    dest_dir = os.path.join(work_dir, leg, case.tag)
    prepare_case(case, dest_dir, data_root)

    command = f"mpirun -n {np or case.np} {exe} {case.cfg_file}"
    log_path = os.path.join(dest_dir, "run.log")
    print(f"[{leg}/{case.tag}] {command}  (cwd={dest_dir})")
    with open(log_path, "w") as log:
        result = subprocess.run(command, shell=True, cwd=dest_dir, stdout=log, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        print(f"[{leg}/{case.tag}] FAILED (exit {result.returncode}), see {log_path}")
        return False
    print(f"[{leg}/{case.tag}] done")
    return True


def find_history_csv(dest_dir, case: Case):
    candidates = sorted(glob.glob(os.path.join(dest_dir, "*.csv")))
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    preferred = [c for c in candidates if os.path.basename(c) in ("history.csv", os.path.splitext(case.cfg_file)[0] + ".csv")]
    return preferred[0] if preferred else candidates[0]


def read_history(path):
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = [h.strip().strip('"') for h in next(reader)]
        rows = []
        for row in reader:
            if not row:
                continue
            rows.append([float(v) for v in row])
    return header, rows


ITER_COLUMNS = {"Time_Iter", "Outer_Iter", "Inner_Iter"}


def compare_histories(before_path, after_path, abs_tol, rel_tol):
    header_b, rows_b = read_history(before_path)
    header_a, rows_a = read_history(after_path)

    if header_b != header_a:
        only_b = [h for h in header_b if h not in header_a]
        only_a = [h for h in header_a if h not in header_b]
        return {
            "ok": False,
            "error": f"column mismatch: only-in-before={only_b} only-in-after={only_a}",
        }

    n = min(len(rows_b), len(rows_a))
    if len(rows_b) != len(rows_a):
        length_note = f"before has {len(rows_b)} rows, after has {len(rows_a)} rows; comparing first {n}"
    else:
        length_note = None

    per_column = []
    for col_idx, name in enumerate(header_b):
        is_iter_col = name in ITER_COLUMNS
        max_abs_diff = 0.0
        first_divergence = None
        for row_idx in range(n):
            a = rows_b[row_idx][col_idx]
            b = rows_a[row_idx][col_idx]
            if math.isnan(a) and math.isnan(b):
                continue
            diff = abs(a - b)
            if is_iter_col:
                if diff > 0 and first_divergence is None:
                    first_divergence = row_idx
                max_abs_diff = max(max_abs_diff, diff)
                continue
            scale = abs_tol + rel_tol * max(abs(a), abs(b))
            if diff > max_abs_diff:
                max_abs_diff = diff
            if diff > scale and first_divergence is None:
                first_divergence = row_idx
        per_column.append(
            {
                "name": name,
                "max_abs_diff": max_abs_diff,
                "first_divergence_row": first_divergence,
                "n_rows": n,
            }
        )

    return {"ok": True, "n_rows": n, "length_note": length_note, "columns": per_column}


def cmd_list(args):
    cases = build_cases(args.iters, args.outer_iters)
    print(f"{'tag':<20} {'numerics':<55} {'cfg':<40}")
    for case in cases.values():
        print(f"{case.tag:<20} {case.numerics:<55} {os.path.join(case.cfg_dir, case.cfg_file):<40}")


def select_cases(cases, only):
    if not only:
        return cases
    tags = set(only.split(","))
    missing = tags - set(cases)
    if missing:
        raise SystemExit(f"Unknown case tag(s): {', '.join(sorted(missing))}. Use 'list' to see valid tags.")
    return {t: c for t, c in cases.items() if t in tags}


def cmd_run(args):
    cases = select_cases(build_cases(args.iters, args.outer_iters), args.only)
    failures = []
    for case in cases.values():
        ok = run_case(case, args.exe, args.work_dir, args.leg, args.data_root, args.np)
        if not ok:
            failures.append(case.tag)
    if failures:
        print(f"\n{len(failures)} case(s) failed to run: {', '.join(failures)}")
        sys.exit(1)


def cmd_run_both(args):
    cases = select_cases(build_cases(args.iters, args.outer_iters), args.only)
    failures = []
    for leg, exe in (("before", args.before_exe), ("after", args.after_exe)):
        for case in cases.values():
            ok = run_case(case, exe, args.work_dir, leg, args.data_root, args.np)
            if not ok:
                failures.append(f"{leg}/{case.tag}")
    if failures:
        print(f"\n{len(failures)} run(s) failed: {', '.join(failures)}")
        sys.exit(1)
    cmd_compare(args)


def cmd_compare(args):
    cases = select_cases(build_cases(args.iters, args.outer_iters), args.only)
    report_lines = []
    any_missing = False
    any_diverged_early = False

    for case in cases.values():
        before_dir = os.path.join(args.work_dir, "before", case.tag)
        after_dir = os.path.join(args.work_dir, "after", case.tag)
        before_csv = find_history_csv(before_dir, case)
        after_csv = find_history_csv(after_dir, case)

        print(f"\n=== {case.tag} ({case.numerics}) ===")
        if not before_csv or not after_csv:
            print(f"  MISSING: before={before_csv} after={after_csv} (run both legs first)")
            any_missing = True
            continue

        result = compare_histories(before_csv, after_csv, args.abs_tol, args.rel_tol)
        if not result["ok"]:
            print(f"  ERROR: {result['error']}")
            any_diverged_early = True
            continue

        if result["length_note"]:
            print(f"  NOTE: {result['length_note']}")

        n_rows = result["n_rows"]
        header_line = f"  {'field':<20}{'max |diff|':>16}{'first divergence':>20}"
        print(header_line)
        report_lines.append(f"### {case.tag} ({case.numerics})\n")
        report_lines.append(f"| field | max |diff| | first divergence (row / {n_rows}) |\n|---|---|---|\n")
        for col in result["columns"]:
            if col["name"] in ITER_COLUMNS:
                continue
            if col["first_divergence_row"] is None:
                div_str = f"matched all {n_rows} rows"
            else:
                div_str = f"row {col['first_divergence_row']} / {n_rows}"
                any_diverged_early = any_diverged_early or col["first_divergence_row"] < n_rows * 0.5
            print(f"  {col['name']:<20}{col['max_abs_diff']:>16.3e}{div_str:>20}")
            report_lines.append(f"| {col['name']} | {col['max_abs_diff']:.3e} | {div_str} |\n")

    if args.out:
        with open(args.out, "w") as f:
            f.write("# Scalar-numerics migration: before/after convergence comparison\n\n")
            f.writelines(report_lines)
        print(f"\nWrote report to {args.out}")

    if any_missing:
        sys.exit(2)
    if any_diverged_early:
        print("\nAt least one field diverged in the first half of its run; inspect the table(s) above.")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work-dir", default=os.path.join(TESTCASES_ROOT, ".numerics_migration_compare"))
    parser.add_argument("--data-root", default=DEFAULT_DATA_ROOT, help="Fallback location for large mesh/data files")
    parser.add_argument("--iters", type=int, default=300, help="ITER override for single-zone cases")
    parser.add_argument("--outer-iters", type=int, default=50, help="OUTER_ITER override for the CHT multizone case")
    parser.add_argument("--only", default=None, help="Comma-separated case tags to restrict to")

    sub = parser.add_subparsers(dest="mode", required=True)

    p_list = sub.add_parser("list", help="Show the case table")
    p_list.set_defaults(func=cmd_list)

    p_run = sub.add_parser("run", help="Run all (or --only) cases with one SU2_CFD binary")
    p_run.add_argument("--exe", required=True, help="Path to the SU2_CFD executable")
    p_run.add_argument("--leg", required=True, choices=["before", "after"])
    p_run.add_argument("--np", type=int, default=None, help="Override mpirun rank count for every case")
    p_run.set_defaults(func=cmd_run)

    p_run_both = sub.add_parser("run-both", help="Run both legs with two binaries, then compare")
    p_run_both.add_argument("--before-exe", required=True)
    p_run_both.add_argument("--after-exe", required=True)
    p_run_both.add_argument("--np", type=int, default=None)
    p_run_both.add_argument("--abs-tol", type=float, default=1e-8)
    p_run_both.add_argument("--rel-tol", type=float, default=1e-6)
    p_run_both.add_argument("--out", default=None, help="Optional path to write a Markdown report")
    p_run_both.set_defaults(func=cmd_run_both)

    p_compare = sub.add_parser("compare", help="Compare history.csv files already produced by 'run'")
    p_compare.add_argument("--abs-tol", type=float, default=1e-8)
    p_compare.add_argument("--rel-tol", type=float, default=1e-6)
    p_compare.add_argument("--out", default=None, help="Optional path to write a Markdown report")
    p_compare.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
