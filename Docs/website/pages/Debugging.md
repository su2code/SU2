---
title: Debugging SU2
permalink: /docs_v7/Debugging/
---

This page shows practical ways to find bugs in SU2: crashes, wrong results, divergence, data races and slow code.
It names the build options, config options and tools that exist for this.

Tips marked **(advice)** are general practice, not features of SU2. Words such as "zone", "halo point" or
"reducer strategy" are explained in the [Glossary](../Glossary/).

---

- [Build for debugging](#build-for-debugging)
- [Check the config file first](#check-the-config-file-first)
- [Reproduce with a small case](#reproduce-with-a-small-case)
- [Use a debugger](#use-a-debugger)
- [Sanitizers](#sanitizers)
- [OpenMP data races](#openmp-data-races)
- [Floating point exceptions](#floating-point-exceptions)
- [Find where a NaN comes from](#find-where-a-nan-comes-from)
- [Algorithmic differentiation](#algorithmic-differentiation)
- [Profiling](#profiling)
- [Valgrind](#valgrind)

---

## Build for debugging

Make a separate build folder for debugging, so your normal build stays fast **(advice)**:

```
./meson.py setup build_debug --buildtype=debug --prefix=$PWD/install_debug
./ninja -C build_debug install
```

The build type changes what the compiler does:

| Build type | Optimization | Debug symbols (`-g`) | `assert` checks |
|---|---|---|---|
| `release` (default) | yes | no | off |
| `debugoptimized` | yes (`-O2`) | yes | off |
| `debug` | no | yes | on |

Notes:

- SU2 sets `b_ndebug=true` in `meson.build`, which turns `assert` off. Only `--buildtype=debug` turns it back on
  (SU2 adds `-UNDEBUG`). To keep the checks with optimization, add the meson option `-Db_ndebug=false` to a
  `debugoptimized` build **(advice)**.
- With `assert` on, SU2 also runs some extra checks that are inside `#ifndef NDEBUG`. For example
  `CSysMatrix` checks the sizes of the vectors in the matrix-vector product, and prints a message when it
  regularizes a small pivot.
- `debugoptimized` is a good default for gdb and profilers: the code is fast enough for real cases, and the stack
  traces have function names and lines **(advice)**. Some variables show as "optimized out"; use `debug` when you
  need them.
- `--warnlevel=3` turns on most compiler warnings. The CI builds with `--warnlevel=3` and `--werror` (warnings
  become errors) for most configurations, so check your change with these flags before you push.
- `-Dcpu-arch` (default `native`) sets `-march`. `-Dcpu-arch=` (empty) removes it. This helps with tools that do
  not know the newest instructions, see [Valgrind](#valgrind).
- SU2 always compiles with `-ffast-math -fno-finite-math-only` (when the compiler accepts these flags). NaN checks
  with `std::isnan` still work. Other effects of fast math are described in
  [Floating point exceptions](#floating-point-exceptions).

Build options in `meson_options.txt` that matter when debugging:

| Option | Use |
|---|---|
| `-Dwith-mpi=disabled` | Serial build. Simpler to debug, and needed for the address sanitizer jobs of the CI. |
| `-Dwith-omp=true` | OpenMP build, for thread problems. |
| `-Denable-tests=true` | Builds the unit tests, see [Running Unit Tests](../Running-Unit-Tests/). |
| `-Denable-mixedprec=true` | Linear algebra in single precision. Try without it if a linear solver fails. |
| `-Denable-autodiff=true`, `-Denable-directdiff=true` | AD builds, see [Algorithmic differentiation](#algorithmic-differentiation). |
| `-Dtracy_enable=true` | Tracy profiler, see [Profiling](#profiling). |
| `-Denable-gprof=true` | gprof profiling, see [Profiling](#profiling). |

To change the options of an existing build folder, use `./meson.py setup build_debug --reconfigure <options>`.
See [Build SU2 on Linux/MacOS](../Build-SU2-Linux-MacOS/) for all build steps.

## Check the config file first

Many problems are config problems. Check the config before you look at the code.

### The dry run

```
SU2_CFD -d case.cfg
```

`-d` (or `--dryrun`) starts `CDummyDriver`. SU2 reads and checks the whole config file, sets up the solvers on an
empty geometry, prints the available volume and history output fields, and stops. It does not read the mesh
points, so it is fast even for large cases. The mesh file must still exist, because SU2 reads the number of
dimensions from it.

### Read the error box

When SU2 stops on an error, it prints a box like this (from `SU2_MPI::Error`):

```
Error in "void CConfig::SetConfig_Parsing(std::istream&)":
-------------------------------------------------------------------------
Line 12 KIND_TURB_MODLE: invalid option name. Check current SU2 options in config_template.cfg.
Did you mean KIND_TURB_MODEL?
------------------------------ Error Exit -------------------------------
```

- The first line names the C++ function that stopped. Search for the message in that function to see the exact
  check.
- The parser collects all errors of the file, up to 30, before it stops. Fix them from the top.
- Errors about the config file come in three kinds: an unknown option name (with a "Did you mean" guess), an
  option that appears twice, and a bad value (`invalid option value`, with the list of allowed values).
- Other errors come later, from the checks in `CConfig::SetPostprocessing` and from the solvers.

The [Error Messages](../Error-Messages/) page lists every `SU2_MPI::Error` message with a link to its code. If an
old config file stops with "invalid option name", look for the option in
[Renamed and Removed Options](../Renamed-Options/). The [Option Reference](../Option-Reference/) lists all
current options.

## Reproduce with a small case

Debugging is faster on a small case that runs in seconds **(advice)**. Good sources are the regression cases in
`TestCases/`. [Running Regression Tests](../Running-Regression-Tests/#run-one-test) explains how to run one test.

When the test script runs a case, it writes a copy of the config file with the number of iterations changed,
named `<config>.cfg.autotest`, and the screen output to `<config>.log`. Both stay in the case folder. To repeat
exactly what the test did, run the `.autotest` file by hand, for example in gdb:

```
cd euler/naca0012
gdb --args SU2_CFD inv_NACA0012_Roe.cfg.autotest
```

Lower `ITER` (or `INNER_ITER`, `OUTER_ITER`, `TIME_ITER`) to the first iteration that shows the problem
**(advice)**.

## Use a debugger

### Serial runs

```
gdb --args SU2_CFD case.cfg
(gdb) run
(gdb) bt
```

After a crash, `bt` prints the call stack. To stop where SU2 reports an error, set a breakpoint before `run`
**(advice)**. `SU2_MPI::Error` ends with `MPI_Abort` in an MPI build, and with `exit` in a build without MPI:

```
(gdb) break MPI_Abort     # build with MPI
(gdb) break exit          # build without MPI
```

To stop at a solver function, use its full name, for example `break CSolver::SetResidual_RMS` **(advice)**.

### Runs with MPI

Some bugs only show in parallel, for example at the halo points. Three ways to use gdb with MPI **(advice)**:

- One terminal per rank (needs a display):

  ```
  mpirun -np 2 xterm -e gdb --args SU2_CFD case.cfg
  ```

- Only one rank in gdb, the others run normally (Open MPI syntax):

  ```
  mpirun -np 1 gdb --args SU2_CFD case.cfg : -np 1 SU2_CFD case.cfg
  ```

- Attach to a running rank: start the run, find the process id with `ps` or `pgrep SU2_CFD`, then
  `gdb -p <pid>`. This is useful for a run that hangs. Run `bt` in each rank to see where they wait.

When only rank 0 prints an error, the bug may be on another rank. `SU2_MPI::Error` prints the message from the
lowest rank that called it.

### The Python wrapper

The Python wrapper (built with `-Denable-pywrapper=true`) runs the same C++ code. Debug it with
`gdb --args python3 your_script.py` **(advice)**.

## Sanitizers

Sanitizers are compiler tools that check the program while it runs:

- the address sanitizer (asan) finds memory errors, such as reading past the end of an array, use after free and
  leaks;
- the thread sanitizer (tsan) finds data races between threads.

The program runs slower, and stops with a report on the first error.

### How the CI uses them

The workflow `.github/workflows/regression.yml` has these jobs:

| Build job | Configurations | Main build flags | Test job |
|---|---|---|---|
| Build SU2 (asan) | `BaseNoMPI-asan`, `ReverseNoMPI-asan` | `--buildtype=debugoptimized -Dwith-mpi=disabled`; the reverse build also has `--optimization=1 -Denable-autodiff=true -Denable-normal=false` | Address Sanitizer Tests: `serial_regression.py`, `serial_regression_AD.py` with `--asan` |
| Build SU2 (tsan) | `BaseOMP-tsan`, `ReverseOMP-tsan` | `--buildtype=debugoptimized -Dwith-omp=true -Denable-mixedprec=true -Denable-tecio=false`; the reverse build also has `-Denable-autodiff=true -Denable-normal=false` | Thread Sanitizer Tests: `hybrid_regression.py`, `hybrid_regression_AD.py` with `--tsan` |

The sanitizer flags themselves are not in these flags. They are set by the build containers
`ghcr.io/su2code/su2/build-su2-asan` and `build-su2-tsan`, made from the
[Docker-Builds](https://github.com/su2code/Docker-Builds) repository. The tsan container has a compiler and
OpenMP library built for the thread sanitizer.

### What `--asan` and `--tsan` do

The test scripts read these arguments with `parse_args()` in `TestCase.py`. With one of them:

- tests with `enabled_with_asan = False` (or `enabled_with_tsan = False`) are skipped;
- the number of iterations is cut to at most 2 after the restart iteration, because sanitized runs are slow;
- a test fails only if SU2 returns an error. The results are not compared with the stored values.

### Run a sanitizer on your computer

The easiest way is the CI container, see [Container Development](../Container-Development/). It also explains how
to work inside the tsan container.

You can also use the meson option `b_sanitize` **(advice)**:

```
./meson.py setup build_asan --buildtype=debugoptimized -Dwith-mpi=disabled -Db_sanitize=address
./meson.py setup build_tsan --buildtype=debugoptimized -Dwith-omp=true -Db_sanitize=thread
```

Use a build without MPI for asan, like the CI: MPI libraries give many reports that are not SU2 bugs
**(advice)**. For tsan with OpenMP, the OpenMP library must also be built for tsan. Otherwise you get false
reports, which is why the CI uses a special container **(advice)**.

## OpenMP data races

A data race happens when two threads write to the same memory, or one writes while another reads, without
synchronization. The result then changes from run to run.

How to look for one:

1. Run the case with threads: `SU2_CFD -t 2 case.cfg` (`-t` sets the number of OpenMP threads per MPI rank). The
   hybrid regression tests use `-t 2`.
2. Run it several times with the same number of threads, and compare the history files. Results that change
   between identical runs point to a race **(advice)**. Small differences between different thread counts can
   also come from a different order of sums, and are not always a bug **(advice)**.
3. Set `EDGE_COLORING_GROUP_SIZE= 0`. This forces the reducer strategy in the flow and scalar solvers: the edge
   fluxes are stored first, then summed for each point, so the edge loop has no shared writes. If the problem
   disappears, look at the edge loops and the edge coloring.
4. Run the case in a tsan build, see [Sanitizers](#sanitizers).

At start-up, SU2 prints either "All ranks use edge coloring." or a warning that some ranks use the fallback (the
reducer strategy), with the coloring efficiency.

Write thread code only with the macros of `Common/include/parallelization/omp_structure.hpp`
(`SU2_OMP_PARALLEL`, `SU2_OMP_FOR_STAT`, `END_SU2_OMP_FOR`, `SU2_OMP_MASTER`, `BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS`, ...).
A block such as `SU2_OMP_PARALLEL` or `SU2_OMP_FOR_STAT` needs its `END_SU2_OMP_...` macro. When you configure a build with both
`-Dwith-omp=true` and `-Denable-autodiff=true`, meson runs the OpDiLib syntax checker (rules in
`su2omp.syntax.json`) and stops if the macros are not used correctly.

See [Parallelism](../Code-Structure/#parallelism-mpi-openmp-and-vectorization) for how the threads share the work.

## Floating point exceptions

SU2 has no config option or build option to stop on floating point exceptions (division by zero, overflow,
invalid operations that make a NaN).

`SU2_CFD/src/SU2_CFD.cpp` has two commented lines for this. To use them, remove the comment marks from both,
and rebuild:

```cpp
#include <fenv.h>                                      // near the top of the file
feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO );  // in main()
```

The program then stops with the signal `SIGFPE` at the first such operation. Run it in gdb to see where.

Limits **(advice)**:

- `feenableexcept` is a GNU extension (Linux with glibc).
- With `-ffast-math` the compiler assumes that math operations do not trap. It may move or add operations, so a
  trap can come from code that is correct. A `debug` build has no optimization, so it gives cleaner results.
- On x86, `-ffast-math` may also link code that sets denormal numbers to zero. Results with very small numbers can
  then differ from a build without fast math.
- Libraries called before or by SU2 (MPI, BLAS) can also trigger the trap.

## Find where a NaN comes from

### Where SU2 stops

SU2 checks for NaN in a few places and stops with `SU2 has diverged (NaN detected).`:

- `CSolver::SetResidual_RMS`, on the root mean square residual of each equation. It also stops when the residual
  is larger than 10^20 (`SU2 has diverged (Residual > 10^20 detected).`).
- `COutput` (`SU2_CFD/src/output/COutput.cpp`), on the history fields used for the convergence check. Here SU2
  stops without saving the solution.
- The DG solver (`CFEM_DG_EulerSolver`) has its own checks.

The NEMO solvers also check each flux and source term with `CNumerics::CheckResidualNaNs`.

The check happens after the iteration is done, so the NaN can come from anywhere in that iteration. The steps
below narrow it down.

### Write the solution before it breaks

Write the volume files at every iteration, and keep one file per iteration:

```
OUTPUT_WRT_FREQ= 1
WRT_VOLUME_OVERWRITE= NO
VOLUME_OUTPUT= SOLUTION, PRIMITIVE, RESIDUAL, LIMITER, MESH_QUALITY, MPI
```

The last file written before the error shows the state just before the NaN. Useful volume output groups:

| Group | Fields | What to look for |
|---|---|---|
| `RESIDUAL` | `RES_DENSITY`, `RES_MOMENTUM-X`, ..., `RES_TKE`, ... | Where the residual grows first. |
| `LIMITER` | `LIMITER_PRESSURE`, ..., `LIMITER_TKE`, ... | Points where the limiter goes to 0 (strong gradients, shocks). |
| `PRIMITIVE` | `PRESSURE`, `TEMPERATURE`, `MACH`, `EDDY_VISCOSITY`, ... | Negative or extreme values. |
| `MESH_QUALITY` | `ORTHOGONALITY`, `ASPECT_RATIO`, `VOLUME_RATIO` | Bad cells near the problem. |
| `MPI` | `RANK` | Problems on a partition border point to halo exchange bugs. |

The groups that exist depend on the solver. Run `SU2_CFD -d case.cfg` to print them, or see the
[Output Fields Reference](../Output-Fields-Reference/).

In a long run, write the files only from a few iterations before the failure: restart from a solution saved
before that **(advice)**.

### History fields to watch

- `RMS_RES` and `MAX_RES` groups: which equation grows first.
- `NONPHYSICAL_POINTS`: the number of points where the flow solver found a non-physical state (for example a
  negative pressure or density). At these points the compressible solver goes back to the solution of the previous
  iteration. A growing number often comes before a NaN. At the end of the run
  SU2 also prints a warning with this number.
- `MIN_CFL`, `MAX_CFL`, `AVG_CFL` (group `CFL_NUMBER`, with `CFL_ADAPT= YES`): a CFL number that keeps growing
  while the residual grows.
- `LINSOL` group (for example `LINSOL_ITER`, `LINSOL_RESIDUAL`): a linear solver that does not converge.

Add them with `HISTORY_OUTPUT` and `SCREEN_OUTPUT`.

### Common checks

These help to tell a bug from a hard case **(advice)**:

- Lower `CFL_NUMBER`, or use `CFL_ADAPT= YES`.
- Run first order (`MUSCL_FLOW= NO`, `MUSCL_TURB= NO`). If that works, look at the gradients and limiters.
- Run serial and in parallel. If only the parallel run fails, look at the halo exchange.
- Run with a `debug` build, so `assert` checks are on.

## Algorithmic differentiation

The discrete adjoint has its own debugging tools: the tape test mode `--tapetests` of the test scripts, the
options `CHECK_TAPE_TYPE` and `CHECK_TAPE_VARIABLES` (for builds with `-Dcodi-tape=Tag`), the tape statistics of
`WRT_AD_STATISTICS`, and the CoDiPack build options. A `debug` or `debugoptimized` AD build also turns on the
internal checks of CoDiPack (`CODI_EnableAssert`).

See [Automatic Differentiation](../Automatic-Differentiation/) for how to use them.

## Profiling

### Built-in timing

`WRT_PERFORMANCE= YES` prints a performance summary at the end of `SU2_CFD`: the time spent in preprocessing,
computation and output, the time per iteration, points per core and ghost points per core. It also prints the
write speed of each output file. The history fields `ITER_TIME` and `WALL_TIME` (group `WALL_TIME`) give the time
per iteration.

### Tracy

[Tracy](../Tracy-Integration/) shows a timeline of the instrumented code sections. The build option is
`-Dtracy_enable=true`; use it with `--buildtype=debugoptimized` (meson prints a warning otherwise). Mark code
sections with `SU2_ZONE_SCOPED` or `SU2_ZONE_SCOPED_N("name")` from `Common/include/tracy_structure.hpp`. These
macros are empty when Tracy is off.

### perf

On Linux, `perf` samples the whole program without changes to the code **(advice)**:

```
perf record -g SU2_CFD case.cfg
perf report
```

Use a `debugoptimized` build, so the report has function names and lines **(advice)**. For an MPI run, profile
one rank, for example with `mpirun -np 1 perf record -g SU2_CFD case.cfg : -np 3 SU2_CFD case.cfg` (Open MPI
syntax) **(advice)**.

### gprof

`-Denable-gprof=true` adds `-pg` to the `SU2_CFD` executables (see
[Build SU2 on Linux/MacOS](../Build-SU2-Linux-MacOS/)). Only the main file `SU2_CFD.cpp` and the link step get the
flag, not the library `SU2core` with the solver code. So the call graph of gprof is mostly empty. Prefer perf or
Tracy **(advice)**.

## Valgrind

Valgrind (memcheck) finds memory errors like the address sanitizer, without a special build, but it is much
slower **(advice)**. Notes:

- The default build uses `-march=native`. Valgrind may not know the newest vector instructions (for example
  AVX-512) and then stops with "unrecognised instruction". Build with `-Dcpu-arch=` (empty) or an older
  architecture for valgrind runs **(advice)**.
- Use a build without MPI, or expect many reports from the MPI library **(advice)**.
- Use a very small case and few iterations **(advice)**.

```
valgrind --track-origins=yes SU2_CFD case.cfg
```

Last checked against SU2 develop `11c177b904`.
