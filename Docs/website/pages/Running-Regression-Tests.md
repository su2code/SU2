---
title: Running Regression Tests
permalink: /docs_v7/Running-Regression-Tests/
---

The regression tests run SU2 on many small cases. For each case, they compare a few numbers of the output with
stored values. This catches changes in the results that nobody intended. Every pull request runs the tests, and a
pull request is merged only when they pass.

This page explains:

- where the tests are;
- how the CI runs them;
- how to run them on your computer;
- how to add a test;
- how to update the stored values.

For the C++ unit tests, see [Running Unit Tests](../Running-Unit-Tests/) and
[Writing Unit Tests](../Writing-Unit-Tests/).

## Where the tests are

The tests use three repositories:

- **SU2** (`su2code/SU2`): the test scripts and the configuration files, in `TestCases/`.
- **TestCases** (`su2code/TestCases`): the meshes, restart files and other large input files. Its folders match
  the folders of `TestCases/` in SU2, and the two sets of files do not overlap. See [Test Cases](../Test-Cases/).
- **Tutorials** (`su2code/Tutorials`): the configuration files and meshes of the tutorials. `tutorials.py` uses them.

The scripts in `TestCases/` are:

| Script | What it runs | Default command |
|---|---|---|
| `serial_regression.py` | Main tests without MPI, plus `SU2_GEO`, `SU2_DEF`, Python scripts and Python wrapper tests | `SU2_CFD` |
| `parallel_regression.py` | Main tests with MPI, plus `SU2_SOL`, `SU2_DEF` and output file tests | `mpirun -n 2 SU2_CFD` |
| `serial_regression_AD.py` | Discrete adjoint and direct differentiation tests without MPI | `SU2_CFD_AD` |
| `parallel_regression_AD.py` | Discrete adjoint tests with MPI | `mpirun -n 2 SU2_CFD_AD` |
| `hybrid_regression.py` | Tests with OpenMP threads | `SU2_CFD -t 2` |
| `hybrid_regression_AD.py` | Discrete adjoint tests with OpenMP threads | `SU2_CFD_AD -t 2` |
| `tutorials.py` | Cases of the Tutorials repository | `mpirun -n 2 SU2_CFD` |
| `vandv.py` | A few verification and validation cases (`TestCases/vandv/`) | `mpirun -n 2 SU2_CFD` |

All scripts use the class `TestCase` from `TestCases/TestCase.py`.

## How a test works

Each test is a `TestCase` object. This entry is copied from `TestCases/serial_regression.py`:

```python
    # NACA0012
    naca0012           = TestCase('naca0012')
    naca0012.cfg_dir   = "euler/naca0012"
    naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
    naca0012.test_iter = 20
    naca0012.test_vals = [-4.031446, -3.524530, 0.273056, 0.009703]
    test_list.append(naca0012)
```

`TestCase.run_test()` does the following:

1. It copies `cfg_file` to `<cfg_file>.autotest` in `cfg_dir`, and sets the number of iterations to `test_iter + 1`.
   The option it changes is:
   - `ITER` by default;
   - `TIME_ITER` if `unsteady` is set;
   - `OUTER_ITER` if `multizone` is set and `unsteady` is not.
2. It runs the command in `cfg_dir`, and writes the screen output to `<cfg_file>.log`
   (here `inv_NACA0012_Roe.cfg.log`).
3. It reads the screen output after the line `Begin Solver`. It finds the first row whose first column is `test_iter`.
4. It takes the **last** `len(test_vals)` columns of that row and compares them with `test_vals`.

The test fails if:

- one absolute difference is larger than `tol`;
- the row is not found;
- the run takes longer than `timeout` seconds;
- the command returns an error.

The configuration file of this test has `SCREEN_OUTPUT = (INNER_ITER, RMS_DENSITY, RMS_ENERGY, LIFT, DRAG)`.
So the four values are the last four columns of the screen output:

```
|  Inner_Iter|    rms[Rho]|   rms[RhoE]|          CL|          CD|
```

### Main attributes of TestCase

| Attribute | Meaning | Default in `TestCase.py` |
|---|---|---|
| `cfg_dir`, `cfg_file` | Folder (relative to `TestCases/`) and configuration file | `"."`, `"default.cfg"` |
| `test_iter` | Iteration whose row is checked | `1` |
| `test_vals` | Stored values, compared with the last columns of the row | `[]` |
| `test_vals_aarch64` | Stored values used instead on ARM64 (`aarch64`, `arm64`) | `[]` |
| `tol` | Allowed absolute difference: one number, or a list with one number per value | `0.0` |
| `tol_aarch64` | Tolerance used instead on ARM64 | `0.0` |
| `timeout` | Maximum run time in seconds | `0` |
| `command` | `TestCase.Command(launch, exec, param)`, e.g. `TestCase.Command("mpirun -n 2", "SU2_CFD")` | empty |
| `unsteady`, `multizone` | Choose which iteration option is changed (see above) | `False` |
| `no_restart` | Writes `RESTART_SOL= NO` in the copied configuration file | `False` |
| `enabled_on_cpu_arch` | CPU types on which the test runs | `["x86_64","amd64","aarch64","arm64"]` |
| `enabled_with_tsan`, `enabled_with_asan` | Run the test in the thread / address sanitizer jobs | `True` |
| `enabled_with_tapetests`, `tapetest_vals` | Run the test in the tape test job, and the expected number of tape errors | `False`, `[]` |
| `enabled_with_regular` | Set to `False` for a test that runs only under a sanitizer | `True` |

In `TestCase.py`, the defaults of `command`, `timeout` and `tol` are empty or zero. Each script fills them in a loop
before it runs its `test_list` (the block `RUN TESTS`). For example, `serial_regression.py` uses `SU2_CFD`,
1600 s and `0.00001`.

Some scripts overwrite values that you set, for all tests in `test_list`:

- `hybrid_regression.py` sets `tol = 1e-4`;
- `hybrid_regression_AD.py` sets the command `SU2_CFD_AD -t 2` and `timeout = 600`;
- `vandv.py` sets `timeout = 300` and `tol = 1e-5`.

### Other kinds of tests

`TestCase` has more methods, for tests that do not check the screen output of `SU2_CFD`:

- `run_filediff()`: runs the command and compares the file `test_file` with `reference_file`.
  - With `tol_file_percent = 0`, the two files must be equal.
  - Otherwise each number may differ by `tol_file_percent` percent. Numbers smaller than `comp_threshold` are not
    compared.
  - On ARM64 it uses `reference_file_aarch64`.
- `run_opt()`: for `shape_optimization.py`. It sets `OPT_ITERATIONS` to `test_iter`.
- `run_geo()`: for `SU2_GEO`.
- `run_def()`: for `SU2_DEF`.

The scripts call these methods one test at a time, after the main loop. The defaults loop has already run by then.
So tests added there must set `command`, `timeout` and `tol` themselves.

## How the CI runs the tests

The workflow is `.github/workflows/regression.yml`. It runs on every push to `develop` and `master`, and on every
pull request to them. `.github/workflows/regression-arm64.yml` runs the same workflow on an ARM64 runner, but only on
pushes to `develop`.

### Builds

The job **Build SU2** compiles SU2 in a container (`ghcr.io/su2code/su2/build-su2`), in several configurations:

- `Base` is the normal build.
- `Reverse` is the discrete adjoint build (`-Denable-autodiff=true`).
- `Forward` is the direct differentiation build (`-Denable-directdiff=true`).

Each one is built with MPI (`MPI`), without MPI (`NoMPI`) and with OpenMP (`OMP`). `ReverseTagNoMPI` is a discrete
adjoint build with `-Dcodi-tape=Tag`.

The jobs **Build SU2 (tsan)** and **Build SU2 (asan)** make the builds for the thread and address sanitizers. These
are compiler tools that detect data races between threads (tsan) and memory errors (asan) while the program runs.

### Test jobs

Each test job runs one script. It uses the `Base`, `Reverse` and `Forward` binaries of one build type (those that
exist):

| Job | Script | Binaries | Extra argument |
|---|---|---|---|
| Regression Tests | `serial_regression.py`, `serial_regression_AD.py` | `NoMPI` | |
| Regression Tests | `parallel_regression.py`, `parallel_regression_AD.py`, `tutorials.py`, `vandv.py` | `MPI` | |
| Regression Tests | `hybrid_regression.py`, `hybrid_regression_AD.py` | `OMP` | |
| Tape Tests | `serial_regression_AD.py` | `TagNoMPI` | `--tapetests` |
| Thread Sanitizer Tests | `hybrid_regression.py`, `hybrid_regression_AD.py` | `OMP-tsan` | `--tsan` |
| Address Sanitizer Tests | `serial_regression.py`, `serial_regression_AD.py` | `NoMPI-asan` | `--asan` |
| Unit Tests | `test_driver`, `test_driver_AD`, `test_driver_DD` | `MPI` | |

Details:

- The Tape Tests and Thread Sanitizer Tests jobs do not run on the ARM64 runner.
- `parse_args()` in `TestCase.py` reads the extra arguments. Only `serial_regression.py`, `serial_regression_AD.py`,
  `hybrid_regression.py` and `hybrid_regression_AD.py` read them.
- With a sanitizer, a test fails only if the run returns an error. The output values are not checked.
- The Unit Tests job with `test_driver` also runs `UnitTests/SU2_CFD/output/check_nemo_convergence_migration.py`.

### Where the meshes come from

The test jobs run the container `ghcr.io/su2code/su2/test-su2`. Its entry script is `test/runTests.sh` in the
repository `su2code/Docker-Builds`. The jobs call it with these arguments:

```
-b <ref of the SU2 commit> -t develop -c develop -s <script>
```

The script:

1. Clones SU2 at the given ref and copies its `TestCases/` folder to `tests/TestCases/`.
2. Clones the `develop` branch of `su2code/TestCases` (option `-c`) and copies it into `tests/TestCases/`.
3. Clones the `develop` branch of `su2code/Tutorials` (option `-t`) to `tests/Tutorials/`.
4. Sets `SU2_RUN`, `PATH` and `PYTHONPATH` to the installed binaries, and sets `SU2_MPI_COMMAND`.
5. Runs `python <script> <arguments>` in `tests/TestCases/`.

This is why `tutorials.py` finds the tutorial cases with paths such as `../Tutorials/compressible_flow/Inviscid_Bump`.

The CI always uses the `develop` branch of TestCases and Tutorials. If your pull request needs new or changed mesh
files, they must be in the `develop` branch of `su2code/TestCases` before the tests can pass. While that TestCases
pull request is open, some developers change `-c develop` in `regression.yml` to their TestCases branch. They change
it back before the SU2 pull request is merged.

## Run the tests on your computer

First, build and install SU2 with the executables that the script uses (see
[Build SU2 on Linux/MacOS](../Build-SU2-Linux-MacOS/)), and set `SU2_RUN`, `PATH` and `PYTHONPATH`. Some scripts
need special builds:

- The AD scripts need a build with `-Denable-autodiff=true`. Some of their tests also need `-Denable-directdiff=true`.
- The hybrid scripts need `-Dwith-omp=true`.
- The Python wrapper tests need `-Denable-pywrapper=true`.

These commands make the same folder layout as the CI, outside the SU2 source tree:

```
mkdir tests
cp -r SU2/TestCases tests/
git clone --depth=1 -b develop https://github.com/su2code/TestCases.git TestData
cp -R TestData/* tests/TestCases/
git clone --depth=1 -b develop https://github.com/su2code/Tutorials.git tests/Tutorials
cd tests/TestCases
python3 serial_regression.py
```

Each test writes its output files in its `cfg_dir`. If you change a script in `SU2/TestCases/`, copy it again to
`tests/TestCases/`. You can also run the containers of the CI on your computer, see
[Container Development](../Container-Development/).

### Run one test

The scripts have no option to choose one test. To run one test, put its entry in a small script in
`tests/TestCases/`. Set `command`, `timeout` and `tol` there, because the defaults of `TestCase` are empty or zero:

```python
import sys
from TestCase import TestCase

naca0012           = TestCase('naca0012')
naca0012.cfg_dir   = "euler/naca0012"
naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
naca0012.test_iter = 20
naca0012.test_vals = [-4.031446, -3.524530, 0.273056, 0.009703]
naca0012.command   = TestCase.Command(exec = "SU2_CFD")
naca0012.timeout   = 1600
naca0012.tol       = 0.00001
sys.exit(0 if naca0012.run_test() else 1)
```

Use the same `command` as the script that has the test, e.g. `TestCase.Command("mpirun -n 2", "SU2_CFD")` for
`parallel_regression.py`. You can also run SU2 by hand in `cfg_dir` and look at the row of iteration `test_iter`.

## Read a failure

To find why a test failed in the CI, open the failed job and the step **Run Tests in Container**. For every test,
the log has a block between `Start Test: <name>` and `End Test: <name>`. The end of the log has a summary, with one
line per test:

```
Summary of the serial tests
  passed - channel
* FAILED - naca0012
```

Search the log for `FAILED`. The block of a failed test prints the full screen output of SU2, the command, the
reason and the values (numbers left out here):

```
naca0012: FAILED
Output for the failed case
...
execution command: SU2_CFD inv_NACA0012_Roe.cfg.autotest > inv_NACA0012_Roe.cfg.log 2>&1
ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=1e-05
CPU architecture=x86_64
test_iter=20
test_vals (stored): ...
sim_vals (computed): ...
delta_vals: ...
```

The reason is one of these messages:

- `exceeded tolerance`: the results changed.
- `Execution timed out`: the test took longer than `timeout`.
- `The iteration number ... could not be found`: SU2 stopped before `test_iter`, or that row is not printed.
- `not able to get to the "Begin solver" section`: SU2 failed during the setup. Read the screen output above.
- `Error in test_vals!`: the number of values does not match the row.

If the summary has no failed test but the job failed, look for an error in the build or in the container step.

## Update the stored values

Update the values only when the change of the results is expected. Say why in the pull request.

1. In the CI log, find the block of each failed test.
2. Copy the numbers of `sim_vals (computed)` into `test_vals` of that test. They are printed with six decimals,
   like the stored values.
3. Search the other scripts for the same case. Some cases are in more than one script with the same values, e.g.
   `naca0012` and `pywrapper_naca0012` in `serial_regression.py`.
4. For a `run_filediff()` test, replace the reference file. Some reference files are in the SU2 repository next to
   the configuration file. Others are in the TestCases repository. `*.ref` is in the SU2 `.gitignore`, so add a
   changed reference file with `git add -f`.
5. Push and check that the CI passes.

Use the values computed by the CI, because your computer can give slightly different numbers.
`test_vals_aarch64` and `reference_file_aarch64` are used only on ARM64. The ARM64 workflow runs after the merge into
`develop`, not on the pull request.

## Add a new test

1. Make the case small. A test should run in seconds or a few minutes. Most tests check iteration 10 or 20.
2. Put the configuration file in a folder of `TestCases/` in SU2, e.g. `TestCases/rans/my_case/`.
3. Put the mesh and other large files in the same folder of the TestCases repository. Open a pull request there
   against `develop` (see [Where the meshes come from](#where-the-meshes-come-from)).
4. In the configuration file, set `SCREEN_OUTPUT` to the values you want to check. Make sure that the row of
   `test_iter` is printed (see `SCREEN_WRT_FREQ_INNER`).
5. Add a `TestCase` entry to the script that matches the build (serial, parallel, AD, hybrid), before its block
   `RUN TESTS`. Add it to `test_list` with `test_list.append(...)`. Set at least `cfg_dir`, `cfg_file`,
   `test_iter` and `test_vals`.
6. To get the first values, set `test_vals` to a list of zeros, with one entry per column you check. Run the test,
   and copy `sim_vals (computed)` from the output. Then check the values with the CI, as explained above.

Last checked against SU2 develop `11c177b904`.
