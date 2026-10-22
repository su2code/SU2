---
title: Automatic Differentiation (AD)
permalink: /docs_v7/Automatic-Differentiation/
---

SU2 computes derivatives of its results with automatic differentiation (AD), also called algorithmic
differentiation. This page explains how AD works in the code, and what you must do so that new code gives correct
derivatives.

This page covers:

- [What AD is in SU2](#what-ad-is-in-su2)
- [Building the AD executables](#building-the-ad-executables)
- [How the discrete adjoint runs](#how-the-discrete-adjoint-runs)
- [Forward mode (direct differentiation)](#forward-mode-direct-differentiation)
- [Writing differentiable code](#writing-differentiable-code)
- [Debugging AD](#debugging-ad)
- [Testing gradients](#testing-gradients)
- [Checklist for a new feature](#checklist-for-a-new-feature)

For a short overview, see [Code Structure](../Code-Structure/#algorithmic-differentiation). For the theory of
expression templates, see [Advanced AD Techniques](../Advanced-AD-Techniques/).

## What AD is in SU2

AD applies the chain rule to every operation that the program runs, so the derivatives are exact up to rounding.
SU2 uses AD in two modes. The forward mode (SU2 calls it "direct differentiation") carries one derivative along
with each value, and gives the derivative of all outputs with respect to one input per run. The reverse mode gives
the derivative of one output with respect to all inputs at once; SU2 uses it for the discrete adjoint.
SU2 gets AD from [CoDiPack](https://github.com/SciCompKL/CoDiPack), an operator-overloading tool: the type
`su2double` replaces `double`, and its operators also record derivative information. In reverse mode this
information goes to the tape: a stored list of all operations and their local derivatives. CoDiPack later replays
the tape backwards to get the adjoints, the derivatives of the output with respect to each recorded value.

The type `su2double` is chosen at compile time in `Common/include/code_config.hpp`:

```cpp
#if defined(CODI_REVERSE_TYPE)  // reverse mode AD
#include "codi.hpp"
// ...
#if defined(HAVE_OMP)
using su2double = codi::RealReverseIndexOpenMPGen<double, double>;
#else
#if defined(CODI_JACOBIAN_LINEAR_TAPE)
using su2double = codi::RealReverse;
// ...
#elif defined(CODI_TAG_TAPE)
using su2double = codi::RealReverseTag;
// ...
#elif defined(CODI_FORWARD_TYPE)  // forward mode AD
#include "codi.hpp"
using su2double = codi::RealForward;
#else  // primal / direct / no AD
using su2double = passivedouble;
#endif
```

`passivedouble` is always a plain `double` (or `float` in a full single-precision build). Use it for values that
must never carry derivatives.

Two more libraries complete the AD support:

- [MeDiPack](https://github.com/SciCompKL/MeDiPack) (`externals/medi`) differentiates MPI communication.
- [OpDiLib](https://github.com/SciCompKL/OpDiLib) (`externals/opdi`) differentiates OpenMP parallel regions.

The `meson.py` script downloads CoDiPack, MeDiPack and OpDiLib at the commits set in `meson_scripts/init.py`.

## Building the AD executables

The general build steps are on the [Build SU2 on Linux/MacOS](../Build-SU2-Linux-MacOS/) page. These meson options
(`meson_options.txt`) control AD:

| Option | Default | Effect |
|---|---|---|
| `-Denable-autodiff=true` | `false` | Builds the reverse mode executables `SU2_CFD_AD` and `SU2_DOT_AD` (and the Python module `pysu2ad` with `-Denable-pywrapper=true`). Defines `CODI_REVERSE_TYPE`. |
| `-Denable-directdiff=true` | `false` | Builds the forward mode executable `SU2_CFD_DIRECTDIFF`. Defines `CODI_FORWARD_TYPE`. |
| `-Denable-normal=false` | `true` | Skips the normal executables (`SU2_CFD`, `SU2_DEF`, `SU2_DOT`, ...). Useful to build only the AD versions. |
| `-Dcodi-tape=<type>` | `JacobianLinear` | The CoDiPack tape of the reverse build, see below. |
| `-Dopdi-backend=auto\|macro\|ompt` | `auto` | How OpDiLib follows the OpenMP regions (reverse build with OpenMP only). |
| `-Dopdi-shared-read-opt=false` | `true` | Turns off the faster "no shared reading" sections, see [OpenMP](#openmp). |

There is no `SU2_DEF_AD` executable. The AD build of `SU2_DEF` is only a library (`SU2core_DEF_AD`), which
`SU2_DOT_AD` and `pysu2ad` link.

A typical reverse build, without the normal executables:

```
./meson.py setup build_ad -Denable-autodiff=true -Denable-normal=false --prefix=$PWD/install
./ninja -C build_ad install
```

`SU2_CFD_AD` runs the discrete adjoint by default: in a reverse build, `MATH_PROBLEM` defaults to
`DISCRETE_ADJOINT`. The normal `SU2_CFD` stops with an error if you ask for `MATH_PROBLEM= DISCRETE_ADJOINT` or
`DIRECT_DIFF`.

### Tape types

The option `-Dcodi-tape` selects the reverse type. Each value defines one macro in `meson.build`, and
`code_config.hpp` turns it into a CoDiPack type:

| `-Dcodi-tape=` | `su2double` is | What it is for |
|---|---|---|
| `JacobianLinear` | `codi::RealReverse` | The default, and the only tape tested regularly. It stores the partial derivatives of each statement. Identifiers (the indices of the adjoint values) grow with every statement. |
| `JacobianReuse` | `codi::RealReverseIndexGen<double, double, codi::ReuseIndexManager<int>>` | Same derivative data, but identifiers of deleted variables are reused. The adjoint vector is smaller. |
| `JacobianMultiUse` | `codi::RealReverseIndex` | Like `JacobianReuse`, and a copy shares the identifier of the original instead of adding a statement. |
| `PrimalLinear`, `PrimalReuse`, `PrimalMultiUse` | `codi::RealReversePrimal`, `codi::RealReversePrimalIndexGen<...>`, `codi::RealReversePrimalIndex` | Primal value tapes. They store the values and the operations, and compute the partial derivatives during the reverse sweep. |
| `Tag` | `codi::RealReverseTag` | Computes no derivatives. It gives each variable a tag, to find errors in the recording. See [Debugging AD](#debugging-ad). |

Notes:

- `meson.build` prints a warning for every tape except `JacobianLinear`.
- With `-Dwith-omp=true`, `meson.build` ignores `-Dcodi-tape`. The type is always the thread-safe
  `codi::RealReverseIndexOpenMPGen<double, double>`, which reuses identifiers. So the `Tag` tape needs a build
  without OpenMP.
- The tapes that reuse identifiers define `CODI_INDEX_REUSE`. With them, `CDiscAdjMultizoneDriver::ComputeAdjoints`
  prints "WARNING: AD types that reuse indices do not support multiple zones."
- With `--buildtype=debug` or `--buildtype=debugoptimized`, meson also defines `CODI_EnableAssert`, which turns on
  the internal checks of CoDiPack.
- With MPI, `SU2_MPI` becomes `CMediMPIWrapper`, built on MeDiPack. With OpenMP, a configure step runs
  `externals/opdi/syntax/check.py` with `su2omp.syntax.json`. It stops the configuration if an `SU2_OMP_*` macro has
  no matching `END_SU2_OMP_*`.

The CI builds `Reverse` and `Forward` versions with MPI, without MPI and with OpenMP, plus `ReverseTagNoMPI`
(`-Dcodi-tape=Tag`). See [Running Regression Tests](../Running-Regression-Tests/#builds).

## How the discrete adjoint runs

The discrete adjoint solves a fixed-point problem. One primal iteration maps the solution `U` to a new solution
`G(U)`. At the converged solution, the adjoint `Ū` satisfies

`Ū = (∂G/∂U)^T Ū + (∂J/∂U)^T`

where `J` is the objective function. SU2 records one primal iteration on the tape. Each evaluation of the tape
gives `(∂G/∂U)^T Ū + (∂J/∂U)^T` for the current `Ū`, so repeated evaluations are a fixed-point iteration for `Ū`.
No hand-written adjoint equations are needed.

The drivers are `CDiscAdjSinglezoneDriver` and `CDiscAdjMultizoneDriver` (`SU2_CFD/src/drivers`). The single-zone
driver is described first.

### Recording the tape

`CDiscAdjSinglezoneDriver::SetRecording` records one iteration. Shortened:

```cpp
void CDiscAdjSinglezoneDriver::SetRecording(RECORDING kind_recording){
  AD::Reset();

  /*--- Prepare for recording by resetting the solution to the initial converged solution. ---*/
  // ... solver->SetRecording(geometry, config) for every adjoint solver ...

  /*---Enable recording and register input of the iteration --- */
  if (kind_recording != RECORDING::CLEAR_INDICES){
    AD::StartRecording();
    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);
  }

  /*--- Set the dependencies of the iteration ---*/
  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0,
                             INST_0, kind_recording);

  /*--- Do one iteration of the direct solver ---*/
  DirectRun(kind_recording);
  // ...
  /*--- Register Output of the iteration ---*/
  iteration->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Extract the objective function and store it --- */
  SetObjFunction();
  // ...
  AD::StopRecording();
}
```

The steps are:

1. `AD::Reset()` deletes the old tape. `CDiscAdjSolver::SetRecording` copies the converged primal solution
   (`GetSolution_Direct`) back into the primal solver. Then it removes the old tape identifiers from the solution
   and from the other registered variables (`AD::ResetInput`).
2. `AD::StartRecording()` switches the tape on, and `RegisterInput` marks the inputs. What is registered depends
   on the `RECORDING` value (`Common/include/option_structure.hpp`):

   | `RECORDING::` | Inputs registered by `CDiscAdjFluidIteration::RegisterInput` |
   |---|---|
   | `CLEAR_INDICES` | Nothing. The tape stays off; this pass only removes old identifiers. |
   | `SOLUTION_VARIABLES` | The solution of each adjoint solver (`RegisterSolution`), the solutions at time `n` and `n-1` for dual time stepping, extra variables (`RegisterSolutionExtra`), and the free-stream or boundary values (`RegisterVariables`). |
   | `MESH_COORDS` | The mesh coordinates (`CGeometry::RegisterCoordinates`). |
   | `MESH_DEFORM` | The undeformed mesh and the boundary displacements, through the adjoint mesh solver `ADJMESH_SOL`. Used when the mesh is deformed (`DEFORM_MESH= YES`). |
   | `SOLUTION_AND_MESH` | Both the solution and the coordinates. |

3. `SetDependencies` recomputes the values that depend on the inputs but are normally computed only once or
   outside the iteration: the geometry (normals, volumes, wall distance) when the coordinates are inputs, and the
   primitive variables, eddy viscosity and halo values. Without this step the tape would miss these dependencies.
4. `DirectRun` runs one primal iteration (`Preprocess`, `Iterate`, `Postprocess` of the direct iteration). Every
   operation on active `su2double` values goes to the tape.
5. `RegisterOutput` marks the new solution as output, and `SetObjFunction` computes `ObjFunc` and registers it
   (only on the master rank).

`MainRecording` records with the solution as input (`RECORDING::SOLUTION_VARIABLES`).
`SecondaryRecording` records with the mesh as input (`MESH_COORDS` or `MESH_DEFORM`). Both first call
`SetRecording(RECORDING::CLEAR_INDICES)`.

The adjoint solver keeps the tape identifiers of the inputs and outputs, so it can later read and write their
adjoints. `CVariable::RegisterContainer` (`SU2_CFD/include/variables/CVariable.hpp`) registers a whole matrix and
stores the identifiers:

```cpp
void RegisterContainer(bool input, su2activematrix& variable, su2matrix<AD::Identifier>* ad_index = nullptr) {
  const auto nPoint = variable.rows();
  SU2_OMP_FOR_STAT(roundUpDiv(nPoint,omp_get_num_threads()))
  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {
    for(unsigned long iVar=0; iVar<variable.cols(); ++iVar) {

      if (input) AD::RegisterInput(variable(iPoint,iVar));
      else AD::RegisterOutput(variable(iPoint,iVar));

      if (ad_index) AD::SetIndex((*ad_index)(iPoint,iVar), variable(iPoint,iVar));
    }
  }
  END_SU2_OMP_FOR
}
```

`CVariable::RegisterSolution(bool input)` calls it with `AD_InputIndex` or `AD_OutputIndex`.

### The fixed-point iteration

`CDiscAdjSinglezoneDriver::Run` evaluates the tape once per adjoint iteration:

```cpp
for (auto Adjoint_Iter = 0ul; Adjoint_Iter < nAdjoint_Iter; Adjoint_Iter++) {
  // ...
  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/
  SetAdjObjFunction();

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/
  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->IterateDiscAdj(geometry_container, solver_container,
                            config_container, ZONE_0, INST_0, false);
  // ... Monitor ...
  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/
  AD::ClearAdjoints();
  // ... Output, quasi-Newton correction ...
}
```

- `InitializeAdjoint` calls `CDiscAdjSolver::SetAdjoint_Output`. It writes the current adjoint solution `Ū` into
  the adjoints of the outputs (`CVariable::SetAdjointSolution`, which uses `AD::SetDerivative` with the stored
  `AD_OutputIndex`).
- `SetAdjObjFunction` seeds the adjoint of `ObjFunc` with 1 (with a time-window weight for unsteady problems).
- `AD::ComputeAdjoint()` evaluates the tape backwards.
- `IterateDiscAdj` calls `CDiscAdjSolver::ExtractAdjoint_Solution`. It reads the adjoints of the inputs
  (`CVariable::GetAdjointSolution`, which uses `AD::GetDerivative` with `AD_InputIndex`), applies the relaxation
  factor and computes the adjoint residuals. `ExtractAdjoint_Variables` reads the sensitivities to the registered
  free-stream values (Mach, angle of attack, pressure, temperature, ...).
- `AD::ClearAdjoints()` sets all adjoints to zero, but keeps the tape.

For a steady problem the tape is recorded once. `Preprocess` records again only when the last recording was of
another kind, for example after a sensitivity recording.

### Extracting sensitivities

After the loop, `Postprocess` calls `SecondaryRecording`. It records again with the mesh as input, seeds the
outputs with the converged adjoint, evaluates the tape once and extracts the mesh sensitivities:

```cpp
if (SecondaryVariables == RECORDING::MESH_COORDS) {
  solver[MainSolver]->SetSensitivity(geometry, config);
}
else { // MESH_DEFORM
  solver[ADJMESH_SOL]->SetSensitivity(geometry, config, solver[MainSolver]);
}
```

`CDiscAdjSolver::SetSensitivity` reads the adjoint of each coordinate and stores it as the sensitivity.
`SU2_DOT_AD` then projects the surface sensitivities onto the design variables.
`CDiscAdjDeformationDriver::SetProjection_AD` (`SU2_DEF/src/drivers`) records the surface deformation with the
design variables as inputs (`AD::RegisterInput(config->GetDV_Value(iDV, iDV_Value))`), seeds the moved surface
coordinates with the sensitivities, and runs `AD::ComputeAdjoint()`.

### Free-stream sensitivities

`CDiscAdjSolver::RegisterVariables` shows how to add a scalar input. It copies the value from the config into a
member, registers the member, and writes it back so that the primal code uses the registered copy:

```cpp
if (!reset) {
  AD::RegisterInput(Mach);
  AD::RegisterInput(Alpha);
  AD::RegisterInput(Temperature);
  AD::RegisterInput(Pressure);
}

/*--- Set the free stream velocity (now using the registered values for Alpha, Beta and Mach). ---*/
if (nDim == 2) {
  config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
  config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
}
// ...
```

`ExtractAdjoint_Variables` then reads `SU2_TYPE::GetDerivative(Mach)` and sums it over the ranks.

### Multizone problems

`CDiscAdjMultizoneDriver` records all zones on one tape and marks positions in it with `AD::Push_TapePosition()`
(`START`, `REGISTERED`, `DEPENDENCIES`, `OBJECTIVE_FUNCTION`, `TRANSFER`, then one pair per zone). It evaluates
only parts of the tape with `AD::ComputeAdjoint(enter, leave)`:

```cpp
const unsigned short enter_izone = iZone*2+1 + ITERATION_READY;
const unsigned short leave_izone = iZone*2 + ITERATION_READY;

AD::ComputeAdjoint(enter_izone, leave_izone);
// ...
if (eval_transfer || config_container[iZone]->GetStructuralProblem())
  AD::ComputeAdjoint(TRANSFER, OBJECTIVE_FUNCTION);
// ...
AD::ComputeAdjoint(DEPENDENCIES, START);
```

It uses two kinds of tape (`Kind_Tape`): `FULL_SOLVER_TAPE` for the iterations, and `OBJECTIVE_FUNCTION_TAPE`,
which contains only the objective function. It solves the coupled adjoint with a block Gauss-Seidel loop over the
zones, with optional quasi-Newton or Krylov acceleration. A single-zone case also runs with this driver when its
configuration sets `MULTIZONE= YES`.

## Forward mode (direct differentiation)

`SU2_CFD_DIRECTDIFF` carries one derivative with each `su2double`. There is no tape. You choose the input with
`DIRECT_DIFF` (`DirectDiff_Var_Map`): `MACH`, `AOA`, `SIDESLIP`, `REYNOLDS`, `TURB2LAM`, `PRESSURE`,
`TEMPERATURE`, `DENSITY`, `VISCOSITY`, `DESIGN_VARIABLES`, and structural values (`YOUNG_MODULUS`,
`POISSON_RATIO`, `STRUCTURAL_DENSITY`, `STRUCTURAL_DEAD_LOAD`, `ELECTRIC_FIELD`).

`CConfig` seeds some of them, for example `SU2_TYPE::SetDerivative(Mach, 1.0)`. The solvers seed the others.
For `DESIGN_VARIABLES`, `CSurfaceMovement::SetSurface_Derivative` seeds every design variable with a non-zero value.
So run one design variable at a time. The history output then has a column `D_<field>` for every coefficient
field, filled with `SU2_TYPE::GetDerivative`. `direct_differentiation.py` runs one `SU2_CFD_DIRECTDIFF` per design
variable.

In forward mode, the `AD::` functions of the reverse mode (recording, preaccumulation, ...) do nothing.

## Writing differentiable code

Most code needs no AD work: if it computes with `su2double`, CoDiPack records it. You must take care of the
cases below.

### The AD API

All functions are in the namespace [`AD`](../../doxygen/namespaceAD.html) (`Common/include/basic_types/ad_structure.hpp`)
or `SU2_TYPE` (`Common/include/basic_types/datatype_structure.hpp`). In a build without reverse AD, the `AD::`
functions are empty, so you can call them in any build.

| Function | Use |
|---|---|
| `AD::StartRecording()`, `AD::StopRecording()`, `AD::TapeActive()` | Switch the tape on and off; ask if it is on. |
| `AD::RegisterInput(x)`, `AD::RegisterOutput(y)` | Mark the inputs and outputs of the recording. |
| `AD::ComputeAdjoint()`, `AD::ComputeAdjoint(enter, leave)` | Evaluate the whole tape, or the part between two positions stored with `AD::Push_TapePosition()`. |
| `AD::ClearAdjoints()`, `AD::Reset()` | Set all adjoints to zero; delete the tape. |
| `SU2_TYPE::SetDerivative(x, v)`, `SU2_TYPE::GetDerivative(x)` | Set or read the adjoint (reverse) or the tangent (forward) of a variable. |
| `AD::SetIndex(index, x)`, `AD::GetDerivative(index)`, `AD::SetDerivative(index, v)` | Store the identifier of `x`, then read or write its adjoint later by identifier. Use `AD::ResizeAdjoints()` first and wrap the calls in `AD::BeginUseAdjoints()` / `AD::EndUseAdjoints()`: these functions do no bounds checks and no locking. |
| `AD::ResetInput(x)` | Make `x` passive again (keeps the value, drops the identifier). |
| `AD::StartPreacc()`, `AD::SetPreaccIn(...)`, `AD::SetPreaccOut(...)`, `AD::EndPreacc()` | Preaccumulation, see below. |
| `AD::PausePreaccumulation()`, `AD::ResumePreaccumulation(wasActive)` | Switch preaccumulation off for a code section. |
| `AD::BeginPassive()`, `AD::EndPassive(wasActive)`; `AD::PauseRecording()`, `AD::ResumeRecording(wasActive)` | Run a section without recording. The two pairs do the same thing. |
| `AD::StartNoSharedReading()`, `AD::EndNoSharedReading()` | Mark an OpenMP section where no two threads read the same active value. |
| `AD::SetExtFuncIn`, `AD::SetExtFuncOut`, `AD::FuncHelper` | External functions: code that is not recorded, with a hand-written reverse function. |
| `AD::PrintStatistics(comm, printingRank)` | Print the tape size. The drivers call it when `WRT_AD_STATISTICS= YES`. |
| `SU2_TYPE::GetValue(x)`, `SU2_TYPE::PassiveValue(x)` | The plain value of `x`, as `passivedouble`. |

The unit test `UnitTests/Common/simple_ad_test.cpp` shows the basic sequence:

```cpp
su2double x = 4.0;

AD::StartRecording();
AD::RegisterInput(x);

su2double y = func(x);

AD::RegisterOutput(y);
AD::StopRecording();
SU2_TYPE::SetDerivative(y, 1.0);
AD::ComputeAdjoint();

CHECK(SU2_TYPE::GetValue(y) == Approx(64));
CHECK(SU2_TYPE::GetDerivative(x) == Approx(48));
```

### Preaccumulation

The tape stores every operation, so it is large: often several times the memory of the primal solver.
Preaccumulation reduces it. For a code section with few inputs and outputs but many operations, CoDiPack computes
the local Jacobian (outputs with respect to inputs) at the end of the section. It then replaces the recorded
operations with this Jacobian. The derivatives do not change.

`CFVMFlowSolverBase::ComputeVorticityAndStrainMag` (`SU2_CFD/include/solvers/CFVMFlowSolverBase.inl`) is a small
example:

```cpp
const su2double vy = nodes->GetVelocity(iPoint, 1);
const su2double y = geometry->nodes->GetCoord(iPoint, 1);
AD::StartPreacc();
AD::SetPreaccIn(VelocityGradient, nDim, nDim);
AD::SetPreaccIn(vy, y);

StrainMag(iPoint) = 0.0;

/*--- Add diagonal part ---*/

for (unsigned long iDim = 0; iDim < nDim; iDim++) {
  StrainMag(iPoint) += pow(VelocityGradient(iDim, iDim), 2);
}
if (config.GetAxisymmetric() && y > EPS) {
  StrainMag(iPoint) += pow(vy / y, 2);
}
// ... off-diagonal terms ...
StrainMag(iPoint) = sqrt(2.0*StrainMag(iPoint));
AD::SetPreaccOut(StrainMag(iPoint));
AD::EndPreacc();
```

Rules:

- Declare **every** active value that the outputs depend on with `AD::SetPreaccIn`. A value that is not declared
  is treated as a constant. The gradient is then wrong, and there is no error message. Inputs can be declared
  anywhere between `StartPreacc` and their first use, as `CSourcePieceWise_TurbSST::ComputeResidual` does for the
  intermittency.
- Declare every output with `AD::SetPreaccOut` and close the section with `AD::EndPreacc()` on **every** return
  path. `CUpwRoeBase_Flow::ComputeResidual` does this before its early return, see
  [Adding a Numerical Scheme](../Adding-a-Numerical-Scheme/#discrete-adjoint).
- `SetPreaccIn` takes lvalues of type `su2double`, single values, 1D arrays (`SetPreaccIn(data, size)`) or 2D
  arrays (`SetPreaccIn(data, size_x, size_y)`). An expression or a temporary does not compile ("rvalues cannot be
  registered"). Store it in a local `su2double` first, as `vy` and `y` above.
- Do not nest sections, and do not call code that starts its own section.
- The section is skipped when the tape is off, or when `PREACC= NO`. So the code must give the same values with
  and without preaccumulation.
- Preaccumulation helps when the section has many operations per input and output. Use as few inputs as possible.
  `CPhysicalGeometry` recomputes intermediate values inside the section for this reason.
- In the vectorized numerics (`SU2_CFD/include/numerics_simd`), `gatherVariables` declares the inputs and
  `stopPreacc` (`SU2_CFD/include/numerics/util.hpp`) ends the section.

### Passive sections

Some values must not be differentiated, for example the Jacobian of the implicit scheme, a monitored maximum, or
a frozen limiter. Wrap them in a passive section, so that they cost no tape memory:

```cpp
/*--- The derivative with respect to strainMax and omegaMax is not required. ---*/
bool wa = AD::PauseRecording();
strainMax = max(strainMax, StrainMag(iPoint));
omegaMax = max(omegaMax, GeometryToolbox::Norm(3, Vorticity));
AD::ResumeRecording(wa);
```

A value computed in a passive section is a constant on the tape. Use this only when that is what you want.

### Values, branches and iteration counters

- Do not use `SU2_TYPE::GetValue` (or a cast to `passivedouble`, `int`, ...) on a value that the result depends on.
  It cuts the dependency. Use it only for output, for logic that has no derivative (indices, counters), or for
  values that must be constant on purpose.
- A branch on an active value is recorded as the branch that was taken. The derivative is the derivative of that
  branch. This is correct away from the switch point. At the switch point (`max`, `min`, `fabs`, limiters) the
  derivative is one-sided. Avoid branches whose outcome changes between the converged primal and the recording,
  because then the tape describes another function.
- The tape is recorded once, at one inner iteration. Code that runs only at some iterations
  (`config->GetInnerIter() == 0`, every n iterations, only when a residual is small, ...) must also run during the
  recording if it affects the result. For example, `CIncEulerSolver::Preprocessing`:

  ```cpp
  /*--- InnerIter is not reset while recording the discrete adjoint tape. ---*/
  if (dual_time && (config->GetInnerIter() == 0 || AD::TapeActive())) {
    RecomputeDensity_time_n(solver_container, config);
  }
  ```

- A value that the iteration uses must be computed from the inputs during the recording. If it is computed
  earlier (in a postprocessing step of the previous iteration, at start-up, ...), the tape sees it as a constant.
  Recompute it in `SetDependencies`, or register it with `RegisterSolutionExtra`, as `CTurbSolver` does for the
  eddy viscosity.
- `static` and global `su2double` variables survive between recordings. Their tape identifiers then point to an old
  tape. Avoid them, or reset them with `AD::ResetInput`.

### MPI

Use the `SU2_MPI` wrapper (`Common/include/parallelization/mpi_structure.hpp`). In AD builds it is
`CMediMPIWrapper`, and `MPI_DOUBLE` means `su2double`: MeDiPack records the communication. Other MPI types
(`MPI_INT`, `MPI_UNSIGNED_LONG`, ...) are passed through. A type that the wrapper does not know stops with
"Conversion not implemented".

Do not send `passivedouble` buffers through `SU2_MPI` with `MPI_DOUBLE`. Use `SelectMPIWrapper<passivedouble>::W`,
which in the reverse build is the plain MPI wrapper:

```cpp
using CPassiveMPI = SelectMPIWrapper<passivedouble>::W;
vector<passivedouble> dirOut(nRecvTotal * nDim, 0.0), dirIn(nSendTotal * nDim, 0.0);
// ...
CPassiveMPI::Sendrecv(dirOut.data() + hp.offR * nDim, hp.nVertexR * nDim, MPI_DOUBLE, hp.receive_from, 1,
                      dirIn.data() + hp.offS * nDim, hp.nVertexS * nDim, MPI_DOUBLE, hp.send_to, 1,
                      SU2_MPI::GetComm(), MPI_STATUS_IGNORE);
```

Halo updates of solver data (`InitiateComms` / `CompleteComms`) are recorded like any other operation. Halo
points are local copies of points owned by other MPI ranks.

### OpenMP

In a reverse build with OpenMP, the `SU2_OMP_*` macros of `omp_structure.hpp` become OpDiLib macros, which record
the parallel structure for the reverse sweep. Rules:

- Use only the `SU2_OMP_*` macros, never `#pragma omp` directly, and close each one with its `END_SU2_OMP_*` macro.
- In the reverse sweep, the adjoint of a value is updated by every thread that read the value in the primal run.
  By default OpDiLib makes these updates atomic, which is safe but slow.
- Between `AD::StartNoSharedReading()` and `AD::EndNoSharedReading()`, OpDiLib uses plain updates. Use these
  sections only where no two threads read the same active value, for example a loop over points that reads only
  data of its own point. `CTurbSASolver::Postprocessing` puts the eddy viscosity loop in such a section.
- `SU2_OMP_SAFE_GLOBAL_ACCESS(...)` runs code on the master thread, with barriers. It also opens a no-shared-reading
  section, because the other threads wait.
- Preaccumulation is thread-local (`PreaccHelper` and `PreaccActive` are `threadprivate`). But in an edge loop
  whose threads share points, it is not safe. This happens with the reducer strategy, which SU2 uses when edge
  coloring is not efficient: the fluxes are stored per edge and summed per point afterwards.
  `CScalarSolver::EdgeFluxResidual` therefore pauses preaccumulation in that case:

  ```cpp
  bool pausePreacc = false;
  if (ReducerStrategy) pausePreacc = AD::PausePreaccumulation();
  else AD::StartNoSharedReading();
  ```

- `SU2_OMP_SIMD_IF_NOT_AD` is a SIMD pragma only in builds without AD. Use it on loops over `su2double`.
- A solver that supports hybrid parallel adjoints returns `true` from `CSolver::GetHasHybridParallel()`. The
  adjoint iterations then open their parallel regions with `SU2_OMP_PARALLEL_(if(...GetHasHybridParallel()))`.

### External functions

Some algorithms are cheaper to differentiate by hand than by recording, mainly linear solvers. `CSysSolve::Solve`
(`Common/src/linear_algebra/CSysSolve.cpp`) declares the right-hand side as input and the solution as output, runs
the primal solve without recording, and adds a reverse function to the tape:

```cpp
BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
  AD::SetExtFuncIn(&LinSysRes[0], LinSysRes.GetLocSize());
  AD::SetExtFuncOut(&LinSysSol[0], LinSysSol.GetLocSize());
  AD::FuncHelper.addUserData(&LinSysRes_tmp);
  // ...
}
END_SU2_OMP_SAFE_GLOBAL_ACCESS
// ...
AD::FuncHelper.callPrimalFuncWithADType(externalFunction);

AD::FuncHelper.addToTape(CSysSolve_b<ScalarType>::Solve_b);
```

`CSysSolve_b::Solve_b` solves the transposed system in the reverse sweep.

### Memory

The tape grows with the number of active operations in one iteration. To keep it small:

- Use `passivedouble` for values that never need derivatives (counters, timers, geometric search structures).
- Preaccumulate small functions with many operations.
- Put code that must not be differentiated in a passive section.
- Check the effect with `WRT_AD_STATISTICS= YES`. The drivers then print the tape statistics after each recording,
  summed over MPI ranks and OpenMP threads.

## Debugging AD

### Tape debug build

A build with `-Dcodi-tape=Tag` (and without OpenMP) uses the tag tape. It computes no derivatives. It gives every
variable a tag and reports when a statement uses a variable with an unexpected tag.

The check runs only in `CDiscAdjMultizoneDriver`, so the configuration needs `MULTIZONE= YES` (one zone is fine).
In a tag build, `CConfig` sets the debug mode, and `CDiscAdjMultizoneDriver::StartSolver` calls `TapeTest` instead
of the solver. `TapeTest` records the tape twice:

1. The **initial** recording gives the inputs the tag 1 (in multizone problems, `10*(iZone+1) + 1`). Errors in
   this recording are mostly preaccumulation errors. Inside a preaccumulation section the current tag is 1337.
2. The **check** recording repeats the recording with the tag 2. A variable that still has the tag 1 was not
   recomputed from the inputs in this recording. It was computed before the inputs were registered, for example
   in a postprocessing step, and the tape treats it as a constant.

Two options select what is checked:

| Option | Values | Meaning |
|---|---|---|
| `CHECK_TAPE_TYPE` | `FULL_SOLVER_TAPE` (default), `OBJECTIVE_FUNCTION_TAPE` | Check the tape of the whole iteration, or only the tape of the objective function. |
| `CHECK_TAPE_VARIABLES` | `SOLVER_VARIABLES` (default), `MESH_COORDINATES` | Register the solution, or the mesh coordinates (plus the solution) as inputs. |

Each rank writes its errors to `debug_run_process<rank>.out`. The screen shows the total, "Total number of tape
inconsistencies: N". A typical line is:

```
Use of variable with bad tag '01', should be '1337'.
```

How to read it:

- `should be '1337'`: a value used inside a preaccumulation section is not declared with `AD::SetPreaccIn`.
- `bad tag '1337'`: a value computed inside a preaccumulation section is used later, but it is not declared with
  `AD::SetPreaccOut`.
- `bad tag '01', should be '02'` (in the check recording): a value is used that was not computed from the inputs
  in this recording.

The file does not say where the error is. Use a debugger with a breakpoint in `AD::tagErrorCallback`
(`ad_structure.hpp`) to get the call stack.

The functions `AD::SetTag`, `AD::ComputeTag`, `AD::SetTapeDebugOption` (with a `AD::TAPE_DEBUG_OPTION` value) and
`AD::SetTagOnVariable` control the tags. `CInterface` uses `TAPE_DEBUG_OPTION::ALLOW_ZONE` to allow values of the
donor zone during the transfer between zones.

### Tape tests in the regression scripts

`TestCase.py` has the argument `--tapetests`. It runs only the tests with `enabled_with_tapetests = True`, and
compares the number of tape inconsistencies with `tapetest_vals`:

```python
discadj_naca0012_via_mz = TestCase('discadj_naca0012_via_mz')
discadj_naca0012_via_mz.cfg_dir = "cont_adj_euler/naca0012"
discadj_naca0012_via_mz.cfg_file = "inv_NACA0012_discadj_multizone.cfg"
# ...
discadj_naca0012_via_mz.enabled_with_tapetests = True
discadj_naca0012_via_mz.tapetest_vals = [0]
```

The CI job "Tape Tests" runs `serial_regression_AD.py --tapetests` with the `ReverseTagNoMPI` build. To run it
yourself, build `SU2_CFD_AD` with `-Dcodi-tape=Tag -Dwith-mpi=disabled` and run
`python3 serial_regression_AD.py --tapetests`. See [Running Regression Tests](../Running-Regression-Tests/).

### Other checks

- **`PREACC= NO`.** Preaccumulation does not change the derivatives. If the gradient changes when you turn it off,
  a preaccumulation section is wrong.
- **Assertions.** A debug build defines `CODI_EnableAssert` and enables the `assert` checks of SU2.
- **Thread sanitizer.** The CI runs `hybrid_regression_AD.py --tsan` with an OpenMP reverse build, to find data
  races in the parallel reverse sweep.

### Common errors

| Symptom | Likely cause |
|---|---|
| Gradient is zero for an input | The input is not registered, or the dependency is cut by `SU2_TYPE::GetValue`, a cast, a passive section or `passivedouble`. |
| Gradient is wrong, no error | A missing `AD::SetPreaccIn`, a value computed outside the recording, or a branch that differs between primal and recording. Run the tape test and compare with finite differences. |
| Adjoint diverges, primal converges | Often the same causes. Also check that the primal is well converged, and try a smaller `CFL_NUMBER` or the quasi-Newton option `QUASI_NEWTON_NUM_SAMPLES`. |
| Compile error "rvalues cannot be registered" | An expression passed to `AD::SetPreaccIn`. Store it in a variable first. |
| "Conversion not implemented" | An MPI type that `CMediMPIWrapper::convertDatatype` does not know. |
| Configure error from `check.py` | An `SU2_OMP_*` macro without its `END_SU2_OMP_*` (OpenMP reverse build). |
| Crash or nonsense in the reverse sweep with OpenMP | A no-shared-reading section where threads do share values, or preaccumulation in a loop whose threads share points. |
| Out of memory | The tape is too big. Check `WRT_AD_STATISTICS= YES`, preaccumulate, use passive sections. |

## Testing gradients

### Compare with finite differences

Compare the AD gradient with finite differences on a small case. The Python scripts in `SU2_PY` compute the
gradient of the objective with respect to the design variables in the configuration (`DEFINITION_DV`):

| Script | Method | Needs |
|---|---|---|
| `discrete_adjoint.py -f <cfg>` | Discrete adjoint: `SU2_CFD`, `SU2_CFD_AD`, then `SU2_DOT_AD`. | Normal and reverse builds. |
| `direct_differentiation.py -f <cfg>` | Forward mode, one `SU2_CFD_DIRECTDIFF` run per design variable. The result is in `DIRECTDIFF/of_grad_directdiff.dat`. | Forward build. `HISTORY_OUTPUT` must contain a derivative group (`D_...`). |
| `finite_differences.py -f <cfg>` | One primal run per design variable, step `FIN_DIFF_STEP`. The result is in the folder `FINDIFF`. | Normal build. |

All three use the same configuration file. The case `TestCases/cont_adj_euler/naca0012` has configurations for
them (`inv_NACA0012.cfg`, `inv_NACA0012_FD.cfg`).

Tips:

- Converge the primal and the adjoint tightly. Most mismatches come from poor convergence.
- Try two or three finite-difference steps. A small step suffers from rounding, a large one from truncation.
- Forward and reverse mode should agree to many digits, because both are exact derivatives of the same discrete
  code. A difference between them points to an AD error, not to convergence or step size.

### Unit tests

For a self-contained piece of code (a fluid model, a numerics class), a Catch2 unit test is faster. The test
driver `test_driver_AD` (`UnitTests/meson.build`, list `su2_cfd_tests_ad`) is built with `-Denable-tests=true`
and `-Denable-autodiff=true`. `UnitTests/SU2_CFD/fluid/CFluidModel_tests_AD.cpp` checks the derivatives that a
fluid model computes by hand against AD:

```cpp
AD::Reset();
AD::StartRecording();
AD::RegisterInput(val_rho);
AD::RegisterInput(val_e);

fluid_model->SetTDState_rhoe(val_rho, val_e);
su2double val_p_pred = fluid_model->GetPressure();
AD::RegisterOutput(val_p_pred);
AD::StopRecording();

SU2_TYPE::SetDerivative(val_p_pred, 1.0);
AD::ComputeAdjoint();

/*--- Extract pressure derivatives from fluid model and AD. ---*/
const su2double dpdrho_e_AD = SU2_TYPE::GetDerivative(val_rho);
// ...
CHECK(SU2_TYPE::GetValue(dpdrho_e_AD) == Approx(SU2_TYPE::GetValue(dpdrho_e_pred)));
```

See [Writing Unit Tests](../Writing-Unit-Tests/).

### Regression tests

The AD regression scripts are `serial_regression_AD.py`, `parallel_regression_AD.py` and
`hybrid_regression_AD.py`. They compare adjoint residuals and sensitivities after a few iterations, and they run
the Python gradient scripts above (`discadj_euler_py`, `directdiff_euler_py`, ...) against stored gradient files.
See [Running Regression Tests](../Running-Regression-Tests/).

When a change breaks the AD tests:

1. Check the primal tests first. If the primal values changed on purpose, the adjoint values change too.
2. If only the AD tests fail, find out why before you change any value. Run the tape test, try `PREACC= NO`, and
   compare with finite differences or forward mode on a small case.
3. If the new values are correct, update them as explained in
   [Update the stored values](../Running-Regression-Tests/#update-the-stored-values), and say in the pull request
   why they changed.

If your feature is meant for the discrete adjoint, add an AD regression test. For a feature with new recorded
code, consider a multizone test with `enabled_with_tapetests = True` and `tapetest_vals = [0]`.

## Checklist for a new feature

- [ ] All values that can carry derivatives are `su2double`. `passivedouble` only for values that never need them.
- [ ] No `SU2_TYPE::GetValue`, casts or passive sections on values that the result depends on.
- [ ] New inputs of the iteration (for example a new boundary value) are registered in `RegisterVariables` and
      extracted in `ExtractAdjoint_Variables`, or registered with `RegisterSolutionExtra`.
- [ ] Values computed outside the iteration are recomputed in `SetDependencies`.
- [ ] Code that runs only at some iterations also runs during the recording (`AD::TapeActive()`), if it affects the
      result.
- [ ] Every preaccumulation section declares all inputs and outputs, and every return path calls `AD::EndPreacc()`.
- [ ] MPI goes through `SU2_MPI` (active data) or `SelectMPIWrapper<passivedouble>::W` (passive data).
- [ ] OpenMP uses the `SU2_OMP_*` macros; no-shared-reading sections contain no shared reads.
- [ ] The feature compiles in the reverse, forward and `Tag` builds.
- [ ] The tape test shows no new inconsistencies.
- [ ] The gradient agrees with finite differences (or forward mode) on a small case.
- [ ] An AD regression test covers the feature.

Last checked against SU2 develop `11c177b904`.
