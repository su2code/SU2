---
title: Adding Output Fields
permalink: /docs_v7/Adding-Output-Fields/
---

This page shows how to add a new field to the output of SU2. A field can go to the screen and history output, or to the volume output (restart and visualization files). The page is for developers. For how users select fields, see [History and Solution Output](../Custom-Output/).

---

- [Terms](#terms)
- [Where the output classes are](#where-the-output-classes-are)
- [How the output works](#how-the-output-works)
- [Add a history field](#add-a-history-field)
- [Add a volume field](#add-a-volume-field)
- [Turbulence, transition and species fields](#turbulence-transition-and-species-fields)
- [How users select the new field](#how-users-select-the-new-field)
- [Checklist](#checklist)

---

## Terms

- **History field**: one value per iteration. It can be printed on screen (`SCREEN_OUTPUT`) and written to the history file (`HISTORY_OUTPUT`).
- **Per-surface history field**: one history value for each marker of a list, for example for each marker in `MARKER_MONITORING`.
- **Volume field**: one value per mesh point. It is written to the restart and visualization files (`VOLUME_OUTPUT`).
- **Group**: a name shared by several fields. Users can request a whole group instead of single fields.

Each field has two labels:

- a **name**, used in the configuration file and in the code, for example `DRAG`;
- a **header**, printed on screen or in the file, for example `CD`.

## Where the output classes are

All output code is in `SU2_CFD/include/output/` and `SU2_CFD/src/output/`. The base class is [`COutput`](../../doxygen/classCOutput.html). `COutputFactory::CreateOutput` picks the child class from the `SOLVER` option:

| `SOLVER` | Output class |
|---|---|
| `EULER`, `NAVIER_STOKES`, `RANS` | `CFlowCompOutput` |
| `INC_EULER`, `INC_NAVIER_STOKES`, `INC_RANS` | `CFlowIncOutput` |
| `NEMO_EULER`, `NEMO_NAVIER_STOKES` | `CNEMOCompOutput` |
| `FEM_EULER`, `FEM_NAVIER_STOKES`, `FEM_RANS`, `FEM_LES` | `CFlowCompFEMOutput` |
| `HEAT_EQUATION` | `CHeatOutput` |
| `ELASTICITY` | `CElasticityOutput` |
| compressible flow with `MATH_PROBLEM= CONTINUOUS_ADJOINT` or `DISCRETE_ADJOINT` | `CAdjFlowCompOutput` |
| incompressible flow with `MATH_PROBLEM= DISCRETE_ADJOINT` | `CAdjFlowIncOutput` |
| `HEAT_EQUATION` with `MATH_PROBLEM= DISCRETE_ADJOINT` | `CAdjHeatOutput` |
| `ELASTICITY` with `MATH_PROBLEM= DISCRETE_ADJOINT` | `CAdjElasticityOutput` |

Each class `CXxxOutput` is in `SU2_CFD/src/output/CXxxOutput.cpp`, with the header `SU2_CFD/include/output/CXxxOutput.hpp`.

Shared parents and other output classes:

- [`CFVMOutput`](../../doxygen/classCFVMOutput.html): parent of the finite volume outputs. It adds the coordinates (`AddCoordinates`, `LoadCoordinates`) and common fields such as the mesh quality (`AddCommonFVMOutputs`, `LoadCommonFVMOutputs`).
- [`CFlowOutput`](../../doxygen/classCFlowOutput.html): parent of the flow outputs (`CFlowCompOutput`, `CFlowIncOutput`, `CNEMOCompOutput`, `CFlowCompFEMOutput`). It holds:
  - the turbulence, transition and species fields;
  - the aerodynamic coefficients;
  - the `MARKER_ANALYZE` fields.

  It also sets the wall fields (skin friction, heat flux, y+) in `LoadSurfaceData`.
- [`CAdjFlowOutput`](../../doxygen/classCAdjFlowOutput.html): parent of the adjoint flow outputs, with the adjoint turbulence and species fields.
- [`CMultizoneOutput`](../../doxygen/classCMultizoneOutput.html): the history of a multizone run. It copies the history fields of every zone and adds the suffix `[<zone index>]` to their names. So you do not need to change it for a new zone field.
- `CBaselineOutput` (used by `SU2_SOL` and `SU2_DOT`) and `CMeshOutput` (used by `SU2_DEF`).

## How the output works

`COutput` calls virtual functions ("hooks") that the child classes implement. To add a field, you change these hooks:

| Hook | When it is called | What it does |
|---|---|---|
| `SetHistoryOutputFields(CConfig*)` | once, at start-up (from `PreprocessHistoryOutput`) | declares the history fields with `AddHistoryOutput` and `AddHistoryOutputPerSurface` |
| `LoadHistoryData(CConfig*, CGeometry*, CSolver**)` | every iteration (from `SetHistoryOutput`) | sets the values with `SetHistoryOutputValue` and `SetHistoryOutputPerSurfaceValue` |
| `SetVolumeOutputFields(CConfig*)` | once, at start-up (from `PreprocessVolumeOutput`) | declares the volume fields with `AddVolumeOutput` |
| `LoadVolumeData(CConfig*, CGeometry*, CSolver**, unsigned long iPoint)` | when files are written, for each point owned by the rank. For time-domain problems, at every call of `SetResultFiles`, to update the time averages. | sets the values with `SetVolumeOutputValue` |
| `LoadSurfaceData(CConfig*, CGeometry*, CSolver**, unsigned long iPoint, unsigned short iMarker, unsigned long iVertex)` | after `LoadVolumeData`, for each vertex of the solid wall markers | sets volume fields that only exist at walls |

Other hooks:

- The DG output uses `LoadVolumeDataFEM` instead of `LoadVolumeData`.
- Multizone output uses `SetMultizoneHistoryOutputFields` and `LoadMultizoneHistoryData`.
- `SetAdditionalScreenOutput` and `WriteAdditionalFiles`.

The functions to declare and set fields are members of `COutput` (`SU2_CFD/include/output/COutput.hpp`):

```cpp
void AddHistoryOutput(string name, string field_name, ScreenOutputFormat format,
                      string groupname, string description,
                      HistoryFieldType field_type = HistoryFieldType::DEFAULT);
void SetHistoryOutputValue(string name, su2double value);

void AddHistoryOutputPerSurface(string name, string field_name, ScreenOutputFormat format,
                                string groupname, const vector<string>& marker_names,
                                HistoryFieldType field_type = HistoryFieldType::DEFAULT);
void SetHistoryOutputPerSurfaceValue(string name, su2double value, unsigned short iMarker);

void AddVolumeOutput(const string& name, const string& field_name,
                     const string& group_name, const string& description);
void SetVolumeOutputValue(const string& name, unsigned long iPoint, su2double value);
```

`ScreenOutputFormat` sets how the value is printed on screen: `INTEGER`, `FIXED`, `SCIENTIFIC` or `PERCENT`.

`HistoryFieldType` says what kind of value the field holds. SU2 uses it to make extra fields and to check convergence:

| `HistoryFieldType` | Use it for | What SU2 adds |
|---|---|---|
| `RESIDUAL` | residuals, stored as `log10` of the value | The field `REL_<name>`: the value minus a reference value (the initial residual, or a larger value seen later). Also `AVG_<group>` for the groups `RMS_RES`, `MAX_RES` and `BGS_RES`. In `CONV_FIELD`, the field converges when it is below `CONV_RESIDUAL_MINVAL`. |
| `COEFFICIENT` | integrated quantities, such as forces or mass flow | `TAVG_<name>` (time average) when `TIME_DOMAIN= YES`. `D_<name>` when `DIRECT_DIFF` is set. `CAUCHY_<name>` when the field is in `CONV_FIELD`. In `CONV_FIELD`, the field converges by the Cauchy criterion (`CONV_CAUCHY_ELEMS`, `CONV_CAUCHY_EPS`). |
| `DEFAULT` | everything else, for example iteration counts | Nothing. A `DEFAULT` field in `CONV_FIELD` never meets the convergence criterion. |
| `AUTO_RESIDUAL`, `AUTO_COEFFICIENT` | fields made by SU2 itself (the `REL_`, `AVG_`, `TAVG_`, `D_`, `CAUCHY_` fields) | Do not use them for new fields. |

These steps are in `COutput::PostprocessHistoryFields`, `COutput::PostprocessHistoryData` and `COutput::ConvergenceMonitoring` (`SU2_CFD/src/output/COutput.cpp`).

## Add a history field

1. Compute the value in the solver. Store it in a member and add a getter. If the output class accesses the solver through a `CSolver*`, add a virtual getter with a default value to `CSolver` (`SU2_CFD/include/solvers/CSolver.hpp`). Override it in your solver.
2. In `SetHistoryOutputFields` of the output class, call `AddHistoryOutput`. Pick a name, a short header, a format, a group, a description and a `HistoryFieldType`.
3. In `LoadHistoryData`, call `SetHistoryOutputValue` with the same name.
4. Use the same condition in both functions. If the field was not added, `SetHistoryOutputValue` stops SU2 with the error "Cannot find output field with name ...".

### Example: the buffet sensor

`CNSSolver::Buffet_Monitoring` computes the buffet sensor and stores it in `Total_Buffet_Metric`. `CSolver` has the virtual getter `GetTotal_Buffet_Metric()` (it returns 0), and `CNSSolver` overrides it.

In `CFlowCompOutput::SetHistoryOutputFields` (`SU2_CFD/src/output/CFlowCompOutput.cpp`):

```cpp
  if (config->GetViscous()) {
    AddHistoryOutput("BUFFET", "Buffet", ScreenOutputFormat::SCIENTIFIC, "AERO_COEFF", "Buffet sensor", HistoryFieldType::COEFFICIENT);
  }
```

In `CFlowCompOutput::LoadHistoryData` (same file):

```cpp
  CSolver* flow_solver = solver[FLOW_SOL];
  // ...
  if (config->GetViscous()) {
    SetHistoryOutputValue("BUFFET", flow_solver->GetTotal_Buffet_Metric());
  }
```

The field joins the group `AERO_COEFF`. Because it is a `COEFFICIENT`:

- it can be used in `CONV_FIELD` with the Cauchy criterion;
- it gets a `TAVG_BUFFET` field in unsteady runs.

### Per-surface fields

`AddHistoryOutputPerSurface` makes one field for each marker in the list. The header of each one is `field_name(marker)`. Set the values with the index of the marker in the list.

This example is from `CFlowOutput::AddAerodynamicCoefficients` and `CFlowOutput::SetAerodynamicCoefficients` (`SU2_CFD/src/output/CFlowOutput.cpp`):

```cpp
  vector<string> Marker_Monitoring;
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++){
    Marker_Monitoring.push_back(config->GetMarker_Monitoring_TagBound(iMarker_Monitoring));
  }
  /// DESCRIPTION: Drag coefficient
  AddHistoryOutputPerSurface("DRAG_ON_SURFACE",       "CD",   ScreenOutputFormat::FIXED, "AERO_COEFF_SURF", Marker_Monitoring, HistoryFieldType::COEFFICIENT);
```

```cpp
  for (unsigned short iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
    SetHistoryOutputPerSurfaceValue("DRAG_ON_SURFACE", flow_solver->GetSurface_CD(iMarker_Monitoring), iMarker_Monitoring);
```

Users can select per-surface fields in `SCREEN_OUTPUT` and `HISTORY_OUTPUT`, but not in `CONV_FIELD`.

### Notes

- `LoadHistoryData` runs every iteration, also when the field is not requested. Keep it cheap. If the value is expensive, compute it in the solver only when needed.
- For residual fields, store `log10` of the residual, as the existing fields do (for example `log10(flow_solver->GetRes_RMS(0))`).
- Always write a description. `SU2_CFD -d` (dry run, `CDummyDriver`) prints the list of fields with their descriptions. History fields with an empty description are left out of that list.
- Python code can read history fields through the wrapper (`GetOutputNames`, `GetOutputValue` of `CDriverBase`).
- Some Python scripts in `SU2_PY` look up the group and type of a history field in `SU2_PY/SU2/io/historyMap.py`. Examples are the objectives of unsteady optimizations in `SU2_PY/SU2/io/config.py`, and `direct_differentiation.py`. If the new field should be usable there:
  1. Write the `AddHistoryOutput(...)` call on one line, at the start of the line. `SU2_PY/updateHistoryMap.py` makes `historyMap.py` by reading the lines that start with `AddHistoryOutput(` in the output source files listed in the script. It does not find calls written in another way.
  2. Run `SU2_PY/updateHistoryMap.py` (it needs `SU2_HOME`).
  3. Commit the new `historyMap.py`.

## Add a volume field

1. Make the value available per point, usually through a getter of the `CVariable` child (`solver->GetNodes()`). If the output class uses a `CVariable*`, add a virtual getter with a default value to `CVariable`. Override it in your variable class.
2. In `SetVolumeOutputFields`, call `AddVolumeOutput` with a name, a header, a group and a description.
3. In `LoadVolumeData`, call `SetVolumeOutputValue(name, iPoint, value)`.
4. Use the same condition in both functions.

**Important:** `LoadVolumeData` must make the same calls, in the same order, for every point. The reason is that `SetVolumeOutputValue` caches the position of each field. It builds the cache at the first point, and reuses it for all other points in the order of the calls. So conditions may depend on the configuration, but not on `iPoint`. The same holds for `LoadSurfaceData`: for each vertex, make either no calls or the full sequence.

### Example: the Roe low-dissipation function

`CVariable` has the virtual getter `GetRoe_Dissipation(iPoint)` (it returns 0), and `CNSVariable` overrides it.

In `CFlowCompOutput::SetVolumeOutputFields` (`SU2_CFD/src/output/CFlowCompOutput.cpp`):

```cpp
  // Roe Low Dissipation
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS) {
    AddVolumeOutput("ROE_DISSIPATION", "Roe_Dissipation", "ROE_DISSIPATION", "Value of the Roe dissipation");
  }
```

In `CFlowCompOutput::LoadVolumeData` (same file):

```cpp
  const auto* Node_Flow = solver[FLOW_SOL]->GetNodes();
  // ...
  if (config->GetKind_RoeLowDiss() != NO_ROELOWDISS){
    SetVolumeOutputValue("ROE_DISSIPATION", iPoint, Node_Flow->GetRoe_Dissipation(iPoint));
  }
```

### Wall fields

Set values that only exist on walls in `LoadSurfaceData`. These values are stored per marker and vertex. `LoadSurfaceData` is called only for the markers where `config->GetSolid_Wall(iMarker)` is true. You still declare the field with `AddVolumeOutput` in `SetVolumeOutputFields`.

From `CFlowOutput::LoadSurfaceData` (`SU2_CFD/src/output/CFlowOutput.cpp`):

```cpp
  if (!config->GetViscous_Wall(iMarker)) return;
  // ...
  SetVolumeOutputValue("Y_PLUS", iPoint, solver[FLOW_SOL]->GetYPlus(iMarker, iVertex));
```

### Groups and restart files

- Common volume groups are `COORDINATES`, `SOLUTION`, `PRIMITIVE`, `RESIDUAL`, `LIMITER`, `GRID_VELOCITY` and `MESH_QUALITY`. You can also make a new group, as `ROE_DISSIPATION` above.
- The restart file holds the groups `COORDINATES`, `SOLUTION`, `SENSITIVITY` and `GRID_VELOCITY` (with the default `WRT_RESTART_COMPACT= YES`).
- The solvers read the restart file **by position**: coordinates first, then the `SOLUTION` fields in the order of the `AddVolumeOutput` calls. So add a field to `SOLUTION` only if it is a new solution variable. In that case, change the restart reading (`LoadRestart` of the solver) at the same time. Put any other new field in another group.
- Surface files (`SURFACE_PARAVIEW`, `SURFACE_CSV`, ...) take their data from the volume fields, on the markers in `MARKER_PLOTTING`. There is no separate list of surface fields.
- For time-averaged fields, use `SetAvgVolumeOutputValue` instead of `SetVolumeOutputValue` (see `CFlowOutput::LoadTimeAveragedData`).

## Turbulence, transition and species fields

Add turbulence, transition and species fields in `CFlowOutput` (`SU2_CFD/src/output/CFlowOutput.cpp`). The compressible, incompressible and NEMO flow outputs share them. The functions below switch on the model (`TurbModelFamily(config->GetKind_Turb_Model())`, `config->GetKind_Trans_Model()`, `config->GetKind_Species_Model()`):

| Function | Content |
|---|---|
| `AddHistoryOutputFields_ScalarRMS_RES`, `AddHistoryOutputFields_ScalarMAX_RES`, `AddHistoryOutputFields_ScalarBGS_RES`, `AddHistoryOutputFieldsScalarLinsol` | declare the history fields |
| `LoadHistoryDataScalar` | set the history values |
| `SetVolumeOutputFieldsScalarSolution`, `SetVolumeOutputFieldsScalarResidual`, `SetVolumeOutputFieldsScalarPrimitive`, `SetVolumeOutputFieldsScalarLimiter`, `SetVolumeOutputFieldsScalarSource`, `SetVolumeOutputFieldsScalarLookup`, `SetVolumeOutputFieldsScalarMisc` | declare the volume fields |
| `LoadVolumeDataScalar` | set the volume values |

`CFlowCompOutput`, `CFlowIncOutput` and `CNEMOCompOutput` call these functions from their own hooks. So a new turbulence or species field added there appears in all three outputs.

For the adjoint outputs, the same pattern is in `CAdjFlowOutput` (`AddHistoryOutputFields_AdjScalarRMS_RES`, `LoadHistoryDataAdjScalar`, `SetVolumeOutputFieldsAdjScalarSolution`, `LoadVolumeDataAdjScalar`, ...).

## How users select the new field

Users request fields or whole groups in the configuration file. SU2 ignores names that do not exist for the current setup, and prints a message at start-up.

```
SCREEN_OUTPUT= INNER_ITER, RMS_DENSITY, LIFT, DRAG, BUFFET
HISTORY_OUTPUT= ITER, RMS_RES, AERO_COEFF
VOLUME_OUTPUT= COORDINATES, SOLUTION, PRIMITIVE, ROE_DISSIPATION
CONV_FIELD= RMS_DENSITY
```

- `SCREEN_OUTPUT` takes field names or groups.
- `HISTORY_OUTPUT` takes groups or field names.
- `VOLUME_OUTPUT` takes field names or groups.
- `CONV_FIELD` takes history field names (not per-surface fields).
- `CUSTOM_OUTPUTS` lets users define new history values from expressions of existing fields, without code changes (see [User Defined Functions](../Custom-Output/#user-defined-functions)).

See [History and Solution Output](../Custom-Output/) for the user documentation, and the [Option Reference](../Option-Reference/) for the options.

## Checklist

1. Run `SU2_CFD -d <case>.cfg`. Check that the field appears in the list, in the right group, with its description.
2. Run a short case with the field in `SCREEN_OUTPUT`, `HISTORY_OUTPUT` or `VOLUME_OUTPUT`, and check the values. For volume fields, open the file in ParaView or Tecplot.
3. For a history field, check the extra fields that its type creates, for example `REL_<name>` or `TAVG_<name>`. If it is meant to work in `CONV_FIELD`, check that too.
4. Check a parallel run. The master rank prints the history values, so they must be global values. Reduce them over all ranks in the solver, for example with `SU2_MPI::Allreduce` (as `CNSSolver::Buffet_Monitoring` does).
5. Add or update a regression test. The test runner compares the last columns of the screen output at iteration `test_iter` with `test_vals` (see `TestCases/TestCase.py`). So:
   - put the new field in the `SCREEN_OUTPUT` of a test case;
   - add its value to `test_vals` in `TestCases/serial_regression.py`, `TestCases/parallel_regression.py` or the `_AD` versions.

   See [Running Regression Tests](../Running-Regression-Tests/).
6. Did you add the field to a group that existing tests print on screen (for example `RMS_RES` or `AERO_COEFF`)? Then the screen columns of those tests change. Run them and update their `test_vals`.
7. Unit tests for output classes are in `UnitTests/SU2_CFD/output/`.
8. If the field is of general interest, document it for users in [History and Solution Output](../Custom-Output/).

Last checked against SU2 develop `11c177b904`.
