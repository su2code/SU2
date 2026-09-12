# Flamelet Solver — Python Wrapper Wall Boundary Conditions

This document explains how to use the Python wrapper (`pysu2`) to set
custom per-vertex boundary conditions for the flamelet (FGM) scalar transport
solver, and how the setup differs between the two `FLAMELET_ENTHALPY_BC` modes.

---

## Background

The flamelet solver transports `nVar` scalar variables. Their order is:

| Index | Variable |
|-------|----------|
| 0 | `PROGVAR` (progress variable) |
| 1 | `ENTH` (total enthalpy) |
| 2 … n\_CV−1 | additional control variables (e.g. `MIXFRAC`) |
| n\_CV … nVar−1 | auxiliary / user-defined species |

Two modes control how the enthalpy boundary condition is derived at thermal
walls (`MARKER_ISOTHERMAL`, `MARKER_HEATFLUX`):

- **`FLOW_MARKERS`** — enthalpy is derived from the flow thermal field (wall
  temperature via `MARKER_ISOTHERMAL`, or wall heat flux via `MARKER_HEATFLUX`).
- **`SPECIES_MARKERS`** (default) — enthalpy and all other scalars are taken
  directly from `MARKER_WALL_SPECIES`, or from the Python wrapper when
  `MARKER_PYTHON_CUSTOM` is active.

---

## Mode 1: `FLAMELET_ENTHALPY_BC = FLOW_MARKERS`

### How it works

The C++ function `BC_HeatFlux_Wall` reads the wall heat flux and applies it as
a Neumann condition on `I_ENTH`.  When the marker is also listed in
`MARKER_PYTHON_CUSTOM`, it instead reads a **per-vertex heat flux** set by the
Python wrapper.

```
driver.SetMarkerCustomNormalHeatFlux(iMarker, iVertex, q_wall)
    → geometry->CustomBoundaryHeatFlux[iMarker][iVertex]
        → BC_HeatFlux_Wall reads geometry->GetCustomBoundaryHeatFlux(...)
```

> `BC_Isothermal_Wall` in `FLOW_MARKERS` mode always converts the wall
> temperature (from `MARKER_ISOTHERMAL`) to enthalpy via `GetEnthFromTemp`.
> Per-vertex customisation of the enthalpy Dirichlet value is not supported in
> this mode — use `SPECIES_MARKERS` instead.

### Limitations in this mode

- Only `I_ENTH` is reachable via the Python wrapper.
- Non-enthalpy scalars (`PROGVAR`, auxiliary species) receive implicit zero
  Neumann at the wall and cannot be set via Python.

### config.cfg

```cfg
FLAMELET_ENTHALPY_BC = FLOW_MARKERS

% Flow solver wall BC — sets KindBC and provides the fallback heat flux value
% used when MARKER_PYTHON_CUSTOM is not active.
MARKER_HEATFLUX = ( burner_wall, 0.0 )

% Enable the per-vertex Python override path.
MARKER_PYTHON_CUSTOM = ( burner_wall )
```

### run.py

```python
import pysu2

driver = pysu2.CSinglezoneDriver("config.cfg", 1, False)

marker_name = "burner_wall"
iMarker     = list(driver.GetMarkerTags()).index(marker_name)
n_vertex    = driver.GetNumberMarkerNodes(iMarker)

def compute_heat_flux(coord):
    """Heat flux [W/m²], positive = into the domain."""
    import math
    r = math.sqrt(coord[0]**2 + coord[1]**2)
    return -5000.0 * math.exp(-r**2 / 0.01)   # Gaussian profile, heat out

for iteration in range(driver.GetNumberIter()):

    for iVertex in range(n_vertex):
        node_id = driver.GetMarkerNode(iMarker, iVertex)
        coord   = driver.GetInitialMeshCoord(node_id)   # [x, y(, z)]
        q_wall  = compute_heat_flux(coord)
        driver.SetMarkerCustomNormalHeatFlux(iMarker, iVertex, q_wall)

    driver.Preprocess(iteration)
    driver.Run()
    driver.Postprocess()
    driver.Monitor(iteration)
    driver.Output(iteration)

driver.Finalize()
```

---

## Mode 2: `FLAMELET_ENTHALPY_BC = SPECIES_MARKERS` (default)

### How it works

Both `BC_HeatFlux_Wall` and `BC_Isothermal_Wall` delegate immediately to
`CSpeciesSolver::BC_Wall_Generic`, which processes **all `nVar` scalars**
independently.  When `MARKER_PYTHON_CUSTOM` is active, the value for each
scalar is taken from the array set by `driver.SetMarkerCustomScalar()`.

```
driver.SetMarkerCustomScalar(iMarker, iVertex, [val_0, val_1, ..., val_nVar-1])
    → CustomBoundaryScalar[iMarker](iVertex, iVar)
        → BC_Wall_Generic reads CustomBoundaryScalar for every iVar
```

The **type** of BC for each scalar (Dirichlet `VALUE` or Neumann `FLUX`) is
read from `MARKER_WALL_SPECIES` in the config and cannot be changed from
Python.  The Python wrapper only overrides the **magnitude** per vertex.

### config.cfg

The example below uses `nVar = 3` (PROGVAR, ENTH, one auxiliary species).
Adjust the number of `BC_TYPE, value` pairs to match the actual number of
scalar variables in your setup.

```cfg
FLAMELET_ENTHALPY_BC = SPECIES_MARKERS

% Flow solver wall BC — still required to set KindBC.
% The choice between MARKER_ISOTHERMAL and MARKER_HEATFLUX only affects
% the flow solver; the species solver uses BC_Wall_Generic in both cases.
MARKER_ISOTHERMAL = ( burner_wall, 300.0 )

% Per-variable BC types for the species solver.
% Format: (marker_name, BC_TYPE, fallback_value, BC_TYPE, fallback_value, ...)
%         One BC_TYPE+value pair per scalar variable, in index order.
% BC_TYPE: FLUX  → Neumann (value is normal flux density [unit/m²])
%          VALUE → Dirichlet (value is the scalar value at the wall)
%
% When MARKER_PYTHON_CUSTOM is active, the fallback_value is ignored and
% the Python-supplied value is used instead.  The BC_TYPE is always from config.
%
% idx 0 PROGVAR : zero Neumann  (no progress variable source at the wall)
% idx 1 ENTH    : Dirichlet     (enthalpy value set per-vertex by Python)
% idx 2 aux_1   : zero Neumann  (no auxiliary species flux at the wall)
MARKER_WALL_SPECIES = ( burner_wall, FLUX, 0.0, VALUE, 0.0, FLUX, 0.0 )

% Enable the per-vertex Python override path.
MARKER_PYTHON_CUSTOM = ( burner_wall )
```

### run.py

```python
import pysu2

driver = pysu2.CSinglezoneDriver("config.cfg", 1, False)

marker_name = "burner_wall"
iMarker     = list(driver.GetMarkerTags()).index(marker_name)
n_vertex    = driver.GetNumberMarkerNodes(iMarker)

def compute_wall_enthalpy(coord):
    """Return wall enthalpy [J/kg] at this vertex."""
    # Example: uniform value corresponding to T=300 K from your LUT.
    return 3.5e5

for iteration in range(driver.GetNumberIter()):

    for iVertex in range(n_vertex):
        node_id = driver.GetMarkerNode(iMarker, iVertex)
        coord   = driver.GetInitialMeshCoord(node_id)

        enth_wall = compute_wall_enthalpy(coord)

        # Pass one value per scalar variable (length = nVar).
        # PROGVAR (idx 0): 0.0 → applied as FLUX → zero Neumann (BC_TYPE from config)
        # ENTH    (idx 1): enth_wall → applied as VALUE → Dirichlet
        # aux_1   (idx 2): 0.0 → applied as FLUX → zero Neumann
        driver.SetMarkerCustomScalar(iMarker, iVertex, [0.0, enth_wall, 0.0])

    driver.Preprocess(iteration)
    driver.Run()
    driver.Postprocess()
    driver.Monitor(iteration)
    driver.Output(iteration)

driver.Finalize()
```

---

## Comparison

| | `FLOW_MARKERS` | `SPECIES_MARKERS` |
|---|---|---|
| **Python API** | `SetMarkerCustomNormalHeatFlux` | `SetMarkerCustomScalar` |
| **Storage** | `geometry->CustomBoundaryHeatFlux` | `CustomBoundaryScalar[iMarker](iVertex, iVar)` |
| **What you set** | Per-vertex heat flux (W/m²) for `I_ENTH` only | Per-vertex values for **all** `nVar` scalars |
| **BC type control** | Hard-coded Neumann in C++ | Per-variable via `MARKER_WALL_SPECIES` in config |
| **`MARKER_ISOTHERMAL` walls** | Enthalpy computed from T via `GetEnthFromTemp`; not overridable by Python | Enthalpy set as Dirichlet `VALUE` from Python |
| **Non-enthalpy scalars** | Not reachable via Python | All controlled via the same `SetMarkerCustomScalar` call |
| **Required config keys** | `MARKER_HEATFLUX`, `MARKER_PYTHON_CUSTOM` | `MARKER_ISOTHERMAL` or `MARKER_HEATFLUX`, `MARKER_WALL_SPECIES`, `MARKER_PYTHON_CUSTOM` |

---

## Notes

- `MARKER_PYTHON_CUSTOM` is a purely additive flag. A marker can appear in
  both `MARKER_HEATFLUX` (or `MARKER_ISOTHERMAL`) and `MARKER_PYTHON_CUSTOM`
  simultaneously.
- The fallback values in `MARKER_WALL_SPECIES` are used when `py_custom` is
  false (e.g., during a plain SU2\_CFD run without the Python driver). They
  allow the same config to be used both ways.
- In `SPECIES_MARKERS` mode, `SetMarkerCustomScalar` must supply a vector of
  exactly `nVar` values. Indices that correspond to `FLUX`-type variables
  should receive `0.0` unless you intentionally want a non-zero flux there.
- The `CHT` (conjugate heat transfer) path (`BC_ConjugateHeat_Interface`)
  always uses `BC_Isothermal_Wall_Generic` and is unaffected by
  `FLAMELET_ENTHALPY_BC`. Python custom BCs are not applicable to CHT markers.
