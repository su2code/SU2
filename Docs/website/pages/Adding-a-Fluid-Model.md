---
title: Adding a Fluid Model
permalink: /docs_v7/Adding-a-Fluid-Model/
---

This page shows how to add a new thermodynamic model to the finite volume flow solvers of SU2. Such a model is an equation of state: the relations between pressure, density, temperature and energy of the fluid. The page also covers, more briefly, new models for viscosity, thermal conductivity and mass diffusivity.

You will learn:

- what the compressible and incompressible solvers ask from a fluid model;
- where the model is created;
- how the model works with the non-dimensionalization, that is, the scaling of all quantities by reference values.

For the overall layout of the code, see [Code Structure](../Code-Structure/). For options in general, see [Adding a Config Option](../Adding-a-Config-Option/).

---

- [Files you will touch](#files-you-will-touch)
- [The FLUID_MODEL option](#the-fluid_model-option)
- [The CFluidModel class](#the-cfluidmodel-class)
- [What the solvers call](#what-the-solvers-call)
- [Existing models](#existing-models)
- [Add a compressible model](#add-a-compressible-model)
- [Add an incompressible model](#add-an-incompressible-model)
- [Non-dimensionalization and free-stream values](#non-dimensionalization-and-free-stream-values)
- [Transport property models](#transport-property-models)
- [NEMO](#nemo)
- [Discrete adjoint](#discrete-adjoint)
- [Tests](#tests)
- [Checklist](#checklist)

---

## Files you will touch

| File | What goes there |
|---|---|
| `Common/include/option_structure.hpp` | A value in `enum ENUM_FLUIDMODEL` and its keyword in `FluidModel_Map` |
| `Common/include/CConfig.hpp`, `Common/src/CConfig.cpp` | Options for the parameters of the model, checks in `CConfig::SetPostprocessing` |
| `config_template.cfg` | The keyword in the comment above `FLUID_MODEL`, and the new options |
| `SU2_CFD/include/fluid/CMyFluid.hpp`, `SU2_CFD/src/fluid/CMyFluid.cpp` | The model class (a header alone is fine for simple models) |
| `SU2_CFD/src/meson.build` | The new `.cpp` in the `fluid/...` list |
| `SU2_CFD/src/solvers/CEulerSolver.cpp` (compressible) or `SU2_CFD/src/solvers/CIncEulerSolver.cpp` (incompressible) | The `new` of the model in `SetNondimensionalization` (two places) |
| `SU2_CFD/src/solvers/CFEM_DG_EulerSolver.cpp` | Optional: the same for the DG solver |
| `SU2_CFD/src/output/CFlowOutput.cpp` | Optional: the description in `CFlowOutput::WriteForcesBreakdown` |
| `UnitTests/SU2_CFD/fluid/CFluidModel_tests.cpp`, `CFluidModel_tests_AD.cpp` | Consistency tests of the model |
| `TestCases/serial_regression.py`, `TestCases/parallel_regression.py` | A regression test |

## The FLUID_MODEL option

This section shows how to register the keyword of your model and its parameters.

The option is registered in `CConfig::SetConfig_Options`:

```cpp
  addEnumOption("FLUID_MODEL", Kind_FluidModel, FluidModel_Map, STANDARD_AIR);
```

The enum and the map are in `Common/include/option_structure.hpp`. `ENUM_FLUIDMODEL` is a plain `enum`, so the values have no prefix: write `IDEAL_GAS`, not `ENUM_FLUIDMODEL::IDEAL_GAS`.

```cpp
enum ENUM_FLUIDMODEL {
  STANDARD_AIR = 0,       /*!< \brief Standard air gas model. */
  IDEAL_GAS = 1,          /*!< \brief Ideal gas model. */
  VW_GAS = 2,             /*!< \brief Van Der Waals gas model. */
  PR_GAS = 3,             /*!< \brief Perfect Real gas model. */
  CONSTANT_DENSITY = 4,   /*!< \brief Constant density gas model. */
  INC_IDEAL_GAS = 5,      /*!< \brief Incompressible ideal gas model. */
  INC_IDEAL_GAS_POLY = 6, /*!< \brief Inc. ideal gas, polynomial gas model. */
  MUTATIONPP = 7,         /*!< \brief Mutation++ gas model for nonequilibrium flow. */
  SU2_NONEQ = 8,          /*!< \brief User defined gas model for nonequilibrium flow. */
  FLUID_MIXTURE = 9,      /*!< \brief Species mixture model. */
  COOLPROP = 10,          /*!< \brief Thermodynamics library. */
  FLUID_FLAMELET = 11,    /*!< \brief lookup table (LUT) method for premixed flamelets. */
  DATADRIVEN_FLUID = 12,  /*!< \brief multi-layer perceptron driven fluid model. */
};
```

To register the model:

1. Add the new value at the end of the enum.
2. Add a `MakePair("MY_FLUID", MY_FLUID)` line to `FluidModel_Map`.
3. Add the keyword to the comment above `FLUID_MODEL` in `config_template.cfg`.

Read the value with `config->GetKind_FluidModel()`.

The model parameters are normal options. The existing ones are next to `FLUID_MODEL` in `CConfig::SetConfig_Options`. Reuse them when they fit:

- `GAMMA_VALUE`, `GAS_CONSTANT`, `CRITICAL_TEMPERATURE`, `CRITICAL_PRESSURE`, `ACENTRIC_FACTOR`;
- `FLUID_NAME` (CoolProp);
- `INTERPOLATION_METHOD` and `FILENAMES_INTERPOLATOR` (data-driven);
- `SPECIFIC_HEAT_CP`, `MOLECULAR_WEIGHT`, `THERMODYNAMIC_PRESSURE` (incompressible);
- ...

## The CFluidModel class

Your model derives from [`CFluidModel`](../../doxygen/classCFluidModel.html) (`SU2_CFD/include/fluid/CFluidModel.hpp`, `SU2_CFD/src/fluid/CFluidModel.cpp`). This class holds the thermodynamic state of **one** point. A solver sets the state with one of the `SetTDState_...` functions, then reads the other quantities with getters.

The state is stored in protected members. A child class must fill them:

- `Density`, `Pressure`, `Temperature`, `StaticEnergy`, `Enthalpy`, `Entropy`, `SoundSpeed2`, `Cp`, `Cv`;
- the derivatives `dPdrho_e`, `dPde_rho`, `dTdrho_e`, `dTde_rho`, used by the compressible solver;
- the derivatives `dhdrho_P`, `dhdP_rho`, `dsdrho_P`, `dsdP_rho`, used by the non-reflecting boundary conditions (NRBC);
- `Mu`, `Kt`, `mass_diffusivity` and their derivatives. The transport models fill these (see [below](#transport-property-models)).

The getters (`GetPressure()`, `GetTemperature()`, `GetSoundSpeed()`, `GetSoundSpeed2()`, `GetdPdrho_e()`, ...) are not virtual. They return these members.

The virtual functions that set the state all do nothing in the base class. Override the ones that your target solver calls:

| Function | Inputs |
|---|---|
| `SetTDState_rhoe(su2double rho, su2double e)` | density, static energy |
| `SetTDState_PT(su2double P, su2double T)` | pressure, temperature |
| `SetTDState_Prho(su2double P, su2double rho)` | pressure, density |
| `SetEnergy_Prho(su2double P, su2double rho)` | pressure, density; only `StaticEnergy` is required |
| `SetTDState_hs(su2double h, su2double s)` | enthalpy, entropy |
| `SetTDState_rhoT(su2double rho, su2double T)` | density, temperature |
| `SetTDState_Ps(su2double P, su2double s)` | pressure, entropy |
| `ComputeDerivativeNRBC_Prho(su2double P, su2double rho)` | pressure, density; sets the state and `dhdrho_P`, `dhdP_rho`, `dsdrho_P`, `dsdP_rho` |
| `SetTDState_T(su2double val_Temperature, const su2double* val_scalars = nullptr)` | temperature (incompressible, with optional species or flamelet scalars) |
| `SetTDState_h(su2double val_enthalpy, const su2double* val_scalars = nullptr)` | enthalpy (incompressible) |

Special models use other virtual functions:

- `SetCpModel` (`CIncIdealGasPolynomial`);
- `EvaluateDataSet` (`CFluidFlamelet`);
- `GetNScalars` (not overridden by any model);
- `GetExtrapolation` (`CDataDrivenFluid`, `CFluidFlamelet`);
- `GetnIter_Newton` (`CDataDrivenFluid`);
- `GetEnthalpyDiffusivity` and `GetGradEnthalpyDiffusivity` (`CFluidScalar`);
- the transport getters `GetLaminarViscosity`, `GetThermalConductivity`, `GetMassDiffusivity`.

A call that is not overridden does nothing and gives no error. So check the [table of callers](#what-the-solvers-call) before you leave a function out.

## What the solvers call

This section lists which functions each solver calls, and so which ones your model must implement.

**Compressible solver** (`CEulerSolver`, `CNSSolver`). The main call is in the variable classes. `CEulerVariable::SetPrimVar` (`SU2_CFD/src/variables/CEulerVariable.cpp`) runs for every point at every iteration:

```cpp
bool CEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool RightVol = true;

  SetVelocity(iPoint);   // Computes velocity and velocity^2
  su2double density      = GetDensity(iPoint);
  su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint);

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());

  // ...
}
```

Other calls in the variable classes:

- `CNSVariable::SetPrimVar` does the same, with the turbulent kinetic energy removed from the energy. It then calls `GetLaminarViscosity()`, `GetThermalConductivity()` and `GetCp()`.
- The `SetSecondaryVar` functions store `GetdPdrho_e()` and `GetdPde_rho()` (Euler). For Navier-Stokes they also store `GetdTdrho_e()`, `GetdTde_rho()`, `Getdmudrho_T()`, `GetdmudT_rho()`, `Getdktdrho_T()`, `GetdktdT_rho()`.
- The solver reaches the model with `GetFluidModel()`, which returns the object of the current OpenMP thread.

So in `SetTDState_rhoe`, a compressible model must set at least:

- `Density`, `StaticEnergy`, `Pressure`, `Temperature`, `SoundSpeed2`;
- `dPdrho_e`, `dPde_rho`, `dTdrho_e`, `dTde_rho`;
- `Entropy`;
- `Cp`, for viscous flows.

The ideal gas sets `Cp`, `Cv` and the temperature derivatives only once, in its constructor.

The other functions are called by:

| Function | Callers (compressible) |
|---|---|
| `SetTDState_rhoe` | the variable classes (above), the `CEulerSolver` constructor, `BC_Riemann`, `BC_TurboRiemann`, `BC_Giles` |
| `SetTDState_PT` | `CEulerSolver::SetNondimensionalization` (free stream), `BC_Riemann`, `BC_TurboRiemann`, `BC_Giles`, `BC_ActDisk` |
| `SetTDState_Prho` | `SetNondimensionalization`, `ComputeConsistentExtrapolation` (MUSCL with non-ideal fluids), `BC_Riemann`, `BC_TurboRiemann`, `BC_Giles`, `PreprocessBC_Giles`, `MixedOut_Average`, `TurboAverageProcess`, `SetFreeStream_TurboSolution`, `CFVMFlowSolverBase::BC_Fluid_Interface`; also `CTurbSASolver` and `CTurbSSTSolver` |
| `SetEnergy_Prho` | `SetNondimensionalization` (non-dimensional free-stream energy), `LowMachPrimitiveCorrection` |
| `SetTDState_rhoT` | `SetNondimensionalization` (`INIT_OPTION= REYNOLDS`) |
| `SetTDState_hs` | `BC_Riemann`, `BC_TurboRiemann` |
| `SetTDState_Ps` | `CTurboOutput` (turbomachinery performance) |
| `ComputeDerivativeNRBC_Prho` | `CEulerSolver::MixedOut_Average`, `CNumerics::ComputeResJacobianGiles` |

Some boundary conditions do not use the fluid model. For example, `CEulerSolver::BC_Inlet` uses ideal-gas relations with `Gamma` and `Gas_Constant`. For non-ideal fluids, the Riemann and Giles boundary conditions go through the fluid model.

**Incompressible solver** (`CIncEulerSolver`, `CIncNSSolver`). The energy variable is the enthalpy.

- `CIncEulerVariable::SetPrimVar` and `CIncNSVariable::SetPrimVar` call `SetTDState_h(Enthalpy, scalar)`. They then read `GetTemperature()`, `GetDensity()`, `GetCp()`, `GetCv()` and `GetEnthalpy()`. For viscous flows they also read `GetLaminarViscosity()` and `GetThermalConductivity()`, after `SetEddyViscosity`.
- Initialization and boundary conditions call `SetTDState_T`.
- The thermodynamic pressure is constant. It is set in the constructor of the model, and `SetPressure` can change it.

## Existing models

Start from an existing model close to yours. All are in `SU2_CFD/include/fluid/` and `SU2_CFD/src/fluid/`.

| Keyword | Class | Solver | Notes |
|---|---|---|---|
| `STANDARD_AIR`, `IDEAL_GAS` | [`CIdealGas`](../../doxygen/classCIdealGas.html) | compressible | Closed-form inverses for every `SetTDState_...` |
| `VW_GAS` | [`CVanDerWaalsGas`](../../doxygen/classCVanDerWaalsGas.html) | compressible | Child of `CIdealGas`; Newton iteration in `SetTDState_PT` |
| `PR_GAS` | [`CPengRobinson`](../../doxygen/classCPengRobinson.html) | compressible | Child of `CIdealGas` |
| `COOLPROP` | [`CCoolProp`](../../doxygen/classCCoolProp.html) | compressible | Calls the CoolProp library (`-Denable-coolprop=true`); dimensional only |
| `DATADRIVEN_FLUID` | [`CDataDrivenFluid`](../../doxygen/classCDataDrivenFluid.html) | compressible | Multi-layer perceptron, a small neural network (`-Denable-mlpcpp=true`), or look-up table; Newton solver for the inverse functions |
| `CONSTANT_DENSITY` | [`CConstantDensity`](../../doxygen/classCConstantDensity.html) | incompressible | Header only |
| `INC_IDEAL_GAS` | [`CIncIdealGas`](../../doxygen/classCIncIdealGas.html) | incompressible | Header only, density from the operating pressure |
| `INC_IDEAL_GAS_POLY` | [`CIncIdealGasPolynomial`](../../doxygen/classCIncIdealGasPolynomial.html) | incompressible | Header only, Cp polynomial in temperature |
| `FLUID_MIXTURE` | [`CFluidScalar`](../../doxygen/classCFluidScalar.html) | incompressible | Species mixture, one transport model per species |
| `FLUID_FLAMELET` | [`CFluidFlamelet`](../../doxygen/classCFluidFlamelet.html) | incompressible | Look-up table of a flamelet manifold (precomputed combustion states) |
| `SU2_NONEQ`, `MUTATIONPP` | `CSU2TCLib`, `CMutationTCLib` (children of [`CNEMOGas`](../../doxygen/classCNEMOGas.html)) | NEMO | See [NEMO](#nemo) |

Good starting points:

- `CVanDerWaalsGas` for a cubic equation of state;
- `CDataDrivenFluid` for a model that inverts its functions with a Newton solver;
- `CIncIdealGas` for an incompressible model.

The state function of the ideal gas, from `SU2_CFD/src/fluid/CIdealGas.cpp`, shows how little a model needs:

```cpp
void CIdealGas::SetTDState_rhoe(su2double rho, su2double e) {
  Density = rho;
  StaticEnergy = e;
  dPdrho_e = Gamma_Minus_One * StaticEnergy;
  dPde_rho = Gamma_Minus_One * Density;
  Pressure = dPde_rho * StaticEnergy;
  Temperature = dTde_rho * StaticEnergy;
  SoundSpeed2 = Gamma * Gas_Constant * Temperature;

  if (ComputeEntropy) Entropy = Cv * log(Temperature) - Gas_Constant * log(Density);
}
```

## Add a compressible model

Follow these steps to add a model for `CEulerSolver` and `CNSSolver`.

1. **Option.** Add the enum value and the keyword ([above](#the-fluid_model-option)). Add options for new parameters: member, getter, `add...Option`, `config_template.cfg`. See [Adding a Config Option](../Adding-a-Config-Option/).
2. **Class.** Write `class CMyFluid final : public CFluidModel` in `SU2_CFD/include/fluid/CMyFluid.hpp`. You can also derive from `CIdealGas` to reuse its members, as `CVanDerWaalsGas` and `CPengRobinson` do. Implement the functions of the [table above](#what-the-solvers-call) that your cases need. The initialization alone needs `SetTDState_rhoe`, `SetTDState_PT`, `SetTDState_Prho` and `SetEnergy_Prho`. It also needs `SetTDState_rhoT` when the free stream is set from the Reynolds number (`INIT_OPTION= REYNOLDS`).
3. **Build.** Add `fluid/CMyFluid.cpp` to the `fluid/...` list in `SU2_CFD/src/meson.build`.
4. **Creation.** In `CEulerSolver::SetNondimensionalization` (`SU2_CFD/src/solvers/CEulerSolver.cpp`), add a `case` to **both** `switch (config->GetKind_FluidModel())` blocks:
   - the first creates `auxFluidModel` with **dimensional** parameters, to compute the free-stream state;
   - the second, inside `SU2_OMP_PARALLEL`, creates one model per thread with **non-dimensional** parameters, and stores it in `FluidModel[thread]`.

   This is the second block for the Van der Waals gas:

   ```cpp
      case VW_GAS:
        FluidModel[thread] = new CVanDerWaalsGas(Gamma, Gas_ConstantND,
                                                 config->GetPressure_Critical() / config->GetPressure_Ref(),
                                                 config->GetTemperature_Critical() / config->GetTemperature_Ref());
        break;
   ```

   Include the header at the top of `CEulerSolver.cpp`. `CEulerSolver::~CEulerSolver` deletes the models. The first `switch` ends with `SU2_MPI::Error("Unknown fluid model.", ...)`, so an incompressible keyword cannot reach the compressible solver.
5. **Optional.** Add the name to the `ModelTable` printed at the end of `SetNondimensionalization`, and to `CFlowOutput::WriteForcesBreakdown`. If the model must also work with the DG solver, repeat step 4 in `CFEM_DG_EulerSolver::SetNondimensionalization`.

When the model is not `STANDARD_AIR` or `IDEAL_GAS`, the solver changes in three ways:

- `CDriver::InitializeNumerics` creates the general numerics: `CUpwGeneralRoe_Flow`, `CUpwGeneralHLLC_Flow`, `CGeneralAvgGrad_Flow`, `CSourceGeneralAxisymmetric_Flow`. They read the secondary variables: `S_i[0]` = dP/drho at constant e, `S_i[1]` = dP/de at constant rho.
- `CEulerSolver::Upwind_Residual` uses the classic path, not the vectorized one. It makes the MUSCL states consistent with `ComputeConsistentExtrapolation`, which calls `SetTDState_Prho`.
- `CConfig::SetPostprocessing` allows only `ROE`, `HLLC` and `JST` (search for "Non-Ideal Compressible Fluids").

## Add an incompressible model

Follow these steps to add a model for `CIncEulerSolver` and `CIncNSSolver`.

1. Add the enum value, the keyword and the options, as for a compressible model.
2. Write the class with `SetTDState_T` and `SetTDState_h`.
   - They must set `Temperature`, `Density` and `Enthalpy`.
   - Set `Cp` and `Cv` too. Constant values can be set in the constructor, as in `CIncIdealGas`.
   - The enthalpy is the sensible enthalpy relative to a reference temperature. The existing models take `STD_REF_TEMP / config->GetTemperature_Ref()` as the last constructor argument.
3. In `CIncEulerSolver::SetNondimensionalization` (`SU2_CFD/src/solvers/CIncEulerSolver.cpp`), add a `case` to both `switch` blocks. The incompressible model headers are included in this file. The two blocks are:
   - the dimensional `auxFluidModel`;
   - the per-thread `FluidModel` loop (`for (auto& fluidModel : FluidModel)`). After the second `switch`, this loop sets the transport models when the flow is viscous.
4. Check `CConfig::SetPostprocessing`:
   - With `INC_DENSITY_MODEL= CONSTANT` or `BOUSSINESQ`, the code sets `Kind_FluidModel = CONSTANT_DENSITY`.
   - With `INC_DENSITY_MODEL= VARIABLE`, it accepts only `INC_IDEAL_GAS`, `INC_IDEAL_GAS_POLY`, `FLUID_MIXTURE` and `FLUID_FLAMELET`. Add the new value there if the density varies.
   - The same function rejects `DATADRIVEN_FLUID` for the incompressible solvers.

## Non-dimensionalization and free-stream values

This section shows how the solver sets the free-stream state and the reference values, and which scales your model must follow. `SetNondimensionalization` creates each fluid model twice:

- a dimensional copy, for the setup;
- non-dimensional copies, one per thread, used during the run.

**Compressible free stream** (`CEulerSolver::SetNondimensionalization`):

1. `Temperature_Ref`, `Viscosity_Ref` and `Conductivity_Ref` are set to 1. So the `...ND()` getters of `CConfig` return dimensional values for the auxiliary model. For example, `GetMu_ConstantND()` returns `Mu_Constant / Viscosity_Ref`.
2. `auxFluidModel` is created.
   - With `FREESTREAM_OPTION= TEMPERATURE_FS`, it gets `SetTDState_PT(Pressure_FreeStream, Temperature_FreeStream)` and gives the density.
   - With `DENSITY_FS`, it gets `SetTDState_Prho` and gives the temperature.
3. `GetSoundSpeed()` turns the Mach number into the free-stream velocity.
4. For viscous flows, `SetLaminarViscosityModel` and `SetThermalConductivityModel` are called on the auxiliary model. With `INIT_OPTION= REYNOLDS`, the density comes from the Reynolds number, and `SetTDState_rhoT` gives the pressure. `GetStaticEnergy()` gives the free-stream energy.
5. The reference values are set from `REF_DIMENSIONALIZATION`: `Pressure_Ref`, `Density_Ref`, `Temperature_Ref`, then `Velocity_Ref = sqrt(Pressure_Ref/Density_Ref)`, `Gas_Constant_Ref = Velocity_Ref^2/Temperature_Ref`, `Viscosity_Ref`, `Conductivity_Ref`, ...
6. The per-thread models are created with non-dimensional parameters. `SetEnergy_Prho(Pressure_FreeStreamND, Density_FreeStreamND)` gives the non-dimensional free-stream energy. `Energy_Ref` is the ratio of the dimensional and the non-dimensional energies.

So the non-dimensional model must be consistent with these scales:

| Quantity | Divided by |
|---|---|
| pressure | `Pressure_Ref` |
| density | `Density_Ref` |
| temperature | `Temperature_Ref` |
| energy and enthalpy | `Velocity_Ref^2` |
| gas constant | `Gas_Constant_Ref` |

The simplest way is to give non-dimensional parameters to the constructor. `VW_GAS` and `PR_GAS` do this with `Gas_ConstantND` and with the critical values divided by the reference values.

Some models cannot be scaled, for example tables or external libraries in physical units. For these, allow only `REF_DIMENSIONALIZATION= DIMENSIONAL`. CoolProp does this in `CConfig::SetPostprocessing`:

```cpp
  /*--- Check if CoolProp is used with non-dimensionalization. ---*/
  if (Kind_FluidModel == COOLPROP && Ref_NonDim != DIMENSIONAL) {
    SU2_MPI::Error("CoolProp can not be used with non-dimensionalization.", CURRENT_FUNCTION);
  }
```

**Incompressible** (`CIncEulerSolver::SetNondimensionalization`):

- `auxFluidModel->SetTDState_T(Temperature_FreeStream, scalar_init)` gives the free-stream density.
- `INIT_OPTION_INC` decides what happens next. Either the thermodynamic pressure is scaled (`SetPressure`) to match the given density, or the density is taken from the model.
- The per-thread models get non-dimensional parameters: `Density_FreeStreamND`, `Specific_Heat_CpND`, `Gas_ConstantND`, `Pressure_ThermodynamicND`. Most of them call `SetTDState_T(Temperature_FreeStreamND)` right away.

## Transport property models

This section shows how to add a model for viscosity, thermal conductivity or mass diffusivity. These are separate objects owned by the fluid model:

```cpp
  unique_ptr<CViscosityModel> LaminarViscosity;       /*!< \brief Laminar Viscosity Model */
  unique_ptr<CConductivityModel> ThermalConductivity; /*!< \brief Thermal Conductivity Model */
  unique_ptr<CDiffusivityModel> MassDiffusivity;      /*!< \brief Mass Diffusivity Model */
```

| Option | Enum and map | Base class | Factory in `CFluidModel.cpp` |
|---|---|---|---|
| `VISCOSITY_MODEL` | `VISCOSITYMODEL`, `ViscosityModel_Map` | [`CViscosityModel`](../../doxygen/classCViscosityModel.html): `SetViscosity(t, rho)` | `MakeLaminarViscosityModel` |
| `CONDUCTIVITY_MODEL`, `TURBULENT_CONDUCTIVITY_MODEL` | `CONDUCTIVITYMODEL`, `ConductivityModel_Map`; `CONDUCTIVITYMODEL_TURB`, `TurbConductivityModel_Map` | [`CConductivityModel`](../../doxygen/classCConductivityModel.html): `SetConductivity(t, rho, mu_lam, mu_turb, cp, dmudrho_t, dmudt_rho)` | `MakeThermalConductivityModel` |
| `DIFFUSIVITY_MODEL` | `DIFFUSIVITYMODEL`, `Diffusivity_Model_Map` | [`CDiffusivityModel`](../../doxygen/classCDiffusivityModel.html): `SetDiffusivity(rho, mu_lam, cp, kt)` | `MakeMassDiffusivityModel` |

Existing children:

- viscosity: `CConstantViscosity`, `CSutherland`, `CPolynomialViscosity`, `CCoolPropViscosity`;
- conductivity: `CConstantConductivity`, `CConstantConductivityRANS`, `CConstantPrandtl`, `CConstantPrandtlRANS`, `CPolynomialConductivity`, `CPolynomialConductivityRANS`, `CCoolPropConductivity`;
- diffusivity: `CConstantDiffusivity`, `CConstantSchmidt`, `CConstantLewisDiffusivity`.

The `FLAMELET` values return `nullptr`, because `CFluidFlamelet` overrides the getters and reads the table.

A model sets its value and its derivatives in protected members. `CSutherland`:

```cpp
  void SetViscosity(su2double t, su2double rho) override {
    const su2double t_ref_inv = 1.0 / t_ref_;
    const su2double t_nondim = t_ref_inv * t;
    mu_ = mu_ref_ * t_nondim * sqrt(t_nondim) * ((t_ref_ + s_) / (t + s_));

    /*--- Set Viscosity Derivatives. ---*/
    const su2double ts_inv = 1.0 / (t + s_);
    dmudrho_t_ = 0.0;
    dmudt_rho_ = mu_ref_ * (t_ref_ + s_) * ts_inv * sqrt(t_nondim) * (1.5 * t_ref_inv - t_nondim * ts_inv);
  }
```

To add a transport model:

1. Add a value to the enum and its map in `option_structure.hpp`, and options for its parameters.
2. Write the class in `SU2_CFD/include/fluid/`, derived from the base class of the table. A header alone is fine.
3. Add a `case` in the factory in `SU2_CFD/src/fluid/CFluidModel.cpp`, and include the header there. Take the parameters from non-dimensional getters of `CConfig`, such as `GetMu_RefND(iSpecies)`. These divide by the reference values, so the same factory serves both the dimensional and the non-dimensional model. The `iSpecies` argument is used by `CFluidScalar`, which creates one model per species.

How the solvers use it:

- In `SetNondimensionalization`, the solvers call `SetLaminarViscosityModel(config)`, `SetThermalConductivityModel(config)` and `SetMassDiffusivityModel(config)` on each fluid model.
- `CFluidModel::GetLaminarViscosity()` then calls `LaminarViscosity->SetViscosity(Temperature, Density)`. So the state must be set before.
- `GetThermalConductivity()` uses `Mu`, `Mu_Turb` and `Cp`. Only the incompressible solver sets `Mu_Turb`, with `SetEddyViscosity`. In the compressible solver it stays zero, and the viscous numerics add the turbulent part of the conductivity with `PRANDTL_TURB`.

## NEMO

Nothing on this page applies to non-equilibrium flows (`NEMO_EULER`, `NEMO_NAVIER_STOKES`). They use a separate interface, [`CNEMOGas`](../../doxygen/classCNEMOGas.html), a child of `CFluidModel`. Its children are `CSU2TCLib` and `CMutationTCLib`. They are created in `CNEMOEulerSolver::SetNondimensionalization` and in `CNEMONumerics`.

## Discrete adjoint

Write your model so that automatic differentiation (AD) works. The discrete adjoint records the primal iteration on the AD tape, the record of every operation it differentiates. The fluid model is called inside that iteration: in `SetPrimVar` for every point, in the boundary conditions, and in the reconstruction. So:

- Write the model with `su2double`. A value converted to `double` has no derivative, whether through `SU2_TYPE::GetValue` or through a library that works on `double`. CoolProp is such a library: the error message of `CCoolProp` says that it cannot be used with `directdiff` or `autodiff`.
- Iterative inversions, such as the Newton loops in `CVanDerWaalsGas::SetTDState_PT` or in `CDataDrivenFluid`, are recorded with all their iterations. This makes the tape larger.
- Avoid hidden state between calls. For example, `CVanDerWaalsGas` keeps the compressibility factor `Zed` of the last call as the first guess of its Newton loop. So the result of `SetTDState_PT` depends on the previous call, within the tolerance.
- The hand-written derivatives (`dPdrho_e`, `dPde_rho`, ...) are recorded and differentiated like any other value. They enter the general Roe and HLLC fluxes (through the secondary variables), the viscous Jacobians and some boundary conditions. Compute them with `su2double` from the inputs, and check them with the unit tests below.

An adjoint regression test with a non-ideal fluid is `transonic_stator` (`TestCases/disc_adj_turbomachinery/transonic_stator_2D/transonic_stator.cfg`, `PR_GAS`). It is in `serial_regression_AD.py` and `parallel_regression_AD.py`.

## Tests

Add unit tests and a regression test for the new model.

**Unit tests.** Add a `TEST_CASE` for the new model to both files below. Build with `-Denable-tests=true`; see [Writing Unit Tests](../Writing-Unit-Tests/).

- `UnitTests/SU2_CFD/fluid/CFluidModel_tests.cpp` has `FluidModelChecks(fluid_model, p, T)`. It checks that `SetTDState_rhoe` inverts `SetTDState_PT`. It also compares `dPdrho_e`, `dPde_rho`, `dTdrho_e`, `dTde_rho` with finite differences.
- `UnitTests/SU2_CFD/fluid/CFluidModel_tests_AD.cpp` has `FluidModelADChecks(fluid_model, p, T)`. It compares the same derivatives, `Cp`, and the NRBC entropy derivatives with AD. It is built as `test_driver_AD`.

**Regression tests.** The non-ideal compressible cases are in `TestCases/nicf/`. They are in the "NICFD" sections of the regression scripts:

| Test | Configuration | Fluid model | Script |
|---|---|---|---|
| `ls89_sa` | `nicf/LS89/turb_SA_PR.cfg` | `PR_GAS` | serial |
| `edge_VW` | `nicf/edge/edge_VW.cfg` | `VW_GAS` | serial, parallel, hybrid |
| `edge_PPR` | `nicf/edge/edge_PPR.cfg` | `PR_GAS` | serial, parallel, hybrid |
| `coolprop_fluidModel`, `coolprop_transportModel` | `nicf/coolprop/fluidModel.cfg`, `transportModel.cfg` | `COOLPROP` (only on `x86_64`) | parallel |
| `datadriven_fluidModel` | `nicf/datadriven/datadriven_nozzle.cfg` | `DATADRIVEN_FLUID` | parallel |

Some `TestCases/turbomachinery/` cases also use `PR_GAS`: `centrifugal_blade`, `centrifugal_stage`, `transonic_stator_2D`.

To add a regression test:

1. Copy one of these cases and switch it to the new model.
2. Add it as a new `TestCase`; see [Running Regression Tests](../Running-Regression-Tests/).
3. Keep the test short (few iterations).
4. If the model is meant for optimization, add an AD test too.

## Checklist

- [ ] Enum value and keyword in `ENUM_FLUIDMODEL` / `FluidModel_Map`, keyword and options in `config_template.cfg`.
- [ ] The `SetTDState_...` functions that your solver calls, and all the state members and derivatives they must set.
- [ ] New `.cpp` in `SU2_CFD/src/meson.build`.
- [ ] Both `switch` blocks in `SetNondimensionalization` (dimensional and non-dimensional).
- [ ] Consistent non-dimensional parameters, or a check that forces `REF_DIMENSIONALIZATION= DIMENSIONAL`.
- [ ] Checks in `CConfig::SetPostprocessing` (density model list for incompressible models).
- [ ] Unit tests with `FluidModelChecks` and `FluidModelADChecks`, and a regression test.

Last checked against SU2 develop `11c177b904`.
