/*!
 * \file CPBIncEulerSolver.cpp
 * \brief Main subroutines for solving pressure based incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/CPBIncEulerSolver.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/fluid/CConstantDensity.hpp"
#include "../../include/fluid/CIncIdealGas.hpp"
#include "../../include/fluid/CIncIdealGasPolynomial.hpp"
#include "../../include/variables/CPBIncNSVariable.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/fluid/CFluidScalar.hpp"
#include "../../include/fluid/CFluidFlamelet.hpp"
#include "../../include/fluid/CFluidModel.hpp"

// Temporarily directly copied from CIncEulerSolver's constructor to adhere to necessary
// initializations, 
CPBIncEulerSolver::CPBIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh,
                                 const bool navier_stokes) :
  CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>(*geometry, *config) {
  SU2_ZONE_SCOPED


  /*--- Based on the navier_stokes boolean, determine if this constructor is
   *    being called by itself, or by its derived class CIncNSSolver. ---*/
  const string description = navier_stokes? "Navier-Stokes" : "Euler";

  unsigned short iMarker;
  ifstream restart_file;
  bool restart = (config->GetRestart() || config->GetRestart_Flow());
  int Unst_RestartIter = 0;
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_MARCHING::TIME_STEPPING;
  bool adjoint = (config->GetContinuous_Adjoint()) || (config->GetDiscrete_Adjoint());
  const bool centered = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  const su2double* scalar_init = nullptr;

  /* A grid is defined as dynamic if there's rigid grid movement or grid deformation AND the problem is time domain */
  dynamic_grid = config->GetDynamic_Grid();

  /*--- Store the multigrid level. ---*/
  MGLevel = iMesh;

  /*--- Check for a restart file to evaluate if there is a change in the angle of attack
   before computing all the non-dimesional quantities. ---*/

  if (restart) {

    if (iMesh == MESH_0) {

      /*--- Modify file name for a dual-time unsteady restart ---*/

      if (dual_time) {
        if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
        else if (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST)
          Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
        else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
      }

      /*--- Modify file name for a time stepping unsteady restart ---*/

      if (time_stepping) {
        if (adjoint) Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
        else Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
      }

    }

    if (config->GetKind_Streamwise_Periodic() == ENUM_STREAMWISE_PERIODIC::MASSFLOW) {
      if (rank==MASTER_NODE) cout << "Setting streamwise periodic pressure drop from restart metadata file." << endl;
    }

    auto filename_ = config->GetFilename("flow", ".meta", Unst_RestartIter);
    Read_SU2_Restart_Metadata(geometry, config, adjoint, filename_);
  }


  /*--- Set the gamma value ---*/

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Define geometry constants in the solver structure.
   * Incompressible flow, primitive variables (P, vx, vy, vz, T, rho, beta, lamMu, EddyMu, Kt_eff, Cp, Cv) ---*/

  nDim = geometry->GetnDim();

  /*--- Make sure to align the sizes with the constructor of CIncEulerVariable. ---*/
  nVar = nDim + 2;
  nPrimVar = nDim + 10;
  /*--- Centered schemes only need gradients for viscous fluxes (T and v, but we need also to include P). ---*/
  nPrimVarGrad = nDim + (centered ? 2 : 4);

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nPrimVarGrad;

  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Store the number of vertices on each marker for deallocation later ---*/

  nVertex.resize(nMarker);
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];

  /*--- Perform the non-dimensionalization for the flow equations using the
   specified reference values. ---*/

  // SetNondimensionalization(config, iMesh);

  /*--- Check if we are executing a verification case. If so, the
   VerificationSolution object will be instantiated for a particular
   option from the available library of verification solutions. Note
   that this is done after SetNondim(), as problem-specific initial
   parameters are needed by the solution constructors. ---*/

  SetVerificationSolution(nDim, nVar, config);

  /*--- Allocate base class members. ---*/

  Allocate(*config);

  /*--- MPI + OpenMP initialization. ---*/

  HybridParallelInitialization(*config, *geometry);

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT) {

    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy);
  }
  else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (" << description << "). MG level: " << iMesh <<"." << endl;
  }

   /*--- Read farfield conditions ---*/

  Density_Inf     = config->GetDensity_FreeStreamND();
  Pressure_Inf    = config->GetPressure_FreeStreamND();
  Velocity_Inf    = config->GetVelocity_FreeStreamND();
  Temperature_Inf = config->GetTemperature_FreeStreamND();
  if (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) scalar_init = config->GetSpecies_Init();
  //GetFluidModel()->SetTDState_T(Temperature_Inf, scalar_init);
  //Enthalpy_Inf = GetFluidModel()->GetEnthalpy();

  /*--- Initialize the secondary values for direct derivative approximations ---*/

  switch (config->GetDirectDiff()) {
    case NO_DERIVATIVE:
      /*--- Default ---*/
      break;
    case D_DENSITY:
      SU2_TYPE::SetDerivative(Density_Inf, 1.0);
      break;
    case D_PRESSURE:
      SU2_TYPE::SetDerivative(Pressure_Inf, 1.0);
      break;
    case D_TEMPERATURE:
      SU2_TYPE::SetDerivative(Temperature_Inf, 1.0);
      break;
    case D_MACH: case D_AOA:
    case D_SIDESLIP: case D_REYNOLDS:
    case D_TURB2LAM: case D_DESIGN:
      /*--- Already done in postprocessing of config ---*/
      break;
    default:
      break;
  }

  SetReferenceValues(*config);

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  if (navier_stokes) {
    nodes = new CPBIncNSVariable(Pressure_Inf, Velocity_Inf, Enthalpy_Inf, nPoint, nDim, nVar, config);
  } else {
    nodes = new CPBIncEulerVariable(Pressure_Inf, Velocity_Inf, Enthalpy_Inf, nPoint, nDim, nVar, config);
  }
  SetBaseClassPointerToNodes();

  if (iMesh == MESH_0) {
    nodes->NonPhysicalEdgeCounter.resize(geometry->GetnEdge()) = 0;
  }

  /*--- Initial comms. ---*/

  CommunicateInitialState(geometry, config);

  /*--- Sizing edge mass flux array ---*/
  if (config->GetBounded_Scalar())
    EdgeMassFluxes.resize(geometry->GetnEdge()) = su2double(0.0);

  /*--- Add the solver name. ---*/
  SolverName = "INC.FLOW";

  /*--- Finally, check that the static arrays will be large enough (keep this
   *    check at the bottom to make sure we consider the "final" values). ---*/
  if((nDim > MAXNDIM) || (nPrimVar > MAXNVAR))
    SU2_MPI::Error("Oops! The CPBIncEulerSolver static array sizes are not large enough.", CURRENT_FUNCTION);

}
