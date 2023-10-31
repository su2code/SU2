/*!
 * \file CIntegration.cpp
 * \brief Implementation of the base class for space and time integration.
 * \author F. Palacios, T. Economon
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/integration/CIntegration.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"


CIntegration::CIntegration() {
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
}

void CIntegration::Space_Integration(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics **numerics,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep,
                                     unsigned short RunTime_EqSystem) {
  unsigned short iMarker, KindBC;

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);
  bool dual_time = ((config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND));

  /*--- Compute inviscid residuals ---*/

  switch (config->GetKind_ConvNumScheme()) {
    case SPACE_CENTERED:
      solver_container[MainSolver]->Centered_Residual(geometry, solver_container, numerics, config, iMesh, iRKStep);
      break;
    case SPACE_UPWIND:
      solver_container[MainSolver]->Upwind_Residual(geometry, solver_container, numerics, config, iMesh);
      break;
  }

  /*--- Compute viscous residuals ---*/
  solver_container[MainSolver]->Viscous_Residual(geometry, solver_container, numerics, config, iMesh, iRKStep);

  /*--- Compute source term residuals ---*/
  solver_container[MainSolver]->Source_Residual(geometry, solver_container, numerics, config, iMesh);

  /*--- Add viscous and convective residuals, and compute the Dual Time Source term ---*/

  if (dual_time)
    solver_container[MainSolver]->SetResidual_DualTime(geometry, solver_container, config, iRKStep, iMesh, RunTime_EqSystem);

  /*--- Pick convective and viscous numerics objects for the current thread. ---*/

  CNumerics* conv_bound_numerics = numerics[CONV_BOUND_TERM + omp_get_thread_num()*MAX_TERMS];
  CNumerics* visc_bound_numerics = numerics[VISC_BOUND_TERM + omp_get_thread_num()*MAX_TERMS];

  /*--- Pause preaccumulation in boundary conditions for hybrid parallel AD. ---*/
  /// TODO: Check if this is really needed.
  //const auto pausePreacc = (omp_get_num_threads() > 1) && AD::PausePreaccumulation();

  /*--- Boundary conditions that depend on other boundaries (they require MPI sincronization)---*/

  solver_container[MainSolver]->BC_Fluid_Interface(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config);

  /*--- Compute Fourier Transformations for markers where NRBC_BOUNDARY is applied---*/

  if (config->GetBoolGiles() && config->GetSpatialFourier()){
    solver_container[MainSolver]->PreprocessBC_Giles(geometry, config, conv_bound_numerics, INFLOW);

    solver_container[MainSolver]->PreprocessBC_Giles(geometry, config, conv_bound_numerics, OUTFLOW);
  }

  /*--- Weak boundary conditions ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    switch (KindBC) {
      case EULER_WALL:
        solver_container[MainSolver]->BC_Euler_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case ACTDISK_INLET:
        solver_container[MainSolver]->BC_ActDisk_Inlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case ENGINE_INFLOW:
        solver_container[MainSolver]->BC_Engine_Inflow(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case INLET_FLOW:
        solver_container[MainSolver]->BC_Inlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case ACTDISK_OUTLET:
        solver_container[MainSolver]->BC_ActDisk_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case ENGINE_EXHAUST:
        solver_container[MainSolver]->BC_Engine_Exhaust(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case SUPERSONIC_INLET:
        solver_container[MainSolver]->BC_Supersonic_Inlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case OUTLET_FLOW:
        solver_container[MainSolver]->BC_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case SUPERSONIC_OUTLET:
        solver_container[MainSolver]->BC_Supersonic_Outlet(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case GILES_BOUNDARY:
        solver_container[MainSolver]->BC_Giles(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case RIEMANN_BOUNDARY:
        if (config->GetBoolTurbomachinery()){
          solver_container[MainSolver]->BC_TurboRiemann(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        }
        else{
          solver_container[MainSolver]->BC_Riemann(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        }
        break;
      case FAR_FIELD:
        solver_container[MainSolver]->BC_Far_Field(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case SYMMETRY_PLANE:
        solver_container[MainSolver]->BC_Sym_Plane(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
    }
  }

  /*--- Modification of the system on the whole domain, like for a strong BC. ---*/
  solver_container[MainSolver]->Impose_Fixed_Values(geometry, config);

  /*--- Strong boundary conditions (Navier-Stokes and Dirichlet type BCs) ---*/

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
    switch (config->GetMarker_All_KindBC(iMarker)) {
      case ISOTHERMAL:
        solver_container[MainSolver]->BC_Isothermal_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case HEAT_FLUX:
        solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case HEAT_TRANSFER:
        solver_container[MainSolver]->BC_HeatTransfer_Wall(geometry, config, iMarker);
        break;
      case CUSTOM_BOUNDARY:
        solver_container[MainSolver]->BC_Custom(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
      case CHT_WALL_INTERFACE:

        if ( MainSolver == FLOW_SOL && (config->GetKind_FluidModel() == FLUID_FLAMELET) ){
          solver_container[MainSolver]->BC_Isothermal_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
          break;
        }
        if ((MainSolver == SPECIES_SOL) || (MainSolver == HEAT_SOL) || ((MainSolver == FLOW_SOL) && ((config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) || config->GetEnergy_Equation()))) {
          solver_container[MainSolver]->BC_ConjugateHeat_Interface(geometry, solver_container, conv_bound_numerics, config, iMarker);
        }
        else {
          solver_container[MainSolver]->BC_HeatFlux_Wall(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        }
        break;
      case SMOLUCHOWSKI_MAXWELL:
        solver_container[MainSolver]->BC_Smoluchowski_Maxwell(geometry, solver_container, conv_bound_numerics, visc_bound_numerics, config, iMarker);
        break;
    }

  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated correctly during the communications. ---*/

  if (config->GetnMarker_Periodic() > 0) {
    solver_container[MainSolver]->BC_Periodic(geometry, solver_container, conv_bound_numerics, config);
  }

  //AD::ResumePreaccumulation(pausePreacc);

}

void CIntegration::Time_Integration(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                                    unsigned short iRKStep, unsigned short RunTime_EqSystem) {

  unsigned short MainSolver = config->GetContainerPosition(RunTime_EqSystem);

  switch (config->GetKind_TimeIntScheme()) {
    case (RUNGE_KUTTA_EXPLICIT):
      solver_container[MainSolver]->ExplicitRK_Iteration(geometry, solver_container, config, iRKStep);
      break;
    case (CLASSICAL_RK4_EXPLICIT):
      solver_container[MainSolver]->ClassicalRK4_Iteration(geometry, solver_container, config, iRKStep);
      break;
    case (EULER_EXPLICIT):
      solver_container[MainSolver]->ExplicitEuler_Iteration(geometry, solver_container, config);
      break;
    case (EULER_IMPLICIT):
      solver_container[MainSolver]->ImplicitEuler_Iteration(geometry, solver_container, config);
      break;
  }

}

void CIntegration::SetDualTime_Geometry(CGeometry *geometry, CSolver *mesh_solver, const CConfig *config, unsigned short iMesh) {

  SU2_OMP_PARALLEL
  {

  geometry->nodes->SetVolume_nM1();
  geometry->nodes->SetVolume_n();

  if (config->GetGrid_Movement()) {
    geometry->nodes->SetCoord_n1();
    geometry->nodes->SetCoord_n();
  }

  if ((iMesh==MESH_0) && config->GetDeform_Mesh()) mesh_solver->SetDualTime_Mesh();

  }
  END_SU2_OMP_PARALLEL
}

void CIntegration::SetDualTime_Solver(const CGeometry *geometry, CSolver *solver, const CConfig *config, unsigned short iMesh) {

  SU2_OMP_PARALLEL
  {
  /*--- Store old solution ---*/
  solver->GetNodes()->Set_Solution_time_n1();
  solver->GetNodes()->Set_Solution_time_n();

  SU2_OMP_SAFE_GLOBAL_ACCESS(solver->ResetCFLAdapt();)

  SU2_OMP_FOR_STAT(roundUpDiv(geometry->GetnPoint(), omp_get_num_threads()))
  for (auto iPoint = 0ul; iPoint < geometry->GetnPoint(); iPoint++) {

    /*--- Initialize the underrelaxation ---*/
    solver->GetNodes()->SetUnderRelaxation(iPoint, 1.0);

    /*--- Initialize the local CFL number ---*/
    solver->GetNodes()->SetLocalCFL(iPoint, config->GetCFL(iMesh));
  }
  END_SU2_OMP_FOR

  }
  END_SU2_OMP_PARALLEL
}
