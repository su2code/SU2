/*!
 * \file CTurbSolver.cpp
 * \brief Main subroutines of CTurbSolver class
 * \author F. Palacios, A. Bueno
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

#include "../../include/solvers/CTurbSolver.hpp"
#include "../../../Common/include/parallelization/omp_structure.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.inl"

/*--- Explicit instantiation of the parent class of CTurbSolver. ---*/
template class CScalarSolver<CTurbVariable>;

CTurbSolver::CTurbSolver(CGeometry* geometry, CConfig *config, bool conservative)
  : CScalarSolver<CTurbVariable>(geometry, config, conservative) {
  /*--- Store if an implicit scheme is used, for use during periodic boundary conditions. ---*/
  SetImplicitPeriodic(config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
}

CTurbSolver::~CTurbSolver() {
  for (auto& mat : SlidingState) {
    for (auto ptr : mat) delete [] ptr;
  }
}

void CTurbSolver::BC_Riemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::BC_TurboRiemann(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Riemann(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT: case STATIC_SUPERSONIC_INFLOW_PT: case STATIC_SUPERSONIC_INFLOW_PD: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case STATIC_PRESSURE:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}


void CTurbSolver::BC_Giles(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  string Marker_Tag         = config->GetMarker_All_TagBound(val_marker);

  switch(config->GetKind_Data_Giles(Marker_Tag))
  {
  case TOTAL_CONDITIONS_PT:case TOTAL_CONDITIONS_PT_1D: case DENSITY_VELOCITY:
    BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  case MIXING_IN:
    if (config->GetBoolTurbMixingPlane()){
      BC_Inlet_MixingPlane(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    else{
      BC_Inlet_Turbo(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    }
    break;

  case STATIC_PRESSURE: case MIXING_OUT: case STATIC_PRESSURE_1D: case RADIAL_EQUILIBRIUM:
    BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
    break;
  }
}

void CTurbSolver::LoadRestart(CGeometry** geometry, CSolver*** solver, CConfig* config, int val_iter,
                              bool val_update_geo) {
  /*--- Restart the solution from file information ---*/

  const string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);

  /*--- To make this routine safe to call in parallel most of it can only be executed by one thread. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

    if (config->GetRead_Binary_Restart()) {
      Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
    } else {
      Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
    }

    /*--- Skip flow variables ---*/

    unsigned short skipVars = nDim + solver[MESH_0][FLOW_SOL]->GetnVar();

    /*--- Adjust the number of solution variables in the incompressible
     restart. We always carry a space in nVar for the energy equation in the
     mean flow solver, but we only write it to the restart if it is active.
     Therefore, we must reduce skipVars here if energy is inactive so that
     the turbulent variables are read correctly. ---*/

    const bool incompressible = (config->GetKind_Regime() == ENUM_REGIME::INCOMPRESSIBLE);
    const bool energy = config->GetEnergy_Equation();
    const bool weakly_coupled_heat = config->GetWeakly_Coupled_Heat();
    const bool flamelet = (config->GetKind_FluidModel() == FLUID_FLAMELET);

    if (incompressible && ((!energy) && (!weakly_coupled_heat) && (!flamelet))) skipVars--;

    /*--- Load data from the restart into correct containers. ---*/

    unsigned long counter = 0;
    for (auto iPoint_Global = 0ul; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      const auto iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        const auto index = counter * Restart_Vars[1] + skipVars;
        for (auto iVar = 0u; iVar < nVar; iVar++) nodes->SetSolution(iPoint_Local, iVar, Restart_Data[index + iVar]);

        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Detect a wrong solution file ---*/

    if (counter != nPointDomain) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" does not match with the mesh file!\n") +
                         string("This can be caused by empty lines at the end of the file."),
                     CURRENT_FUNCTION);
    }

  }  // end safe global access, pre and postprocessing are thread-safe.
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  /*--- MPI solution and compute the eddy viscosity ---*/

  solver[MESH_0][TURB_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][TURB_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);

  /*--- For turbulent+species simulations the solver Pre-/Postprocessing is done by the species solver. ---*/
  if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE && config->GetKind_Trans_Model() == TURB_TRANS_MODEL::NONE) {
    solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER,
                                            RUNTIME_FLOW_SYS, false);
    solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  }

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (auto iMesh = 1u; iMesh <= config->GetnMGLevels(); iMesh++) {
    MultigridRestriction(*geometry[iMesh - 1], solver[iMesh - 1][TURB_SOL]->GetNodes()->GetSolution(),
                         *geometry[iMesh], solver[iMesh][TURB_SOL]->GetNodes()->GetSolution());
    solver[iMesh][TURB_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][TURB_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);

    if (config->GetKind_Species_Model() == SPECIES_MODEL::NONE) {
      solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS,
                                            false);
      solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
    }
  }

  /*--- Go back to single threaded execution. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    /*--- Delete the class memory that is used to load the restart. ---*/

    delete[] Restart_Vars;
    Restart_Vars = nullptr;
    delete[] Restart_Data;
    Restart_Data = nullptr;
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS
}

void CTurbSolver::Impose_Fixed_Values(const CGeometry *geometry, const CConfig *config){
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Check whether turbulence quantities are fixed to far-field values on a half-plane. ---*/
  if(config->GetTurb_Fixed_Values()){

    /*--- Form normalized far-field velocity ---*/
    const su2double* velocity_inf = config->GetVelocity_FreeStreamND();
    su2double velmag_inf = GeometryToolbox::Norm(nDim, velocity_inf);
    if(velmag_inf==0)
      SU2_MPI::Error("Far-field velocity is zero, cannot fix turbulence quantities to inflow values.", CURRENT_FUNCTION);
    su2double unit_velocity_inf[MAXNDIM];
    for(unsigned short iDim=0; iDim<nDim; iDim++)
      unit_velocity_inf[iDim] = velocity_inf[iDim] / velmag_inf;

    SU2_OMP_FOR_DYN(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if( GeometryToolbox::DotProduct(nDim, geometry->nodes->GetCoord(iPoint), unit_velocity_inf)
        < config->GetTurb_Fixed_Values_MaxScalarProd() ) {
        /*--- Set the solution values and zero the residual ---*/
        nodes->SetSolution_Old(iPoint, Solution_Inf);
        nodes->SetSolution(iPoint, Solution_Inf);
        LinSysRes.SetBlock_Zero(iPoint);
        if (implicit) {
          /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
          for(unsigned long iVar=0; iVar<nVar; iVar++)
            Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
        }
      }
    }
    END_SU2_OMP_FOR
  }

}
