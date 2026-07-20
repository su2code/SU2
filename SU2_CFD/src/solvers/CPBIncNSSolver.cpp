/*!
 * \file CPBIncNSSolver.cpp
 * \brief Main subroutines for solving pressure based Navier-Stokes incompressible flow.
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

#include "../../include/solvers/CPBIncNSSolver.hpp"
#include "../../include/variables/CPBIncNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

//TODO: almost everything in this file is a direct or indirect copy of the code in 
// the CIncNSSolver class. this will need to be refactored.

/*--- Explicit instantiation of the parent class of CPBIncEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CPBIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>;


CPBIncNSSolver::CPBIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
  CPBIncEulerSolver(geometry, config, iMesh, true) {
  SU2_ZONE_SCOPED

  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Tke_Inf         = config->GetTke_FreeStreamND();

  /*--- Initialize the secondary values for direct derivative approximations ---*/

  switch (config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      break;
  }

  /*--- Set the initial Streamwise periodic pressure drop value. ---*/

  if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
    // Note during restarts, the flow.meta is read first. But that sets the cfg-value so we are good here.
    SPvals.Streamwise_Periodic_PressureDrop = config->GetStreamwise_Periodic_PressureDrop();

}

void CPBIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED

  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);
  const bool wall_functions = config->GetWall_Functions();

  /*--- Common preprocessing steps (implemented by CPBIncEulerSolver) ---*/

  CPBIncEulerSolver::Preprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Compute gradient for MUSCL reconstruction ---*/

  if (config->GetReconstructionGradientRequired() && muscl && !center) {
    switch (config->GetKind_Gradient_Method_Recon()) {
      case GREEN_GAUSS:
        SetPrimitive_Gradient_GG(geometry, config, true); break;
      case LEAST_SQUARES:
      case WEIGHTED_LEAST_SQUARES:
        SetPrimitive_Gradient_LS(geometry, config, true); break;
      default: break;
    }
  }

  /*--- Compute gradient of the primitive variables ---*/

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetPrimitive_Gradient_GG(geometry, config);
  }
  else if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetPrimitive_Gradient_LS(geometry, config);
  }

  /*--- Compute the limiters ---*/

  if (muscl && !center && limiter && !van_albada && !Output) {
    SetPrimitive_Limiter(geometry, config);
  }

  ComputeVorticityAndStrainMag(*config, geometry, iMesh);

  // TODO: Currently I haven ot yet needed these functions, these are however similar if not exactly the same as the ones implemented by CIncNSSolver, therefore this has to be reused.
  /*--- Compute the TauWall from the wall functions ---*/
  // if (wall_functions) {
  //   SU2_OMP_SAFE_GLOBAL_ACCESS(SetTau_Wall_WF(geometry, solver_container, config);)
  // }

  // /*--- Compute recovered pressure and temperature for streamwise periodic flow ---*/
  // if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
  //   Compute_Streamwise_Periodic_Recovered_Values(config, geometry, iMesh);

}

unsigned long CPBIncNSSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {
  SU2_ZONE_SCOPED

  unsigned long iPoint, nonPhysicalPoints = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0, LES_Mode = 0.0;
  const su2double* scalar = nullptr;
  const TURB_MODEL turb_model = config->GetKind_Turb_Model();
  const SPECIES_MODEL species_model = config->GetKind_Species_Model();

  bool tkeNeeded = (turb_model == TURB_MODEL::SST);

  AD::StartNoSharedReading();

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Retrieve the value of the kinetic energy (if needed) ---*/

    if (turb_model != TURB_MODEL::NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
        DES_LengthScale = solver_container[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
        LES_Mode = solver_container[TURB_SOL]->GetNodes()->GetLES_Mode(iPoint); 
      }
    }

    /*--- Retrieve scalar values (if needed) ---*/
    if (species_model != SPECIES_MODEL::NONE && solver_container[SPECIES_SOL] != nullptr) {
      scalar = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint);
    }

    /*--- Incompressible flow, primitive variables --- */

    bool physical = static_cast<CPBIncNSVariable*>(nodes)->SetPrimVar(iPoint,eddy_visc, turb_ke, GetFluidModel(), scalar);

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Set the DES length scale ---*/

    nodes->SetDES_LengthScale(iPoint,DES_LengthScale);

    /*--- Set the LES sensor ---*/

    nodes->SetLES_Mode(iPoint, LES_Mode);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  return nonPhysicalPoints;

}

/* TODO: The current function does the same for all types of walls as the PB formulation is decoupled from
  the energy equation. This function is very close to the version used by CIncNSSolver and will likely benefit from 
  a base class to avoid code duplication. */
void CPBIncNSSolver::BC_Wall_Generic(const CGeometry *geometry, const CConfig *config,
                                   unsigned short val_marker, unsigned short kind_boundary) {
  SU2_ZONE_SCOPED

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool energy = config->GetEnergy_Equation();
  const bool py_custom = config->GetMarker_All_PyCustom(val_marker);

  /*--- Variables for streamwise periodicity ---*/
  const bool streamwise_periodic = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
  const bool streamwise_periodic_temperature = config->GetStreamwise_Periodic_Temperature();

  /*--- Identify the boundary by string name ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get wall function treatment from config. ---*/

  //const auto Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  // nijso: we do not have a special treatment yet for heated walls
  // the wall function model is written for heat flux, we have to implement isothermal wall conditions
  //if (Wall_Function != WALL_FUNCTIONS::NONE)
  //  SU2_MPI::Error("Wall function treatment not implemented yet.", CURRENT_FUNCTION);

  /*--- Loop over all of the vertices on this boundary marker ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Compute dual-grid area and boundary normal ---*/

    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

    const su2double Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Impose the value of the velocity as a strong boundary
     condition (Dirichlet). Fix the velocity and remove any
     contribution to the residual at this node. ---*/

    if (dynamic_grid) {
      nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
    } else {
      su2double zero[MAXNDIM] = {0.0};
      nodes->SetVelocity_Old(iPoint, zero);
    }

    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      LinSysRes(iPoint, iDim) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

    if (implicit) {
      for (unsigned short iVar = 0; iVar < nDim; iVar++)
        Jacobian.DeleteValsRowi(iPoint, iVar);
    }
  }
  END_SU2_OMP_FOR
}

void CPBIncNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_Wall_Generic(geometry, config, val_marker, HEAT_FLUX);
}

void CPBIncNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_Wall_Generic(geometry, config, val_marker, ISOTHERMAL);
}

void CPBIncNSSolver::BC_HeatTransfer_Wall(const CGeometry *geometry, const CConfig *config, const unsigned short val_marker) {
  SU2_ZONE_SCOPED

  BC_Wall_Generic(geometry, config, val_marker, HEAT_TRANSFER);
}

void CPBIncNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *numerics, CConfig *config) {
  SU2_ZONE_SCOPED
                                      
  Viscous_Residual_impl(iEdge, geometry, solver_container, numerics, config);

}