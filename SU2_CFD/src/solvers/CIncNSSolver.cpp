/*!
 * \file CIncNSSolver.cpp
 * \brief Main subroutines for solving Navier-Stokes incompressible flow.
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

#include "../../include/solvers/CIncNSSolver.hpp"
#include "../../include/variables/CIncNSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

/*--- Explicit instantiation of the parent class of CIncEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CIncEulerVariable, ENUM_REGIME::INCOMPRESSIBLE>;


CIncNSSolver::CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
  CIncEulerSolver(geometry, config, iMesh, true) {

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

void CIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE);
  const bool wall_functions = config->GetWall_Functions();

  /*--- Common preprocessing steps (implemented by CEulerSolver) ---*/

  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

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

  /*--- Compute the TauWall from the wall functions ---*/

  if (wall_functions) {
    SU2_OMP_SAFE_GLOBAL_ACCESS(SetTau_Wall_WF(geometry, solver_container, config);)
    // nijso: we have to set this as well??
    // seteddyviscfirstpoint
  }

  /*--- Compute recovered pressure and temperature for streamwise periodic flow ---*/
  if (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE)
    Compute_Streamwise_Periodic_Recovered_Values(config, geometry, iMesh);
}

void CIncNSSolver::GetStreamwise_Periodic_Properties(const CGeometry *geometry,
                                                     CConfig *config,
                                                     const unsigned short iMesh) {

  /*---------------------------------------------------------------------------------------------*/
  // 1. Evaluate massflow, area avg density & Temperature and Area at streamwise periodic outlet.
  // 2. Update delta_p is target massflow is chosen.
  // 3. Loop Heatflux markers and integrate heat across the boundary. Only if energy equation is on.
  /*---------------------------------------------------------------------------------------------*/

  /*-------------------------------------------------------------------------------------------------*/
  /*--- 1. Evaluate Massflow [kg/s], area-averaged density [kg/m^3] and Area [m^2] at the         ---*/
  /*---    (there can be only one) streamwise periodic outlet/donor marker. Massflow is obviously ---*/
  /*---    needed for prescribed massflow but also for the additional source and heatflux         ---*/
  /*---    boundary terms of the energy equation. Area and the avg-density are used for the       ---*/
  /*---    Pressure-Drop update in case of a prescribed massflow.                                 ---*/
  /*-------------------------------------------------------------------------------------------------*/

  su2double Area_Local            = 0.0,
            MassFlow_Local        = 0.0,
            Average_Density_Local = 0.0,
            Temperature_Local     = 0.0;

  for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- Only "outlet"/donor periodic marker. GetnMarker_Periodic counts from 1: First all the inlet marker from
          the periodic marker pairs and then the all the outlets. I.e. in order to get the outlet of the first pair
          we need to divide the number of periodic markers by 2 and add 1 (because count starts at 1). ---*/
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY &&
        config->GetMarker_All_PerBound(iMarker) == (config->GetnMarker_Periodic()/2 + 1)) {

      for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        if (geometry->nodes->GetDomain(iPoint)) {

          /*--- A = dot_prod(n_A*n_A), with n_A beeing the area-normal. ---*/

          const auto AreaNormal = geometry->vertex[iMarker][iVertex]->GetNormal();

          const su2double FaceArea = GeometryToolbox::Norm(nDim, AreaNormal);

          /*--- m_dot = dot_prod(n*v) * A * rho, with n beeing unit normal. ---*/
          MassFlow_Local += nodes->GetProjVel(iPoint, AreaNormal) * nodes->GetDensity(iPoint);

          Area_Local += FaceArea;

          Average_Density_Local += FaceArea * nodes->GetDensity(iPoint);

          /*--- Due to periodicity, temperatures are equal one the inlet(1) and outlet(2) ---*/
          Temperature_Local += FaceArea * nodes->GetTemperature(iPoint);

        } // if domain
      } // loop vertices
    } // loop periodic boundaries
  } // loop MarkerAll

  // MPI Communication: Sum Area, Sum rho*A & T*A and divide by AreaGlobal, sum massflow
  su2double Area_Global(0), Average_Density_Global(0), MassFlow_Global(0), Temperature_Global(0);
  SU2_MPI::Allreduce(&Area_Local,            &Area_Global,            1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&Average_Density_Local, &Average_Density_Global, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&MassFlow_Local,        &MassFlow_Global,        1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());
  SU2_MPI::Allreduce(&Temperature_Local,     &Temperature_Global,     1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

  Average_Density_Global /= Area_Global;
  Temperature_Global /= Area_Global;

  /*--- Set solver variables ---*/
  SPvals.Streamwise_Periodic_MassFlow = MassFlow_Global;
  SPvals.Streamwise_Periodic_InletTemperature = Temperature_Global;
  SPvals.Streamwise_Periodic_BoundaryArea = Area_Global;
  SPvals.Streamwise_Periodic_AvgDensity = Average_Density_Global;

  if (config->GetEnergy_Equation()) {
    /*---------------------------------------------------------------------------------------------*/
    /*--- 3. Compute the integrated Heatflow [W] for the energy equation source term, heatflux  ---*/
    /*---    boundary term and recovered Temperature. The computation is not completely clear.  ---*/
    /*---    Here the Heatflux from all Boundary markers in the config-file is used.            ---*/
    /*---------------------------------------------------------------------------------------------*/

    su2double HeatFlow_Local = 0.0, HeatFlow_Global = 0.0;

    /*--- Loop over all heatflux Markers ---*/
    for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

      if (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {

        /*--- Identify the boundary by string name and retrive heatflux from config ---*/
        const auto Marker_StringTag = config->GetMarker_All_TagBound(iMarker);
        const su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_StringTag);

        for (auto iVertex = 0ul; iVertex < geometry->nVertex[iMarker]; iVertex++) {

          auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (!geometry->nodes->GetDomain(iPoint)) continue;

          const auto AreaNormal = geometry->vertex[iMarker][iVertex]->GetNormal();

          const su2double FaceArea = GeometryToolbox::Norm(nDim, AreaNormal);

          HeatFlow_Local += FaceArea * (-1.0) * Wall_HeatFlux/config->GetHeat_Flux_Ref();
        } // loop Vertices
      } // loop Heatflux marker
    } // loop AllMarker

    /*--- MPI Communication sum up integrated Heatflux from all processes ---*/
    SU2_MPI::Allreduce(&HeatFlow_Local, &HeatFlow_Global, 1, MPI_DOUBLE, MPI_SUM, SU2_MPI::GetComm());

    /*--- Set the solver variable Integrated Heatflux ---*/
    SPvals.Streamwise_Periodic_IntegratedHeatFlow = HeatFlow_Global;
  } // if energy
}


void CIncNSSolver::Compute_Streamwise_Periodic_Recovered_Values(CConfig *config, const CGeometry *geometry,
                                                                const unsigned short iMesh) {

  const bool energy = (config->GetEnergy_Equation() && config->GetStreamwise_Periodic_Temperature());
  const auto InnerIter = config->GetInnerIter();

  /*--- Reference node on inlet periodic marker to compute relative distance along periodic translation vector. ---*/
  const auto ReferenceNode = geometry->GetStreamwise_Periodic_RefNode();

  /*--- Compute square of the distance between the 2 periodic surfaces. ---*/
  const su2double norm2_translation = GeometryToolbox::SquaredNorm(nDim, config->GetPeriodic_Translation(0));

  /*--- Compute recoverd pressure and temperature for all points ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {

    /*--- First, compute helping terms based on relative distance (0,l) between periodic markers ---*/
    su2double dot_product = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      dot_product += fabs( (geometry->nodes->GetCoord(iPoint,iDim) - ReferenceNode[iDim]) * config->GetPeriodic_Translation(0)[iDim]);

    /*--- Second, substract/add correction from reduced pressure/temperature to get recoverd pressure/temperature ---*/
    const su2double Pressure_Recovered = nodes->GetPressure(iPoint) - SPvals.Streamwise_Periodic_PressureDrop /
                                         norm2_translation * dot_product;
    nodes->SetStreamwise_Periodic_RecoveredPressure(iPoint, Pressure_Recovered);

    /*--- InnerIter > 0 as otherwise MassFlow in the denominator would be zero ---*/
    if (energy && InnerIter > 0) {
      su2double Temperature_Recovered = nodes->GetTemperature(iPoint);
      Temperature_Recovered += SPvals.Streamwise_Periodic_IntegratedHeatFlow /
                              (SPvals.Streamwise_Periodic_MassFlow * nodes->GetSpecificHeatCp(iPoint) * norm2_translation) * dot_product;
      nodes->SetStreamwise_Periodic_RecoveredTemperature(iPoint, Temperature_Recovered);
    }
  } // for iPoint
  END_SU2_OMP_FOR

  /*--- Compute the integrated Heatflux Q into the domain, and massflow over periodic markers ---*/
  SU2_OMP_SAFE_GLOBAL_ACCESS(GetStreamwise_Periodic_Properties(geometry, config, iMesh);)
}

void CIncNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *numerics, CConfig *config) {

  Viscous_Residual_impl(iEdge, geometry, solver_container, numerics, config);
}

unsigned long CIncNSSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0;
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
      }
    }

    /*--- Retrieve scalar values (if needed) ---*/
    if (species_model != SPECIES_MODEL::NONE && solver_container[SPECIES_SOL] != nullptr) {
      scalar = solver_container[SPECIES_SOL]->GetNodes()->GetSolution(iPoint);
    }

    /*--- Incompressible flow, primitive variables --- */

    bool physical = static_cast<CIncNSVariable*>(nodes)->SetPrimVar(iPoint,eddy_visc, turb_ke, GetFluidModel(), scalar);

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Set the DES length scale ---*/

    nodes->SetDES_LengthScale(iPoint,DES_LengthScale);

  }
  END_SU2_OMP_FOR

  AD::EndNoSharedReading();

  return nonPhysicalPoints;

}

void CIncNSSolver::BC_Wall_Generic(const CGeometry *geometry, const CConfig *config,
                                   unsigned short val_marker, unsigned short kind_boundary) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool energy = config->GetEnergy_Equation();
  const bool py_custom = config->GetMarker_All_PyCustom(val_marker);

  /*--- Variables for streamwise periodicity ---*/
  const bool streamwise_periodic = (config->GetKind_Streamwise_Periodic() != ENUM_STREAMWISE_PERIODIC::NONE);
  const bool streamwise_periodic_temperature = config->GetStreamwise_Periodic_Temperature();
  su2double Cp, thermal_conductivity, dot_product, scalar_factor;

  /*--- Identify the boundary by string name ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux, temperature or heat transfer coefficient from config ---*/

  su2double Wall_HeatFlux = 0.0, Twall = 0.0, Tinfinity = 0.0, Transfer_Coefficient = 0.0;

  switch (kind_boundary) {
    case HEAT_FLUX:
      Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
      if (config->GetIntegrated_HeatFlux()) {
        Wall_HeatFlux /= geometry->GetSurfaceArea(config, val_marker);
      }
      break;
    case ISOTHERMAL:
      Twall = config->GetIsothermal_Temperature(Marker_Tag) / config->GetTemperature_Ref();
      break;
    case HEAT_TRANSFER:
      Transfer_Coefficient = config->GetWall_HeatTransfer_Coefficient(Marker_Tag) * config->GetTemperature_Ref() /
                             config->GetHeat_Flux_Ref();
      Tinfinity = config->GetWall_HeatTransfer_Temperature(Marker_Tag) / config->GetTemperature_Ref();
      break;
    default:
      SU2_MPI::Error("Unknown type of boundary condition.", CURRENT_FUNCTION);
      break;
  }

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
      LinSysRes(iPoint, iDim+1) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

    if (implicit) {
      for (unsigned short iVar = 1; iVar <= nDim; iVar++)
        Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
    }

    if (!energy) continue;

    switch(kind_boundary) {
    case HEAT_FLUX:

      /*--- Apply a weak boundary condition for the energy equation.
      Compute the residual due to the prescribed heat flux. ---*/

      if (py_custom) {
        Wall_HeatFlux = geometry->GetCustomBoundaryHeatFlux(val_marker, iVertex) / config->GetHeat_Flux_Ref();
      }
      LinSysRes(iPoint, nDim+1) -= Wall_HeatFlux*Area;

      /*--- With streamwise periodic flow and heatflux walls an additional term is introduced in the boundary formulation ---*/
      if (streamwise_periodic && streamwise_periodic_temperature) {

        Cp = nodes->GetSpecificHeatCp(iPoint);
        thermal_conductivity = nodes->GetThermalConductivity(iPoint);

        /*--- Scalar factor of the residual contribution ---*/
        const su2double norm2_translation = GeometryToolbox::SquaredNorm(nDim, config->GetPeriodic_Translation(0));
        scalar_factor = SPvals.Streamwise_Periodic_IntegratedHeatFlow*thermal_conductivity / (SPvals.Streamwise_Periodic_MassFlow * Cp * norm2_translation);

        /*--- Dot product ---*/
        dot_product = GeometryToolbox::DotProduct(nDim, config->GetPeriodic_Translation(0), Normal);

        LinSysRes(iPoint, nDim+1) += scalar_factor*dot_product;
      } // if streamwise_periodic
      break;

    case HEAT_TRANSFER:
      Twall = nodes->GetTemperature(iPoint);
      Wall_HeatFlux = Transfer_Coefficient * (Tinfinity - Twall);

      /*--- Apply a weak boundary condition for the energy equation.
      Compute the residual due to the prescribed temperature and transfer coefficient.
      Note that for the Heat_Transfer wall, basically a heatflux wall that depends on the local Temperature is applied. ---*/
      LinSysRes(iPoint, nDim+1) -= Wall_HeatFlux*Area;

      if (implicit) {
        Jacobian.AddVal2Diag(iPoint, nDim+1, Transfer_Coefficient*Area);
      }
      break;

    case ISOTHERMAL:
      if (py_custom) {
        Twall = geometry->GetCustomBoundaryTemperature(val_marker, iVertex) / config->GetTemperature_Ref();
      }
      const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      const auto Coord_i = geometry->nodes->GetCoord(iPoint);
      const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
      su2double Edge_Vector[MAXNDIM];
      GeometryToolbox::Distance(nDim, Coord_j, Coord_i, Edge_Vector);
      su2double dist_ij_2 = GeometryToolbox::SquaredNorm(nDim, Edge_Vector);
      su2double dist_ij = sqrt(dist_ij_2);

      /*--- Compute the normal gradient in temperature using Twall ---*/

      su2double dTdn = -(nodes->GetTemperature(Point_Normal) - Twall)/dist_ij;

      /*--- Get thermal conductivity ---*/

      su2double thermal_conductivity = nodes->GetThermalConductivity(iPoint);

      /*--- Apply a weak boundary condition for the energy equation.
      Compute the residual due to the prescribed heat flux. ---*/

      LinSysRes(iPoint, nDim+1) -= thermal_conductivity*dTdn*Area;

      /*--- Jacobian contribution for temperature equation. ---*/

      if (implicit) {
        su2double proj_vector_ij = 0.0;
        if (dist_ij_2 > 0.0)
          proj_vector_ij = GeometryToolbox::DotProduct(nDim, Edge_Vector, Normal) / dist_ij_2;
        Jacobian.AddVal2Diag(iPoint, nDim+1, thermal_conductivity*proj_vector_ij);
      }
      break;
    } // switch
  }
  END_SU2_OMP_FOR
}

void CIncNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {

  BC_Wall_Generic(geometry, config, val_marker, HEAT_FLUX);
}

void CIncNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {

  BC_Wall_Generic(geometry, config, val_marker, ISOTHERMAL);
}

void CIncNSSolver::BC_HeatTransfer_Wall(const CGeometry *geometry, const CConfig *config, const unsigned short val_marker) {

  BC_Wall_Generic(geometry, config, val_marker, HEAT_TRANSFER);
}

void CIncNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                              CConfig *config, unsigned short val_marker) {

  const su2double Temperature_Ref = config->GetTemperature_Ref();
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool energy = config->GetEnergy_Equation();

  /*--- Identify the boundary ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall function treatment.---*/

  if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::NONE) {
    SU2_MPI::Error("Wall function treatment not implemented yet.", CURRENT_FUNCTION);
  }

  /*--- Loop over boundary points ---*/

  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

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
      LinSysRes(iPoint, iDim+1) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/

    if (implicit) {
      for (unsigned short iVar = 1; iVar <= nDim; iVar++)
        Jacobian.DeleteValsRowi(iPoint*nVar+iVar);
      if (energy) Jacobian.DeleteValsRowi(iPoint*nVar+nDim+1);
    }

    if (!energy) continue;

    su2double Tconjugate = GetConjugateHeatVariable(val_marker, iVertex, 0) / Temperature_Ref;
    su2double Twall = 0.0;

    if ((config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX) ||
        (config->GetKind_CHT_Coupling() == CHT_COUPLING::AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

      /*--- Compute closest normal neighbor ---*/

      const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      const auto Coord_i = geometry->nodes->GetCoord(iPoint);
      const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
      su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_j, Coord_i);

      /*--- Compute wall temperature from both temperatures ---*/

      su2double thermal_conductivity = nodes->GetThermalConductivity(iPoint);
      su2double There = nodes->GetTemperature(Point_Normal);
      su2double HF_FactorHere = thermal_conductivity*config->GetViscosity_Ref()/dist_ij;
      su2double HF_FactorConjugate = GetConjugateHeatVariable(val_marker, iVertex, 2);

      Twall = (There*HF_FactorHere + Tconjugate*HF_FactorConjugate)/(HF_FactorHere + HF_FactorConjugate);
    }
    else if ((config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_NEUMANN_HEATFLUX) ||
             (config->GetKind_CHT_Coupling() == CHT_COUPLING::DIRECT_TEMPERATURE_ROBIN_HEATFLUX)) {

      /*--- (Directly) Set wall temperature to conjugate temperature. ---*/

      Twall = Tconjugate;
    }
    else {
      SU2_MPI::Error("Unknown CHT coupling method.", CURRENT_FUNCTION);
    }

    /*--- Strong imposition of the temperature on the fluid zone. ---*/

    LinSysRes(iPoint, nDim+1) = 0.0;
    nodes->SetSolution_Old(iPoint, nDim+1, Twall);
    nodes->SetEnergy_ResTruncError_Zero(iPoint);
  }
  END_SU2_OMP_FOR
}

void CIncNSSolver::SetTau_Wall_WF(CGeometry *geometry, CSolver **solver_container, const CConfig *config) {
  /*--- The wall function implemented herein is based on Nichols and Nelson, AIAA J. v32 n6 2004. ---*/

  unsigned long notConvergedCounter = 0;  /*--- Counts the number of wall cells that are not converged ---*/
  unsigned long smallYPlusCounter = 0;    /*--- Counts the number of wall cells where y+ < 5 ---*/

  const auto max_iter = config->GetwallModel_MaxIter();
  const su2double relax = config->GetwallModel_RelFac();

  /*--- Typical constants from boundary layer theory ---*/

  const su2double kappa = config->GetwallModel_Kappa();
  const su2double B = config->GetwallModel_B();

  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    if (!config->GetViscous_Wall(iMarker)) continue;

    /*--- Identify the boundary by string name ---*/

    const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);

    /*--- Jump to another BC if it is not wall function ---*/

    if (config->GetWallFunction_Treatment(Marker_Tag) != WALL_FUNCTIONS::STANDARD_FUNCTION)
      continue;

    /*--- Loop over all of the vertices on this boundary marker ---*/

    SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
    for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

      const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
      /*--- On the finest mesh compute also on halo nodes to avoid communication of tau wall. ---*/
      if ((!geometry->nodes->GetDomain(iPoint)) && !(MGLevel==MESH_0)) continue;

      /*--- Get coordinates of the current vertex and nearest normal point ---*/

      const auto Coord = geometry->nodes->GetCoord(iPoint);
      const auto Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

      /*--- Compute dual-grid area and boundary normal ---*/

      const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

      const su2double Area = GeometryToolbox::Norm(nDim, Normal);

      su2double UnitNormal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Get the velocity, pressure, and temperature at the nearest
       (normal) interior point. ---*/

      su2double Vel[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        Vel[iDim] = nodes->GetVelocity(Point_Normal,iDim);

      /*--- Compute the wall-parallel velocity at first point off the wall ---*/

      const su2double VelNormal = GeometryToolbox::DotProduct(int(MAXNDIM), Vel, UnitNormal);

      su2double VelTang[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

      const su2double VelTangMod = GeometryToolbox::Norm(int(MAXNDIM), VelTang);

      /*--- Compute normal distance of the interior point from the wall ---*/

      su2double WallDist[MAXNDIM] = {0.0};
      GeometryToolbox::Distance(nDim, Coord, Coord_Normal, WallDist);

      su2double WallDistMod = GeometryToolbox::Norm(int(MAXNDIM), WallDist);

      su2double Density_Wall = nodes->GetDensity(iPoint);

      /*--- Compute the shear stress at the wall in the regular fashion
       *    by using the stress tensor on the surface ---*/

      su2double tau[MAXNDIM][MAXNDIM] = {{0.0}};
      const su2double Lam_Visc_Wall = nodes->GetLaminarViscosity(iPoint);
      su2double Eddy_Visc_Wall = nodes->GetEddyViscosity(iPoint);

      CNumerics::ComputeStressTensor(nDim, tau, nodes->GetVelocityGradient(iPoint), Lam_Visc_Wall);

      su2double TauTangent[MAXNDIM] = {0.0};
      GeometryToolbox::TangentProjection(nDim, tau, UnitNormal, TauTangent);

      const su2double WallShearStress = GeometryToolbox::Norm(int(MAXNDIM), TauTangent);

      /*--- Calculate the quantities from boundary layer theory and
       *    iteratively solve for a new wall shear stress. Use the current wall
       *    shear stress as a starting guess for the wall function. ---*/

      unsigned long counter = 0;
      su2double diff = 1.0;
      su2double U_Tau = max(1.0e-6,sqrt(WallShearStress/Density_Wall));
      /*--- Use minimum y+ as defined in the config, in case the routine below for computing y+ does not converge ---*/
      su2double Y_Plus = 0.99*config->GetwallModel_MinYPlus();

      const su2double Y_Plus_Start = Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall;

      /*--- Automatic switch off when y+ < "limit" according to Nichols & Nelson (2004) ---*/

      if (Y_Plus_Start < config->GetwallModel_MinYPlus()) {
        smallYPlusCounter++;
        continue;
      }

      /*--- Convergence criterium for the Newton solver, note that 1e-10 is too large ---*/
      const su2double tol = 1e-12;
      while (fabs(diff) > tol) {

        /*--- Friction velocity and u+ ---*/

        const su2double U_Plus = VelTangMod / U_Tau;

        /*--- Y+ defined by White & Christoph ---*/

        const su2double kUp = kappa * U_Plus;

        // incompressible adiabatic result
        const su2double Y_Plus_White = exp(kUp) * exp(-kappa * B);

        /*--- Spalding's universal form for the BL velocity with the
         *    outer velocity form of White & Christoph above. ---*/
        Y_Plus = U_Plus + Y_Plus_White + (exp(-kappa * B)* (1.0 - kUp - 0.5 * kUp * kUp - kUp * kUp * kUp / 6.0));

        /*--- incompressible formulation ---*/
        Eddy_Visc_Wall = Lam_Visc_Wall * kappa*exp(-kappa*B) * (exp(kUp) -1.0 - kUp - kUp * kUp / 2.0);

        Eddy_Visc_Wall = max(1.0e-6, Eddy_Visc_Wall);

        /* --- Define function for Newton method to zero --- */

        diff = (Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall) - Y_Plus;

        /* --- Gradient of function defined above wrt U_Tau --- */

        const su2double dyp_dup = 1.0 + exp(-kappa * B) * (kappa * exp(kUp) - kappa - kUp - 0.5 * kUp * kUp);
        const su2double dup_dutau = - U_Plus / U_Tau;
        const su2double grad_diff = Density_Wall * WallDistMod / Lam_Visc_Wall - dyp_dup * dup_dutau;

        /* --- Newton Step --- */

        U_Tau = U_Tau - relax*(diff / grad_diff);

        counter++;
        if (counter > max_iter) {
          notConvergedCounter++;
          // use some safe values for convergence
          Y_Plus = 30.0;
          Eddy_Visc_Wall = 1.0;
          U_Tau = 1.0;
          break;
        }
      }

      /*--- Calculate an updated value for the wall shear stress
       *    using the y+ value, the definition of y+, and the definition of
       *    the friction velocity. ---*/
      YPlus[iMarker][iVertex] = Y_Plus;
      EddyViscWall[iMarker][iVertex] = Eddy_Visc_Wall;
      UTau[iMarker][iVertex] = U_Tau;

      const su2double Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

      /*--- Store this value for the wall shear stress at the node.  ---*/
      nodes->SetTau_Wall(iPoint, Tau_Wall);

    }
    END_SU2_OMP_FOR
  }

  if (config->GetComm_Level() == COMM_FULL) {
    static unsigned long globalCounter1, globalCounter2;

    ompMasterAssignBarrier(globalCounter1,0, globalCounter2,0);

    SU2_OMP_ATOMIC
    globalCounter1 += notConvergedCounter;

    SU2_OMP_ATOMIC
    globalCounter2 += smallYPlusCounter;

    BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
      SU2_MPI::Allreduce(&globalCounter1, &notConvergedCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());
      SU2_MPI::Allreduce(&globalCounter2, &smallYPlusCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

      if (rank == MASTER_NODE) {
        if (notConvergedCounter)
          cout << "Warning: Computation of wall coefficients (y+) did not converge in "
               << notConvergedCounter << " points." << endl;

        if (smallYPlusCounter)
          cout << "Warning: y+ < " << config->GetwallModel_MinYPlus() << " in " << smallYPlusCounter
               << " points, for which the wall model is not active." << endl;
      }
    }
    END_SU2_OMP_SAFE_GLOBAL_ACCESS
  }

}
