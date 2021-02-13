/*!
 * \file CIncNSSolver.cpp
 * \brief Main subroutines for solving Navier-Stokes incompressible flow.
 * \author F. Palacios, T. Economon
 * \version 7.1.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
template class CFVMFlowSolverBase<CIncEulerVariable, INCOMPRESSIBLE>;


CIncNSSolver::CIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
  CIncEulerSolver(geometry, config, iMesh, true) {

  /*--- Read farfield conditions from config ---*/

  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Tke_Inf         = config->GetTke_FreeStreamND();

  /*--- Allocates a 2D array with variable "outer" sizes and init to 0. ---*/

  auto Alloc2D = [](unsigned long M, const unsigned long* N, su2double**& X) {
    X = new su2double* [M];
    for(unsigned long i = 0; i < M; ++i)
      X[i] = new su2double [N[i]] ();
  };

  /*--- Allocates a 3D array with variable "middle" sizes and init to 0. ---*/

  auto Alloc3D = [](unsigned long M, const unsigned long* N, unsigned long P, su2double***& X) {
    X = new su2double** [M];
    for(unsigned long i = 0; i < M; ++i) {
      X[i] = new su2double* [N[i]];
      for(unsigned long j = 0; j < N[i]; ++j)
        X[i][j] = new su2double [P] ();
    }
  };
    
  /*--- Set the SGS model in case an LES simulation is carried out.
   Make a distinction between the SGS models used and set SGSModel and
  SGSModelUsed accordingly. ---*/

  SGSModel = NULL;
  SGSModelUsed = false;

  switch( config->GetKind_SGS_Model() ) {
    /* No LES, so no SGS model needed.
     Set the pointer to NULL and the boolean to false. */
    case NO_SGS_MODEL: case IMPLICIT_LES:
      SGSModel = NULL;
      SGSModelUsed = false;
      break;
    case SMAGORINSKY:
      SGSModel     = new CSmagorinskyModel;
      SGSModelUsed = true;
      break;
    case WALE:
      SGSModel     = new CWALEModel;
      SGSModelUsed = true;
      break;
    case VREMAN:
      SGSModel     = new CVremanModel;
      SGSModelUsed = true;
      break;
    case SIGMA:
      SGSModel     = new CSIGMAModel;
      SGSModelUsed = true;
      break;
    default:
      SU2_MPI::Error("Unknown SGS model encountered", CURRENT_FUNCTION);
  }

  /*--- Set the wall model to NULL ---*/

  WallModel = NULL;

  if (config->GetWall_Models()){

    /*--- Set the WMLES class  ---*/
    /*--- First allocate the auxiliary variables ---*/

    Alloc2D(nMarker, nVertex, TauWall_WMLES);
    Alloc2D(nMarker, nVertex, HeatFlux_WMLES);
    Alloc3D(nMarker, nVertex, nDim, FlowDirTan_WMLES);
    Alloc3D(nMarker, nVertex, nDim, VelTimeFilter_WMLES);

    /*--- Check if the Wall models or Wall functions are unique. ---*/
    /*--- Note: At this moment, all markers must use the same wall model/function ---*/

    vector<unsigned short> WallFunctions_;
    vector<string> WallFunctionsMarker_;
    for(unsigned short iMarker=0; iMarker<nMarker; ++iMarker) {
      switch (config->GetMarker_All_KindBC(iMarker)) {
        case ISOTHERMAL:
        case HEAT_FLUX: {
          string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          if(config->GetWallFunction_Treatment(Marker_Tag) != NO_WALL_FUNCTION){
            WallFunctions_.push_back(config->GetWallFunction_Treatment(Marker_Tag));
            WallFunctionsMarker_.push_back(Marker_Tag);
          }
          break;
        }
        default:  /* Just to avoid a compiler warning. */
          break;
      }
    }

    if (!WallFunctions_.empty()){
      sort(WallFunctions_.begin(), WallFunctions_.end());
      vector<unsigned short>::iterator it = std::unique( WallFunctions_.begin(), WallFunctions_.end() );
      WallFunctions_.erase(it, WallFunctions_.end());

      if(WallFunctions_.size() == 1) {
        switch (config->GetWallFunction_Treatment(WallFunctionsMarker_[0])) {
          case EQUILIBRIUM_WALL_MODEL:
            WallModel = new CWallModel1DEQ(config,WallFunctionsMarker_[0]);
            break;
          case LOGARITHMIC_WALL_MODEL:
            WallModel = new CWallModelLogLaw(config,WallFunctionsMarker_[0]);
            break;
          case ALGEBRAIC_WALL_MODEL:
            WallModel = new CWallModelAlgebraic(config,WallFunctionsMarker_[0]);
            break;
          case APGLL_WALL_MODEL:
            WallModel = new CWallModelAPGLL(config,WallFunctionsMarker_[0]);
            break;
          case TEMPLATE_WALL_MODEL:
            WallModel = new CWallModelTemplate(config,WallFunctionsMarker_[0]);
          break;

          default:
            break;
        }
      }
      else{
        SU2_MPI::Error("Wall function/model type must be the same in all wall BCs", CURRENT_FUNCTION);
      }
    }
  }

  /*--- Initialize the secondary values for direct derivative approxiations ---*/

  switch (config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      break;
  }
}

CIncNSSolver::~CIncNSSolver(void) {

  unsigned short iMarker;

  if (TauWall_WMLES != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      delete [] TauWall_WMLES[iMarker];
    }
    delete [] TauWall_WMLES;
  }

  if (HeatFlux_WMLES != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      delete [] HeatFlux_WMLES[iMarker];
    }
    delete [] HeatFlux_WMLES;
  }

  if (FlowDirTan_WMLES != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (FlowDirTan_WMLES[iMarker] != nullptr) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) delete[] FlowDirTan_WMLES[iMarker][iVertex];
        delete[] FlowDirTan_WMLES[iMarker];
      }
    }
    delete[] FlowDirTan_WMLES;
  }

  if (VelTimeFilter_WMLES != nullptr) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (VelTimeFilter_WMLES[iMarker] != nullptr) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) delete[] VelTimeFilter_WMLES[iMarker][iVertex];
        delete[] VelTimeFilter_WMLES[iMarker];
      }
    }
    delete[] VelTimeFilter_WMLES;
  }

  if( SGSModel  != NULL) delete SGSModel;
  if( WallModel != NULL) delete WallModel;

}

void CIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE);

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

  ComputeVorticityAndStrainMag<1>(*config, iMesh);
  
  /*--- Calculate the eddy viscosity using a SGS model ---*/

  if (SGSModelUsed){
    Setmut_LES(geometry, solver_container, config);
  }
  
  /*--- Compute the wall shear stress from the wall model ---*/

   if (config->GetWall_Models() && (iRKStep==0)){
     SetTauWallHeatFlux_WMLES1stPoint(geometry, solver_container, config, iRKStep);
   }
}

void CIncNSSolver::Viscous_Residual(unsigned long iEdge, CGeometry *geometry, CSolver **solver_container,
                                    CNumerics *numerics, CConfig *config) {

  Viscous_Residual_impl(iEdge, geometry, solver_container, numerics, config);
}

unsigned long CIncNSSolver::SetPrimitive_Variables(CSolver **solver_container, const CConfig *config) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  su2double eddy_visc = 0.0, turb_ke = 0.0, DES_LengthScale = 0.0;
  unsigned short turb_model = config->GetKind_Turb_Model();

  bool tkeNeeded = ((turb_model == SST) || (turb_model == SST_SUST));

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Retrieve the value of the kinetic energy (if needed) ---*/

    if (turb_model != NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

      if (config->GetKind_HybridRANSLES() != NO_HYBRIDRANSLES){
        DES_LengthScale = solver_container[TURB_SOL]->GetNodes()->GetDES_LengthScale(iPoint);
      }
    }
    
    if (turb_model == NONE && SGSModelUsed)
        eddy_visc = solver_container[FLOW_SOL]->GetNodes()->GetEddyViscosity(iPoint);

    /*--- Incompressible flow, primitive variables --- */

    bool physical = static_cast<CIncNSVariable*>(nodes)->SetPrimVar(iPoint,eddy_visc, turb_ke, GetFluidModel());

    /* Check for non-realizable states for reporting. */

    if (!physical) nonPhysicalPoints++;

    /*--- Set the DES length scale ---*/

    nodes->SetDES_LengthScale(iPoint,DES_LengthScale);

  }

  return nonPhysicalPoints;

}

void CIncNSSolver::BC_Wall_Generic(const CGeometry *geometry, const CConfig *config,
                                   unsigned short val_marker, unsigned short kind_boundary) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool energy = config->GetEnergy_Equation();

  /*--- Identify the boundary by string name ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux or temperature from config ---*/

  su2double Wall_HeatFlux = 0.0, Twall = 0.0;

  if (kind_boundary == HEAT_FLUX)
    Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag)/config->GetHeat_Flux_Ref();
  else if (kind_boundary == ISOTHERMAL)
    Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();
  else
    SU2_MPI::Error("Unknown type of boundary condition", CURRENT_FUNCTION);

  /*--- Get wall function treatment from config. ---*/

  const auto Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if (Wall_Function != NO_WALL_FUNCTION)
    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);

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

    if (kind_boundary == HEAT_FLUX) {

      /*--- Apply a weak boundary condition for the energy equation.
      Compute the residual due to the prescribed heat flux. ---*/

      LinSysRes(iPoint, nDim+1) -= Wall_HeatFlux*Area;
    }
    else { // ISOTHERMAL

      auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      auto Coord_i = geometry->nodes->GetCoord(iPoint);
      auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
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
    }
  }
}

void CIncNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {

  BC_Wall_Generic(geometry, config, val_marker, HEAT_FLUX);
}

void CIncNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver**, CNumerics*,
                                    CNumerics*, CConfig *config, unsigned short val_marker) {

  BC_Wall_Generic(geometry, config, val_marker, ISOTHERMAL);
}

void CIncNSSolver::BC_ConjugateHeat_Interface(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                              CConfig *config, unsigned short val_marker) {

  const su2double Temperature_Ref = config->GetTemperature_Ref();
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool energy = config->GetEnergy_Equation();

  /*--- Identify the boundary ---*/

  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall function treatment.---*/

  const auto Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  if (Wall_Function != NO_WALL_FUNCTION) {
    SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);
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

    if ((config->GetKind_CHT_Coupling() == AVERAGED_TEMPERATURE_NEUMANN_HEATFLUX) ||
        (config->GetKind_CHT_Coupling() == AVERAGED_TEMPERATURE_ROBIN_HEATFLUX)) {

      /*--- Compute closest normal neighbor ---*/

      auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Get coordinates of i & nearest normal and compute distance ---*/

      auto Coord_i = geometry->nodes->GetCoord(iPoint);
      auto Coord_j = geometry->nodes->GetCoord(Point_Normal);
      su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_j, Coord_i);

      /*--- Compute wall temperature from both temperatures ---*/

      su2double thermal_conductivity = nodes->GetThermalConductivity(iPoint);
      su2double There = nodes->GetTemperature(Point_Normal);
      su2double HF_FactorHere = thermal_conductivity*config->GetViscosity_Ref()/dist_ij;
      su2double HF_FactorConjugate = GetConjugateHeatVariable(val_marker, iVertex, 2);

      Twall = (There*HF_FactorHere + Tconjugate*HF_FactorConjugate)/(HF_FactorHere + HF_FactorConjugate);
    }
    else if ((config->GetKind_CHT_Coupling() == DIRECT_TEMPERATURE_NEUMANN_HEATFLUX) ||
             (config->GetKind_CHT_Coupling() == DIRECT_TEMPERATURE_ROBIN_HEATFLUX)) {

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
}

void CIncNSSolver::BC_WallModel(CGeometry      *geometry,
                                CSolver        **solver_container,
                                CNumerics      *conv_numerics,
                                CNumerics      *visc_numerics,
                                CConfig        *config,
                                unsigned short val_marker) {

  unsigned short iDim, iVar;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool HeatFlux_Prescribed = false;
  bool energy = config->GetEnergy_Equation();

  /*--- Allocation of variables necessary for convective fluxes. ---*/
  su2double ProjVelocity_i, Wall_HeatFlux;
  su2double *V_reflected, *V_domain;
  
  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  if(config->GetMarker_All_KindBC(val_marker) == HEAT_FLUX) {
    HeatFlux_Prescribed = true;
  }
  
  /*--- Jacobian, initialized to zero if needed. ---*/
  su2double **Jacobian_i = nullptr;
  if (implicit) {
    Jacobian_i = new su2double* [nVar];
    for (auto iVar = 0u; iVar < nVar; iVar++)
      Jacobian_i[iVar] = new su2double [nVar] ();
  }

  /*--- Loop over boundary points ---*/

  for (auto iVertex = 0u; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*-------------------------------------------------------------------------------*/
      /*--- Step 1: For the convective fluxes, create a reflected state of the      ---*/
      /*---         Primitive variables by copying all interior values to the       ---*/
      /*---         reflected. Only the velocity is mirrored for the wall model     ---*/
      /*---         and negative for wall functions (weakly impose v = 0)           ---*/
      /*---         axis. Based on the Upwind_Residual routine.                     ---*/
      /*-------------------------------------------------------------------------------*/

      su2double Normal[MAXNDIM] = {0.0};
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);

      su2double Area = GeometryToolbox::Norm(nDim, Normal);
      
      su2double UnitNormal[MAXNDIM] = {0.0};
      for (auto iDim = 0u; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Allocate the reflected state at the symmetry boundary. ---*/
      V_reflected = GetCharacPrimVar(val_marker, iVertex);

      /*--- Grid movement ---*/
      if (config->GetGrid_Movement())
        conv_numerics->SetGridVel(geometry->nodes->GetGridVel(iPoint), geometry->nodes->GetGridVel(iPoint));

      /*--- Normal vector for this vertex (negate for outward convention). ---*/
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Get current solution at this boundary node ---*/
      V_domain = nodes->GetPrimitive(iPoint);

      /*--- Set the reflected state based on the boundary node. ---*/
      for(iVar = 0; iVar < nPrimVar; iVar++)
        V_reflected[iVar] = nodes->GetPrimitive(iPoint,iVar);

      /*--- Compute velocity in normal direction (ProjVelcity_i=(v*n)) and substract from
       velocity in normal direction: v_r = v - (v*n)n ---*/
      ProjVelocity_i = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        ProjVelocity_i += nodes->GetVelocity(iPoint,iDim)*UnitNormal[iDim];


      if (nodes->GetTauWall_Flag(iPoint)){
        /*--- Scalars are copied and the velocity is mirrored along the wall boundary,
         i.e. the velocity in normal direction is substracted twice. ---*/

        /*--- Force the velocity to be tangential ---*/
        for (iDim = 0; iDim < nDim; iDim++)
          V_reflected[iDim+1] = nodes->GetVelocity(iPoint,iDim) - 2.0 * ProjVelocity_i*UnitNormal[iDim];

        /*--- Set Primitive and Secondary for numerics class. ---*/

        conv_numerics->SetPrimitive(V_domain, V_reflected);
        conv_numerics->SetSecondary(nodes->GetSecondary(iPoint), nodes->GetSecondary(iPoint));

        /*--- Compute the residual using an upwind scheme. ---*/

        auto residual = conv_numerics->ComputeResidual(config);

        LinSysRes.AddBlock(iPoint, residual);

        /*--- Jacobian contribution for implicit integration. ---*/
        if (implicit)
          Jacobian.AddBlock2Diag(iPoint, residual.jacobian_i);
      }
      else{

        /*--- Apply a weak boundary condition for the energy equation.
         Compute the residual due to the prescribed heat flux.
         The convective part will be zero if the grid is not moving. ---*/

        su2double Res_Conv = 0.0;
        su2double Res_Visc = 0.0;
        
        if (energy){
          if (HeatFlux_Prescribed){
            Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag)/config->GetHeat_Flux_Ref();
          }
          else{
            
            /*--- Compute closest normal neighbor ---*/

            const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

            /*--- Get coordinates of i & nearest normal and compute distance ---*/

            const auto Coord_i = geometry->nodes->GetCoord(iPoint);
            const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);

            su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);
            
            /*--- Get transport coefficients ---*/
            
            const su2double Temperature_Ref = config->GetTemperature_Ref();
            const su2double Prandtl_Lam = config->GetPrandtl_Lam();
            const su2double Prandtl_Turb = config->GetPrandtl_Turb();
            const su2double Gas_Constant = config->GetGas_ConstantND();
            const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
            su2double laminar_viscosity    = nodes->GetLaminarViscosity(iPoint);
            su2double eddy_viscosity       = nodes->GetEddyViscosity(iPoint);
            su2double thermal_conductivity = Cp * (laminar_viscosity/Prandtl_Lam + eddy_viscosity/Prandtl_Turb);

            /*--- Get temperatures for heat flux calculation  ---*/
            
            const su2double Twall = config->GetIsothermal_Temperature(Marker_Tag) / Temperature_Ref;
            const su2double There = nodes->GetTemperature(Point_Normal);

            /*--- Compute the normal temperature gradient  ---*/

            su2double dTdn = -(There - Twall)/dist_ij;
            
            Wall_HeatFlux = thermal_conductivity * dTdn;
            
          }
          Res_Visc = Wall_HeatFlux * Area;
        }
        
        /*--- Impose the value of the velocity as a strong boundary
         condition (Dirichlet). Fix the velocity and remove any
         contribution to the residual at this node. ---*/

        if (dynamic_grid) {
          nodes->SetVelocity_Old(iPoint, geometry->nodes->GetGridVel(iPoint));
        }
        else {
          su2double zero[MAXNDIM] = {0.0};
          nodes->SetVelocity_Old(iPoint, zero);
        }

        for (auto iDim = 0u; iDim < nDim; iDim++)
          LinSysRes(iPoint, iDim+1) = 0.0;
        nodes->SetVel_ResTruncError_Zero(iPoint);

        /*--- If the wall is moving, there are additional residual contributions
         due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

        if (dynamic_grid) {
          if (implicit) {
            for (auto iVar = 0u; iVar < nVar; ++iVar)
              Jacobian_i[nDim+1][iVar] = 0.0;
          }

          const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

//          AddDynamicGridResidualContribution(iPoint, Point_Normal, geometry, UnitNormal,
//                                             Area, geometry->nodes->GetGridVel(iPoint),
//                                             Jacobian_i, Res_Conv, Res_Visc);
        }

        /*--- Convective and viscous contributions to the residual at the wall ---*/

        LinSysRes(iPoint, nDim+1) += Res_Conv - Res_Visc;

        /*--- Enforce the no-slip boundary condition in a strong way by
         modifying the velocity-rows of the Jacobian (1 on the diagonal).
         And add the contributions to the Jacobian due to energy. ---*/

        if (implicit) {
          if (dynamic_grid) {
            Jacobian.AddBlock2Diag(iPoint, Jacobian_i);
          }

          for (auto iVar = 1u; iVar <= nDim; iVar++) {
            auto total_index = iPoint*nVar+iVar;
            Jacobian.DeleteValsRowi(total_index);
          }
        }
      }


      /*-------------------------------------------------------*/
      /*-------------------------------------------------------*/
      /*--- Viscous residual contribution of the wall model ---*/
      /*--- TODO: Build the jacobian contribution of the WM ---*/
      /*-------------------------------------------------------*/
      /*-------------------------------------------------------*/

      su2double Res_Viscous[5]={0.0,0.0,0.0,0.0,0.0};
      
      if (nodes->GetTauWall_Flag(iPoint)){

        /*--- Weakly enforce the WM heat flux for the energy equation---*/
        su2double velWall_tan = 0.;
        su2double DirTanWM[MAXNDIM] = {0.};

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          DirTanWM[iDim] = GetFlowDirTan_WMLES(val_marker,iVertex,iDim);

        const su2double TauWall = GetTauWall_WMLES(val_marker,iVertex);
        const su2double Wall_HeatFlux = GetHeatFlux_WMLES(val_marker, iVertex);

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          velWall_tan +=  nodes->GetVelocity(iPoint,iDim) * DirTanWM[iDim];

        Res_Viscous[0] = 0.0;
        Res_Viscous[nDim+1] = 0.0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Res_Viscous[iDim+1] = 0.0;

        for (unsigned short iDim = 0; iDim < nDim; iDim++)
          Res_Viscous[iDim+1] = - TauWall * DirTanWM[iDim] * Area;
        
        if (energy)
          Res_Viscous[nDim+1] = (Wall_HeatFlux - TauWall * velWall_tan) * Area;
      }
      LinSysRes.SubtractBlock(iPoint, Res_Viscous);

    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  if (Jacobian_i)
    for (auto iVar = 0u; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
  delete [] Jacobian_i;
}

void CIncNSSolver::Setmut_LES(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long iPoint;
  su2double Grad_Vel[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  su2double lenScale, muTurb, rho;

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /* Get Density */
    rho = nodes->GetSolution(iPoint, 0);

    /* Velocity Gradients */
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      for (unsigned short jDim = 0 ; jDim < nDim; jDim++)
        Grad_Vel[iDim][jDim] = nodes->GetGradient_Primitive(iPoint, iDim+1, jDim);

    /* Distance to the wall. */
    su2double dist = geometry->nodes->GetWall_Distance(iPoint); // Is the distance to the wall used in any SGS calculation?

    /* Length Scale for the SGS model: Cubic root of the volume. */
    su2double Vol = geometry->nodes->GetVolume(iPoint) + geometry->nodes->GetPeriodicVolume(iPoint);
    lenScale = pow(Vol,1./3.);

    /* Compute the eddy viscosity. */
    if (nDim == 2){
      muTurb = SGSModel->ComputeEddyViscosity_2D(rho, Grad_Vel[0][0], Grad_Vel[1][0],
                                                 Grad_Vel[0][1], Grad_Vel[1][1],
                                                 lenScale, dist);
    }
    else{
      muTurb = SGSModel->ComputeEddyViscosity_3D(rho, Grad_Vel[0][0], Grad_Vel[1][0], Grad_Vel[2][0],
                                               Grad_Vel[0][1], Grad_Vel[1][1], Grad_Vel[2][1],
                                               Grad_Vel[0][2], Grad_Vel[1][2], Grad_Vel[2][2],
                                               lenScale, dist);
    }
    /* Set eddy viscosity. */
    nodes->SetEddyViscosity(iPoint, muTurb);
  }
}

void CIncNSSolver::SetTauWallHeatFlux_WMLES1stPoint(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep) {

  unsigned short iDim, iMarker;
  bool CalculateWallModel = false;

  su2double Vel[3], VelNormal, VelTang[3], VelTangMod, WallDist[3], WallDistMod;
  su2double T_Normal, P_Normal, mu_Normal;
  su2double TimeFilter = config->GetDelta_UnstTimeND()/ (config->GetTimeFilter_WMLES() / config->GetTime_Ref());

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
       (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

     /*--- Identify the boundary by string name ---*/
     string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

     /*--- Identify if this marker is a wall model one---*/
     switch (config->GetWallFunction_Treatment(Marker_Tag)) {
       case EQUILIBRIUM_WALL_MODEL:
       case LOGARITHMIC_WALL_MODEL:
       case ALGEBRAIC_WALL_MODEL:
       case APGLL_WALL_MODEL:
       case TEMPLATE_WALL_MODEL:
         CalculateWallModel = true;
         break;

       case NO_WALL_FUNCTION:
       case STANDARD_WALL_FUNCTION:
         CalculateWallModel = false;
       default:
         break;
     }

     /*--- If not just continue to the next---*/
     if (!CalculateWallModel) continue;

     /*--- Determine the prescribed heat flux or prescribed temperature. ---*/
     bool HeatFlux_Prescribed = false, Temperature_Prescribed = false;
     su2double Wall_HeatFlux = 0.0, Wall_Temperature = 0.0;

     if(config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) {
       HeatFlux_Prescribed = true;
       Wall_HeatFlux       = config->GetWall_HeatFlux(Marker_Tag);
     }
     else {
       Temperature_Prescribed = true;
       Wall_Temperature       = config->GetIsothermal_Temperature(Marker_Tag);
     }

     /*--- Loop over all of the vertices on this boundary marker ---*/
      
     for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {
    
       const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
       const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

       /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
      
       if (!geometry->nodes->GetDomain(iPoint)) continue;

       /*--- Check if the normal point lies on a surface. If yes, skip the wall model calculation.
        For complex geometries, it is possible to have only one element connecting 2 boundaries,
        i.e. very small gaps. This will avoid unphysical shear stress and possible divergence. --*/
       
       if (geometry->nodes->GetPhysicalBoundary(Point_Normal)){
         nodes->SetTauWall_Flag(iPoint,false);
         continue;
       }
      
       /*--- Get coordinates of the current vertex and nearest normal point ---*/

       const auto Coord = geometry->nodes->GetCoord(iPoint);
       const auto Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

       /*--- Compute dual-grid area and boundary normal ---*/

       const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

       su2double Area = GeometryToolbox::Norm(nDim, Normal);

       su2double UnitNormal[MAXNDIM] = {0.0};
       for (auto iDim = 0u; iDim < nDim; iDim++)
         UnitNormal[iDim] = -Normal[iDim]/Area;

       /*--- Compute normal distance of the interior point from the wall ---*/

       for (iDim = 0; iDim < nDim; iDim++)
         WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);

       WallDistMod = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         WallDistMod += WallDist[iDim]*WallDist[iDim];
       WallDistMod = sqrt(WallDistMod);

       /*--- Get the velocity, pressure, and temperature at the nearest
        (normal) interior point. ---*/

       for (iDim = 0; iDim < nDim; iDim++){
         Vel[iDim]   = nodes->GetVelocity(Point_Normal,iDim);
       }

       P_Normal  = nodes->GetPressure(Point_Normal);
       T_Normal  = nodes->GetTemperature(Point_Normal);
       mu_Normal = nodes->GetLaminarViscosity(Point_Normal);

       /*--- Filter the input LES velocity ---*/

       long curAbsTimeIter = (config->GetTimeIter() - config->GetRestart_Iter());
       if (curAbsTimeIter > 0){

         /*--- Old input LES velocity and GradP---*/
         su2double Vel_old[MAXNDIM]   = {0.};
         for (iDim = 0; iDim < nDim; iDim++){
           Vel_old[iDim]   = VelTimeFilter_WMLES[iMarker][iVertex][iDim];
         }
         /*--- Now filter the LES velocity ---*/
         for (iDim = 0; iDim < nDim; iDim++){
           Vel[iDim] = (1.0 - TimeFilter) * Vel_old[iDim] + TimeFilter * Vel[iDim];
         }
       }

       /*--- Update input LES velocity if it is the 1st inner iteration---*/
       if (config->GetInnerIter() == 0){
         for (iDim = 0; iDim < nDim; iDim++){
           VelTimeFilter_WMLES[iMarker][iVertex][iDim] = Vel[iDim];
         }
       }

       /*--- Compute dimensional variables before calling the Wall Model ---*/
       for (iDim = 0; iDim < nDim; iDim++ ){
         Vel[iDim] *= config->GetVelocity_Ref();
       }
       P_Normal *= config->GetPressure_Ref();
       T_Normal *= config->GetTemperature_Ref();
       mu_Normal *= (config->GetPressure_Ref()/config->GetVelocity_Ref());

       /*--- Compute the wall-parallel velocity ---*/

       VelNormal = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         VelNormal += Vel[iDim] * UnitNormal[iDim];
       for (iDim = 0; iDim < nDim; iDim++)
         VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

       VelTangMod = 0.0;
       for (iDim = 0; iDim < nDim; iDim++)
         VelTangMod += VelTang[iDim]*VelTang[iDim];
       VelTangMod = sqrt(VelTangMod);
       VelTangMod = max(VelTangMod,1.e-25);

       su2double dirTan[MAXNDIM] = {0.0};
       for(iDim = 0; iDim<nDim; iDim++) dirTan[iDim] = VelTang[iDim]/VelTangMod;
      
       su2double Rho_Normal = nodes->GetDensity(iPoint) * config->GetDensity_Ref();
       
       P_Normal = Rho_Normal * config->GetGas_Constant() * T_Normal;
       
       /* Compute the wall shear stress and heat flux vector using
        the wall model. */
       su2double tauWall, qWall, ViscosityWall, kOverCvWall;
       WallModel->WallShearStressAndHeatFlux(WallDistMod, T_Normal, VelTangMod, mu_Normal, P_Normal, 0.0,
                                             Wall_HeatFlux, HeatFlux_Prescribed,
                                             Wall_Temperature, Temperature_Prescribed,
                                             GetFluidModel(), tauWall, qWall, ViscosityWall,
                                             kOverCvWall);

       /*--- Compute the non-dimensional values if necessary. ---*/
       tauWall /= config->GetPressure_Ref();
       qWall   /= (config->GetPressure_Ref() * config->GetVelocity_Ref());
       ViscosityWall /= (config->GetPressure_Ref()/config->GetVelocity_Ref());
       nodes->SetLaminarViscosity(iPoint, ViscosityWall);

       /*--- Set tau wall value and flag for flux computation---*/
       nodes->SetTauWall_Flag(iPoint,true);
       nodes->SetTauWall(iPoint, tauWall);
       
       /*--- Set tau wall projected to the flow direction for pos-processing only---*/
       for(iDim = 0; iDim<nDim; iDim++)
         nodes->SetTauWallDir(iPoint, iDim, tauWall*dirTan[iDim]);

       /*--- Set tau wall value and heat flux for boundary conditions---*/
       TauWall_WMLES[iMarker][iVertex] = tauWall;
       HeatFlux_WMLES[iMarker][iVertex] = qWall;
       for (iDim = 0; iDim < nDim; iDim++)
         FlowDirTan_WMLES[iMarker][iVertex][iDim] = dirTan[iDim];

     }
   }
 }
}
