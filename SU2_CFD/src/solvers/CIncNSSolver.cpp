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

  /*--- Initialize the secondary values for direct derivative approxiations ---*/

  switch (config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      break;
  }
}

void CIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                 unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const auto InnerIter = config->GetInnerIter();
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED);
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = (config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE);
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

  ComputeVorticityAndStrainMag<1>(*config, iMesh);


  /*--- Compute the TauWall from the wall functions ---*/

  if (wall_functions) {
    SU2_OMP_MASTER
    SetTauWall_WF(geometry, solver_container, config);
    // nijso: we have to set this as well??
    // seteddyviscfirstpoint
    SU2_OMP_BARRIER
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

  //const auto Wall_Function = config->GetWallFunction_Treatment(Marker_Tag);
  // nijso: we do not have a special treatment yet for heated walls
  // the wall function model is written for heat flux, we have to implement isothermal wall conditions
  //if (Wall_Function != NO_WALL_FUNCTION)
  //  SU2_MPI::Error("Wall function treament not implemented yet", CURRENT_FUNCTION);

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




void CIncNSSolver::SetTauWall_WF(CGeometry *geometry, CSolver **solver_container, const CConfig *config) {
  /*---  The wall function implemented herein is based on Nichols and Nelson AIAAJ v32 n6 2004.
   At this moment, the wall function is only available for adiabatic flows.
   ---*/
  su2double RefDensity,RefVel2; 
  if ((config->GetRef_Inc_NonDim() == DIMENSIONAL) || (config->GetRef_Inc_NonDim() == INITIAL_VALUES)) {
    RefDensity = Density_Inf;
    RefVel2 = 0.0;
    for (auto iDim = 0u; iDim < nDim; iDim++) 
      RefVel2 += Velocity_Inf[iDim] * Velocity_Inf[iDim];
  } else if (config->GetRef_Inc_NonDim() == REFERENCE_VALUES) {
    RefDensity = config->GetInc_Density_Ref();
    RefVel2 = config->GetInc_Velocity_Ref() * config->GetInc_Velocity_Ref();
  }

  unsigned long notConvergedCounter=0;

  su2double grad_diff;
  su2double U_Tau, Y_Plus;
  const su2double Gas_Constant = config->GetGas_ConstantND();
  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  su2double Eddy_Visc;
  constexpr unsigned short max_iter = 50;
  const su2double tol = 1e-12; // 1e-10 is too large
  const su2double relax = 0.5; 
  
    /*--- Compute the recovery factor ---*/
  // Molecular (Laminar) Prandtl number (see Nichols & Nelson, nomenclature )
  const su2double Recovery = pow(config->GetPrandtl_Lam(), (1.0/3.0));

  /*--- Typical constants from boundary layer theory ---*/

  const su2double kappa = 0.41; // put model constants in config file
  const su2double B = 5.5;


  for (auto iMarker = 0u; iMarker < config->GetnMarker_All(); iMarker++) {

    if (!config->GetViscous_Wall(iMarker)) continue;

    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

      /*--- Identify the boundary by string name ---*/

      const auto Marker_Tag = config->GetMarker_All_TagBound(iMarker);

      /*--- Jump to another BC if it is not wall function ---*/

      if (config->GetWallFunction_Treatment(Marker_Tag) != STANDARD_WALL_FUNCTION)
        continue;

      /*--- Get the specified wall heat flux from config ---*/
      // note that we can get the heat flux from the temperature gradient
      su2double q_w = config->GetWall_HeatFlux(Marker_Tag);
      // heat flux from temperature: q_w = h*(T_wall - T_fluid)

      /*--- Loop over all of the vertices on this boundary marker ---*/

      SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
      for (auto iVertex = 0u; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        const auto iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        const auto Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();

        /*--- Check if the node belongs to the domain (i.e, not a halo node)
         and the neighbor is not part of the physical boundary ---*/

        if (!geometry->nodes->GetDomain(iPoint)) continue;

          /*--- Get coordinates of the current vertex and nearest normal point ---*/

          const auto Coord = geometry->nodes->GetCoord(iPoint);
          const auto Coord_Normal = geometry->nodes->GetCoord(Point_Normal);

          /*--- Compute dual-grid area and boundary normal ---*/

          const auto Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

          su2double Area = GeometryToolbox::Norm(nDim, Normal);

          su2double UnitNormal[MAXNDIM] = {0.0};
          for (auto iDim = 0u; iDim < nDim; iDim++)
            UnitNormal[iDim] = -Normal[iDim]/Area;

          /*--- Get the velocity, pressure, and temperature at the nearest
           (normal) interior point. ---*/

          su2double Vel[MAXNDIM] = {0.0};
          for (auto iDim = 0u; iDim < nDim; iDim++)
            Vel[iDim] = nodes->GetVelocity(Point_Normal,iDim);

         // su2double P_Normal = nodes->GetPressure(Point_Normal);
         // su2double T_Normal = nodes->GetTemperature(Point_Normal);

          /*--- Compute the wall-parallel velocity at first point off the wall ---*/
             
          su2double VelNormal = GeometryToolbox::DotProduct(int(MAXNDIM), Vel, UnitNormal);

          su2double VelTang[MAXNDIM] = {0.0};
          for (auto iDim = 0u; iDim < nDim; iDim++)
            VelTang[iDim] = Vel[iDim] - VelNormal*UnitNormal[iDim];

          su2double VelTangMod = GeometryToolbox::Norm(int(MAXNDIM), VelTang);

          /*--- Compute normal distance of the interior point from the wall ---*/

          su2double WallDist[MAXNDIM] = {0.0};
          GeometryToolbox::Distance(nDim, Coord, Coord_Normal, WallDist);
  
          su2double WallDistMod = GeometryToolbox::Norm(int(MAXNDIM), WallDist);

          /*--- Compute mach number ---*/

          // M_Normal = VelTangMod / sqrt(Gamma * Gas_Constant * T_Normal);

          /*--- Compute the wall temperature using the Crocco-Buseman equation ---*/

          //T_Normal = T_Wall * (1.0 + 0.5*Gamma_Minus_One*Recovery*u_normal*u_normal);
          // this means that T_Wall = T_Normal/(1.0+0.5*Gamma_Minus_One*Recovery*u_normal*u_normal)
          //T_Wall = T_Normal/(1.0+0.5*Gamma_Minus_One*Recovery*VelTangMod*VelTangMod);
          // in incompressible flows, we can assume that there is no velocity-related temperature change
          // Prandtl: T+ = Pr*y+
          su2double T_Wall = nodes->GetTemperature(iPoint); 

          /*--- Extrapolate the pressure from the interior & compute the
           wall density using the equation of state ---*/

          /*--- incompressible formulation ---*/
          //su2double P_Wall = nodes->GetPressure(iPoint);
          su2double Density_Wall = nodes->GetDensity(iPoint);
          su2double Conductivity_Wall = nodes->GetThermalConductivity(iPoint);
          //su2double Density_Normal = nodes->GetDensity(Point_Normal);
          su2double Lam_Visc_Normal = nodes->GetLaminarViscosity(Point_Normal);

          /*--- Compute the shear stress at the wall in the regular fashion
           by using the stress tensor on the surface ---*/

          su2double tau[MAXNDIM][MAXNDIM] = {{0.0}}, TauElem[MAXNDIM] = {0.0};
          su2double Lam_Visc_Wall = nodes->GetLaminarViscosity(iPoint);
          CNumerics::ComputeStressTensor(nDim, tau, nodes->GetGradient_Primitive(iPoint)+1, Lam_Visc_Wall);

          for (auto iDim = 0u; iDim < nDim; iDim++) {
            TauElem[iDim] = GeometryToolbox::DotProduct(nDim, tau[iDim], UnitNormal);
          }

          /*--- Compute wall shear stress as the magnitude of the wall-tangential
           component of the shear stress tensor---*/

          su2double TauNormal = GeometryToolbox::DotProduct(nDim, TauElem, UnitNormal);

          su2double TauTangent[MAXNDIM] = {0.0};
          for (auto iDim = 0u; iDim < nDim; iDim++)
            TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];

          su2double WallShearStress = GeometryToolbox::Norm(int(MAXNDIM), TauTangent);


          /*--- Calculate the quantities from boundary layer theory and
           iteratively solve for a new wall shear stress. Use the current wall
           shear stress as a starting guess for the wall function. ---*/

          unsigned long counter = 0; su2double diff = 1.0;
          U_Tau = sqrt(WallShearStress/Density_Wall);
          Y_Plus = 0.0; // to avoid warning

          su2double Y_Plus_Start = Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall;

          /*--- Automatic switch off when y+ < 5 according to Nichols & Nelson (2004) ---*/

          if (Y_Plus_Start < 5.0) {
            continue;
          }

          while (fabs(diff) > tol) {

            /*--- Friction velocity and u+ ---*/

            su2double U_Plus = VelTangMod/U_Tau;

            /*--- Gamma, Beta, Q, and Phi, defined by Nichols & Nelson (2004) page 1110 ---*/

            su2double Gam  = Recovery*U_Tau*U_Tau/(2.0*Cp*T_Wall);
            /*--- nijso: heated wall needs validation testcase! ---*/
            su2double Beta = q_w*Lam_Visc_Wall/(Density_Wall*T_Wall*Conductivity_Wall*U_Tau); // TODO: nonzero heatflux needs validation case
            su2double Q    = sqrt(Beta*Beta + 4.0*Gam);
            su2double Phi  = asin(-1.0*Beta/Q);

            /*--- Y+ defined by White & Christoph (compressibility and heat transfer) negative value for (2.0*Gam*U_Plus - Beta)/Q ---*/

            su2double Y_Plus_White = exp((kappa/sqrt(Gam))*(asin((2.0*Gam*U_Plus - Beta)/Q) - Phi))*exp(-1.0*kappa*B);

            /*--- Spalding's universal form for the BL velocity with the
             outer velocity form of White & Christoph above. ---*/

            su2double kUp = kappa*U_Plus;
            su2double Y_Plus = U_Plus + Y_Plus_White - (exp(-1.0*kappa*B)* (1.0 + kUp + 0.5*kUp*kUp + kUp*kUp*kUp/6.0));

            su2double dypw_dyp = 2.0*Y_Plus_White*(kappa*sqrt(Gam)/Q)*sqrt(1.0 - pow(2.0*Gam*U_Plus - Beta,2.0)/(Q*Q));
 
            Eddy_Visc = Lam_Visc_Wall*(1.0 + dypw_dyp - kappa*exp(-1.0*kappa*B)*
                                             (1.0 + kappa*U_Plus + kappa*kappa*U_Plus*U_Plus/2.0)
                                             - Lam_Visc_Normal/Lam_Visc_Wall);
            Eddy_Visc = max(1.0e-6, Eddy_Visc);

            /* --- Define function for Newton method to zero --- */

            diff = (Density_Wall * U_Tau * WallDistMod / Lam_Visc_Wall) - Y_Plus;

            /* --- Gradient of function defined above --- */

            grad_diff = Density_Wall * WallDistMod / Lam_Visc_Wall + VelTangMod / (U_Tau * U_Tau) +
                      kappa /(U_Tau * sqrt(Gam)) * asin(U_Plus * sqrt(Gam)) * Y_Plus_White -
                      exp(-1.0 * B * kappa) * (0.5 * pow(VelTangMod * kappa / U_Tau, 3) +
                      pow(VelTangMod * kappa / U_Tau, 2) + VelTangMod * kappa / U_Tau) / U_Tau;

            /* --- Newton Step --- */

            U_Tau = U_Tau - relax*(diff / grad_diff);

            counter++;


            if (counter > max_iter) {
              notConvergedCounter++;
              //cout << "Warning: Y+ did not converge within the max number of iterations!" << endl;
              //cout << "diff = " << endl;
              // do not break, use some safe values for convergence
              //break;
              Y_Plus = 30.0;
              Eddy_Visc = 1.0;
              U_Tau = 1.0;
            }

          }

          /*--- Calculate an updated value for the wall shear stress
            using the y+ value, the definition of y+, and the definition of
            the friction velocity. ---*/
          
          YPlus[iMarker][iVertex] = Y_Plus;
          EddyViscWall[iMarker][iVertex] = Eddy_Visc;
          UTau[iMarker][iVertex] = U_Tau;

          // wall model value
          su2double Tau_Wall = (1.0/Density_Wall)*pow(Y_Plus*Lam_Visc_Wall/WallDistMod,2.0);

          for (auto iDim = 0u; iDim < nDim; iDim++)
            CSkinFriction[iMarker][iVertex][iDim] = (Tau_Wall/WallShearStress)*TauTangent[iDim] / (0.5 * RefDensity * RefVel2);


          nodes->SetTauWall(iPoint, Tau_Wall);
          // for compressible flow:
          //nodes->SetTemperature(iPoint,T_Wall);
          //nodes->SetSolution(iPoint, 0, Density_Wall);
          //nodes->SetPrimitive(iPoint, nDim + 1, P_Wall);
          // for incompressible flow:
          // ...? 

      }
    }
  }

  if (notConvergedCounter>0) {
    cout << "Warning: computation of wall coefficients (y+) did not converge in " << notConvergedCounter<< " points"<<endl;
  }

}
