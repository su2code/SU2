/*!
 * \file CEulerSolver.hpp
 * \brief Headers of the CEulerSolver class
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
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

#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../include/solvers/CNEMONSSolver.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

/*--- Explicit instantiation of the parent class of CEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CNEMOEulerVariable, COMPRESSIBLE>;


CNEMONSSolver::CNEMONSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) :
               CNEMOEulerSolver(geometry, config, iMesh, true) {

  Viscosity_Inf      = config->GetViscosity_FreeStreamND();
  Prandtl_Lam        = config->GetPrandtl_Lam();
  Prandtl_Turb       = config->GetPrandtl_Turb();

  /*--- Initialize the secondary values for direct derivative approxiations ---*/
  switch(config->GetDirectDiff()) {
    case D_VISCOSITY:
      SU2_TYPE::SetDerivative(Viscosity_Inf, 1.0);
      break;
    default:
      /*--- Already done upstream. ---*/
      break;
  }

  /* Auxiliary vector for storing primitives for gradient computation in viscous flow */
  /* V = [Y1, ... , Yn, T, Tve, ... ] */
  primitives_aux = new su2double[nPrimVar];

}

CNEMONSSolver::~CNEMONSSolver(void) {

  if (primitives_aux != nullptr) delete [] primitives_aux;

}

void CNEMONSSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) {

  unsigned long iPoint, iVar;
  unsigned short iSpecies, RHO_INDEX, RHOS_INDEX;
  
  auto& gradient = nodes->GetGradient_Primitive();

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Modify species density to mass concentration ---*/
  for ( iPoint = 0; iPoint < nPoint; iPoint++){
    for( iVar = 0; iVar < nPrimVar; iVar++) {
      primitives_aux[iVar] = nodes->GetPrimitive(iPoint, iVar);
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      primitives_aux[RHOS_INDEX+iSpecies] = primitives_aux[RHOS_INDEX+iSpecies]/primitives_aux[RHO_INDEX];
    for( iVar = 0; iVar < nPrimVar; iVar++)
      nodes->SetPrimitive_Aux(iPoint, iVar, primitives_aux[iVar] );
  }

  const auto& primitives = nodes->GetPrimitive_Aux();

  computeGradientsGreenGauss(this, PRIMITIVE_GRADIENT, PERIODIC_PRIM_GG, *geometry,
                             *config, primitives, 0, nPrimVarGrad, gradient);
}

void CNEMONSSolver::SetPrimitive_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) {

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;

  if (reconstruction)
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
  else
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

  
  unsigned long iPoint, iVar;
  unsigned short iSpecies, RHO_INDEX, RHOS_INDEX;

  auto& rmatrix = nodes->GetRmatrix();
  auto& gradient = nodes->GetGradient_Primitive();


  PERIODIC_QUANTITIES kindPeriodicComm = weighted? PERIODIC_PRIM_LS : PERIODIC_PRIM_ULS;

  const auto& primitives = nodes->GetPrimitive();
  
  computeGradientsLeastSquares(this, PRIMITIVE_GRADIENT, kindPeriodicComm, *geometry, *config,
                               weighted, primitives, 0, nPrimVarGrad, gradient, rmatrix);
}

void CNEMONSSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep) {

  bool implicit, err;
  unsigned short iVar, jVar;
  unsigned long iPoint, jPoint, iEdge;

  /*--- Determine time integration scheme ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  CNumerics* numerics = numerics_container[VISC_TERM];

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points, coordinates and normal vector in edge ---*/
    iPoint = geometry->edges->GetNode(iEdge, 0);
    jPoint = geometry->edges->GetNode(iEdge, 1);
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                       geometry->nodes->GetCoord(jPoint) );
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Primitive variables, and gradient ---*/
    numerics->SetConservative   (nodes->GetSolution(iPoint),
                                 nodes->GetSolution(jPoint) );
    numerics->SetConsVarGradient(nodes->GetGradient(iPoint),
                                 nodes->GetGradient(jPoint) );
    numerics->SetPrimitive      (nodes->GetPrimitive(iPoint),
                                 nodes->GetPrimitive(jPoint) );
    numerics->SetPrimVarGradient(nodes->GetGradient_Primitive(iPoint),
                                 nodes->GetGradient_Primitive(jPoint) );

    /*--- Pass supplementary information to CNumerics ---*/
    numerics->SetdPdU  (nodes->GetdPdU(iPoint),   nodes->GetdPdU(jPoint));
    numerics->SetdTdU  (nodes->GetdTdU(iPoint),   nodes->GetdTdU(jPoint));
    numerics->SetdTvedU(nodes->GetdTvedU(iPoint), nodes->GetdTvedU(jPoint));
    numerics->SetEve   (nodes->GetEve(iPoint),    nodes->GetEve(jPoint));
    numerics->SetCvve  (nodes->GetCvve(iPoint),   nodes->GetCvve(jPoint));

    /*--- Species diffusion coefficients ---*/
    numerics->SetDiffusionCoeff(nodes->GetDiffusionCoeff(iPoint),
                                nodes->GetDiffusionCoeff(jPoint) );

    /*--- Laminar viscosity ---*/
    numerics->SetLaminarViscosity(nodes->GetLaminarViscosity(iPoint),
                                  nodes->GetLaminarViscosity(jPoint) );

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
                                     nodes->GetThermalConductivity(jPoint));

    /*--- Vib-el. thermal conductivity ---*/
    numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                        nodes->GetThermalConductivity_ve(jPoint) );

    /*--- Compute and update residual ---*/
    numerics->ComputeResidual(Res_Visc, Jacobian_i, Jacobian_j, config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (Res_Visc[iVar] != Res_Visc[iVar]) err = true;
    if (implicit)
      for (iVar = 0; iVar < nVar; iVar++)
        for (jVar = 0; jVar < nVar; jVar++)
          if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
              (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
            err = true;

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
      LinSysRes.AddBlock(jPoint, Res_Visc);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
        Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
        Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
      }
    }
  } //iEdge
}

void CNEMONSSolver::Friction_Forces(CGeometry *geometry, CConfig *config) {

//  unsigned short Boundary, Monitoring, iMarker, iMarker_Monitoring, iDim, jDim;
//  unsigned short VEL_INDEX, T_INDEX, TVE_INDEX;
//  unsigned long iVertex, iPoint, iPointNormal;
//  su2double **Grad_PrimVar, Delta, Viscosity, ThermalCond, ThermalCond_ve, TauNormal,
//  FrictionVel, *Normal, *Coord, *Coord_Normal, Area, Force[3], MomentDist[3],
//  RefDensity, Density, div_vel, RefVel2, dTn, dTven, pnorm, Alpha, Beta, RefLength,
//  RefArea, *Origin, factor, MaxNorm = 8.0, WallShearStress, WallDistMod, WallDist[3],
//  Vel[3], Velocity_Inf[3], UnitNormal[3], TauElem[3], TauTangent[3], Tau[3][3],
//  MomentX_Force[3], MomentY_Force[3], MomentZ_Force[3];
//  string Marker_Tag, Monitoring_Tag;
//
//#ifdef HAVE_MPI
//  su2double MyAllBound_CD_Visc, MyAllBound_CL_Visc, MyAllBound_CSF_Visc, MyAllBound_CMx_Visc, MyAllBound_CMy_Visc, MyAllBound_CMz_Visc, MyAllBound_CoPx_Visc, MyAllBound_CoPy_Visc, MyAllBound_CoPz_Visc, MyAllBound_CFx_Visc, MyAllBound_CFy_Visc, MyAllBound_CFz_Visc, MyAllBound_CT_Visc, MyAllBound_CQ_Visc, MyAllBound_HF_Visc, MyAllBound_MaxHF_Visc, *MySurface_CL_Visc = nullptr, *MySurface_CD_Visc = nullptr, *MySurface_CSF_Visc = nullptr, *MySurface_CEff_Visc = nullptr, *MySurface_CFx_Visc = nullptr, *MySurface_CFy_Visc = nullptr, *MySurface_CFz_Visc = nullptr, *MySurface_CMx_Visc = nullptr, *MySurface_CMy_Visc = nullptr, *MySurface_CMz_Visc = nullptr, *MySurface_HF_Visc = nullptr, *MySurface_MaxHF_Visc;
//#endif
//
//  /*--- Retrieve index information from CVariable ---*/
//  VEL_INDEX = nodes->GetVelIndex();
//  T_INDEX   = nodes->GetTIndex();
//  TVE_INDEX = nodes->GetTveIndex();
//
//  /*--- Retrieve data from CConfig ---*/
//  pnorm = 0.0; //cat: getbackpnorm config->GetPnormHeat();
//
//  /*--- Calculate angle of attack & sideslip ---*/
//  Alpha = config->GetAoA()*PI_NUMBER/180.0;
//  Beta  = config->GetAoS()*PI_NUMBER/180.0;
//
//  /*--- Determine reference geometrical parameters ---*/
//  RefArea    = config->GetRefArea();
//  RefLength  = config->GetRefLength();
//  Origin     = config->GetRefOriginMoment(0);
//
//  /*--- Get reference values from the freestream node. ---*/
//  RefVel2 = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    Velocity_Inf[iDim] = node_infty->GetVelocity(0,iDim);
//    RefVel2 += Velocity_Inf[iDim]*Velocity_Inf[iDim];
//  }
//  RefDensity  = node_infty->GetDensity(0);
//
//  factor = 1.0 / (0.5*RefDensity*RefArea*RefVel2);
//
//  /*-- Initialization --*/
//  AllBound_CMx_Visc  = 0.0; AllBound_CMy_Visc   = 0.0; AllBound_CMz_Visc = 0.0;
//  AllBound_CFx_Visc  = 0.0; AllBound_CFy_Visc   = 0.0; AllBound_CFz_Visc = 0.0;
//  AllBound_CD_Visc   = 0.0; AllBound_CL_Visc    = 0.0;
//  AllBound_HF_Visc   = 0.0; AllBound_MaxHF_Visc = 0.0;
//  AllBound_CEff_Visc = 0.0;
//
//  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
//    Surface_CL_Visc[iMarker_Monitoring]  = 0.0; Surface_CD_Visc[iMarker_Monitoring]    = 0.0;
//    Surface_CSF_Visc[iMarker_Monitoring] = 0.0; Surface_CEff_Visc[iMarker_Monitoring]  = 0.0;
//    Surface_CFx_Visc[iMarker_Monitoring] = 0.0; Surface_CFy_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CFz_Visc[iMarker_Monitoring] = 0.0; Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CMy_Visc[iMarker_Monitoring] = 0.0; Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_HF_Visc[iMarker_Monitoring]  = 0.0; Surface_MaxHF_Visc[iMarker_Monitoring] = 0.0;
//  }
//
//  /*--- Loop over the Navier-Stokes markers ---*/
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    /*--- Identify boundary information ---*/
//    Boundary   = config->GetMarker_All_KindBC(iMarker);
//    Monitoring = config->GetMarker_All_Monitoring(iMarker);
//
//    /*--- Obtain the origin for the moment computation for a particular marker ---*/
//    if (Monitoring == YES) {
//      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
//        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
//        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
//        if (Marker_Tag == Monitoring_Tag)
//          Origin = config->GetRefOriginMoment(iMarker_Monitoring);
//      }
//    }
//
//    /*--- Solve for the forces ---*/
//    if ((Boundary == HEAT_FLUX              ) ||
//        (Boundary == HEAT_FLUX_CATALYTIC    ) ||
//        (Boundary == HEAT_FLUX_NONCATALYTIC ) ||
//        (Boundary == ISOTHERMAL             ) ||
//        (Boundary == ISOTHERMAL_CATALYTIC   ) ||
//        (Boundary == ISOTHERMAL_NONCATALYTIC) ||
//        (Boundary == SMOLUCHOWSKI_MAXWELL)) {
//
//      /*--- Forces initialization at each Marker ---*/
//      CD_Visc[iMarker]   = 0.0; CL_Visc[iMarker]    = 0.0; CSF_Visc[iMarker]    = 0.0;
//      CFx_Visc[iMarker]  = 0.0; CFy_Visc[iMarker]   = 0.0; CFz_Visc[iMarker]    = 0.0;
//      CMx_Visc[iMarker]  = 0.0; CMy_Visc[iMarker]   = 0.0; CMz_Visc[iMarker]    = 0.0;
//      CoPx_Visc[iMarker] = 0.0; CoPy_Visc[iMarker]  = 0.0; CoPz_Visc[iMarker]   = 0.0;
//      CT_Visc[iMarker]   = 0.0; CQ_Visc[iMarker]    = 0.0; CMerit_Visc[iMarker] = 0.0;
//      HF_Visc[iMarker]   = 0.0; MaxHF_Visc[iMarker] = 0.0; CEff_Visc[iMarker]   = 0.0;
//
//      for (iDim = 0; iDim < nDim; iDim++) {
//        ForceViscous[iDim]  = 0.0; MomentViscous[iDim] = 0.0;
//        MomentX_Force[iDim] = 0.0; MomentY_Force[iDim] = 0.0; MomentZ_Force[iDim] = 0.0;
//      }
//
//      /*--- Loop over the boundary points ---*/
//      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
//
//        /*--- Acquire & calculate geometric parameters ---*/
//        iPoint       = geometry->vertex[iMarker][iVertex]->GetNode();
//        iPointNormal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
//        Coord        = geometry->nodes->GetCoord(iPoint);
//        Coord_Normal = geometry->nodes->GetCoord(iPointNormal);
//        Normal       = geometry->vertex[iMarker][iVertex]->GetNormal();
//
//        Area         = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          Area += Normal[iDim]*Normal[iDim];
//        Area = sqrt(Area);
//
//        for (iDim = 0; iDim < nDim; iDim++) {
//          UnitNormal[iDim] = Normal[iDim]/Area;
//          MomentDist[iDim] = Coord[iDim] - Origin[iDim];
//        }
//
//        /*--- Get vertex flow parameters ---*/
//        Grad_PrimVar   = nodes->GetGradient_Primitive(iPoint);
//        Viscosity      = nodes->GetLaminarViscosity(iPoint);
//        ThermalCond    = nodes->GetThermalConductivity(iPoint);
//        ThermalCond_ve = nodes->GetThermalConductivity_ve(iPoint);
//        Density        = nodes->GetDensity(iPoint);
//
//        /*--- Calculate the viscous stress tensor ---*/
//        div_vel = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          div_vel += Grad_PrimVar[VEL_INDEX+iDim][iDim];
//
//        for (iDim = 0; iDim < nDim; iDim++) {
//          for (jDim = 0 ; jDim < nDim; jDim++) {
//            Delta = 0.0; if (iDim == jDim) Delta = 1.0;
//            Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[VEL_INDEX+jDim][iDim] +
//                Grad_PrimVar[VEL_INDEX+iDim][jDim]  )
//                - TWO3*Viscosity*div_vel*Delta;
//          }
//          TauElem[iDim] = 0.0;
//          for (jDim = 0; jDim < nDim; jDim++)
//            TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
//        }
//
//        /*--- Compute wall shear stress (using the stress tensor) ---*/
//        TauNormal = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          TauNormal += TauElem[iDim] * UnitNormal[iDim];
//        for (iDim = 0; iDim < nDim; iDim++)
//          TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
//
//        WallShearStress = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          WallShearStress += TauTangent[iDim]*TauTangent[iDim];
//        WallShearStress = sqrt(WallShearStress);
//
//        for (iDim = 0; iDim < nDim; iDim++)
//          WallDist[iDim] = (Coord[iDim] - Coord_Normal[iDim]);
//        WallDistMod = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++)
//          WallDistMod += WallDist[iDim]*WallDist[iDim];
//        WallDistMod = sqrt(WallDistMod);
//
//        /*--- Compute wall skin friction coefficient, and heat flux on the wall ---*/
//        for (iDim = 0; iDim < nDim; iDim++)
//          CSkinFriction[iMarker][iDim][iVertex] = TauTangent[iDim] / (0.5*RefDensity*RefVel2);
//
//        /*--- Compute y+ and non-dimensional velocity ---*/
//        FrictionVel = sqrt(fabs(WallShearStress)/Density);
//        YPlus[iMarker][iVertex] = WallDistMod*FrictionVel/(Viscosity/Density);
//
//        /*--- Compute heat flux on the wall ---*/
//        dTn = 0.0; dTven = 0.0;
//        for (iDim = 0; iDim < nDim; iDim++) {
//          dTn   += Grad_PrimVar[T_INDEX][iDim]*UnitNormal[iDim];
//          dTven += Grad_PrimVar[TVE_INDEX][iDim]*UnitNormal[iDim];
//        }
//
//        HeatFlux[iMarker][iVertex] = ThermalCond*dTn + ThermalCond_ve*dTven;
//        HF_Visc[iMarker] += HeatFlux[iMarker][iVertex]*Area;
//        MaxHF_Visc[iMarker] += pow(HeatFlux[iMarker][iVertex], pnorm)*Area;
//
//        /*--- Compute viscous forces, and moment using the stress tensor ---*/
//        if ((geometry->nodes->GetDomain(iPoint)) && (Monitoring == YES)) {
//
//          for (iDim = 0; iDim < nDim; iDim++) {
//            Force[iDim] = TauElem[iDim]*Area*factor;
//            ForceViscous[iDim] += Force[iDim];
//          }
//
//          if (iDim == 3) {
//            MomentViscous[0] += (Force[2]*MomentDist[1] - Force[1]*MomentDist[2])/RefLength;
//            MomentX_Force[1] += (-Force[1]*Coord[2]);
//            MomentX_Force[2] += (Force[2]*Coord[1]);
//
//            MomentViscous[1] += (Force[0]*MomentDist[2] - Force[2]*MomentDist[0])/RefLength;
//            MomentY_Force[2] += (-Force[2]*Coord[0]);
//            MomentY_Force[0] += (Force[0]*Coord[2]);
//          }
//          MomentViscous[2] += (Force[1]*MomentDist[0] - Force[0]*MomentDist[1])/RefLength;
//          MomentZ_Force[0] += (-Force[0]*Coord[1]);
//          MomentZ_Force[1] += (Force[1]*Coord[0]);
//        }
//      }
//
//      /*--- Transform ForceInviscid into CLift and CDrag ---*/
//      if  (Monitoring == YES) {
//
//        if (nDim == 2) {
//          CD_Visc[iMarker]     =  ForceViscous[0]*cos(Alpha) + ForceViscous[1]*sin(Alpha);
//          CL_Visc[iMarker]     = -ForceViscous[0]*sin(Alpha) + ForceViscous[1]*cos(Alpha);
//          CEff_Visc[iMarker]   = CL_Visc[iMarker] / (CD_Visc[iMarker]+EPS);
//          CFx_Visc[iMarker]    = ForceViscous[0];
//          CFy_Visc[iMarker]    = ForceViscous[1];
//          CMz_Visc[iMarker]    = MomentViscous[2];
//          CoPx_Visc[iMarker]   = MomentZ_Force[1];
//          CoPy_Visc[iMarker]   = -MomentZ_Force[0];
//          CT_Visc[iMarker]     = -CFx_Visc[iMarker];
//          CQ_Visc[iMarker]     = -CMz_Visc[iMarker];
//          CMerit_Visc[iMarker] = CT_Visc[iMarker] / (CQ_Visc[iMarker]+EPS);
//          MaxHF_Visc[iMarker]  = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
//        }
//
//        if (nDim == 3) {
//          CD_Visc[iMarker]     =  ForceViscous[0]*cos(Alpha)*cos(Beta) + ForceViscous[1]*sin(Beta) + ForceViscous[2]*sin(Alpha)*cos(Beta);
//          CL_Visc[iMarker]     = -ForceViscous[0]*sin(Alpha) + ForceViscous[2]*cos(Alpha);
//          CSF_Visc[iMarker]    = -ForceViscous[0]*sin(Beta)*cos(Alpha) + ForceViscous[1]*cos(Beta) - ForceViscous[2]*sin(Beta)*sin(Alpha);
//          CEff_Visc[iMarker]   = CL_Visc[iMarker]/(CD_Visc[iMarker] + EPS);
//          CFx_Visc[iMarker]    = ForceViscous[0];
//          CFy_Visc[iMarker]    = ForceViscous[1];
//          CFz_Visc[iMarker]    = ForceViscous[2];
//          CMx_Visc[iMarker]    = MomentViscous[0];
//          CMy_Visc[iMarker]    = MomentViscous[1];
//          CMz_Visc[iMarker]    = MomentViscous[2];
//          CoPx_Visc[iMarker]   =  -MomentY_Force[0];
//          CoPz_Visc[iMarker]   = MomentY_Force[2];
//          CT_Visc[iMarker]     = -CFz_Visc[iMarker];
//          CQ_Visc[iMarker]     = -CMz_Visc[iMarker];
//          CMerit_Visc[iMarker] = CT_Visc[iMarker] / (CQ_Visc[iMarker] + EPS);
//          MaxHF_Visc[iMarker]  = pow(MaxHF_Visc[iMarker], 1.0/MaxNorm);
//        }
//
//        AllBound_CD_Visc    += CD_Visc[iMarker];
//        AllBound_CL_Visc    += CL_Visc[iMarker];
//        AllBound_CSF_Visc   += CSF_Visc[iMarker];
//        AllBound_CFx_Visc   += CFx_Visc[iMarker];
//        AllBound_CFy_Visc   += CFy_Visc[iMarker];
//        AllBound_CFz_Visc   += CFz_Visc[iMarker];
//        AllBound_CMx_Visc   += CMx_Visc[iMarker];
//        AllBound_CMy_Visc   += CMy_Visc[iMarker];
//        AllBound_CMz_Visc   += CMz_Visc[iMarker];
//        AllBound_CoPx_Visc  += CoPx_Visc[iMarker];
//        AllBound_CoPy_Visc  += CoPy_Visc[iMarker];
//        AllBound_CoPz_Visc  += CoPz_Visc[iMarker];
//        AllBound_CT_Visc    += CT_Visc[iMarker];
//        AllBound_CQ_Visc    += CQ_Visc[iMarker];
//        AllBound_HF_Visc    += HF_Visc[iMarker];
//        AllBound_MaxHF_Visc += pow(MaxHF_Visc[iMarker], MaxNorm);
//
//        /*--- Compute the coefficients per surface ---*/
//        for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
//          Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
//          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
//          if (Marker_Tag == Monitoring_Tag) {
//            Surface_CL_Visc[iMarker_Monitoring]    += CL_Visc[iMarker];
//            Surface_CD_Visc[iMarker_Monitoring]    += CD_Visc[iMarker];
//            Surface_CSF_Visc[iMarker_Monitoring]   += CSF_Visc[iMarker];
//            Surface_CEff_Visc[iMarker_Monitoring]  += CEff_Visc[iMarker];
//            Surface_CFx_Visc[iMarker_Monitoring]   += CFx_Visc[iMarker];
//            Surface_CFy_Visc[iMarker_Monitoring]   += CFy_Visc[iMarker];
//            Surface_CFz_Visc[iMarker_Monitoring]   += CFz_Visc[iMarker];
//            Surface_CMx_Visc[iMarker_Monitoring]   += CMx_Visc[iMarker];
//            Surface_CMy_Visc[iMarker_Monitoring]   += CMy_Visc[iMarker];
//            Surface_CMz_Visc[iMarker_Monitoring]   += CMz_Visc[iMarker];
//            Surface_HF_Visc[iMarker_Monitoring]    += HF_Visc[iMarker];
//            Surface_MaxHF_Visc[iMarker_Monitoring] += pow(MaxHF_Visc[iMarker],MaxNorm);
//          }
//        }
//      }
//    }
//  }
//
//  /*--- Update some global coeffients ---*/
//  AllBound_CEff_Visc   = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
//  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
//  AllBound_MaxHF_Visc  = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
//
//#ifdef HAVE_MPI
//
//  /*--- Add AllBound information using all the nodes ---*/
//  MyAllBound_CD_Visc    = AllBound_CD_Visc;   AllBound_CD_Visc = 0.0;
//  MyAllBound_CL_Visc    = AllBound_CL_Visc;   AllBound_CL_Visc = 0.0;
//  MyAllBound_CSF_Visc   = AllBound_CSF_Visc;  AllBound_CSF_Visc = 0.0;
//  AllBound_CEff_Visc    = 0.0;
//  MyAllBound_CMx_Visc   = AllBound_CMx_Visc;  AllBound_CMx_Visc = 0.0;
//  MyAllBound_CMy_Visc   = AllBound_CMy_Visc;  AllBound_CMy_Visc = 0.0;
//  MyAllBound_CMz_Visc   = AllBound_CMz_Visc;  AllBound_CMz_Visc = 0.0;
//  MyAllBound_CoPx_Visc  = AllBound_CoPx_Visc; AllBound_CoPx_Visc = 0.0;
//  MyAllBound_CoPy_Visc  = AllBound_CoPy_Visc; AllBound_CoPy_Visc = 0.0;
//  MyAllBound_CoPz_Visc  = AllBound_CoPz_Visc; AllBound_CoPz_Visc = 0.0;
//  MyAllBound_CFx_Visc   = AllBound_CFx_Visc;  AllBound_CFx_Visc = 0.0;
//  MyAllBound_CFy_Visc   = AllBound_CFy_Visc;  AllBound_CFy_Visc = 0.0;
//  MyAllBound_CFz_Visc   = AllBound_CFz_Visc;  AllBound_CFz_Visc = 0.0;
//  MyAllBound_CT_Visc    = AllBound_CT_Visc;   AllBound_CT_Visc = 0.0;
//  MyAllBound_CQ_Visc    = AllBound_CQ_Visc;   AllBound_CQ_Visc = 0.0;
//  AllBound_CMerit_Visc  = 0.0;
//  MyAllBound_HF_Visc    = AllBound_HF_Visc;   AllBound_HF_Visc = 0.0;
//  MyAllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, MaxNorm);
//  AllBound_MaxHF_Visc = 0.0;
//
//  SU2_MPI::Allreduce(&MyAllBound_CD_Visc, &AllBound_CD_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CL_Visc, &AllBound_CL_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CSF_Visc, &AllBound_CSF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  AllBound_CEff_Visc = AllBound_CL_Visc / (AllBound_CD_Visc + EPS);
//  SU2_MPI::Allreduce(&MyAllBound_CMx_Visc, &AllBound_CMx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CMy_Visc, &AllBound_CMy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CMz_Visc, &AllBound_CMz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CFx_Visc, &AllBound_CFx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CFy_Visc, &AllBound_CFy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CFz_Visc, &AllBound_CFz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CoPx_Visc, &AllBound_CoPx_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CoPy_Visc, &AllBound_CoPy_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CoPz_Visc, &AllBound_CoPz_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CT_Visc, &AllBound_CT_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_CQ_Visc, &AllBound_CQ_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  AllBound_CMerit_Visc = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
//  SU2_MPI::Allreduce(&MyAllBound_HF_Visc, &AllBound_HF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(&MyAllBound_MaxHF_Visc, &AllBound_MaxHF_Visc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  AllBound_MaxHF_Visc = pow(AllBound_MaxHF_Visc, 1.0/MaxNorm);
//
//  /*--- Add the forces on the surfaces using all the nodes ---*/
//
//  MySurface_CL_Visc    = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CD_Visc    = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CSF_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CEff_Visc  = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CFx_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CFy_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CFz_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CMx_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CMy_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_CMz_Visc   = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_HF_Visc    = new su2double[config->GetnMarker_Monitoring()];
//  MySurface_MaxHF_Visc = new su2double[config->GetnMarker_Monitoring()];
//
//  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
//
//    MySurface_CL_Visc[iMarker_Monitoring]    = Surface_CL_Visc[iMarker_Monitoring];
//    MySurface_CD_Visc[iMarker_Monitoring]    = Surface_CD_Visc[iMarker_Monitoring];
//    MySurface_CSF_Visc[iMarker_Monitoring]   = Surface_CSF_Visc[iMarker_Monitoring];
//    MySurface_CEff_Visc[iMarker_Monitoring]  = Surface_CEff_Visc[iMarker_Monitoring];
//    MySurface_CFx_Visc[iMarker_Monitoring]   = Surface_CFx_Visc[iMarker_Monitoring];
//    MySurface_CFy_Visc[iMarker_Monitoring]   = Surface_CFy_Visc[iMarker_Monitoring];
//    MySurface_CFz_Visc[iMarker_Monitoring]   = Surface_CFz_Visc[iMarker_Monitoring];
//    MySurface_CMx_Visc[iMarker_Monitoring]   = Surface_CMx_Visc[iMarker_Monitoring];
//    MySurface_CMy_Visc[iMarker_Monitoring]   = Surface_CMy_Visc[iMarker_Monitoring];
//    MySurface_CMz_Visc[iMarker_Monitoring]   = Surface_CMz_Visc[iMarker_Monitoring];
//    MySurface_HF_Visc[iMarker_Monitoring]    = Surface_HF_Visc[iMarker_Monitoring];
//    MySurface_MaxHF_Visc[iMarker_Monitoring] = Surface_MaxHF_Visc[iMarker_Monitoring];
//
//    Surface_CL_Visc[iMarker_Monitoring]    = 0.0;
//    Surface_CD_Visc[iMarker_Monitoring]    = 0.0;
//    Surface_CSF_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CEff_Visc[iMarker_Monitoring]  = 0.0;
//    Surface_CFx_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CFy_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CFz_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CMx_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CMy_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_CMz_Visc[iMarker_Monitoring]   = 0.0;
//    Surface_HF_Visc[iMarker_Monitoring]    = 0.0;
//    Surface_MaxHF_Visc[iMarker_Monitoring] = 0.0;
//  }
//
//  SU2_MPI::Allreduce(MySurface_CL_Visc, Surface_CL_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CD_Visc, Surface_CD_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CSF_Visc, Surface_CSF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++)
//    Surface_CEff_Visc[iMarker_Monitoring] = Surface_CL_Visc[iMarker_Monitoring] / (Surface_CD_Visc[iMarker_Monitoring] + EPS);
//  SU2_MPI::Allreduce(MySurface_CFx_Visc, Surface_CFx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CFy_Visc, Surface_CFy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CFz_Visc, Surface_CFz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CMx_Visc, Surface_CMx_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CMy_Visc, Surface_CMy_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_CMz_Visc, Surface_CMz_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_HF_Visc, Surface_HF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//  SU2_MPI::Allreduce(MySurface_MaxHF_Visc, Surface_MaxHF_Visc, config->GetnMarker_Monitoring(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//
//  delete [] MySurface_CL_Visc;   delete [] MySurface_CD_Visc;  delete [] MySurface_CSF_Visc;
//  delete [] MySurface_CEff_Visc; delete [] MySurface_CFx_Visc; delete [] MySurface_CFy_Visc;
//  delete [] MySurface_CFz_Visc;  delete [] MySurface_CMx_Visc; delete [] MySurface_CMy_Visc;
//  delete [] MySurface_CMz_Visc;  delete [] MySurface_HF_Visc;  delete [] MySurface_MaxHF_Visc;
//
//#endif
//
//  /*--- Update the total coefficients (note that all the nodes have the same value)---*/
//
//  Total_CD      += AllBound_CD_Visc;
//  Total_CL      += AllBound_CL_Visc;
//  Total_CSF     += AllBound_CSF_Visc;
//  Total_CEff     = Total_CL / (Total_CD + EPS);
//  Total_CFx     += AllBound_CFx_Visc;
//  Total_CFy     += AllBound_CFy_Visc;
//  Total_CFz     += AllBound_CFz_Visc;
//  Total_CMx     += AllBound_CMx_Visc;
//  Total_CMy     += AllBound_CMy_Visc;
//  Total_CMz     += AllBound_CMz_Visc;
//  Total_CoPx    += AllBound_CoPx_Visc;
//  Total_CoPy    += AllBound_CoPy_Visc;
//  Total_CoPz    += AllBound_CoPz_Visc;
//  Total_CT      += AllBound_CT_Visc;
//  Total_CQ      += AllBound_CQ_Visc;
//  Total_CMerit   = AllBound_CT_Visc / (AllBound_CQ_Visc + EPS);
//  Total_Heat     = AllBound_HF_Visc;
//  Total_MaxHeat  = AllBound_MaxHF_Visc;
//
//  /*--- Update the total coefficients per surface (note that all the nodes have the same value)---*/
//  for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
//    Surface_CL[iMarker_Monitoring]   += Surface_CL_Visc[iMarker_Monitoring];
//    Surface_CD[iMarker_Monitoring]   += Surface_CD_Visc[iMarker_Monitoring];
//    Surface_CSF[iMarker_Monitoring]  += Surface_CSF_Visc[iMarker_Monitoring];
//    Surface_CEff[iMarker_Monitoring]  = Surface_CL[iMarker_Monitoring] / (Surface_CD[iMarker_Monitoring] + EPS);
//    Surface_CFx[iMarker_Monitoring]  += Surface_CFx_Visc[iMarker_Monitoring];
//    Surface_CFy[iMarker_Monitoring]  += Surface_CFy_Visc[iMarker_Monitoring];
//    Surface_CFz[iMarker_Monitoring]  += Surface_CFz_Visc[iMarker_Monitoring];
//    Surface_CMx[iMarker_Monitoring]  += Surface_CMx_Visc[iMarker_Monitoring];
//    Surface_CMy[iMarker_Monitoring]  += Surface_CMy_Visc[iMarker_Monitoring];
//    Surface_CMz[iMarker_Monitoring]  += Surface_CMz_Visc[iMarker_Monitoring];
//  }
}

void CNEMONSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *sour_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {

  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iVar;
  unsigned short T_INDEX, TVE_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double *Normal, Area;
  su2double **GradV;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 1.0;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Get the locations of the primitive variables ---*/
  T_INDEX    = nodes->GetTIndex();
  TVE_INDEX  = nodes->GetTveIndex();

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual on the boundary with the specified heat flux ---*/
      // Note: Contributions from qtr and qve are used for proportional control
      //       to drive the solution toward the specified heatflux more quickly.
      GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }
      ktr = nodes->GetThermalConductivity(iPoint);
      //      kve = nodes->GetThermalConductivity_ve();
      //			Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
      //                                   Wall_HeatFlux*Area;
      //      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
      //                                   Wall_HeatFlux*Area;
      //
      //			/*--- Apply viscous residual to the linear system ---*/
      //      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Apply the no-slip condition in a strong way ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }
      if (implicit) {
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
    }
  }
}

void CNEMONSSolver::BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                                 CSolver **solution_container,
                                                 CNumerics *conv_numerics,
                                                 CNumerics *sour_numerics,
                                                 CConfig *config,
                                                 unsigned short val_marker) {

  /*--- Call standard "HeatFlux" wall to apply no-slip & energy b.c.'s ---*/
  BC_HeatFlux_Wall(geometry, solution_container, conv_numerics,
                   sour_numerics, config, val_marker);

  //	/*--- Local variables ---*/
  //  bool implicit;
  //	unsigned short iDim, iSpecies, iVar;
  //  unsigned short RHOS_INDEX, RHO_INDEX, T_INDEX, TVE_INDEX;
  //	unsigned long iVertex, iPoint;
  //	double pcontrol;
  //  su2double rho, Ys, eves, hs;
  //	double *Normal, Area;
  //  su2double *Ds, *V, *dYdn, SdYdn;
  //  su2double **GradV, **GradY;
  //
  //  /*--- Assign booleans ---*/
  //	implicit = (config->GetKind_TimeIntScheme_NEMO() == EULER_IMPLICIT);
  //
  //  /*--- Set "Proportional control" coefficient ---*/
  //  pcontrol = 0.6;
  //
  //  /*--- Get the locations of the primitive variables ---*/
  //  RHOS_INDEX = nodes->GetRhosIndex();
  //  RHO_INDEX  = nodes->GetRhoIndex();
  //  T_INDEX    = nodes->GetTIndex();
  //  TVE_INDEX  = nodes->GetTveIndex();
  //
  //  /*--- Allocate arrays ---*/
  //  dYdn = new su2double[nSpecies];
  //  GradY = new su2double*[nSpecies];
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    GradY[iSpecies] = new su2double[nDim];
  //
  //	/*--- Loop over all of the vertices on this boundary marker ---*/
  //	for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
  //		iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
  //
  //		/*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
  //		if (geometry->nodes->GetDomain()) {
  //
  //			/*--- Compute dual-grid area and boundary normal ---*/
  //			Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
  //			Area = 0.0;
  //			for (iDim = 0; iDim < nDim; iDim++)
  //				Area += Normal[iDim]*Normal[iDim];
  //			Area = sqrt (Area);
  //
  //			/*--- Initialize the convective & viscous residuals to zero ---*/
  //			for (iVar = 0; iVar < nVar; iVar++)
  //				Res_Visc[iVar] = 0.0;
  //
  //      /*--- Get temperature gradient information ---*/
  //      V = nodes->GetPrimVar();
  //      GradV  = nodes->GetGradient_Primitive();
  //
  //      /*--- Rename for convenience ---*/
  //      rho = V[RHO_INDEX];
  //      Ds  = nodes->GetDiffusionCoeff();
  //
  //      /*--- Calculate normal derivative of mass fraction ---*/
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //        Ys = V[RHOS_INDEX+iSpecies]/rho;
  //        dYdn[iSpecies] = 0.0;
  //        for (iDim = 0; iDim < nDim; iDim++)
  //          dYdn[iSpecies] += 1.0/rho * (GradV[RHOS_INDEX+iSpecies][iDim] -
  //                                       Ys*GradV[RHO_INDEX][iDim])*Normal[iDim];
  //      }
  //
  //      /*--- Calculate supplementary quantities ---*/
  //      SdYdn = 0.0;
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //        SdYdn += rho*Ds[iSpecies]*dYdn[iSpecies];
  //
  //      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //        Ys   = V[RHOS_INDEX+iSpecies]/rho;
  //        eves = nodes->CalcEve(config, V[TVE_INDEX], iSpecies);
  //        hs   = nodes->CalcHs(config, V[T_INDEX], eves, iSpecies);
  //        Res_Visc[iSpecies] = rho*Ds[iSpecies]*dYdn[iSpecies] - Ys*SdYdn;
  //        Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
  //        Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
  //      }
  //
  //			/*--- Viscous contribution to the residual at the wall ---*/
  //      LinSysRes.SubtractBlock(iPoint, Res_Visc);
  //		}
  //	}
  //
  //  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
  //    delete [] GradY[iSpecies];
  //  delete [] GradY;
  //  delete [] dYdn;
}

void CNEMONSSolver::BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                              CSolver **solution_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *sour_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {

  /*--- Local variables ---*/
  bool implicit, catalytic;
  unsigned short iDim, iSpecies, iVar;
  unsigned short T_INDEX, TVE_INDEX, RHOS_INDEX, RHO_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double rho, Ys, eves, hs;
  su2double *Normal, Area;
  su2double *Ds, *V, *dYdn, SdYdn;
  su2double **GradV, **GradY;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  catalytic = false;

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 0.6;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Get the locations of the primitive variables ---*/
  T_INDEX    = nodes->GetTIndex();
  TVE_INDEX  = nodes->GetTveIndex();
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Allocate arrays ---*/
  dYdn = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        Res_Sour[iVar] = 0.0;
      }

      /*--- Assign wall velocity to "Vector" array ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;

      /*--- Set the residual, truncation error, and velocity value ---*/
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Get temperature gradient information ---*/
      V = nodes->GetPrimitive(iPoint);
      GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }

      if (catalytic) {
        cout << "NEED TO IMPLEMENT CATALYTIC BOUNDARIES IN HEATFLUX!!!" << endl;
        exit(1);
      }
      else {

        /*--- Rename for convenience ---*/
        rho = V[RHO_INDEX];
        Ds  = nodes->GetDiffusionCoeff(iPoint);

        /*--- Calculate normal derivative of mass fraction ---*/
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          Ys = V[RHOS_INDEX+iSpecies]/rho;
          dYdn[iSpecies] = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            dYdn[iSpecies] += 1.0/rho * (GradV[RHOS_INDEX+iSpecies][iDim] -
                Ys*GradV[RHO_INDEX][iDim])*Normal[iDim];
        }

        /*--- Calculate supplementary quantities ---*/
        SdYdn = 0.0;
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
          SdYdn += rho*Ds[iSpecies]*dYdn[iSpecies];

        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          Ys   = V[RHOS_INDEX+iSpecies]/rho;
          //eves = nodes->CalcEve(config, V[TVE_INDEX], iSpecies);
          //hs   = nodes->CalcHs(config, V[T_INDEX], eves, iSpecies);
          //          Res_Visc[iSpecies] = -rho*Ds[iSpecies]*dYdn[iSpecies] + Ys*SdYdn;
          //          Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
          //          Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
        }
      }

      /*--- Get node thermal conductivity ---*/
      ktr = nodes->GetThermalConductivity(iPoint);
      kve = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Set the residual on the boundary with the specified heat flux ---*/
      // Note: Contributions from qtr and qve are used for proportional control
      //       to drive the solution toward the specified heatflux more quickly.
      Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
          Wall_HeatFlux*Area;
      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
          Wall_HeatFlux*Area;

      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      //      /*--- Apply the non-catalytic wall boundary ---*/
      //      // Note: We are re-calculating the chemistry residual and adding it to
      //      //       the linear system to eliminate the contribution from the solution
      //      //       (convention is to subtract sources)
      //      sour_numerics->SetConservative(nodes->GetSolution(),
      //                                     nodes->GetSolution() );
      //      sour_numerics->SetPrimitive   (nodes->GetPrimVar() ,
      //                                     nodes->GetPrimVar()  );
      //      sour_numerics->SetdPdU        (nodes->GetdPdU()    ,
      //                                     nodes->GetdPdU()     );
      //      sour_numerics->SetdTdU        (nodes->GetdTdU()    ,
      //                                     nodes->GetdTdU()     );
      //      sour_numerics->SetdTvedU      (nodes->GetdTvedU()  ,
      //                                     nodes->GetdTvedU()   );
      //      sour_numerics->SetVolume      (geometry->nodes->GetVolume());
      //      sour_numerics->ComputeChemistry(Res_Sour, Jacobian_i, config);
      //      LinSysRes.AddBlock(iPoint, Res_Sour);
      //      if (implicit)
      //        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      /*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal)/
       Note that we need to add a contribution for moving walls to the Jacobian. ---*/
      if (implicit) {
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for (iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  delete [] dYdn;
}

void CNEMONSSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                       CSolver **solution_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *sour_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {

  unsigned short iDim, iVar, jVar;
  unsigned short RHOS_INDEX, T_INDEX, TVE_INDEX, RHOCVTR_INDEX, RHOCVVE_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve, Ti, Tvei, Tj, Tvej, *dTdU, *dTvedU, Twall, dij, theta,
  Area, *Normal, UnitNormal[3], *Coord_i, *Coord_j, C;

  bool implicit   = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool ionization = config->GetIonization();

  if (ionization) {
    cout << "BC_ISOTHERMAL: NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Define 'proportional control' constant ---*/
  C = 5;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag);

  /*--- Extract necessary indices ---*/
  RHOS_INDEX    = nodes->GetRhosIndex();
  T_INDEX       = nodes->GetTIndex();
  TVE_INDEX     = nodes->GetTveIndex();
  RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();
  RHOCVVE_INDEX = nodes->GetRhoCvveIndex();

  /*--- Loop over boundary points to calculate energy flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);

      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = sqrt(dij);

      /*--- Calculate geometrical parameters ---*/
      theta = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        theta += UnitNormal[iDim]*UnitNormal[iDim];
      }

      /*--- Initialize viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar ++)
        Res_Visc[iVar] = 0.0;

      /*--- Store the corrected velocity at the wall which will
       be zero (v = 0), unless there is grid motion (v = u_wall)---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Calculate the gradient of temperature ---*/
      Ti   = nodes->GetTemperature(iPoint);
      Tj   = nodes->GetTemperature(jPoint);
      Tvei = nodes->GetTemperature_ve(iPoint);
      Tvej = nodes->GetTemperature_ve(jPoint);

      /*--- Rename variables for convenience ---*/
      ktr     = nodes->GetThermalConductivity(iPoint);
      kve     = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Twall-Ti) + kve*(Twall-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Twall-Tvei) *C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;

        dTdU   = nodes->GetdTdU(iPoint);
        dTvedU = nodes->GetdTvedU(iPoint);
        for (iVar = 0; iVar < nVar; iVar++) {
          Jacobian_i[nSpecies+nDim][iVar]   = -(ktr*theta/dij*dTdU[iVar] +
                                                kve*theta/dij*dTvedU[iVar])*Area;
          Jacobian_i[nSpecies+nDim+1][iVar] = - kve*theta/dij*dTvedU[iVar]*Area;
        }
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      } // implicit
    }
  } 
}

void CNEMONSSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                                   CSolver **solution_container,
                                                   CNumerics *conv_numerics,
                                                   CNumerics *sour_numerics,
                                                   CConfig *config,
                                                   unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_Isothermal_Wall(geometry, solution_container, conv_numerics,
                     sour_numerics, config, val_marker);

}

void CNEMONSSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                                CSolver **solution_container,
                                                CNumerics *conv_numerics,
                                                CNumerics *sour_numerics,
                                                CConfig *config,
                                                unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_Isothermal_Wall(geometry, solution_container, conv_numerics,
                     sour_numerics, config, val_marker);

  ///////////// FINITE DIFFERENCE METHOD ///////////////
  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar, RHOS_INDEX, RHO_INDEX, T_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double pcontrol, rho, *eves, *dTdU, *dTvedU, *Cvtr, *Cvve, *Normal, Area, Ru, RuSI,
  dij, UnitNormal[3], *Di, *Dj, *Vi, *Vj, *Yj, *dYdn, SdYdn, **GradY, **dVdU;
  const su2double *Yst;
  vector<su2double> hs, Cvtrs, Ms;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 0.6;

  /*--- Get species mass fractions at the wall ---*/
  Yst = config->GetWall_Catalycity();

  /*--- Get universal information ---*/
  RuSI = UNIVERSAL_GAS_CONSTANT;
  Ru   = 1000.0*RuSI;
  Ms   = FluidModel->GetMolarMass();
  
  /*--- Get the locations of the primitive variables ---*/
  RHOS_INDEX    = nodes->GetRhosIndex();
  RHO_INDEX     = nodes->GetRhoIndex();
  T_INDEX       = nodes ->GetTIndex();

  /*--- Allocate arrays ---*/
  Yj    = new su2double[nSpecies];
  dYdn  = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];
  dVdU = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    dVdU[iVar] = new su2double[nVar];

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dij += (geometry->nodes->GetCoord(jPoint, iDim) -
                geometry->nodes->GetCoord(iPoint, iDim))
            * (geometry->nodes->GetCoord(jPoint, iDim) -
               geometry->nodes->GetCoord(iPoint, iDim));
      }
      dij = sqrt(dij);


      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;


      /*--- Initialize the viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Res_Visc[iVar] = 0.0;

      /*--- Get primitive information ---*/
      Vi   = nodes->GetPrimitive(iPoint);
      Vj   = nodes->GetPrimitive(jPoint);
      Di   = nodes->GetDiffusionCoeff(iPoint);
      Dj   = nodes->GetDiffusionCoeff(jPoint);
      eves = nodes->GetEve(iPoint);
      hs   = FluidModel->GetSpeciesEnthalpy(Vi[T_INDEX], eves);
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)      
        Yj[iSpecies] = Vj[RHOS_INDEX+iSpecies]/Vj[RHO_INDEX];
      rho    = Vi[RHO_INDEX];
      dTdU   = nodes->GetdTdU(iPoint);
      dTvedU = nodes->GetdTvedU(iPoint);

      /*--- Calculate normal derivative of mass fraction ---*/
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dYdn[iSpecies] = (Yst[iSpecies]-Yj[iSpecies])/dij;

      /*--- Calculate supplementary quantities ---*/
      SdYdn = 0.0;
      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        SdYdn += rho*Di[iSpecies]*dYdn[iSpecies];

      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
        Res_Visc[iSpecies]         = -(-rho*Di[iSpecies]*dYdn[iSpecies]
                                       +Yst[iSpecies]*SdYdn            )*Area;
        Res_Visc[nSpecies+nDim]   += (Res_Visc[iSpecies]*hs[iSpecies]  )*Area;
        Res_Visc[nSpecies+nDim+1] += (Res_Visc[iSpecies]*eves[iSpecies])*Area;
      }

      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      if (implicit) {

        /*--- Initialize the transformation matrix ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++) {
            dVdU[iVar][jVar] = 0.0;
            Jacobian_j[iVar][jVar] = 0.0;
            Jacobian_i[iVar][jVar] = 0.0;
          }

        /*--- Populate transformation matrix ---*/
        // dYsdrk
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            dVdU[iSpecies][jSpecies] += -1.0/rho*Yst[iSpecies];
          dVdU[iSpecies][iSpecies] += 1.0/rho;
        }
        for (iVar = 0; iVar < nVar; iVar++) {
          dVdU[nSpecies+nDim][iVar]   = dTdU[iVar];
          dVdU[nSpecies+nDim+1][iVar] = dTvedU[iVar];
        }

        /*--- Calculate supplementary quantities ---*/
        Cvtrs = FluidModel->GetSpeciesCvTraRot();
        Cvve = nodes->GetCvve(iPoint);

        /*--- Take the primitive var. Jacobian & store in Jac. jj ---*/
        // Species mass fraction
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            Jacobian_j[iSpecies][jSpecies] += -Yst[iSpecies]*rho*Di[jSpecies]/dij;
          Jacobian_j[iSpecies][iSpecies] += rho*Di[iSpecies]/dij - SdYdn;
        }

        // Temperature
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
            Jacobian_j[nSpecies+nDim][iSpecies] += Jacobian_j[jSpecies][iSpecies]*hs[iSpecies];
          }
          Jacobian_j[nSpecies+nDim][nSpecies+nDim] += Res_Visc[iSpecies]/Area*(Ru/Ms[iSpecies] +
                                                                               Cvtrs[iSpecies]  );
          Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
        }

        // Vib.-El. Temperature
        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
          for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
            Jacobian_j[nSpecies+nDim+1][iSpecies] += Jacobian_j[jSpecies][iSpecies]*eves[iSpecies];
          Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
        }

        /*--- Multiply by the transformation matrix and store in Jac. ii ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (jVar = 0; jVar < nVar; jVar++)
            for (kVar = 0; kVar < nVar; kVar++)
              Jacobian_i[iVar][jVar] += Jacobian_j[iVar][kVar]*dVdU[kVar][jVar]*Area;

        /*--- Apply to the linear system ---*/
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  delete [] dYdn;
  delete [] Yj;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] dVdU[iVar];
  delete [] dVdU;
}

void CNEMONSSolver::BC_Smoluchowski_Maxwell(CGeometry *geometry,
                                CSolver **solution_container,
                                CNumerics *conv_numerics,
                                CNumerics *visc_numerics,
                                CConfig *config,
                                unsigned short val_marker) {


  unsigned short iDim, jDim, iVar, jVar, iSpecies;
  unsigned short T_INDEX, TVE_INDEX, VEL_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve;
  su2double Ti, Tvei, Tj, Tvej, *dTdU, *dTvedU;
  su2double Twall, Tslip, dij, theta;
  su2double Pi;
  su2double Area, *Normal, UnitNormal[3];
  su2double *Coord_i, *Coord_j;
  su2double C;

  su2double TMAC, TAC;
  su2double Viscosity, Lambda;
  su2double Density, GasConstant, Ru;

  su2double **Grad_PrimVar;
  su2double Vector_Tangent_dT[3], Vector_Tangent_dTve[3], Vector_Tangent_HF[3];
  su2double dTn, dTven;

  su2double TauElem[3], TauTangent[3];
  su2double Tau[3][3];
  su2double TauNormal;
  su2double div_vel=0, Delta;

  vector<su2double> Ms;
  
  bool implicit   = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool ionization = config->GetIonization();

  if (ionization) {
    cout << "BC_SMOLUCHOWSKI_MAXWELL: NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION" << endl;
    exit(1);
  }

  /*--- Define 'proportional control' constant ---*/
  C = 5;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature and accomodation coefficients---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag);
  TMAC  = 1.0;
  TAC   = 1.0;

  /*--- Extract necessary indices ---*/
  T_INDEX       = nodes->GetTIndex();
  VEL_INDEX     = nodes->GetVelIndex();
  TVE_INDEX     = nodes->GetTveIndex();

  /*--- Loop over boundary points to calculate energy flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = -Normal[iDim]/Area;

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);

      dij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dij += (Coord_j[iDim] - Coord_i[iDim])*(Coord_j[iDim] - Coord_i[iDim]);
      dij = sqrt(dij);

      /*--- Calculate geometrical parameters ---*/
      /*theta = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        theta += UnitNormal[iDim]*UnitNormal[iDim];
      }*/

      /*--- Calculate Pressure ---*/
      Pi   = nodes->GetPressure(iPoint);

      /*--- Calculate the gradient of temperature ---*/
      Ti   = nodes->GetTemperature(iPoint);
      Tj   = nodes->GetTemperature(jPoint);
      Tvei = nodes->GetTemperature_ve(iPoint);
      Tvej = nodes->GetTemperature_ve(jPoint);

      /*--- Rename variables for convenience ---*/
      ktr  = nodes->GetThermalConductivity(iPoint);
      kve  = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Retrieve Flow Data ---*/
      Viscosity = nodes->GetLaminarViscosity(iPoint);
      Density = nodes->GetDensity(iPoint);

      Ms = FluidModel->GetMolarMass();

      /*--- Calculate specific gas constant --- */
      GasConstant=0;
      for(iSpecies=0;iSpecies<nSpecies;iSpecies++)
        GasConstant+=UNIVERSAL_GAS_CONSTANT*1000.0/Ms[iSpecies]*nodes->GetMassFraction(iPoint,iSpecies);
      
      /*--- Calculate temperature gradients normal to surface---*/ //Doubt about minus sign
      dTn   = - (Ti-Tj)/dij;
      dTven = - (Tvei-Tvej)/dij;

      /*--- Calculate molecular mean free path ---*/
      Lambda = Viscosity/Density*sqrt(PI_NUMBER/(2*GasConstant*Ti));

      /*--- Calculate Temperature Slip ---*/
      Tslip = ((2-TAC)/TAC)*2*Gamma/(Gamma+1)/Prandtl_Lam*Lambda*dTn+Twall;

      /*--- Retrieve Primitive Gradients ---*/
      Grad_PrimVar = nodes->GetGradient_Primitive(iPoint);

      /*--- Calculate temperature gradients tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_Tangent_dT[iDim]   = Grad_PrimVar[T_INDEX][iDim] - dTn * UnitNormal[iDim];
        Vector_Tangent_dTve[iDim] = Grad_PrimVar[TVE_INDEX][iDim] - dTven * UnitNormal[iDim];
      }

      /*--- Calculate Heatflux tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++) 
        Vector_Tangent_HF[iDim] = ktr*Vector_Tangent_dT[iDim]+kve*Vector_Tangent_dTve[iDim];
      
      /*--- Initialize viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar ++)
        Res_Visc[iVar] = 0.0;

      div_vel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        div_vel += Grad_PrimVar[VEL_INDEX+iDim][iDim];

      for (iDim = 0; iDim < nDim; iDim++) {
        for (jDim = 0 ; jDim < nDim; jDim++) {
          Delta = 0.0; if (iDim == jDim) Delta = 1.0;
          Tau[iDim][jDim] = Viscosity*(Grad_PrimVar[VEL_INDEX+jDim][iDim] +
              Grad_PrimVar[VEL_INDEX+iDim][jDim]  )
              - TWO3*Viscosity*div_vel*Delta;
        }
        TauElem[iDim] = 0.0;
        for (jDim = 0; jDim < nDim; jDim++)
          TauElem[iDim] += Tau[iDim][jDim]*UnitNormal[jDim];
      }

      /*--- Compute wall shear stress (using the stress tensor) ---*/
      TauNormal = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        TauNormal += TauElem[iDim] * UnitNormal[iDim];
      for (iDim = 0; iDim < nDim; iDim++) {
        TauTangent[iDim] = TauElem[iDim] - TauNormal * UnitNormal[iDim];
      }

      /*--- Store the Slip Velocity at the wall */
      for (iDim = 0; iDim < nDim; iDim++)
        Vector[iDim] = - Lambda/Viscosity*(2-TMAC)/TMAC*(TauTangent[iDim])-3/4*(Gamma-1)/Gamma*Prandtl_Lam/Pi*Vector_Tangent_HF[iDim];

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes.SetBlock_Zero(iPoint, nSpecies+iDim);
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Tslip-Ti) + kve*(Tslip-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Tslip-Tvei)*C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }
  } 
}