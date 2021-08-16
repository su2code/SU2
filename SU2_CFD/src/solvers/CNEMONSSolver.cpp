/*!
 * \file CNEMONSSolver.cpp
 * \brief Headers of the CNEMONSSolver class
 * \author S. R. Copeland, F. Palacios, W. Maier.
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/solvers/CNEMONSSolver.hpp"
#include "../../include/variables/CNEMONSVariable.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CFVMFlowSolverBase.inl"

/*--- Explicit instantiation of the parent class of CNEMOEulerSolver,
 *    to spread the compilation over two cpp files. ---*/
template class CFVMFlowSolverBase<CNEMOEulerVariable, ENUM_REGIME::COMPRESSIBLE>;

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

}

void CNEMONSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                              unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  const unsigned long InnerIter = config->GetInnerIter();
  const bool limiter    = (config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = config->GetKind_SlopeLimit_Flow() == VAN_ALBADA_EDGE;
  const bool muscl      = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center     = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;

  /*--- Common preprocessing steps (implemented by CNEMOEulerSolver) ---*/

  CommonPreprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  /*--- Compute gradient for MUSCL reconstruction. ---*/

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

  /*--- Evaluate the vorticity and strain rate magnitude ---*/
  //TODO: THIS NEEDS TO BE UPDATED ASAP!
  StrainMag_Max = 0.0;
  Omega_Max = 0.0;
  //nodes->SetVorticity_StrainMag();

  //su2double strainMax = 0.0;
  su2double omegaMax = 0.0;

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    //su2double StrainMag = nodes->GetStrainMag(iPoint);
    const su2double* Vorticity = nodes->GetVorticity(iPoint);
    su2double Omega = sqrt(Vorticity[0]*Vorticity[0]+ Vorticity[1]*Vorticity[1]+ Vorticity[2]*Vorticity[2]);

    //strainMax = max(strainMax, StrainMag);
    omegaMax = max(omegaMax, Omega);

  }

  //StrainMag_Max = max(StrainMag_Max, strainMax);
  Omega_Max = max(Omega_Max, omegaMax);

  if ((iMesh == MESH_0) && (config->GetComm_Level() == COMM_FULL)) {

    su2double MyOmega_Max = Omega_Max;
    //su2double MyStrainMag_Max = StrainMag_Max;
    //SU2_MPI::Allreduce(&MyStrainMag_Max, &StrainMag_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());
    SU2_MPI::Allreduce(&MyOmega_Max, &Omega_Max, 1, MPI_DOUBLE, MPI_MAX, SU2_MPI::GetComm());

  }
}

void CNEMONSSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, const CConfig *config, bool reconstruction) {

  unsigned long iPoint, iVar;
  unsigned short iSpecies, RHO_INDEX, RHOS_INDEX;
  auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
  const auto comm = reconstruction? PRIMITIVE_GRAD_REC : PRIMITIVE_GRADIENT;
  const auto commPer = reconstruction? PERIODIC_PRIM_GG_R : PERIODIC_PRIM_GG;

  /*--- Get indices of species & mixture density ---*/
  RHOS_INDEX = nodes->GetRhosIndex();
  RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Modify species density to mass concentration ---*/
  for ( iPoint = 0; iPoint < nPoint; iPoint++){
    su2double primitives_aux[MAXNVAR] = {0.0};
    for( iVar = 0; iVar < nPrimVar; iVar++)
      primitives_aux[iVar] = nodes->GetPrimitive(iPoint, iVar);
    for ( iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      primitives_aux[RHOS_INDEX+iSpecies] = primitives_aux[RHOS_INDEX+iSpecies]/primitives_aux[RHO_INDEX];
    for( iVar = 0; iVar < nPrimVar; iVar++)
      nodes->SetPrimitive_Aux(iPoint, iVar, primitives_aux[iVar] );
  }

  const auto& primitives = nodes->GetPrimitive_Aux();

  computeGradientsGreenGauss(this, comm, commPer, *geometry, *config, primitives, 0, nPrimVarGrad, gradient);
}

unsigned long CNEMONSSolver::SetPrimitive_Variables(CSolver **solver_container,CConfig *config, bool Output) {

  unsigned long iPoint, nonPhysicalPoints = 0;
  const unsigned short turb_model = config->GetKind_Turb_Model();
  //const bool tkeNeeded = (turb_model == SST) || (turb_model == SST_SUST);

  bool nonphysical = true;

  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Retrieve the value of the kinetic energy (if needed). ---*/

    su2double eddy_visc = 0.0;
    //su2double turb_ke = 0.0;

    if (turb_model != NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      //if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

       nodes->SetEddyViscosity(iPoint, eddy_visc);
    }

    /*--- Incompressible flow, primitive variables ---*/

    nonphysical = nodes->SetPrimVar(iPoint,FluidModel);

    /* Check for non-realizable states for reporting. */

    if (nonphysical) nonPhysicalPoints++;

  }

  return nonPhysicalPoints;
}

void CNEMONSSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solution_container,
                                     CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep) {

  bool err;
  unsigned short iVar;
  unsigned long iPoint, jPoint, iEdge;

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

    /*--- Eddy viscosity ---*/
    numerics->SetEddyViscosity(nodes->GetEddyViscosity(iPoint),
                               nodes->GetEddyViscosity(jPoint) );

    /*--- Thermal conductivity ---*/
    numerics->SetThermalConductivity(nodes->GetThermalConductivity(iPoint),
                                     nodes->GetThermalConductivity(jPoint));

    /*--- Vib-el. thermal conductivity ---*/
    numerics->SetThermalConductivity_ve(nodes->GetThermalConductivity_ve(iPoint),
                                        nodes->GetThermalConductivity_ve(jPoint) );

    /*--- Compute and update residual ---*/
    auto residual = numerics->ComputeResidual(config);

    /*--- Check for NaNs before applying the residual to the linear system ---*/
    err = false;
    for (iVar = 0; iVar < nVar; iVar++)
      if (residual[iVar] != residual[iVar]) err = true;
    //if (implicit)
    //  for (iVar = 0; iVar < nVar; iVar++)
    //    for (jVar = 0; jVar < nVar; jVar++)
    //      if ((Jacobian_i[iVar][jVar] != Jacobian_i[iVar][jVar]) ||
    //          (Jacobian_j[iVar][jVar] != Jacobian_j[iVar][jVar])   )
    //        err = true;

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, residual);
      LinSysRes.AddBlock(jPoint, residual);
      //if (implicit) {
      //  Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      //  Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
      //  Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
      //  Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
      //}
    }
  } //iEdge
}

void CNEMONSSolver::BC_HeatFluxNonCatalytic_Wall(CGeometry *geometry,
                                                 CSolver **solver_container,
                                                 CNumerics *conv_numerics,
                                                 CNumerics *sour_numerics,
                                                 CConfig *config,
                                                 unsigned short val_marker) {

  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iVar;
  unsigned short T_INDEX, TVE_INDEX, RHOCVTR_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double dTdn, dTvedn, ktr, kve, pcontrol, Wall_HeatFlux;
  su2double *Normal, Area, *V;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  pcontrol = 1.0;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag);

  /*--- Get the locations of the primitive variables ---*/
  T_INDEX       = nodes->GetTIndex();
  TVE_INDEX     = nodes->GetTveIndex();
  RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
      }

      /*--- Set the residual on the boundary with the specified heat flux ---*/
      // Note: Contributions from qtr and qve are used for proportional control
      //       to drive the solution toward the specified heatflux more quickly.
      V      = nodes->GetPrimitive(iPoint);
      const auto GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }
      ktr = nodes->GetThermalConductivity(iPoint);
      kve = nodes->GetThermalConductivity_ve(iPoint);

      /*--- Scale thermal conductivity with turb ---*/
      // TODO: Need to determine proper way to incorporate eddy viscosity
      // This is only scaling Kve by same factor as ktr
      su2double Mass = 0.0;
      auto&     Ms   = FluidModel->GetSpeciesMolarMass();
      su2double tmp1, scl, Cptr;
      su2double Ru=1000.0*UNIVERSAL_GAS_CONSTANT;
      su2double eddy_viscosity = nodes->GetEddyViscosity(iPoint);
      for (unsigned short iSpecies=0; iSpecies<nSpecies; iSpecies++)
        Mass += V[iSpecies]*Ms[iSpecies];
      Cptr = V[RHOCVTR_INDEX]+Ru/Mass;
      tmp1 = Cptr*(eddy_viscosity/Prandtl_Turb);
      scl  = tmp1/ktr;
      ktr += Cptr*(eddy_viscosity/Prandtl_Turb);
      kve  = kve*(1.0+scl);
      //Cpve = V[RHOCVVE_INDEX]+Ru/Mass;
      //kve += Cpve*(val_eddy_viscosity/Prandtl_Turb);

      /*--- Compute residual ---*/
      Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
                                      Wall_HeatFlux*Area;
      Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
                                      Wall_HeatFlux*Area;
      /*--- Apply viscous residual to the linear system ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Apply the no-slip condition in a strong way ---*/
      for (iDim = 0; iDim < nDim; iDim++) Vector[iDim] = 0.0;
      nodes->SetVelocity_Old(iPoint,Vector);
      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes(iPoint, nSpecies+iDim) = 0.0;
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

void CNEMONSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *sour_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {

  string Marker_Tag, Catalytic_Tag;
  unsigned short iMarker_Catalytic;
  bool catalytic;

  catalytic = false;
  iMarker_Catalytic = 0;

  Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  while( iMarker_Catalytic < config->GetnWall_Catalytic()){

    Catalytic_Tag = config->GetWall_Catalytic_TagBound(iMarker_Catalytic);

    if (Marker_Tag == Catalytic_Tag){

      catalytic = true;
      BC_HeatFluxCatalytic_Wall(geometry, solver_container, conv_numerics,
                                sour_numerics, config, val_marker);
      break;

    } else {

      iMarker_Catalytic++;

    }
  }

  if(!catalytic) BC_HeatFluxNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                              sour_numerics, config, val_marker);


}

void CNEMONSSolver::BC_HeatFluxCatalytic_Wall(CGeometry *geometry,
                                              CSolver **solver_container,
                                              CNumerics *conv_numerics,
                                              CNumerics *sour_numerics,
                                              CConfig *config,
                                              unsigned short val_marker) {

  SU2_MPI::Error("BC_HEATFLUX with catalytic wall: Not operational in NEMO.", CURRENT_FUNCTION);
  //TODO: SCALE WITH EDDY VISC
  /*--- Local variables ---*/
  bool implicit, catalytic;
  unsigned short iDim, iSpecies, iVar;
  unsigned short T_INDEX, TVE_INDEX, RHOS_INDEX, RHO_INDEX;
  unsigned long iVertex, iPoint, total_index;
  su2double Wall_HeatFlux, dTdn, dTvedn, ktr, kve, pcontrol;
  su2double rho, Ys;
  su2double *Normal, Area;
  su2double *Ds, *V, *dYdn, SdYdn;
  su2double **GradY;

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
      Area = GeometryToolbox::Norm(nDim, Normal);

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
        LinSysRes(iPoint, nSpecies+iDim) = 0.0;
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Get temperature gradient information ---*/
      V = nodes->GetPrimitive(iPoint);
      const auto GradV  = nodes->GetGradient_Primitive(iPoint);
      dTdn   = 0.0;
      dTvedn = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
        dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
      }

      if (catalytic) {
        SU2_MPI::Error("NEED TO IMPLEMENT CATALYTIC BOUNDARIES IN HEATFLUX!!!",CURRENT_FUNCTION);
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
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *sour_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {

  string Marker_Tag, Catalytic_Tag;
  unsigned short iMarker_Catalytic;
  bool catalytic;

  catalytic = false;
  iMarker_Catalytic = 0;

  Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  while( iMarker_Catalytic < config->GetnWall_Catalytic()){

    Catalytic_Tag = config->GetWall_Catalytic_TagBound(iMarker_Catalytic);

    if (Marker_Tag == Catalytic_Tag){

      catalytic = true;
      BC_IsothermalCatalytic_Wall(geometry, solver_container, conv_numerics,
                                  sour_numerics, config, val_marker);
      break;
    } else {
      iMarker_Catalytic++;
    }
  }

  if(!catalytic) BC_IsothermalNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                                sour_numerics, config, val_marker);
}

void CNEMONSSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                                   CSolver **solver_container,
                                                   CNumerics *conv_numerics,
                                                   CNumerics *sour_numerics,
                                                   CConfig *config,
                                                   unsigned short val_marker) {

  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve, Ti, Tvei, Tj, Tvej, Twall, dij, theta,
  Area, *Normal, UnitNormal[3], *Coord_i, *Coord_j, C;
  su2double *V;
  bool ionization = config->GetIonization();

  if (ionization) {
    SU2_MPI::Error("NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION",CURRENT_FUNCTION);
  }

  /*--- Extract required indices ---*/
  unsigned short RHOCVTR_INDEX = nodes->GetRhoCvtrIndex();

  /*--- Define 'proportional control' constant ---*/
  C = 5;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag);

  /*--- Loop over boundary points to calculate energy flux ---*/
  for(iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);
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
        LinSysRes(iPoint, nSpecies+iDim) = 0.0;
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

      /*--- Scale thermal conductivity with turb ---*/
      // TODO: Need to determine proper way to incorporate eddy viscosity
      // This is only scaling Kve by same factor as ktr
      V = nodes->GetPrimitive(iPoint);
      su2double Mass = 0.0;
      auto&     Ms   = FluidModel->GetSpeciesMolarMass();
      su2double tmp1, scl, Cptr;
      su2double Ru=1000.0*UNIVERSAL_GAS_CONSTANT;
      su2double eddy_viscosity=nodes->GetEddyViscosity(iPoint);
      for (unsigned short iSpecies=0; iSpecies<nSpecies; iSpecies++)
        Mass += V[iSpecies]*Ms[iSpecies];
      Cptr = V[RHOCVTR_INDEX]+Ru/Mass;
      tmp1 = Cptr*(eddy_viscosity/Prandtl_Turb);
      scl  = tmp1/ktr;
      ktr += Cptr*(eddy_viscosity/Prandtl_Turb);
      kve  = kve*(1.0+scl);
      //Cpve = V[RHOCVVE_INDEX]+Ru/Mass;
      //kve += Cpve*(val_eddy_viscosity/Prandtl_Turb);

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Twall-Ti) + kve*(Twall-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Twall-Tvei) *C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      //if (implicit) {
      //  for (iVar = 0; iVar < nVar; iVar++)
      //    for (jVar = 0; jVar < nVar; jVar++)
      //      Jacobian_i[iVar][jVar] = 0.0;
      //
      //  dTdU   = nodes->GetdTdU(iPoint);
      //  dTvedU = nodes->GetdTvedU(iPoint);
      //  for (iVar = 0; iVar < nVar; iVar++) {
      //    Jacobian_i[nSpecies+nDim][iVar]   = -(ktr*theta/dij*dTdU[iVar] +
      //                                          kve*theta/dij*dTvedU[iVar])*Area;
      //    Jacobian_i[nSpecies+nDim+1][iVar] = - kve*theta/dij*dTvedU[iVar]*Area;
      //  }
      //  Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      //} // implicit
    }
  }
}

void CNEMONSSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                                CSolver **solver_container,
                                                CNumerics *conv_numerics,
                                                CNumerics *sour_numerics,
                                                CConfig *config,
                                                unsigned short val_marker) {

  SU2_MPI::Error("BC_ISOTHERMAL with catalytic wall: Not operational in NEMO.", CURRENT_FUNCTION);

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_IsothermalNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                 sour_numerics, config, val_marker);

  ///////////// FINITE DIFFERENCE METHOD ///////////////
  /*--- Local variables ---*/
  bool implicit;
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar,
      RHOS_INDEX, RHO_INDEX, T_INDEX, TVE_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double rho, *eves, *dTdU, *dTvedU, *Cvve, *Normal, Area, Ru, RuSI,
  dij, *Di, *Vi, *Vj, *Yj, *dYdn, SdYdn, **GradY, **dVdU;
  const su2double *Yst;
  vector<su2double> hs, Cvtrs;

  /*--- Assign booleans ---*/
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  //su2double pcontrol = 0.6;

  /*--- Get species mass fractions at the wall ---*/
  Yst = config->GetWall_Catalycity();

  /*--- Get universal information ---*/
  RuSI     = UNIVERSAL_GAS_CONSTANT;
  Ru       = 1000.0*RuSI;
  auto& Ms = FluidModel->GetSpeciesMolarMass();

  /*--- Get the locations of the primitive variables ---*/
  RHOS_INDEX  = nodes->GetRhosIndex();
  RHO_INDEX   = nodes->GetRhoIndex();
  T_INDEX     = nodes->GetTIndex();
  TVE_INDEX   = nodes->GetTveIndex();

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
      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Initialize the viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Res_Visc[iVar] = 0.0;

      /*--- Get primitive information ---*/
      Vi   = nodes->GetPrimitive(iPoint);
      Vj   = nodes->GetPrimitive(jPoint);
      Di   = nodes->GetDiffusionCoeff(iPoint);
      eves = nodes->GetEve(iPoint);
      hs   = FluidModel->ComputeSpeciesEnthalpy(Vi[T_INDEX], Vi[TVE_INDEX], eves);
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
                                            CSolver **solver_container,
                                            CNumerics *conv_numerics,
                                            CNumerics *visc_numerics,
                                            CConfig *config,
                                            unsigned short val_marker) {

  unsigned short iDim, iVar, iSpecies;
  unsigned short T_INDEX, TVE_INDEX, VEL_INDEX;
  unsigned long iVertex, iPoint, jPoint;
  su2double ktr, kve, Mass = 0.0;
  su2double Ti, Tvei, Tj, Tvej;
  su2double Twall, Tslip, Tslip_ve, dij;
  su2double Pi;
  su2double Area, UnitNormal[MAXNDIM];
  su2double C, alpha_V, alpha_T;

  su2double TMAC, TAC;
  su2double Viscosity, Eddy_Visc, Lambda;
  su2double Density, GasConstant;

  su2double Vector_Tangent_dT[MAXNDIM] = {0.0}, Vector_Tangent_dTve[MAXNDIM] = {0.0}, Vector_Tangent_HF[MAXNDIM] = {0.0};
  su2double dTn, dTven;
  su2double rhoCvtr, rhoCvve;

  su2double TauTangent[MAXNDIM] = {0.0};
  su2double Tau[MAXNDIM][MAXNDIM] = {{0.0}};

  bool ionization = config->GetIonization();

  if (ionization) {
    SU2_MPI::Error("NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION", CURRENT_FUNCTION);
  }

  /*--- Define 'proportional control' constant ---*/
  C = 1;

  /*---Define under-relaxation factors --- */
  alpha_V = 0.1;
  alpha_T = 1.0;

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
      const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);

      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim] = Normal[iDim]/Area;

      /*--- Compute closest normal neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      const auto Coord_i = geometry->nodes->GetCoord(iPoint);
      const auto Coord_j = geometry->nodes->GetCoord(jPoint);

      dij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

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

      /*--- Retrieve Cv*density ---*/
      rhoCvtr = nodes->GetRhoCv_tr(iPoint);
      rhoCvve = nodes->GetRhoCv_ve(iPoint);

      /*--- Retrieve Flow Data ---*/
      Viscosity = nodes->GetLaminarViscosity(iPoint);
      Eddy_Visc = nodes->GetEddyViscosity(iPoint);
      Density   = nodes->GetDensity(iPoint);
      Gamma     = nodes->GetGamma(iPoint);

      /*--- Incorporate turbulence effects ---*/
      const auto& Ms = FluidModel->GetSpeciesMolarMass();
      su2double  Ru = 1000.0*UNIVERSAL_GAS_CONSTANT;
      su2double  tmp1, scl, Cptr;
      su2double *Vi = nodes->GetPrimitive(iPoint);

      for (iSpecies=0; iSpecies<nSpecies; iSpecies++)
        Mass += Vi[iSpecies]*Ms[iSpecies];
      Cptr = rhoCvtr + Ru/Mass;
      tmp1 = Cptr*(Eddy_Visc/Prandtl_Turb);
      scl  = tmp1/ktr;
      ktr += Cptr*(Eddy_Visc/Prandtl_Turb);
      kve  = kve*(1.0+scl);

      /*--- Retrieve Primitive Gradients ---*/
      const auto Grad_PrimVar = nodes->GetGradient_Primitive(iPoint);

      /*--- Calculate specific gas constant --- */
      //TODO: Move to fluidmodel?
      GasConstant = 0.0;
      for(iSpecies = 0; iSpecies<nSpecies; iSpecies++)
        GasConstant+=UNIVERSAL_GAS_CONSTANT*1000.0/Ms[iSpecies]*nodes->GetMassFraction(iPoint,iSpecies);

      /*--- Calculate temperature gradients normal to surface---*/ //Doubt about minus sign
      dTn = GeometryToolbox::DotProduct(nDim, Grad_PrimVar[T_INDEX], UnitNormal);
      dTven = GeometryToolbox::DotProduct(nDim, Grad_PrimVar[TVE_INDEX], UnitNormal);

      /*--- Calculate molecular mean free path ---*/
      Lambda = Viscosity/Density*sqrt(PI_NUMBER/(2.0*GasConstant*Ti));

      /*--- Calculate Temperature Slip ---*/
      Tslip    = ((2.0-TAC)/TAC)*2.0*Gamma/(Gamma+1.0)/Prandtl_Lam*Lambda*dTn+Twall;

      if (dTven==0) Tslip_ve = Twall;
      else Tslip_ve = (Tslip-Twall)*(kve*rhoCvtr/dTn)/(ktr*rhoCvve/dTven)+Twall;

      /*--- Calculate temperature gradients tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_Tangent_dT[iDim]   = Grad_PrimVar[T_INDEX][iDim] - dTn * UnitNormal[iDim];
        Vector_Tangent_dTve[iDim] = Grad_PrimVar[TVE_INDEX][iDim] - dTven * UnitNormal[iDim];
      }

      /*--- Calculate Heatflux tangent to surface ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Vector_Tangent_HF[iDim] = -ktr*Vector_Tangent_dT[iDim]-kve*Vector_Tangent_dTve[iDim];

      /*--- Initialize viscous residual to zero ---*/
      for (iVar = 0; iVar < nVar; iVar ++)
        Res_Visc[iVar] = 0.0;

      CNumerics::ComputeStressTensor(nDim, Tau, Grad_PrimVar+VEL_INDEX, Viscosity);

      GeometryToolbox::TangentProjection(nDim, Tau, UnitNormal, TauTangent);

      /*--- Store the Slip Velocity at the wall */
      for (iDim = 0; iDim < nDim; iDim++)
        Vector[iDim] = Lambda/Viscosity*(2.0-TMAC)/TMAC*(TauTangent[iDim])-3.0/4.0*(Gamma-1.0)/Gamma*Prandtl_Lam/Pi*Vector_Tangent_HF[iDim];

      /*--- Apply under-relaxation ---*/
      Tslip    = (1.0-alpha_T)*nodes->GetTemperature(iPoint)+(alpha_T)*Tslip;
      Tslip_ve = (1.0-alpha_T)*nodes->GetTemperature_ve(iPoint)+(alpha_T)*Tslip_ve;

      for (iDim = 0; iDim < nDim; iDim++){
        Vector[iDim] = (1.0-alpha_V)*nodes->GetVelocity(iPoint,iDim)+(alpha_V)*Vector[iDim];
      }

      nodes->SetVelocity_Old(iPoint,Vector);

      for (iDim = 0; iDim < nDim; iDim++) {
        LinSysRes(iPoint, nSpecies+iDim) = 0.0;
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Apply to the linear system ---*/
      Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                   (ktr*(Tslip-Ti) + kve*(Tslip_ve-Tvei))*C)*Area/dij;
      Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Tslip_ve-Tvei)*C)*Area/dij;

      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }
  }
}
