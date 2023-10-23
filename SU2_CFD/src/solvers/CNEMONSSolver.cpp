/*!
 * \file CNEMONSSolver.cpp
 * \brief Headers of the CNEMONSSolver class
 * \author S. R. Copeland, F. Palacios, W. Maier.
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

  const auto InnerIter = config->GetInnerIter();
  const bool limiter = (config->GetKind_SlopeLimit_Flow() != LIMITER::NONE) && (InnerIter <= config->GetLimiterIter());
  const bool van_albada = config->GetKind_SlopeLimit_Flow() == LIMITER::VAN_ALBADA_EDGE;
  const bool muscl = config->GetMUSCL_Flow() && (iMesh == MESH_0);
  const bool center = config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED;
  const bool wall_functions = config->GetWall_Functions();

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

  /*--- Compute vorticity and strain mag. ---*/

  ComputeVorticityAndStrainMag(*config, geometry, iMesh);

  /*--- Compute the TauWall from the wall functions ---*/

  if (wall_functions) {
    SetTau_Wall_WF(geometry, solver_container, config);
  }
}

unsigned long CNEMONSSolver::SetPrimitive_Variables(CSolver **solver_container,CConfig *config, bool Output) {

  /*--- Number of non-physical points, local to the thread, needs
   *    further reduction if function is called in parallel ---*/
  unsigned long nonPhysicalPoints = 0;

  const TURB_MODEL turb_model = config->GetKind_Turb_Model();
  //const bool tkeNeeded = (turb_model == TURB_MODEL::SST);

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint ++) {

    /*--- Retrieve the value of the kinetic energy (if needed). ---*/

    su2double eddy_visc = 0.0; //su2double turb_ke = 0.0;

    if (turb_model != TURB_MODEL::NONE && solver_container[TURB_SOL] != nullptr) {
      eddy_visc = solver_container[TURB_SOL]->GetNodes()->GetmuT(iPoint);
      //if (tkeNeeded) turb_ke = solver_container[TURB_SOL]->GetNodes()->GetSolution(iPoint,0);

      nodes->SetEddyViscosity(iPoint, eddy_visc);
    }

    /*--- Compressible flow, primitive variables. ---*/

    bool nonphysical = nodes->SetPrimVar(iPoint,FluidModel);

    /* Check for non-realizable states for reporting. */

    nonPhysicalPoints += nonphysical;

  }
  END_SU2_OMP_FOR
  return nonPhysicalPoints;
}

void CNEMONSSolver::SetPrimitive_Gradient_GG(CGeometry *geometry, const CConfig *config, bool reconstruction) {

  auto& gradient = reconstruction ? nodes->GetGradient_Reconstruction() : nodes->GetGradient_Primitive();
  const auto comm = reconstruction? PRIMITIVE_GRAD_REC : PRIMITIVE_GRADIENT;
  const auto commPer = reconstruction? PERIODIC_PRIM_GG_R : PERIODIC_PRIM_GG;

  /*--- Get indices of species & mixture density ---*/
  const unsigned short RHOS_INDEX = nodes->GetRhosIndex();
  const unsigned short RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Modify species density to mass concentration ---*/
  for (auto iPoint = 0ul; iPoint < nPoint; iPoint++) {
    su2double primitives_aux[MAXNVAR] = {0.0};
    for (auto iVar = 0ul; iVar < nPrimVar; iVar++)
      primitives_aux[iVar] = nodes->GetPrimitive(iPoint, iVar);
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      primitives_aux[RHOS_INDEX+iSpecies] = primitives_aux[RHOS_INDEX+iSpecies]/primitives_aux[RHO_INDEX];
    for (auto iVar = 0ul; iVar < nPrimVar; iVar++)
      nodes->SetPrimitive_Aux(iPoint, iVar, primitives_aux[iVar] );
  }

  const auto& primitives = nodes->GetPrimitive_Aux();

  computeGradientsGreenGauss(this, comm, commPer, *geometry, *config, primitives, 0, nPrimVarGrad, gradient);
}

void CNEMONSSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics **numerics_container,
                                     CConfig *config, unsigned short iMesh,
                                     unsigned short iRKStep) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  CNumerics* numerics = numerics_container[VISC_TERM];

  for (auto iEdge = 0ul; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points, coordinates and normal vector in edge ---*/
    const auto iPoint = geometry->edges->GetNode(iEdge, 0);
    const auto jPoint = geometry->edges->GetNode(iEdge, 1);
    numerics->SetCoord(geometry->nodes->GetCoord(iPoint),
                       geometry->nodes->GetCoord(jPoint) );
    numerics->SetNormal(geometry->edges->GetNormal(iEdge));

    /*--- Primitive variables, and gradient ---*/
    numerics->SetConservative   (nodes->GetSolution(iPoint),
                                 nodes->GetSolution(jPoint) );
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
    bool err = CNumerics::CheckResidualNaNs(implicit, nVar, residual);

    /*--- Update the residual and Jacobian ---*/
    if (!err) {
      LinSysRes.SubtractBlock(iPoint, residual);
      LinSysRes.AddBlock(jPoint, residual);
      if (implicit) {
        Jacobian.UpdateBlocks(iEdge, iPoint, jPoint, residual.jacobian_i, residual.jacobian_j);
      }
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
  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
  if (config->GetIntegrated_HeatFlux()) {
    Wall_HeatFlux /= geometry->GetSurfaceArea(config, val_marker);
  }
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);

  /*--- Set "Proportional control" coefficient ---*/
  const su2double pcontrol = 1.0;

  /*--- Get the locations of the primitive variables ---*/
  const unsigned short T_INDEX = nodes->GetTIndex();
  const unsigned short TVE_INDEX = nodes->GetTveIndex();

  /*--- Loop over all of the vertices on this boundary marker ---*/
  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Compute dual-grid area and boundary normal ---*/
    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
    const su2double Area = GeometryToolbox::Norm(nDim, Normal);

    /*--- Initialize the convective & viscous residuals to zero ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar++) {Res_Visc[iVar] = 0.0;}

    /*--- Set the residual on the boundary with the specified heat flux ---*/
    // TODO: Look into this!
    // Note: Contributions from qtr and qve are used for proportional control
    //       to drive the solution toward the specified heatflux more quickly.
    const auto GradV  = nodes->GetGradient_Primitive(iPoint);
    su2double dTdn   = 0.0;
    su2double dTvedn = 0.0;
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      dTdn   += GradV[T_INDEX][iDim]*Normal[iDim];
      dTvedn += GradV[TVE_INDEX][iDim]*Normal[iDim];
    }
    const su2double ktr = nodes->GetThermalConductivity(iPoint);
    const su2double kve = nodes->GetThermalConductivity_ve(iPoint);

    /*--- Compute residual ---*/
    Res_Visc[nSpecies+nDim]   += pcontrol*(ktr*dTdn+kve*dTvedn) +
                                    Wall_HeatFlux*Area;
    Res_Visc[nSpecies+nDim+1] += pcontrol*(kve*dTvedn) +
                                    Wall_HeatFlux*Area;

    /*--- Impose the value of the velocity as a strong boundary
    condition (Dirichlet). Fix the velocity and remove any
    contribution to the residual at this node. ---*/
    su2double zero[MAXNDIM] = {0.0};
    nodes->SetVelocity_Old(iPoint, zero);

    for (auto iDim = 0ul; iDim < nDim; iDim++)
      LinSysRes(iPoint, nSpecies+iDim) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Apply viscous residual to the linear system ---*/
    LinSysRes.SubtractBlock(iPoint, Res_Visc);

    if (implicit) {
      /*--- Enforce the no-slip boundary condition in a strong way ---*/
      for (int iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
        auto total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }

  }
  END_SU2_OMP_FOR
}

void CNEMONSSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *sour_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {

  bool catalytic = config->GetCatalytic_Wall(val_marker);

  if (catalytic) { BC_HeatFluxCatalytic_Wall(geometry, solver_container, conv_numerics,
                                             sour_numerics, config, val_marker);

  } else { BC_HeatFluxNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                        sour_numerics, config, val_marker);
  }
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
  unsigned long iPoint, total_index;
  su2double rho, Ys;
  su2double *Normal, Area;
  su2double *Ds, *dYdn, SdYdn;
  su2double **GradY;

  /*--- Assign booleans ---*/
  bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  bool catalytic = false;

  /*--- Set "Proportional control" coefficient ---*/
  su2double pcontrol = 0.6;

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall heat flux from config ---*/
  su2double Wall_HeatFlux = config->GetWall_HeatFlux(Marker_Tag) / config->GetHeat_Flux_Ref();
  if (config->GetIntegrated_HeatFlux()) {
    Wall_HeatFlux /= geometry->GetSurfaceArea(config, val_marker);
  }

  /*--- Get the locations of the primitive variables ---*/
  const unsigned short T_INDEX    = nodes->GetTIndex();
  const unsigned short TVE_INDEX  = nodes->GetTveIndex();
  const unsigned short RHOS_INDEX = nodes->GetRhosIndex();
  const unsigned short RHO_INDEX  = nodes->GetRhoIndex();

  /*--- Allocate arrays ---*/
  dYdn = new su2double[nSpecies];
  GradY = new su2double*[nSpecies];
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (auto iVar = 0ul; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        Res_Sour[iVar] = 0.0;
      }

      /*--- Assign wall velocity to "Vector" array ---*/
      for (auto iDim = 0ul; iDim < nDim; iDim++) Vector[iDim] = 0.0;

      /*--- Set the residual, truncation error, and velocity value ---*/
      nodes->SetVelocity_Old(iPoint,Vector);
      for (auto iDim = 0ul; iDim < nDim; iDim++) {
        LinSysRes(iPoint, nSpecies+iDim) = 0.0;
        nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
      }

      /*--- Get temperature gradient information ---*/
      const auto V = nodes->GetPrimitive(iPoint);
      const auto GradV  = nodes->GetGradient_Primitive(iPoint);
      su2double dTdn   = 0.0;
      su2double dTvedn = 0.0;
      for (auto iDim = 0ul; iDim < nDim; iDim++) {
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
        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
          Ys = V[RHOS_INDEX+iSpecies]/rho;
          dYdn[iSpecies] = 0.0;
          for (auto iDim = 0ul; iDim < nDim; iDim++)
            dYdn[iSpecies] += 1.0/rho * (GradV[RHOS_INDEX+iSpecies][iDim] -
                Ys*GradV[RHO_INDEX][iDim])*Normal[iDim];
        }

        /*--- Calculate supplementary quantities ---*/
        SdYdn = 0.0;
        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
          SdYdn += rho*Ds[iSpecies]*dYdn[iSpecies];

        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
          Ys   = V[RHOS_INDEX+iSpecies]/rho;
          //eves = nodes->CalcEve(config, V[TVE_INDEX], iSpecies);
          //hs   = nodes->CalcHs(config, V[T_INDEX], eves, iSpecies);
          //          Res_Visc[iSpecies] = -rho*Ds[iSpecies]*dYdn[iSpecies] + Ys*SdYdn;
          //          Res_Visc[nSpecies+nDim]   += Res_Visc[iSpecies]*hs;
          //          Res_Visc[nSpecies+nDim+1] += Res_Visc[iSpecies]*eves;
        }
      }

      /*--- Get node thermal conductivity ---*/
      su2double ktr = nodes->GetThermalConductivity(iPoint);
      su2double kve = nodes->GetThermalConductivity_ve(iPoint);

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
      //        Jacobian.AddBlock2Diag(iPoint, iPoint, Jacobian_i);

      /*--- Only change velocity-rows of the Jacobian (includes 1 in the diagonal)/
       Note that we need to add a contribution for moving walls to the Jacobian. ---*/
      if (implicit) {
        /*--- Enforce the no-slip boundary condition in a strong way ---*/
        for (auto iVar = nSpecies; iVar < nSpecies+nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }

    }
  }

  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  delete [] dYdn;
}

void CNEMONSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container,
                                       CNumerics *conv_numerics, CNumerics *sour_numerics,
                                       CConfig *config, unsigned short val_marker) {

  bool catalytic = config->GetCatalytic_Wall(val_marker);

  if (catalytic) { BC_IsothermalCatalytic_Wall(geometry, solver_container, conv_numerics,
                                             sour_numerics, config, val_marker);

  } else { BC_IsothermalNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                          sour_numerics, config, val_marker);
  }
}

void CNEMONSSolver::BC_IsothermalNonCatalytic_Wall(CGeometry *geometry,
                                                   CSolver **solver_container,
                                                   CNumerics *conv_numerics,
                                                   CNumerics *sour_numerics,
                                                   CConfig *config,
                                                   unsigned short val_marker) {

  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool ionization = config->GetIonization();
  su2double UnitNormal[MAXNDIM] = {0.0};

  if (ionization) {
    SU2_MPI::Error("NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION",CURRENT_FUNCTION);
  }

  /*--- Define 'proportional control' constant ---*/
  const su2double C = 5;

  /*--- Identify the boundary ---*/
  const auto Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/
  const su2double Twall = config->GetIsothermal_Temperature(Marker_Tag);

  su2double **Jacobian_i = nullptr;
  if (implicit) {
    Jacobian_i = new su2double* [nVar];
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      Jacobian_i[iVar] = new su2double [nVar] ();
  }

  /*--- Loop over boundary points to calculate energy flux ---*/
  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Compute dual-grid area and boundary normal ---*/
    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
    const su2double Area = GeometryToolbox::Norm(nDim, Normal);
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    /*--- Compute closest normal neighbor ---*/
    const auto Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Compute distance between wall & normal neighbor ---*/
    const auto Coord_i = geometry->nodes->GetCoord(iPoint);
    const auto Coord_j = geometry->nodes->GetCoord(Point_Normal);

    const su2double dist_ij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

    /*--- Store the corrected velocity at the wall which will
     be zero (v = 0), unless there is grid motion (v = u_wall)---*/
    su2double zero[MAXNDIM] = {0.0};
    nodes->SetVelocity_Old(iPoint, zero);

    /*--- Initialize viscous residual to zero ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar ++) {Res_Visc[iVar] = 0.0;}

    for (auto iDim = 0ul; iDim < nDim; iDim++)
      LinSysRes(iPoint, nSpecies+iDim) = 0.0;
    nodes->SetVel_ResTruncError_Zero(iPoint);

    /*--- Calculate the gradient of temperature ---*/
    const su2double Ti   = nodes->GetTemperature(iPoint);
    const su2double Tj   = nodes->GetTemperature(Point_Normal);
    const su2double Tvei = nodes->GetTemperature_ve(iPoint);
    const su2double Tvej = nodes->GetTemperature_ve(Point_Normal);

    /*--- Rename variables for convenience ---*/
    const su2double ktr = nodes->GetThermalConductivity(iPoint);
    const su2double kve = nodes->GetThermalConductivity_ve(iPoint);

    /*--- Apply to the linear system ---*/
    Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                 (ktr*(Twall-Ti) + kve*(Twall-Tvei))*C)*Area/dist_ij;
    Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Twall-Tvei) *C)*Area/dist_ij;

    /*--- Calculate Jacobian for implicit time stepping ---*/
    if (implicit) {

      /*--- Initialize Jacobian to zero ---*/
      for (auto iVar = 0ul; iVar < nVar; iVar++)
        for (auto jVar = 0ul; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;

      /*--- Add contributions to the Jacobian from the weak enforcement of the energy equations. ---*/
      const auto dTdU   = nodes->GetdTdU(iPoint);
      const auto dTvedU = nodes->GetdTvedU(iPoint);
      const su2double theta = GeometryToolbox::SquaredNorm(nDim, UnitNormal);

      for (auto iVar = 0ul; iVar < nVar; iVar++) {
        Jacobian_i[nSpecies+nDim][iVar]   = -(ktr*theta/dist_ij*dTdU[iVar] +
                                              kve*theta/dist_ij*dTvedU[iVar])*Area;
        Jacobian_i[nSpecies+nDim+1][iVar] = - kve*theta/dist_ij*dTvedU[iVar]*Area;
      }
    } // implicit

    /*--- Convective and viscous contributions to the residual at the wall ---*/
    LinSysRes.SubtractBlock(iPoint, Res_Visc);

    /*--- Enforce the no-slip boundary condition in a strong way by
     modifying the velocity-rows of the Jacobian (1 on the diagonal).
     And add the contributions to the Jacobian due to energy. ---*/
    if (implicit) {
      Jacobian.SubtractBlock2Diag(iPoint, Jacobian_i);

      for (auto iVar = 1u; iVar <= nDim; iVar++) {
        auto total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }
    }
  }
  END_SU2_OMP_FOR

  if (Jacobian_i)
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
  delete [] Jacobian_i;

}

void CNEMONSSolver::BC_IsothermalCatalytic_Wall(CGeometry *geometry,
                                                CSolver **solver_container,
                                                CNumerics *conv_numerics,
                                                CNumerics *sour_numerics,
                                                CConfig *config,
                                                unsigned short val_marker) {

  /*--- Call standard isothermal BC to apply no-slip and energy b.c.'s ---*/
  BC_IsothermalNonCatalytic_Wall(geometry, solver_container, conv_numerics,
                                 sour_numerics, config, val_marker);

  /*--- Local variables ---*/
  su2double **GradY, **dVdU;

  /*--- Assign booleans ---*/
  const bool implicit = (config->GetKind_TimeIntScheme() == EULER_IMPLICIT);
  const bool Supercatalytic_Wall = config->GetSupercatalytic_Wall();

  /*--- Get universal information ---*/
  const su2double RuSI = UNIVERSAL_GAS_CONSTANT;
  const su2double Ru = 1000.0*RuSI;
  const auto& Ms = FluidModel->GetSpeciesMolarMass();

  /*--- Get the locations of the primitive variables ---*/
  const unsigned short RHOS_INDEX  = nodes->GetRhosIndex();
  const unsigned short RHO_INDEX   = nodes->GetRhoIndex();
  const unsigned short T_INDEX     = nodes->GetTIndex();
  const unsigned short TVE_INDEX   = nodes->GetTveIndex();

  /*--- Allocate arrays ---*/
  GradY = new su2double*[nSpecies];
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    GradY[iSpecies] = new su2double[nDim];
  dVdU = new su2double*[nVar];
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    dVdU[iVar] = new su2double[nVar];

  /*--- Loop over all of the vertices on this boundary marker ---*/
  for (auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->nodes->GetDomain(iPoint)) {

      /*--- Compute closest normal neighbor ---*/
      const auto jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Compute distance between wall & normal neighbor ---*/
      const auto Coord_i = geometry->nodes->GetCoord(iPoint);
      const auto Coord_j = geometry->nodes->GetCoord(jPoint);
      const su2double dij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

      /*--- Compute dual-grid area and boundary normal ---*/
      const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      const su2double Area = GeometryToolbox::Norm(nDim, Normal);

      /*--- Initialize the viscous residual to zero ---*/
      for (auto iVar = 0ul; iVar < nVar; iVar++) Res_Visc[iVar] = 0.0;

      /*--- Get primitive information ---*/
      const auto& Vi = nodes->GetPrimitive(iPoint);
      const auto& Vj = nodes->GetPrimitive(jPoint);
      const auto& Di = nodes->GetDiffusionCoeff(iPoint);
      const auto& eves = nodes->GetEve(iPoint);
      const auto& hs = FluidModel->ComputeSpeciesEnthalpy(Vi[T_INDEX], Vi[TVE_INDEX], eves);
      const su2double rho = Vi[RHO_INDEX];
      const auto& dTdU = nodes->GetdTdU(iPoint);
      const auto& dTvedU = nodes->GetdTvedU(iPoint);

      if (Supercatalytic_Wall) {

        const auto& Yst = config->GetSupercatalytic_Wall_Composition();

        /*--- Calculate supplementary quantities ---*/
        su2double dYdn = 0.0, SdYdn = 0.0;
        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
          dYdn = (Yst[iSpecies]-Vj[RHOS_INDEX+iSpecies]/Vj[RHO_INDEX])/dij;
          SdYdn += rho*Di[iSpecies]*dYdn;
        }

        /*--- Calculate species residual at wall ---*/
        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
          dYdn = (Yst[iSpecies]-Vj[RHOS_INDEX+iSpecies]/Vj[RHO_INDEX])/dij;
          Res_Visc[iSpecies]  = -(-rho*Di[iSpecies]*dYdn+Yst[iSpecies]*SdYdn)*Area;
        }

        if (implicit) {
          /*--- Initialize the transformation matrix ---*/
          for (auto iVar = 0ul; iVar < nVar; iVar++)
            for (auto jVar = 0ul; jVar < nVar; jVar++) {
              dVdU[iVar][jVar] = 0.0;
              Jacobian_j[iVar][jVar] = 0.0;
              Jacobian_i[iVar][jVar] = 0.0;
            }

          /*--- Populate transformation matrix ---*/
          // dYsdrk
          for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
            for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++)
              dVdU[iSpecies][jSpecies] += -1.0/rho*Yst[iSpecies];
            dVdU[iSpecies][iSpecies] += 1.0/rho;
          }
          for (auto iVar = 0ul; iVar < nVar; iVar++) {
            dVdU[nSpecies+nDim][iVar]   = dTdU[iVar];
            dVdU[nSpecies+nDim+1][iVar] = dTvedU[iVar];
          }

          /*--- Calculate supplementary quantities ---*/
          const auto& Cvtrs = FluidModel->GetSpeciesCvTraRot();
          const auto& Cvve = nodes->GetCvve(iPoint);

          /*--- Take the primitive var. Jacobian & store in Jac. jj ---*/
          // Species mass fraction
          for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
            for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++)
              Jacobian_j[iSpecies][jSpecies] += -Yst[iSpecies]*rho*Di[jSpecies]/dij;
            Jacobian_j[iSpecies][iSpecies] += rho*Di[iSpecies]/dij - SdYdn;
          }

          // Temperature
          for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
            for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
              Jacobian_j[nSpecies+nDim][iSpecies] += Jacobian_j[jSpecies][iSpecies]*hs[iSpecies];
            }
            Jacobian_j[nSpecies+nDim][nSpecies+nDim] += Res_Visc[iSpecies]/Area*(Ru/Ms[iSpecies] +
                                                                                 Cvtrs[iSpecies]  );
            Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
          }

          // Vib.-El. Temperature
          for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
            for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++)
              Jacobian_j[nSpecies+nDim+1][iSpecies] += Jacobian_j[jSpecies][iSpecies]*eves[iSpecies];
            Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += Res_Visc[iSpecies]/Area*Cvve[iSpecies];
          }

          /*--- Multiply by the transformation matrix and store in Jac. ii ---*/
          for (auto iVar = 0ul; iVar < nVar; iVar++)
            for (auto jVar = 0ul; jVar < nVar; jVar++)
              for (auto kVar = 0ul; kVar < nVar; kVar++)
                Jacobian_i[iVar][jVar] += Jacobian_j[iVar][kVar]*dVdU[kVar][jVar]*Area;

          /*--- Apply to the linear system ---*/
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
        }

      } else {

        if (implicit) {
          SU2_MPI::Error("Implicit not currently available for partially catalytic wall.",CURRENT_FUNCTION);
        }

        /*--- Identify the boundary ---*/
        string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

        /*--- Get isothermal wall temp ----*/
        const su2double Tw = config->GetIsothermal_Temperature(Marker_Tag);

        /*--- Get wall catalytic efficiency ----*/
        const su2double gam = config->GetCatalytic_Efficiency();

        /*--- Get cataltyic reaction map ---*/
        const auto& RxnTable = FluidModel->GetCatalyticRecombination();

        /*--- Common catalytic flux factor ---*/
        const su2double factor = gam*rho*sqrt(RuSI*Tw/2/PI_NUMBER)*Area;

        /*--- Compute catalytic recombination flux ---*/
        // Ref: 10.2514/6.2022-1636
        // ws = gam_s*Ys*rho_wall*sqrt(Ru*Tw/(2*Pi*M_combine)*Area
        for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
          int Index = SU2_TYPE::Int(RxnTable(iSpecies,1));
          Res_Visc[iSpecies] = RxnTable(iSpecies,0)*factor*Vi[Index]/Vi[RHO_INDEX]*sqrt(1/Ms[Index]);
        }
      }

      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
        Res_Visc[nSpecies+nDim]   += (Res_Visc[iSpecies]*hs[iSpecies])*Area;
        Res_Visc[nSpecies+nDim+1] += (Res_Visc[iSpecies]*eves[iSpecies])*Area;
      }

      /*--- Viscous contribution to the residual at the wall ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);
    }
  }

  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    delete [] GradY[iSpecies];
  delete [] GradY;
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    delete [] dVdU[iVar];
  delete [] dVdU;
}

void CNEMONSSolver::BC_Smoluchowski_Maxwell(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CNumerics *conv_numerics,
                                            CNumerics *visc_numerics,
                                            CConfig *config,
                                            unsigned short val_marker) {


  const bool ionization = config->GetIonization();

  if (ionization) {
    SU2_MPI::Error("NEED TO TAKE A CLOSER LOOK AT THE JACOBIAN W/ IONIZATION", CURRENT_FUNCTION);
  }

  /*--- Define 'proportional control' constant ---*/
  const su2double C = 1;

  /*---Define under-relaxation factors --- */
  const su2double alpha_V = 0.1;
  const su2double alpha_T = 1.0;

  /*--- Identify the boundary ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature and accomodation coefficients---*/
  const su2double Twall = config->GetIsothermal_Temperature(Marker_Tag);
  const su2double TMAC  = 1.0;
  const su2double TAC   = 1.0;

  /*--- Extract necessary indices ---*/
  const unsigned short T_INDEX   = nodes->GetTIndex();
  const unsigned short VEL_INDEX = nodes->GetVelIndex();
  const unsigned short TVE_INDEX = nodes->GetTveIndex();

  /*--- Loop over boundary points to calculate energy flux ---*/
  SU2_OMP_FOR_DYN(OMP_MIN_SIZE)
  for(auto iVertex = 0ul; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    const auto iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (!geometry->nodes->GetDomain(iPoint)) continue;

    /*--- Compute dual-grid area and boundary normal ---*/
    const auto Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
    const su2double Area = GeometryToolbox::Norm(nDim, Normal);

    su2double UnitNormal[MAXNDIM] = {0.0};
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    /*--- Compute closest normal neighbor ---*/
    const auto jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

    /*--- Compute distance between wall & normal neighbor ---*/
    const auto Coord_i = geometry->nodes->GetCoord(iPoint);
    const auto Coord_j = geometry->nodes->GetCoord(jPoint);

    const su2double dij = GeometryToolbox::Distance(nDim, Coord_i, Coord_j);

    /*--- Calculate Pressure ---*/
    const su2double Pi   = nodes->GetPressure(iPoint);

    /*--- Calculate the gradient of temperature ---*/
    const su2double Ti   = nodes->GetTemperature(iPoint);
    const su2double Tj   = nodes->GetTemperature(jPoint);
    const su2double Tvei = nodes->GetTemperature_ve(iPoint);
    const su2double Tvej = nodes->GetTemperature_ve(jPoint);

    /*--- Rename variables for convenience ---*/
    const su2double ktr  = nodes->GetThermalConductivity(iPoint);
    const su2double kve  = nodes->GetThermalConductivity_ve(iPoint);

    /*--- Retrieve Cv*density ---*/
    const su2double rhoCvtr = nodes->GetRhoCv_tr(iPoint);
    const su2double rhoCvve = nodes->GetRhoCv_ve(iPoint);

    /*--- Retrieve Flow Data ---*/
    const  su2double Viscosity = nodes->GetLaminarViscosity(iPoint);
    const su2double Density   = nodes->GetDensity(iPoint);
    const su2double Gamma     = nodes->GetGamma(iPoint);

    /*--- Retrieve Primitive Gradients ---*/
    const auto Grad_PrimVar = nodes->GetGradient_Primitive(iPoint);

    /*--- Calculate specific gas constant --- */
    const su2double GasConstant = FluidModel->ComputeGasConstant();

    /*--- Calculate temperature gradients normal to surface---*/ //Doubt about minus sign
    const su2double dTn = GeometryToolbox::DotProduct(nDim, Grad_PrimVar[T_INDEX], UnitNormal);
    const su2double dTven = GeometryToolbox::DotProduct(nDim, Grad_PrimVar[TVE_INDEX], UnitNormal);

    /*--- Calculate molecular mean free path ---*/
    const su2double Lambda = Viscosity/Density*sqrt(PI_NUMBER/(2.0*GasConstant*Ti));

    /*--- Calculate Temperature Slip ---*/
    su2double Tslip = ((2.0-TAC)/TAC)*2.0*Gamma/(Gamma+1.0)/Prandtl_Lam*Lambda*dTn+Twall;

    su2double Tslip_ve = Twall;
    if (dTven !=0) Tslip_ve = (Tslip-Twall)*(kve*rhoCvtr/dTn)/(ktr*rhoCvve/dTven)+Twall;

    /*--- Calculate temperature gradients tangent to surface ---*/
    su2double Vector_Tangent_dT[MAXNDIM] = {0.0}, Vector_Tangent_dTve[MAXNDIM] = {0.0};
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      Vector_Tangent_dT[iDim]   = Grad_PrimVar[T_INDEX][iDim] - dTn * UnitNormal[iDim];
      Vector_Tangent_dTve[iDim] = Grad_PrimVar[TVE_INDEX][iDim] - dTven * UnitNormal[iDim];
    }

    /*--- Calculate Heatflux tangent to surface ---*/
    su2double Vector_Tangent_HF[MAXNDIM] = {0.0};
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      Vector_Tangent_HF[iDim] = -ktr*Vector_Tangent_dT[iDim]-kve*Vector_Tangent_dTve[iDim];

    /*--- Initialize viscous residual to zero ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar ++) Res_Visc[iVar] = 0.0;

    /*--- Compute wall shear stress (using the stress tensor) ---*/
    su2double Tau[MAXNDIM][MAXNDIM] = {{0.0}};
    CNumerics::ComputeStressTensor(nDim, Tau, Grad_PrimVar+VEL_INDEX, Viscosity);

    su2double TauTangent[MAXNDIM] = {0.0};
    GeometryToolbox::TangentProjection(nDim,Tau,UnitNormal,TauTangent);

    /*--- Store the Slip Velocity at the wall */
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      Vector[iDim] = Lambda/Viscosity*(2.0-TMAC)/TMAC*(TauTangent[iDim])-3.0/4.0*(Gamma-1.0)/Gamma*Prandtl_Lam/Pi*Vector_Tangent_HF[iDim];

    /*--- Apply under-relaxation ---*/
    Tslip    = (1.0-alpha_T)*nodes->GetTemperature(iPoint)+(alpha_T)*Tslip;
    Tslip_ve = (1.0-alpha_T)*nodes->GetTemperature_ve(iPoint)+(alpha_T)*Tslip_ve;

    for (auto iDim = 0ul; iDim < nDim; iDim++){
      Vector[iDim] = (1.0-alpha_V)*nodes->GetVelocity(iPoint,iDim)+(alpha_V)*Vector[iDim];
    }

    nodes->SetVelocity_Old(iPoint,Vector);

    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      LinSysRes(iPoint, nSpecies+iDim) = 0.0;
      nodes->SetVal_ResTruncError_Zero(iPoint,nSpecies+iDim);
    }

    /*--- Apply to the linear system ---*/
    Res_Visc[nSpecies+nDim]   = ((ktr*(Ti-Tj)    + kve*(Tvei-Tvej)) +
                                 (ktr*(Tslip-Ti) + kve*(Tslip_ve-Tvei))*C)*Area/dij;
    Res_Visc[nSpecies+nDim+1] = (kve*(Tvei-Tvej) + kve*(Tslip_ve-Tvei)*C)*Area/dij;

    LinSysRes.SubtractBlock(iPoint, Res_Visc);
  }
  END_SU2_OMP_FOR
}

void CNEMONSSolver::SetTau_Wall_WF(CGeometry *geometry, CSolver **solver_container, const CConfig *config) {
    SU2_MPI::Error("Wall Functions not yet operational in NEMO.", CURRENT_FUNCTION);
}

