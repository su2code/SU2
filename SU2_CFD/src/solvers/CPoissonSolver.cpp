/*!
 * \file CPoissonSolver.cpp
 * \brief Main subroutines for solving the Poisson equation
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

#include "../../include/solvers/CPoissonSolver.hpp"
#include <cstddef>
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../include/solvers/CScalarSolver.inl"

/*--- Explicit instantiation of the parent class of CPoissonSolver. ---*/
template class CScalarSolver<CPoissonVariable>;

CPoissonSolver::CPoissonSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh)
  : CScalarSolver<CPoissonVariable>(geometry, config, false, LINEAR_SOLVER_MODE::POISSON) {
  SU2_ZONE_SCOPED

  /*--- Dimension of the problem --> pressure deviation is the only conservative variable ---*/

  nVar = 1;
  nPrimVar = 1;
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Define some structures for locating max residuals ---*/

  Residual_RMS.resize(nVar,0.0);
  Residual_Max.resize(nVar,0.0);
  Point_Max.resize(nVar,0);
  Point_Max_Coord.resize(nVar,nDim) = su2double(0.0);


  /*--- Initialization of the structure of the whole Jacobian ---*/

  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (poisson equation) MG level: " << iMesh << "." << endl;
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config, ReducerStrategy, false, false, true);
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  if (ReducerStrategy) EdgeFluxes.Initialize(geometry->GetnEdge(), geometry->GetnEdge(), nVar, nullptr);

  if (config->GetExtraOutput()) {
    if (nDim == 2) { nOutputVariables = 13; }
    else if (nDim == 3) { nOutputVariables = 19; }
    OutputVariables.Initialize(nPoint, nPointDomain, nOutputVariables, 0.0);
    OutputHeadingNames = new string[nOutputVariables];
  }

  /*--- Initialize the nodes vector. ---*/

  nodes = new CPoissonVariable(0.0, nPoint, nDim, nVar, config);

  SetBaseClassPointerToNodes();

  /*--- Communicate and store volume and the number of neighbors for any dual CVs that lie on on periodic markers. ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_VOLUME);
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_NEIGHBORS);
  }

  /*--- MPI solution ---*/

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- Add the solver name. ---*/

  SolverName = "POISSON";

}


void CPoissonSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)
                               
  /*--- Reset pressure corrections to zero for next iteration. ---*/
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->SetSolution(iPoint,0,0.0);

  /*--- Communicate updated Poisson solution (which should now be zero everywhere) ---*/
  solver_container[POISSON_SOL]->InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  solver_container[POISSON_SOL]->CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    solver_container[POISSON_SOL]->InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    solver_container[POISSON_SOL]->CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- Compute the gradients only after the solution has been reset to zero ---*/  
  CommonPreprocessing(geometry, config, Output);

  /*--- Reset flag for strong BCs. ---*/
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++)
    nodes->ResetStrongBC(iPoint);

  /*--- Need to clear EdgeFluxes and Jacobian. ---*/
  if (!Output) {
    if (ReducerStrategy) EdgeFluxes.SetValZero();
    Jacobian.SetValZero();
  }

  
}

void CPoissonSolver::Postprocessing(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CConfig *config,
                                    unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- The postprocessing function is responsible for applying the pressure and 
  velocity corrections based on the solution of the poisson equation. ---*/

  /*--- Compute gradients of the pressure correction p' so we can use it to find the velocity corrections ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) 
    SetSolution_Gradient_GG(geometry, config,false);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) 
    SetSolution_Gradient_LS(geometry, config,false);

  /*--- Start of computing the corrections ---*/
  unsigned long iEdge, iPoint, jPoint, iMarker, iVertex;
  unsigned short iDim, iVar, KindBC;
  su2double Vel, Current_Pressure, factor, PCorr_Ref, Vol, delT, Density;
  string Marker_Tag;
  su2activevector Pressure_Correc, alpha_p;
  su2activematrix vel_corr;
  long Pref_local;
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool piso = (config->GetPISO_corrections() > 1);

  CSolver* flow_solution = solver_container[FLOW_SOL];
  CVariable* flow_nodes = flow_solution->GetNodes();

  /*--- Allocate corrections and relaxation ---*/
  Pressure_Correc.resize(nPointDomain);
  vel_corr.resize(nPointDomain,nDim);
  alpha_p.resize(nPointDomain);

  /*--- Combine all pressure corrections into a vector for easy access ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    Pressure_Correc[iPoint] = nodes->GetSolution(iPoint,0);

  /*--- Define a reference pressure ---*/
  // TODO: look at this, currently copied (but working?) logic from old solver.
  unsigned long PRef_Point = 1;
  Pref_local = geometry->GetGlobal_to_Local_Point(PRef_Point);
  PCorr_Ref = 0.0;
  if (Pref_local >= 0)
    if(geometry->nodes->GetDomain(Pref_local))
      PCorr_Ref = 0.0;//Pressure_Correc[Pref_local];

  /*--- Compute Velocity Corrections and under relaxation factor for the pressure. ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    factor = 0.0;
    Vol = geometry->nodes->GetVolume(iPoint);
    delT = flow_nodes->GetDelta_Time(iPoint);
    const auto view = flow_solution->Jacobian.GetBlockView(iPoint, iPoint);
    for (iDim = 0; iDim < nDim; iDim++) {
      vel_corr[iPoint][iDim] = nodes->GetMomCoeff(iPoint)*(nodes->GetGradient(iPoint,0,iDim));
      if (implicit) factor += view(iDim, iDim);
    }

    /*--- The PISO algorithm should not underrelax pressure, for explicit iterations we simply do not have the necessary jacobian information to compute this ---*/
    alpha_p[iPoint] = config->GetRelaxation_Factor_PBFlow();
    if (implicit && !piso)
      alpha_p[iPoint] *= (Vol/delT) / (factor+(Vol/delT));      
  }

  /*--- Reassign strong boundary conditions ---*/
  /*--- For now I only have velocity inlet and fully developed outlet. Will need to add other types of inlet/outlet conditions
   *  where different treatment of pressure might be needed. Symmetry and Euler wall are weak BCs. ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag  = config->GetMarker_All_TagBound(iMarker);
    switch (KindBC) {
      case EULER_WALL: case SYMMETRY_PLANE:
        break;

      /*--- Nothing at MPI boundaries ---*/
      case SEND_RECEIVE:
        break;

      /*--- Only a fully developed outlet is implemented. For pressure, a dirichlet
            BC has to be applied and no correction is necessary. Velocity has a neumann BC. ---*/
      case OUTLET_FLOW:{
        auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);
        switch (Kind_Outlet) {
          case INC_OUTLET_TYPE::PRESSURE_OUTLET:{
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if (geometry->nodes->GetDomain(iPoint))
                Pressure_Correc[iPoint] = PCorr_Ref;
            }
            break;
          }
          //TODO: other outlet types
          default: 
            SU2_MPI::Error("The requested outflow boundary condition has not yet been implemented for the pressure based poisson solver", CURRENT_FUNCTION);
            break;
        }
        break;
      }

      /*--- Only a fixed velocity inlet is implemented now. Along with the wall boundaries,
        * the velocity is known and thus no correction is necessary.---*/
      case ISOTHERMAL: case HEAT_FLUX: case INLET_FLOW: {
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) {
            for (iDim = 0; iDim < nDim; iDim++)
              vel_corr[iPoint][iDim] = 0.0;
            alpha_p[iPoint] = 1.0;
            }
        }
        break;
      }

      /*--- Farfield is treated as a fully developed flow for pressure and a fixed pressure is
      * used, thus no correction is necessary. The treatment for velocity depends on whether the
      * flow is into the domain or out. If flow is in, a dirichlet bc is applied and no correction
      * is made, otherwise a Neumann BC is used and velocity is adjusted. ---*/

      case FAR_FIELD:
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->nodes->GetDomain(iPoint)) {
            // Check if the boundary condition is an inlet or not
            if (nodes->GetStrongBC(iPoint)) {
              for (iDim = 0; iDim < nDim; iDim++)
                vel_corr[iPoint][iDim] = 0.0;
            }
            Pressure_Correc[iPoint] = PCorr_Ref;
          }
        }
        
        break;

      default: 
        SU2_MPI::Error("The requested boundary condition has not yet been implemented for the pressure based poisson solver", CURRENT_FUNCTION);
        break;
    }
  }
  

  /*--- Apply corrections for new solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    /*--- Velocity corrections ---*/
    for (iVar = 0; iVar < nDim; iVar++) {
      Vel = flow_nodes->GetVelocity(iPoint,iVar);
      Vel -= vel_corr[iPoint][iVar];
      Density = flow_nodes->GetDensity(iPoint);
      flow_nodes->SetSolution(iPoint,iVar,Density*Vel);
    }
    flow_nodes->SetVelocity(iPoint);
    /*--- Pressure corrections ---*/
    Current_Pressure = flow_nodes->GetPressure(iPoint);
    Current_Pressure += alpha_p[iPoint]*(Pressure_Correc[iPoint] - PCorr_Ref);
    flow_nodes->SetPrimitive(iPoint,0,Current_Pressure);
  }

  /*--- periodic communication for both the momentum and the poisson equations as both are now updated ---*/
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    solver_container[FLOW_SOL]->InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    solver_container[FLOW_SOL]->CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);

    solver_container[FLOW_SOL]->InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
    solver_container[FLOW_SOL]->CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_PRESSURE);
  }

  /*--- Communicate updated velocities and pressure ---*/
  solver_container[FLOW_SOL]->InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  solver_container[FLOW_SOL]->CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  solver_container[FLOW_SOL]->InitiateComms(geometry, config, MPI_QUANTITIES::PRESSURE_VAR);
  solver_container[FLOW_SOL]->CompleteComms(geometry, config, MPI_QUANTITIES::PRESSURE_VAR);

}


void CPoissonSolver::SetMomCoeff(CGeometry *geometry, CSolver **solver_container, CConfig *config, bool periodic, unsigned short iMesh) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iPoint, jPoint, iNeigh;
  su2double Mom_Coeff, Mom_Coeff_nb, Vol, delT;
  bool simplec = (config->GetKind_PBIter() == ENUM_PBITER::SIMPLEC);
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  const CSolver* flow_solution = solver_container[FLOW_SOL];
  const CVariable* flow_nodes = flow_solution->GetNodes();
  

  if (!periodic)  {

    if (implicit) {
      /* First sum up the momentum coefficient using the jacobian from given point and it's neighbors. ---*/
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        /*--- Self contribution of the coefficient a_p. Note that this coefficient should be the same for all variable directions. therefore just the x-momentum coefficient is taken. ---*/
        Mom_Coeff = flow_solution->Jacobian.GetBlockView(iPoint, iPoint)(0,0);
      
        Mom_Coeff_nb = 0.0;

        if (simplec) {
          for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
            jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
            Mom_Coeff_nb += flow_solution->Jacobian.GetBlockView(iPoint, iPoint)(0,0);
          }
        }

        Vol = geometry->nodes->GetVolume(iPoint); delT = flow_nodes->GetDelta_Time(iPoint);

        /*--- Add simplec neighbour contributions and optional time dependent term. ---*/
        Mom_Coeff = Mom_Coeff - Mom_Coeff_nb - config->GetRCFactor()*(Vol/delT);

        /*--- Invert the momentum coefficient to 1/a_p and scale by the volume so it can be used as diffusion coefficient in the poisson eq ---*/
        Mom_Coeff = Vol/Mom_Coeff;

        nodes->SetMomCoeff(iPoint, Mom_Coeff);
      }
    }
    else {
      for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

        Vol = geometry->nodes->GetVolume(iPoint); delT = flow_nodes->GetDelta_Time(iPoint);

        Mom_Coeff = delT;

        nodes->SetMomCoeff(iPoint, Mom_Coeff);
      }
    }

    /*--- Insert MPI call here. ---*/
    InitiateComms(geometry, config, MPI_QUANTITIES::MOM_COEFF);
    CompleteComms(geometry, config, MPI_QUANTITIES::MOM_COEFF); 
  }
}


void CPoissonSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  SU2_ZONE_SCOPED
 
  CNumerics* numerics = numerics_container[VISC_TERM + omp_get_thread_num() * MAX_TERMS];

  bool pausePreacc = false;
  if (ReducerStrategy)
    pausePreacc = AD::PausePreaccumulation();
  else
    AD::StartNoSharedReading();

  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];
      Viscous_Residual(iEdge, geometry, solver_container, numerics, config);
    }
    END_SU2_OMP_FOR
  }

  /*--- Restore preaccumulation and adjoint evaluation state. ---*/
  AD::ResumePreaccumulation(pausePreacc);
  if (!ReducerStrategy) AD::EndNoSharedReading();

  if (ReducerStrategy) {
    SumEdgeFluxes(geometry);
    Jacobian.SetDiagonalAsColumnSum();
  }
}

void CPoissonSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics_container,
                                  CConfig *config, unsigned short iMesh) {
  SU2_ZONE_SCOPED

  unsigned short iVar, iDim, KindBC;
  unsigned long iPoint, jPoint, iEdge, iMarker, iVertex;
  su2double Edge_Vector[MAXNDIM], dist_ij_2;
  su2double *Coord_i, *Coord_j;
  su2double MassFlux_Part, MassFlux_Avg, Mom_Coeff[MAXNDIM], *Normal,Vel_Avg, Grad_Avg;
  su2double Area, MeanDensity, Vol , TimeStep;
  su2double GradP_f[MAXNDIM], GradP_in[MAXNDIM], GradP_proj, RhieChowInterp, Coeff_Mom, PsCorr[MAXNDIM], PsCorrFace;
  su2double *Flow_Dir, Flow_Dir_Mag,*GridVel_i,*GridVel_j;
  string Marker_Tag;
  Normal = new su2double [MAXNDIM];

  bool unsteady = (config->GetTime_Marching() != TIME_MARCHING::STEADY);

  const CSolver* flow_solution = solver_container[FLOW_SOL];
  const CVariable* flow_nodes = flow_solution->GetNodes();

  /*--- Mass flux is computed over all edges ---*/
  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      iPoint = geometry->edges->GetNode(iEdge,0); jPoint = geometry->edges->GetNode(iEdge,1);
      geometry->edges->GetNormal(iEdge, Normal);

      Area = GeometryToolbox::Norm(nDim, Normal);

      MeanDensity = 0.5*(flow_nodes->GetDensity(iPoint) + flow_nodes->GetDensity(jPoint));
      
      if (dynamic_grid) {
        GridVel_i = geometry->nodes->GetGridVel(iPoint);
        GridVel_j = geometry->nodes->GetGridVel(jPoint);
      }

      /*--- Face average mass flux. ---*/
      MassFlux_Avg = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Vel_Avg = 0.5*(flow_nodes->GetVelocity(iPoint,iDim) + flow_nodes->GetVelocity(jPoint,iDim));
        MassFlux_Avg += MeanDensity*Vel_Avg*Normal[iDim];
      }
      
      if (dynamic_grid)
        for (iDim = 0; iDim < nDim; iDim++) {
          Vel_Avg = 0.5*(GridVel_i[iDim]+GridVel_j[iDim]);
          MassFlux_Avg -= MeanDensity*Vel_Avg*Normal[iDim];
        }

      /*--- Rhie Chow interpolation ---*/
      Coord_i = geometry->nodes->GetCoord(iPoint);
      Coord_j = geometry->nodes->GetCoord(jPoint);
      dist_ij_2 = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
        dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
      }
      /*--- 1. Interpolate the pressure gradient based on node values ---*/
      for (iDim = 0; iDim < nDim; iDim++) {
        Grad_Avg = 0.5*(flow_nodes->GetGradient_Primitive(iPoint,0,iDim) + flow_nodes->GetGradient_Primitive(jPoint,0,iDim));
        GradP_in[iDim] = Grad_Avg;
      }

      /*--- 2. Compute pressure gradient at the face ---*
          Eq 15.62 F Moukalled, L Mangani M. Darwish OpenFOAM and uFVM book. ---*/
      GradP_proj = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        GradP_proj += GradP_in[iDim]*Edge_Vector[iDim];
      }
      if (dist_ij_2 != 0.0) {
        for (iDim = 0; iDim < nDim; iDim++) {
          GradP_f[iDim] = GradP_in[iDim] - (GradP_proj - (flow_nodes->GetPressure(jPoint) - flow_nodes->GetPressure(iPoint)))*Edge_Vector[iDim]/ dist_ij_2;
        }
      }

      /*--- Correct the massflux by adding the pressure terms.
      * --- GradP_f is the gradient computed directly at the face and GradP_in is the
      * --- gradient linearly interpolated based on node values. This effectively adds a third
      * --- order derivative of pressure to remove odd-even decoupling of pressure and velocities.
      * --- GradP_f = (p_F^n - p_P^n)/ds , GradP_in = 0.5*(GradP_P^n + GradP_F^n)---*/
      RhieChowInterp = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {

        /*--- Linearly interpolated coefficient. ---*/

        Coeff_Mom = 0.5*(nodes->GetMomCoeff(iPoint) + nodes->GetMomCoeff(jPoint));

        /*--- Difference of pressure gradients. ---*/

        RhieChowInterp += Coeff_Mom*(GradP_f[iDim] - GradP_in[iDim])*Normal[iDim]*MeanDensity;

        /*--- Save the pressure gradient contribution for the correction term used in the next iteration. ---*/
        PsCorr[iDim] = -Coeff_Mom*(GradP_f[iDim] - GradP_in[iDim]);
      }

      /*--- Rhie Chow correction for time step must go here ---*/
      // It is currently not included.

      /*--- Calculate the mass flux at the face including the linearly interpolated velocities, pressure 
      *    gradient difference contribution. The correction for time stepping (both pseudo and dual time) is currently not 
      *    included as testing has shown it to have either no effect or worsen convergence. ---*/ 
      MassFlux_Part = MassFlux_Avg - RhieChowInterp;
      
      auto residual = CNumerics::ResidualType<>(&MassFlux_Part, nullptr, nullptr);

      if (geometry->nodes->GetDomain(iPoint)) LinSysRes.AddBlock(iPoint, residual);
      if (geometry->nodes->GetDomain(jPoint)) LinSysRes.SubtractBlock(jPoint, residual);

    }
    END_SU2_OMP_FOR
  }

  /*--- Now add corrections to the previously computed mass fluxes for boundary conditions which alter the mass flux ---*/

  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    KindBC = config->GetMarker_All_KindBC(iMarker);
    Marker_Tag  = config->GetMarker_All_TagBound(iMarker);

    switch (KindBC) {
      /*--- Wall boundaries have zero mass flux (irrespective of grid movement) ---*/
      case EULER_WALL: case ISOTHERMAL: case HEAT_FLUX: case SYMMETRY_PLANE:
      break;

      /*--- Nothing has to happen at MPI boundaries*/
      case SEND_RECEIVE:
        break;

      case INLET_FLOW:
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (!geometry->nodes->GetDomain(iPoint)) continue;

          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);            
              
          MassFlux_Part = 0.0;
          if (dynamic_grid) {
            GridVel_i = geometry->nodes->GetGridVel(iPoint);
            for (iDim = 0; iDim < nDim; iDim++)
              MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
          } 
          else
            for (iDim = 0; iDim < nDim; iDim++)
            MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim))*Normal[iDim];

          auto residual = CNumerics::ResidualType<>(&MassFlux_Part, nullptr, nullptr);

          if (geometry->nodes->GetDomain(iPoint)) LinSysRes.AddBlock(iPoint, residual);

        }
        break;

      case FAR_FIELD:
        /*--- Treat the farfield as a fully developed outlet for pressure. ---*/
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (geometry->nodes->GetDomain(iPoint)) {
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

            if (dynamic_grid)
              GridVel_i = geometry->nodes->GetGridVel(iPoint);
                
            MassFlux_Part = 0.0;
            if (dynamic_grid)
              for (iDim = 0; iDim < nDim; iDim++)
                MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
            else
             for (iDim = 0; iDim < nDim; iDim++)
              MassFlux_Part -= nodes->GetDensity(iPoint)*(nodes->GetVelocity(iPoint, iDim))*Normal[iDim];
  
            auto residual = CNumerics::ResidualType<>(&MassFlux_Part, nullptr, nullptr);
            LinSysRes.AddBlock(iPoint, residual);    

          }
        }
        break;


      case OUTLET_FLOW:{
      /*--- Note I am assuming a fully developed outlet, thus the pressure value is prescribed
       * -- and a dirichlet bc has to be applied along outlet faces. The Massflux, which forms the RHS
       * -- of the equation, is set to zero to enforce the dirichlet bc. ---*/

        auto Kind_Outlet = config->GetKind_Inc_Outlet(Marker_Tag);

        switch (Kind_Outlet) {
          case INC_OUTLET_TYPE::PRESSURE_OUTLET:
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

              if (geometry->nodes->GetDomain(iPoint)) {
                geometry->vertex[iMarker][iVertex]->GetNormal(Normal);

                if (dynamic_grid)
                  GridVel_i = geometry->nodes->GetGridVel(iPoint);
                
                MassFlux_Part = 0.0;
                if (dynamic_grid)
                  for (iDim = 0; iDim < nDim; iDim++)
                    MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
                else
                  for (iDim = 0; iDim < nDim; iDim++)
                    MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim))*Normal[iDim];

                auto residual = CNumerics::ResidualType<>(&MassFlux_Part, nullptr, nullptr);

                if (geometry->nodes->GetDomain(iPoint)) LinSysRes.AddBlock(iPoint, residual);
              }
            }
            break;
          default:
            SU2_MPI::Error("Requested type of outlet boundary condition not available", CURRENT_FUNCTION);
            break;
        }
        break;
      }

      default:
        SU2_MPI::Error("Invalid boundary condition for mass flux correction", CURRENT_FUNCTION);
        break;
        
    }
  }

}

void CPoissonSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  SU2_ZONE_SCOPED

  /*--- No actual time integration is done here. The routine is used as a means to solve the linear equation 
   * resulting from the poisson equation. The linear system is solved using the jacobian matrix in a way
   * consistent with the rest of the code. The time step is set to zero and no under-relaxation is applied to the
   * jacobian matrix. ---*/
  
  unsigned long iPoint, total_index, IterLinSol = 0;;
  unsigned short iVar;
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;

  SetResToZero();

  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/
    local_Res_TruncError = nodes->GetResTruncError(iPoint);

    /*--- Read the volume ---*/
    Vol = geometry->nodes->GetVolume(iPoint);

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar] );
      LinSysSol[total_index] = 0.0;
      Residual_RMS[iVar] += LinSysRes[total_index]*LinSysRes[total_index];
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->nodes->GetGlobalIndex(iPoint), geometry->nodes->GetCoord(iPoint));
    }
  }

  /*--- Initialize residual and solution at the ghost points ---*/
  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }
  }

  /*--- Solve or smooth the linear system ---*/
  IterLinSol = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  /*--- Store the value of the residual. ---*/
  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    SetIterLinSolver(IterLinSol);
    SetResLinSolver(System.GetResidual());
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint, iVar, LinSysSol[iPoint*nVar+iVar]);
     }
  }

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic() / 2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);

}

void CPoissonSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  /*--- Zero flux (Neumann) BC on pressure ---*/
}

// void CPoissonSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//   /*--- Zero flux (Neumann) BC on pressure ---*/
// }


void CPoissonSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
                                CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar;
  su2double pressureDeviation = 0.0;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (!geometry->nodes->GetDomain(iPoint)) continue;
    /*--- The farfield boundary is considered as an inlet-outlet boundary, where flow 
      * can either enter or leave. For pressure, it is treated as a fully developed flow
      * and a dirichlet BC is applied. For velocity, based on the sign of massflux, either 
      * a dirichlet or a neumann BC is applied (in source routine). ---*/		

    LinSysRes.SetBlock_Zero(iPoint);

    nodes->SetSolution(iPoint, &pressureDeviation);
    nodes->SetSolution_Old(iPoint,&pressureDeviation);
    nodes->SetStrongBC(iPoint);
    Jacobian.DeleteValsRowi(iPoint, 0);
  }
}

void CPoissonSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
	/*--- Zero flux (Neumann) BC on pressure ---*/
}

void CPoissonSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iVertex, iPoint;
  unsigned short iDim, iVar;
  su2double pressureDeviation = 0.0;

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (!geometry->nodes->GetDomain(iPoint)) continue;
    
    /*--- apply a dirichlet boundary condition as pressure is prescribed*/

    LinSysRes.SetBlock_Zero(iPoint);

    nodes->SetSolution(iPoint, &pressureDeviation);
    nodes->SetSolution_Old(iPoint,&pressureDeviation);
    // nodes->SetStrongBC(iPoint);
    Jacobian.DeleteValsRowi(iPoint, 0);
  }
}