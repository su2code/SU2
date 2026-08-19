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
  : CScalarSolver<CPoissonVariable>(geometry, config, false, false, LINEAR_SOLVER_MODE::POISSON) {
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

  // PseudoTimeCorr.resize(geometry->GetnEdge(),nDim) = su2double(0.0);

}


void CPoissonSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh,
                                unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  SU2_ZONE_SCOPED
  SU2_OMP_SAFE_GLOBAL_ACCESS(config->SetGlobalParam(config->GetKind_Solver(), RunTime_EqSystem);)
                               
  /*--- Reset pressure corrections to zero for next iteration. ---*/
  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    nodes->SetSolution(iPoint,0,0.0);
  }
  END_SU2_OMP_FOR

  /*--- Communicate updated Poisson solution (which should now be zero everywhere) ---*/
  solver_container[POISSON_SOL]->InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  solver_container[POISSON_SOL]->CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    solver_container[POISSON_SOL]->InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    solver_container[POISSON_SOL]->CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }

  /*--- Compute the gradients only after the solution has been reset to zero ---*/  
  CommonPreprocessing(geometry, config, Output);

  /*--- Need to clear EdgeFluxes and Jacobian. ---*/
  if (!Output) {
    LinSysRes.SetValZero();
    if (ReducerStrategy) EdgeFluxes.SetValZero();
    Jacobian.SetValZero();
  }

  
}

void CPoissonSolver::Postprocessing(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CConfig *config,
                                    unsigned short iMesh) {
  SU2_ZONE_SCOPED

  /*--- Compute gradients of the pressure correction p' so we can use it to find the velocity corrections ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) 
    SetSolution_Gradient_GG(geometry, config,false);

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) 
    SetSolution_Gradient_LS(geometry, config,false);

}


void CPoissonSolver::SetMomCoeff(CGeometry *geometry, CSolver **solver_container, CConfig *config, bool periodic, unsigned short iMesh) {

  bool simplec = (config->GetKind_PBIter() == ENUM_PBITER::SIMPLEC);
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  const CSolver* flow_solution = solver_container[FLOW_SOL];
  const CVariable* flow_nodes = flow_solution->GetNodes();
  
  if (implicit) {
    /* First sum up the momentum coefficient using the jacobian from given point and it's neighbors. ---*/
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      /*--- Self contribution of the coefficient a_p. Note that this coefficient should be the same for all variable directions. therefore just the x-momentum coefficient is taken. ---*/
      su2double Mom_Coeff = flow_solution->Jacobian.GetBlockView(iPoint, iPoint)(1,1);
    
      su2double Mom_Coeff_nb = 0.0;

      if (simplec) {
        for (unsigned long iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
          auto jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
          Mom_Coeff_nb += flow_solution->Jacobian.GetBlockView(iPoint, jPoint)(1,1);
        }
      }

      su2double Vol = geometry->nodes->GetVolume(iPoint); 
      su2double delT = flow_nodes->GetDelta_Time(iPoint);

      /*--- Add simplec neighbour contributions and optional time dependent term. ---*/
      Mom_Coeff = Mom_Coeff - Mom_Coeff_nb - config->GetRCFactor()*(Vol/delT);

      /*--- Invert the momentum coefficient to 1/a_p and scale by the volume so it can be used as diffusion coefficient in the poisson eq ---*/
      Mom_Coeff = Vol/Mom_Coeff;

      nodes->SetMomCoeff(iPoint, Mom_Coeff);
    }
    END_SU2_OMP_FOR
  }
  else {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

      su2double delT = flow_nodes->GetDelta_Time(iPoint);

      su2double Mom_Coeff = delT;

      nodes->SetMomCoeff(iPoint, Mom_Coeff);
    }
    END_SU2_OMP_FOR
  }

  /*--- Insert MPI call here. ---*/
  InitiateComms(geometry, config, MPI_QUANTITIES::MOM_COEFF);
  CompleteComms(geometry, config, MPI_QUANTITIES::MOM_COEFF); 
}


void CPoissonSolver::ComputeHbyA(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

  unsigned short iVar, jVar, iDim, jDim;
  unsigned long iPoint, jPoint, iNeigh;
  su2double Mom_Coeff, Mom_Coeff_nb, Vol, delT;
  bool simplec = (config->GetKind_PBIter() == ENUM_PBITER::SIMPLEC);
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  const CSolver* flow_solution = solver_container[FLOW_SOL];
  const CVariable* flow_nodes = flow_solution->GetNodes();

  /*--- First exchange momentum correction which is required to compute H. ---*/
  InitiateComms(geometry, config, MPI_QUANTITIES::VEL_CORRECTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::VEL_CORRECTION); 
  
  if (implicit) {
    SU2_OMP_FOR_STAT(omp_chunk_size)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      for (iDim = 0; iDim < nDim; ++iDim) {
        su2double H = 0.0;
        su2double A_P = flow_solution->Jacobian.GetBlockView(iPoint, iPoint)(1,1);
        for (iNeigh = 0; iNeigh < geometry->nodes->GetnPoint(iPoint); iNeigh++) {
          jPoint = geometry->nodes->GetPoint(iPoint,iNeigh);
          H -= flow_solution->Jacobian.GetBlockView(iPoint, jPoint)(1,1) * nodes->GetVelocityCorrection(jPoint, iDim);
        }
        nodes->SetHbyACorrection(iPoint, iDim, H/A_P);
      }
    }
    END_SU2_OMP_FOR
  }
  else {
    SU2_MPI::Error("HbyA is currently not supported for an explicit momentum solver.", CURRENT_FUNCTION);
  }

  /*--- Exchange HbyA with MPI call. ---*/
  InitiateComms(geometry, config, MPI_QUANTITIES::HBYA_CORRECTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::HBYA_CORRECTION); 
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

  su2double *GridVel_i;
  
  const CSolver* flow_solver = solver_container[FLOW_SOL];
  const CVariable* flow_nodes = flow_solver->GetNodes();

  const auto& edgeVelocities = *(flow_solver->GetEdgeVelocity());

  /*--- Mass flux is computed over all edges ---*/
  for (auto color : EdgeColoring) {
    SU2_OMP_FOR_DYN(nextMultiple(OMP_MIN_SIZE, color.groupSize))
    for (auto k = 0ul; k < color.size; ++k) {
      auto iEdge = color.indices[k];

      auto iPoint = geometry->edges->GetNode(iEdge,0); auto jPoint = geometry->edges->GetNode(iEdge,1);
      su2double Normal[MAXNDIM] = {0.0};
      geometry->edges->GetNormal(iEdge, Normal);

      su2double MeanDensity = 0.5*(flow_nodes->GetDensity(iPoint) + flow_nodes->GetDensity(jPoint));

      su2double MassFlux_Part = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim) 
        MassFlux_Part += edgeVelocities[iEdge][iDim] * Normal[iDim] * MeanDensity;

      /*--- Add the mass flux to the source term for the poisson equation ---*/
      auto residual = CNumerics::ResidualType<>(&MassFlux_Part, nullptr, nullptr);

      if (geometry->nodes->GetDomain(iPoint)) LinSysRes.AddBlock(iPoint, residual);
      if (geometry->nodes->GetDomain(jPoint)) LinSysRes.SubtractBlock(jPoint, residual);

      /*--- Only for the second pressure correction in the case PISO is used, we need the additional HbyA(u') term ---*/
      // TODO: currently its just set to zero and does not contribute for the first piso correctin but would be nice if this entire block would be skipped otherwise.
      su2double MeanHbyA = 0.0;
      for (unsigned short iDim = 0; iDim < nDim; ++iDim)
        MeanHbyA += 0.5 * (nodes->GetHbyACorrection(iPoint, iDim) + nodes->GetHbyACorrection(jPoint, iDim)) * Normal[iDim]; 

      auto residualHbyA = CNumerics::ResidualType<>(&MeanHbyA, nullptr, nullptr);
      if (geometry->nodes->GetDomain(iPoint)) LinSysRes.AddBlock(iPoint, residualHbyA);
      if (geometry->nodes->GetDomain(jPoint)) LinSysRes.SubtractBlock(jPoint, residualHbyA);

    }
    END_SU2_OMP_FOR
  }

  /*--- Now add corrections to the previously computed mass fluxes for boundary conditions which alter the mass flux ---*/
  unsigned short iDim, KindBC;
  unsigned long  iMarker, iVertex, iPoint, jPoint;
  string Marker_Tag;
  su2double MassFlux_Part = 0.0, Normal[MAXNDIM];

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
                MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim)-GridVel_i[iDim])*Normal[iDim];
            else
             for (iDim = 0; iDim < nDim; iDim++)
              MassFlux_Part -= flow_nodes->GetDensity(iPoint)*(flow_nodes->GetVelocity(iPoint, iDim))*Normal[iDim];
  
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
  
  unsigned long total_index;

  /*--- Local residual variables for current thread ---*/
  su2double resMax[MAXNVAR] = {0.0}, resRMS[MAXNVAR] = {0.0};
  unsigned long idxMax[MAXNVAR] = {0};

  SetResToZero();

  /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
  SU2_OMP_FOR_(schedule(static,omp_chunk_size) SU2_NOWAIT)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Multigrid contribution to residual. ---*/
    su2double *local_Res_TruncError = nodes->GetResTruncError(iPoint);

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      LinSysRes(iPoint, iVar) = - (LinSysRes(iPoint, iVar) + local_Res_TruncError[iVar] );
      LinSysSol(iPoint, iVar) = 0.0;

      /*--- "Add" residual at (iPoint,iVar) to local residual variables. ---*/
      ResidualReductions_PerThread(iPoint, iVar, LinSysRes(iPoint, iVar), resRMS, resMax, idxMax);
    }
  }
  END_SU2_OMP_FOR

  /*--- "Add" residuals from all threads to global residual variables. ---*/
  ResidualReductions_FromAllThreads(geometry, config, resRMS, resMax, idxMax);

  /*--- Solve or smooth the linear system. ---*/

  SU2_OMP_FOR_(schedule(static,OMP_MIN_SIZE) SU2_NOWAIT)
  for (unsigned long iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }
  END_SU2_OMP_FOR

  auto iter = System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  BEGIN_SU2_OMP_SAFE_GLOBAL_ACCESS {
    SetIterLinSolver(iter);
    SetResLinSolver(System.GetResidual());
  }
  END_SU2_OMP_SAFE_GLOBAL_ACCESS

  SU2_OMP_FOR_STAT(omp_chunk_size)
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      nodes->AddSolution(iPoint, iVar, LinSysSol(iPoint,iVar));
     }
  }
  END_SU2_OMP_FOR

  InitiateComms(geometry, config, MPI_QUANTITIES::SOLUTION);
  CompleteComms(geometry, config, MPI_QUANTITIES::SOLUTION);

}

void CPoissonSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  /*--- Zero flux (Neumann) BC on pressure ---*/
}

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
    Jacobian.DeleteValsRowi(iPoint, 0);
  }
}