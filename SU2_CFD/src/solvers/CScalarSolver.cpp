/*!
 * \file solution_direct_scalar.cpp
 * \brief Main subroutines for solving scalar transport equations.
 * \author T. Economon
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/solvers/CScalarSolver.hpp"

CScalarSolver::CScalarSolver(void) : CSolver() {
  
  FlowPrimVar_i    = NULL;
  FlowPrimVar_j    = NULL;
  lowerlimit       = NULL;
  upperlimit       = NULL;
  nVertex          = NULL;
  nMarker          = 0;
  Inlet_ScalarVars = NULL;
  Scalar_Inf       = NULL;
  
  /*--- The turbulence models are always solved implicitly, so set the
   implicit flag in case we have periodic BCs. ---*/
  
  SetImplicitPeriodic(true);
  
}

CScalarSolver::CScalarSolver(CGeometry* geometry, CConfig *config) : CSolver() {
  
  FlowPrimVar_i    = NULL;
  FlowPrimVar_j    = NULL;
  lowerlimit       = NULL;
  upperlimit       = NULL;
  nMarker          = config->GetnMarker_All();
  Inlet_ScalarVars = NULL;
  Scalar_Inf       = NULL;
  
  /*--- Store the number of vertices on each marker for deallocation later ---*/
  
  nVertex = new unsigned long[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++)
    nVertex[iMarker] = geometry->nVertex[iMarker];
  
  /*--- The turbulence models are always solved implicitly, so set the
   implicit flag in case we have periodic BCs. ---*/
  
  SetImplicitPeriodic(true);
  
}

CScalarSolver::~CScalarSolver(void) {
  
  if (Inlet_ScalarVars != NULL) {
    for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Inlet_ScalarVars[iMarker] != NULL) {
        for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          delete [] Inlet_ScalarVars[iMarker][iVertex];
        }
        delete [] Inlet_ScalarVars[iMarker];
      }
    }
    delete [] Inlet_ScalarVars;
  }
  
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  if (lowerlimit != NULL)    delete [] lowerlimit;
  if (upperlimit != NULL)    delete [] upperlimit;
  if (nVertex != NULL)       delete [] nVertex;
  if (Scalar_Inf != NULL)    delete [] Scalar_Inf;
  
}

void CScalarSolver::Upwind_Residual(CGeometry *geometry,
                                    CSolver **solver_container,
                                    CNumerics *numerics,
                                    CConfig *config,
                                    unsigned short iMesh) {
  
  su2double *Scalar_i, *Scalar_j, *Limiter_i = NULL, *Limiter_j = NULL;
  su2double *V_i, *V_j, **Gradient_i, **Gradient_j;
  su2double Project_Grad_i, Project_Grad_j;
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool muscl         = (config->GetMUSCL_Scalar() && (iMesh == MESH_0));
  bool limiter       = (config->GetKind_SlopeLimit_Scalar() != NO_LIMITER);
  bool grid_movement = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
    V_j = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(jPoint);
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Scalar variables w/o reconstruction ---*/
    
    Scalar_i = nodes->GetSolution(iPoint);
    Scalar_j = nodes->GetSolution(jPoint);
    numerics->SetScalarVar(Scalar_i, Scalar_j);
    
    /*--- Grid Movement ---*/
    
    if (grid_movement)
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Mean flow primitive variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(iPoint);
      Gradient_j = solver_container[FLOW_SOL]->GetNodes()->GetGradient_Primitive(jPoint);
      if (limiter) {
        Limiter_i = solver_container[FLOW_SOL]->GetNodes()->GetLimiter_Primitive(iPoint);
        Limiter_j = solver_container[FLOW_SOL]->GetNodes()->GetLimiter_Primitive(jPoint);
      }
      
      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnPrimVarGrad(); iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          FlowPrimVar_i[iVar] = V_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          FlowPrimVar_i[iVar] = V_i[iVar] + Project_Grad_i;
          FlowPrimVar_j[iVar] = V_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetPrimitive(FlowPrimVar_i, FlowPrimVar_j);
      
      /*--- Scalar variables using gradient reconstruction and limiters ---*/
      
      Gradient_i = nodes->GetGradient(iPoint);
      Gradient_j = nodes->GetGradient(jPoint);
      if (limiter) {
        Limiter_i = nodes->GetLimiter(iPoint);
        Limiter_j = nodes->GetLimiter(jPoint);
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0.0; Project_Grad_j = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
        }
        if (limiter) {
          Solution_i[iVar] = Scalar_i[iVar] + Limiter_i[iVar]*Project_Grad_i;
          Solution_j[iVar] = Scalar_j[iVar] + Limiter_j[iVar]*Project_Grad_j;
        }
        else {
          Solution_i[iVar] = Scalar_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Scalar_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetScalarVar(Solution_i, Solution_j);
      
    }
    
    /*--- Add and subtract residual ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    LinSysRes.AddBlock(iPoint, Residual);
    LinSysRes.SubtractBlock(jPoint, Residual);
    
    /*--- Implicit part ---*/
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CScalarSolver::Viscous_Residual(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *numerics,
                                     CConfig *config,
                                     unsigned short iMesh,
                                     unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and normal vector ---*/
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint),
                           solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(jPoint));
    
    /*--- Scalar variables w/o reconstruction. ---*/
    
    numerics->SetScalarVar(nodes->GetSolution(iPoint),
                           nodes->GetSolution(jPoint));
    
    /*--- Scalar variable gradients. ---*/
    
    numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                   nodes->GetGradient(jPoint));
    
    /*--- Mass diffusivity coefficients. ---*/
    
    numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                nodes->GetDiffusivity(jPoint));
    
    /*--- Compute residuals and Jacobians ---*/
    
    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
    
    /*--- Add/subtract residual and update Jacobians ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);
    
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);
    
  }
  
}

void CScalarSolver::ImplicitEuler_Iteration(CGeometry *geometry,
                                            CSolver **solver_container,
                                            CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, *local_Res_TruncError, Vol;
  
  bool scalar_clipping = config->GetScalar_Clipping();
  su2double scalar_clipping_min = config->GetScalar_Clipping_Min();
  su2double scalar_clipping_max = config->GetScalar_Clipping_Max();
  
  bool adjoint = ( config->GetContinuous_Adjoint() ||
                  (config->GetDiscrete_Adjoint() && config->GetFrozen_Visc_Disc()));
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  if (incompressible) SetPreconditioner(geometry, solver_container, config);
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Read the residual ---*/
    
    local_Res_TruncError = nodes->GetResTruncError(iPoint);
    
    /*--- Read the volume ---*/
    
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta = Vol / (config->GetCFLRedCoeff_Scalar()*solver_container[FLOW_SOL]->GetNodes()->GetDelta_Time(iPoint));
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]),
                 geometry->node[iPoint]->GetGlobalIndex(),
                 geometry->node[iPoint]->GetCoord());
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
  
  System.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  if (!adjoint) {
    
    /*--- Update the scalar solution. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (scalar_clipping) {
        nodes->AddClippedSolution(iPoint, 0, config->GetRelaxation_Factor_Scalar()*LinSysSol[iPoint],
                                  scalar_clipping_min, scalar_clipping_max);
      }
      else {
        nodes->AddSolution(iPoint, 0, config->GetRelaxation_Factor_Scalar()*LinSysSol[iPoint]);
      }
    }
    
  }
  
  /*--- Communicate the solution at the periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_IMPLICIT);
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CScalarSolver::BC_Sym_Plane(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {
  
  /*--- Convective fluxes across symmetry plane are equal to zero. ---*/
  
}

void CScalarSolver::BC_Euler_Wall(CGeometry *geometry,
                                  CSolver **solver_container,
                                  CNumerics *numerics,
                                  CConfig *config,
                                  unsigned short val_marker) {
  
  /*--- Convective fluxes across euler wall are equal to zero. ---*/
  
}

void CScalarSolver::BC_HeatFlux_Wall(CGeometry *geometry,
                                     CSolver **solver_container,
                                     CNumerics *conv_numerics,
                                     CNumerics *visc_numerics,
                                     CConfig *config,
                                     unsigned short val_marker) {
  
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}

void CScalarSolver::BC_Isothermal_Wall(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker) {
  
  /*--- Convective fluxes across viscous walls are equal to zero. ---*/
  
}

void CScalarSolver::BC_Far_Field(CGeometry *geometry,
                                 CSolver **solver_container,
                                 CNumerics *conv_numerics,
                                 CNumerics *visc_numerics,
                                 CConfig *config,
                                 unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *Normal, *V_infty, *V_domain;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Set scalar variable at the wall, and at infinity ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Solution_j[iVar] = Scalar_Inf[iVar];
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set Normal (it is necessary to change the sign) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Compute residuals and Jacobians ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Add residuals and Jacobians ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  delete [] Normal;
  
}

void CScalarSolver::BC_Inlet(CGeometry *geometry,
                             CSolver **solver_container,
                             CNumerics *conv_numerics,
                             CNumerics *visc_numerics,
                             CConfig *config,
                             unsigned short val_marker) {
  
  unsigned short iDim, iVar;
  unsigned long iVertex, iPoint;
  su2double *V_inlet, *V_domain, *Normal;
  
  Normal = new su2double[nDim];
  
  bool grid_movement  = config->GetGrid_Movement();
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Set the scalar variable states (prescribed for an inflow) ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Solution_j[iVar] = Inlet_ScalarVars[val_marker][iVertex][iVar];
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set various other quantities in the conv_numerics class ---*/
      
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      
      unsigned long Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_inlet);
      
      /*--- Scalar variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetScalarVar(Solution_i, Solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                          nodes->GetGradient(iPoint));
      
      /*--- Mass diffusivity coefficients. ---*/
      
      visc_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                       nodes->GetDiffusivity(iPoint));
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  delete[] Normal;
  
}

void CScalarSolver::BC_Outlet(CGeometry *geometry,
                              CSolver **solver_container,
                              CNumerics *conv_numerics,
                              CNumerics *visc_numerics,
                              CConfig *config,
                              unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;
  
  bool grid_movement  = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Allocate the value at the outlet ---*/
      
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->GetNodes()->GetPrimitive(iPoint);
      
      /*--- Set various quantities in the solver class ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Set the scalar variables. Here we use a Neumann BC such
       that the scalar variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> ScalarVar_internal,
       Solution_j --> ScalarVar_outlet ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = nodes->GetSolution(iPoint,iVar);
        Solution_j[iVar] = nodes->GetSolution(iPoint,iVar);
      }
      conv_numerics->SetScalarVar(Solution_i, Solution_j);
      
      /*--- Set Normal (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
        Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Jacobian contribution for implicit integration ---*/
      
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
      /*--- Viscous contribution, commented out because serious convergence problems ---*/
      
      unsigned long Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
      visc_numerics->SetNormal(Normal);
      
      /*--- Conservative variables w/o reconstruction ---*/
      
      visc_numerics->SetPrimitive(V_domain, V_outlet);
      
      /*--- Scalar variables w/o reconstruction, and its gradients ---*/
      
      visc_numerics->SetScalarVar(Solution_i, Solution_j);
      visc_numerics->SetScalarVarGradient(nodes->GetGradient(iPoint),
                                          nodes->GetGradient(iPoint));
      
      /*--- Mass diffusivity coefficients. ---*/
      
      visc_numerics->SetDiffusionCoeff(nodes->GetDiffusivity(iPoint),
                                       nodes->GetDiffusivity(iPoint));
      
      /*--- Compute residual, and Jacobians ---*/
      
      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      
      /*--- Subtract residual, and update Jacobians ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual);
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  
}

void CScalarSolver::BC_Periodic(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config) {
  
  /*--- Complete residuals for periodic boundary conditions. We loop over
   the periodic BCs in matching pairs so that, in the event that there are
   adjacent periodic markers, the repeated points will have their residuals
   accumulated corectly during the communications. For implicit calculations
   the Jacobians and linear system are also correctly adjusted here. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_RESIDUAL);
  }
  
}

void CScalarSolver::SetInletAtVertex(su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex) {
  unsigned short iVar;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Inlet_ScalarVars[iMarker][iVertex][iVar] = val_inlet[Inlet_Position+iVar];
  }
  
}

su2double CScalarSolver::GetInletAtVertex(su2double *val_inlet,
                                          unsigned long val_inlet_point,
                                          unsigned short val_kind_marker,
                                          string val_marker,
                                          CGeometry *geometry,
                                          CConfig *config) {
  
  /*--- Local variables ---*/
  
  unsigned short iMarker, iDim, iVar;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};
  
  /*--- Alias positions within inlet file for readability ---*/
  
  if (val_kind_marker == INLET_FLOW) {
    
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          
          if (iPoint == val_inlet_point) {
            
            /*-- Compute boundary face area for this vertex. ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);
            
            /*--- Access and store the inlet variables for this vertex. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
              val_inlet[Inlet_Position+iVar] = Inlet_ScalarVars[iMarker][iVertex][iVar];
            
            /*--- Exit once we find the point. ---*/
            
            return Area;
            
          }
        }
      }
    }
    
  }
  
  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/
  
  return Area;
  
}

void CScalarSolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {
  
  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Inlet_ScalarVars[iMarker][iVertex][iVar] = Scalar_Inf[iVar];
  }
  
}



void CScalarSolver::SetResidual_DualTime(CGeometry *geometry,
                                         CSolver **solver_container,
                                         CConfig *config,
                                         unsigned short iRKStep,
                                         unsigned short iMesh,
                                         unsigned short RunTime_EqSystem) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, jVar, iMarker, iDim;
  unsigned long iPoint, jPoint, iEdge, iVertex;
  
  su2double *U_time_nM1, *U_time_n, *U_time_nP1;
  su2double Volume_nM1, Volume_nP1, TimeStep;
  su2double Density_nM1, Density_n, Density_nP1;
  su2double *Normal = NULL, *GridVel_i = NULL, *GridVel_j = NULL, Residual_GCL;
  
  bool implicit      = (config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Store the physical time step ---*/
  
  TimeStep = config->GetDelta_UnstTimeND();
  
  /*--- Compute the dual time-stepping source term for static meshes ---*/
  
  if (!grid_movement) {
    
    /*--- Loop over all nodes (excluding halos) ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);
      
      /*--- CV volume at time n+1. As we are on a static mesh, the volume
       of the CV will remained fixed for all time steps. ---*/
      
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      /*--- Compute the dual time-stepping source term based on the chosen
       time discretization scheme (1st- or 2nd-order).---*/
      
      /*--- Get the density to compute the conservative variables. ---*/
      
      if (incompressible){
        /*--- This is temporary and only valid for constant-density problems:
         density could also be temperature dependent, but as it is not a part
         of the solution vector it's neither stored for previous time steps
         nor updated with the solution at the end of each iteration. */
        Density_nM1 = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        Density_n   = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        Density_nP1 = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      }
      else{
        Density_nM1 = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint)[0];
        Density_n   = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint)[0];
        Density_nP1 = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint)[0];
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = ( Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*Volume_nP1 / TimeStep;
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = ( 3.0*Density_nP1*U_time_nP1[iVar] - 4.0*Density_n*U_time_n[iVar]
                            +1.0*Density_nM1*U_time_nM1[iVar])*Volume_nP1 / (2.0*TimeStep);
      }
      
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetTime_Marching() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
          if (config->GetTime_Marching() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
    
  } else {
    
    /*--- For unsteady flows on dynamic meshes (rigidly transforming or
     dynamically deforming), the Geometric Conservation Law (GCL) should be
     satisfied in conjunction with the ALE formulation of the governing
     equations. The GCL prevents accuracy issues caused by grid motion, i.e.
     a uniform free-stream should be preserved through a moving grid. First,
     we will loop over the edges and boundaries to compute the GCL component
     of the dual time source term that depends on grid velocities. ---*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Get indices for nodes i & j plus the face normal ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      Normal = geometry->edge[iEdge]->GetNormal();
      
      /*--- Grid velocities stored at nodes i & j ---*/
      
      GridVel_i = geometry->node[iPoint]->GetGridVel();
      GridVel_j = geometry->node[jPoint]->GetGridVel();
      
      /*--- Compute the GCL term by averaging the grid velocities at the
       edge mid-point and dotting with the face normal. ---*/
      
      Residual_GCL = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual_GCL += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
      
      /*--- Compute the GCL component of the source term for node i ---*/
      
      U_time_n = nodes->GetSolution_time_n(iPoint);
      
      /*--- Multiply by density at node i  ---*/
      
      if (incompressible) Density_n = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint); // Temporary fix
      else Density_n = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint)[0];
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Compute the GCL component of the source term for node j ---*/
      
      U_time_n = nodes->GetSolution_time_n(jPoint);
      
      /*--- Multiply by density at node j ---*/
      
      if (incompressible) Density_n = solver_container[FLOW_SOL]->GetNodes()->GetDensity(jPoint); // Temporary fix
      else Density_n = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n(jPoint)[0];
      for (iVar = 0; iVar < nVar; iVar++)
        Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
      
      LinSysRes.SubtractBlock(jPoint, Residual);
      
    }
    
    /*---  Loop over the boundary edges ---*/
    
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
          
          /*--- Get the index for node i plus the boundary face normal ---*/
          
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          
          /*--- Grid velocities stored at boundary node i ---*/
          
          GridVel_i = geometry->node[iPoint]->GetGridVel();
          
          /*--- Compute the GCL term by dotting the grid velocity with the face
           normal. The normal is negated to match the boundary convention. ---*/
          
          Residual_GCL = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            Residual_GCL -= 0.5*(GridVel_i[iDim]+GridVel_i[iDim])*Normal[iDim];
          
          /*--- Compute the GCL component of the source term for node i ---*/
          
          U_time_n = nodes->GetSolution_time_n(iPoint);
          
          /*--- Multiply by density at node i  ---*/
          
          if (incompressible) Density_n = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint); // Temporary fix
          else Density_n = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint)[0];
          for (iVar = 0; iVar < nVar; iVar++)
            Residual[iVar] = Density_n*U_time_n[iVar]*Residual_GCL;
          
          LinSysRes.AddBlock(iPoint, Residual);
        }
    }
    
    /*--- Loop over all nodes (excluding halos) to compute the remainder
     of the dual time-stepping source term. ---*/
    
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Retrieve the solution at time levels n-1, n, and n+1. Note that
       we are currently iterating on U^n+1 and that U^n & U^n-1 are fixed,
       previous solutions that are stored in memory. ---*/
      
      U_time_nM1 = nodes->GetSolution_time_n1(iPoint);
      U_time_n   = nodes->GetSolution_time_n(iPoint);
      U_time_nP1 = nodes->GetSolution(iPoint);
      
      /*--- CV volume at time n-1 and n+1. In the case of dynamically deforming
       grids, the volumes will change. On rigidly transforming grids, the
       volumes will remain constant. ---*/
      
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
      
      if (incompressible){
        /*--- This is temporary and only valid for constant-density problems:
         density could also be temperature dependent, but as it is not a part
         of the solution vector it's neither stored for previous time steps
         nor updated with the solution at the end of each iteration. */
        Density_nM1 = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        Density_n   = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
        Density_nP1 = solver_container[FLOW_SOL]->GetNodes()->GetDensity(iPoint);
      }
      else{
        Density_nM1 = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n1(iPoint)[0];
        Density_n   = solver_container[FLOW_SOL]->GetNodes()->GetSolution_time_n(iPoint)[0];
        Density_nP1 = solver_container[FLOW_SOL]->GetNodes()->GetSolution(iPoint)[0];
      }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        if (config->GetTime_Marching() == DT_STEPPING_1ST)
          Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(Volume_nP1/TimeStep);
        if (config->GetTime_Marching() == DT_STEPPING_2ND)
          Residual[iVar] = (Density_nP1*U_time_nP1[iVar] - Density_n*U_time_n[iVar])*(3.0*Volume_nP1/(2.0*TimeStep))
          + (Density_nM1*U_time_nM1[iVar] - Density_n*U_time_n[iVar])*(Volume_nM1/(2.0*TimeStep));
      }
      
      /*--- Store the residual and compute the Jacobian contribution due
       to the dual time source term. ---*/
      
      LinSysRes.AddBlock(iPoint, Residual);
      if (implicit) {  // TDE density in Jacobian
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) Jacobian_i[iVar][jVar] = 0.0;
          if (config->GetTime_Marching() == DT_STEPPING_1ST)
            Jacobian_i[iVar][iVar] = Volume_nP1/TimeStep;
          if (config->GetTime_Marching() == DT_STEPPING_2ND)
            Jacobian_i[iVar][iVar] = (3.0*Volume_nP1)/(2.0*TimeStep);
        }
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }
  
}

void CScalarSolver::LoadRestart(CGeometry **geometry,
                                CSolver ***solver,
                                CConfig *config,
                                int val_iter,
                                bool val_update_geo) {
  
  /*--- Restart the solution from file information ---*/
  
  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Solution_Fine;
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetTime_Marching() == TIME_STEPPING);
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();
  
  string UnstExt, text_line;
  ifstream restart_file;
  
  string restart_filename = config->GetFilename(config->GetSolution_FileName(), "", val_iter);
  
  bool turbulent = ((config->GetKind_Solver() == RANS) ||
                    (config->GetKind_Solver() == DISC_ADJ_RANS));
  
  unsigned short turbSkip = 0;
  if (turbulent) turbSkip = solver[MESH_0][TURB_SOL]->GetnVar();
  
  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
  
  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }
  
  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  
  /*--- Skip flow variables ---*/
  
  unsigned short skipVars = 0;
  
  if (nDim == 2) skipVars += 6;
  if (nDim == 3) skipVars += 8;
  
  /*--- Skip turbulent variables if necessary variables ---*/
  
  if (turbulent) skipVars += turbSkip;
  
  /*--- Load data from the restart into correct containers. ---*/
  
  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/
    
    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);
    
    if (iPoint_Local > -1) {
      
      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/
      
      index = counter*Restart_Vars[1] + skipVars;
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
      iPoint_Global_Local++;
      
      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
    
  }
  
  /*--- Detect a wrong solution file ---*/
  
  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
  
#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1,
                     MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }
  
  /*--- MPI solution ---*/
  
  solver[MESH_0][SCALAR_SOL]->InitiateComms(geometry[MESH_0], config, SOLUTION);
  solver[MESH_0][SCALAR_SOL]->CompleteComms(geometry[MESH_0], config, SOLUTION);
  
  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][SCALAR_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
  
  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
  
  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][SCALAR_SOL]->GetNodes()->GetSolution(Point_Fine);
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][SCALAR_SOL]->GetNodes()->SetSolution(iPoint,Solution);
    }
    solver[iMesh][SCALAR_SOL]->InitiateComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][SCALAR_SOL]->CompleteComms(geometry[iMesh], config, SOLUTION);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][SCALAR_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }
  
  /*--- Delete the class memory that is used to load the restart. ---*/
  
  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;
  
}
