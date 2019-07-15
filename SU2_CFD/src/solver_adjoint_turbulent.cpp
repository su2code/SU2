/*!
 * \file solution_adjoint_turbulent.cpp
 * \brief Main subrotuines for solving adjoint problems (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, A. Bueno, T. Economon
 * \version 6.2.0 "Falcon"
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/solver_structure.hpp"
#include "../include/variables/CAdjTurbVariable.hpp"

CAdjTurbSolver::CAdjTurbSolver(void) : CSolver() {}

CAdjTurbSolver::CAdjTurbSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  unsigned long iPoint;
  unsigned short iDim, iVar, nLineLets;

  nDim = geometry->GetnDim();
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Dimension of the problem  ---*/
  switch (config->GetKind_Turb_Model()) {
    case SA :     nVar = 1; break;
    case SA_NEG : nVar = 1; break;
    case SST :    nVar = 2; break;
  }
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar+1;
  
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  Residual   = new su2double [nVar]; Residual_RMS = new su2double[nVar];
  Residual_i = new su2double [nVar]; Residual_j = new su2double [nVar];
  Residual_Max = new su2double [nVar]; Point_Max = new unsigned long[nVar];
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  
  Solution   = new su2double [nVar];
  Solution_i = new su2double [nVar];
  Solution_j = new su2double [nVar];
  
  /*--- Define some auxiliar vector related with the geometry ---*/
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Define some auxiliar vector related with the flow solution ---*/
  FlowSolution_i = new su2double [nDim+2]; FlowSolution_j = new su2double [nDim+2];
  
  /*--- Point to point Jacobians ---*/
  Jacobian_ii = new su2double* [nVar];
  Jacobian_ij = new su2double* [nVar];
  Jacobian_ji = new su2double* [nVar];
  Jacobian_jj = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_ii[iVar] = new su2double [nVar];
    Jacobian_ij[iVar] = new su2double [nVar];
    Jacobian_ji[iVar] = new su2double [nVar];
    Jacobian_jj[iVar] = new su2double [nVar];
  }
  
  /*--- Initialization of the structure of the whole Jacobian ---*/
  Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
  
  if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
      (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
    nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
    if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
  }
  
  Jacobian.SetValZero();
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar+1];
    for (iVar = 0; iVar < nVar+1; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Far-Field values and initizalization ---*/
  node = new CVariable* [nPoint];
  bool restart = config->GetRestart();
  
  if (!restart || (iMesh != MESH_0)) {
    PsiNu_Inf = 0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CAdjTurbVariable(PsiNu_Inf, nDim, nVar, config);
    }
  }
  else {
    unsigned long index;
    su2double dull_val;
    string filename, AdjExt, text_line;
    ifstream restart_file;
    
    /*--- Restart the solution from file information ---*/
    string mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);

    restart_file.open(filename.data(), ios::in);
    
    /*--- In case there is no file ---*/
    if (restart_file.fail()) {
      SU2_MPI::Error(string("There is no adjoint restart file ") + filename, CURRENT_FUNCTION);
    }
    
    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0; unsigned long iPoint_Global_Local = 0;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/
    
    getline (restart_file, text_line);
    
    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file, text_line);
      
      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        
        if (nDim == 2) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        if (nDim == 3) point_line >> index >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> dull_val >> Solution[0];
        node[iPoint_Local] = new CAdjTurbVariable(Solution[0], nDim, nVar, config);
        iPoint_Global_Local++;
      }

    }
    
    /*--- Detect a wrong solution file ---*/
    
    if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
    
#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
        SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                       string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }
    
    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CAdjTurbVariable(Solution[0], nDim, nVar, config);
    }
    
    /*--- Close the restart file ---*/
    restart_file.close();
    
  }
  
  /*--- MPI solution ---*/
  
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
}

CAdjTurbSolver::~CAdjTurbSolver(void) {
}

void CAdjTurbSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  
  for (iVertex = 0; iVertex<geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      Solution[0] = 0.0;
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint);
      
    }
  }
  
}

void CAdjTurbSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iPoint, iVertex;
  
  for (iVertex = 0; iVertex<geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {
      
      Solution[0] = 0.0;
      
      /*--- Set the solution values and zero the residual ---*/
      node[iPoint]->SetSolution_Old(Solution);
      LinSysRes.SetBlock_Zero(iPoint);
      
      /*--- Change rows of the Jacobian (includes 1 in the diagonal) ---*/
      Jacobian.DeleteValsRowi(iPoint);
      
    }
  }
  
}

void CAdjTurbSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iPoint, iVertex;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Set Normal ---*/
    conv_numerics->SetNormal(geometry->vertex[val_marker][iVertex]->GetNormal());

    /*--- Set Conservative variables (for convection) ---*/
    su2double* U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
    conv_numerics->SetConservative(U_i, NULL);
    
    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    su2double* TurbPsi_i = node[iPoint]->GetSolution();
    conv_numerics->SetTurbAdjointVar(TurbPsi_i, NULL);
    
    /*--- Add Residuals and Jacobians ---*/
    conv_numerics->ComputeResidual(Residual, Jacobian_ii, NULL, config);
    LinSysRes.AddBlock(iPoint, Residual);
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
    
  }
  
}

void CAdjTurbSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  unsigned long iPoint;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
  /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  
    /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();
  
  /*--- Gradient of the adjoint turbulent variables ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Gradient of the turbulent variables ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[TURB_SOL]->SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[TURB_SOL]->SetSolution_Gradient_LS(geometry, config);
  
}

void CAdjTurbSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
  unsigned long iEdge, iPoint, jPoint;
  su2double *U_i, *U_j, *TurbPsi_i, *TurbPsi_j, **TurbVar_Grad_i, **TurbVar_Grad_j;
//  su2double *Limiter_i = NULL, *Limiter_j = NULL, **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j;
//  unsigned short iDim, iVar;
  
  bool muscl   = config->GetMUSCL_AdjTurb();
  bool limiter = (config->GetKind_SlopeLimit_AdjTurb() != NO_LIMITER);
  
  if (muscl) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    if (limiter) SetSolution_Limiter(geometry, config);
  }
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
    U_j = solver_container[FLOW_SOL]->node[jPoint]->GetSolution();
    numerics->SetConservative(U_i, U_j);
    
    /*--- Set normal vectors and length ---*/
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    TurbPsi_i = node[iPoint]->GetSolution();
    TurbPsi_j = node[jPoint]->GetSolution();
    numerics->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);
    
    /*--- Gradient of turbulent variables w/o reconstruction ---*/
    TurbVar_Grad_i = solver_container[TURB_SOL]->node[iPoint]->GetGradient();
    TurbVar_Grad_j = solver_container[TURB_SOL]->node[jPoint]->GetGradient();
    numerics->SetTurbVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);
    
//    if (muscl) {
//      
//      /*--- Conservative solution using gradient reconstruction ---*/
//      for (iDim = 0; iDim < nDim; iDim++) {
//        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
//        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
//      }
//      Gradient_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient();
//      Gradient_j = solver_container[FLOW_SOL]->node[jPoint]->GetGradient();
//      for (iVar = 0; iVar < solver_container[FLOW_SOL]->GetnVar(); iVar++) {
//        Project_Grad_i = 0; Project_Grad_j = 0;
//        for (iDim = 0; iDim < nDim; iDim++) {
//          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
//          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
//        }
//        FlowSolution_i[iVar] = U_i[iVar] + Project_Grad_i;
//        FlowSolution_j[iVar] = U_j[iVar] + Project_Grad_j;
//      }
//      numerics->SetConservative(FlowSolution_i, FlowSolution_j);
//      
//      /*--- Adjoint turbulent variables using gradient reconstruction ---*/
//      Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
//      if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
//      for (iVar = 0; iVar < nVar; iVar++) {
//        Project_Grad_i = 0; Project_Grad_j = 0;
//        for (iDim = 0; iDim < nDim; iDim++) {
//          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim];
//          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim];
//        }
//        if (limiter) {
//          Solution_i[iVar] = TurbPsi_i[iVar] + Project_Grad_i*Limiter_i[iVar];
//          Solution_j[iVar] = TurbPsi_j[iVar] + Project_Grad_j*Limiter_j[iVar];
//        }
//        else {
//          Solution_i[iVar] = TurbPsi_i[iVar] + Project_Grad_i;
//          Solution_j[iVar] = TurbPsi_j[iVar] + Project_Grad_j;
//        }
//      }
//      numerics->SetTurbVar(Solution_i, Solution_j);
//      
//    }
    
    /*--- Set normal vectors and length ---*/
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Add and Subtract Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
    
  }
  
}

void CAdjTurbSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config,
                                      unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;
  su2double *Coord_i, *Coord_j;
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge ---*/
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    /*--- Points coordinates, and set normal vectors and length ---*/
    Coord_i = geometry->node[iPoint]->GetCoord();
    Coord_j = geometry->node[jPoint]->GetCoord();
    numerics->SetCoord(Coord_i, Coord_j);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Conservative variables w/o reconstruction, turbulent variables w/o reconstruction,
     and turbulent adjoint variables w/o reconstruction ---*/
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), solver_container[FLOW_SOL]->node[jPoint]->GetSolution());
    numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), solver_container[TURB_SOL]->node[jPoint]->GetSolution());
    numerics->SetTurbAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
    /*--- Viscosity ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                  solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());
    
    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    numerics->SetTurbAdjointGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Compute residual in a non-conservative way, and update ---*/
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Update adjoint viscous residual ---*/
    LinSysRes.AddBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
    Jacobian.AddBlock(iPoint, jPoint, Jacobian_ij);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
    
  }
  
}

void CAdjTurbSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics, CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;
  su2double *U_i, **GradPrimVar_i, *TurbVar_i;
  su2double **TurbVar_Grad_i, *TurbPsi_i, **PsiVar_Grad_i; // Gradients
  
  /*--- Piecewise source term ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
        
    /*--- Conservative variables w/o reconstruction ---*/
    U_i = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
    numerics->SetConservative(U_i, NULL);
    
    /*--- Gradient of primitive variables w/o reconstruction ---*/
    GradPrimVar_i = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
    numerics->SetPrimVarGradient(GradPrimVar_i, NULL);
    
    /*--- Laminar viscosity of the fluid ---*/
    numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(), 0.0);
    
    /*--- Turbulent variables w/o reconstruction ---*/
    TurbVar_i = solver_container[TURB_SOL]->node[iPoint]->GetSolution();
    numerics->SetTurbVar(TurbVar_i, NULL);
    
    /*--- Gradient of Turbulent Variables w/o reconstruction ---*/
    TurbVar_Grad_i = solver_container[TURB_SOL]->node[iPoint]->GetGradient();
    numerics->SetTurbVarGradient(TurbVar_Grad_i, NULL);
    
    /*--- Turbulent adjoint variables w/o reconstruction ---*/
    TurbPsi_i = node[iPoint]->GetSolution();
    numerics->SetTurbAdjointVar(TurbPsi_i, NULL);
    
    /*--- Gradient of Adjoint flow variables w/o reconstruction
     (for non-conservative terms depending on gradients of flow adjoint vars.) ---*/
    PsiVar_Grad_i = solver_container[ADJFLOW_SOL]->node[iPoint]->GetGradient();
    numerics->SetAdjointVarGradient(PsiVar_Grad_i, NULL);

    /*--- Set volume and distances to the surface ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
    
    /*--- Add and Subtract Residual ---*/
    numerics->ComputeResidual(Residual, Jacobian_ii, NULL, config);
    LinSysRes.AddBlock(iPoint, Residual);
    Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
    
  }
  
//  /*--- Conservative Source Term ---*/
//  su2double **TurbVar_Grad_j;
//  unsigned long jPoint, iEdge;
//
//  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
//    
//    /*--- Points in edge ---*/
//    iPoint = geometry->edge[iEdge]->GetNode(0);
//    jPoint = geometry->edge[iEdge]->GetNode(1);
//    
//    /*--- Gradient of turbulent variables w/o reconstruction ---*/
//    TurbVar_Grad_i = solver_container[TURB_SOL]->node[iPoint]->GetGradient();
//    TurbVar_Grad_j = solver_container[TURB_SOL]->node[jPoint]->GetGradient();
//    second_numerics->SetTurbVarGradient(TurbVar_Grad_i, TurbVar_Grad_j);
//    
//    /*--- Turbulent adjoint variables w/o reconstruction ---*/
//    TurbPsi_i = node[iPoint]->GetSolution();
//    TurbPsi_j = node[jPoint]->GetSolution();
//    second_numerics->SetTurbAdjointVar(TurbPsi_i, TurbPsi_j);
//    
//    /*--- Set normal vectors and length ---*/
//    second_numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
//    
//    /*--- Add and Subtract Residual ---*/
//    second_numerics->ComputeResidual(Residual, Jacobian_ii, Jacobian_jj, config);
//    LinSysRes.AddBlock(iPoint, Residual);
//    LinSysRes.SubtractBlock(jPoint, Residual);
//    Jacobian.AddBlock(iPoint, iPoint, Jacobian_ii);
//    Jacobian.AddBlock(iPoint, jPoint, Jacobian_jj);
//    Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ii);
//    Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
//    
//  }
  
}

void CAdjTurbSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol;
  
  /*--- Set maximum residual to zero ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Build implicit system ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    /*--- Read the volume ---*/
    
    Vol = geometry->node[iPoint]->GetVolume();
    
    /*--- Modify matrix diagonal to assure diagonal dominance ---*/
    
    Delta = Vol / (config->GetCFLRedCoeff_AdjTurb()*solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());
    
    Jacobian.AddVal2Diag(iPoint, Delta);
    
    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = -LinSysRes[total_index];
      LinSysSol[total_index] = 0.0;
      AddRes_RMS(iVar, LinSysRes[total_index]*LinSysRes[total_index]);
      AddRes_Max(iVar, fabs(LinSysRes[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
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
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++)
      node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
  }
  
  /*--- MPI solution ---*/
    
  InitiateComms(geometry, config, SOLUTION);
  CompleteComms(geometry, config, SOLUTION);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}
