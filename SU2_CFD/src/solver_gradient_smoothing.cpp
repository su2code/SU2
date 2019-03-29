/*!
 * \file solver_gradient_smoothing.cpp
 * \brief Main subroutines for the gradient smoothing problem.
 * \author T. Dick
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
#include <algorithm>

CGradientSmoothingSolver::CGradientSmoothingSolver(void) : CSolver() {

  nElement = 0;
  nDim = 0;
  nMarker = 0;

  GradN_X = NULL;
  GradN_x = NULL;

  nPoint = 0;
  nPointDomain = 0;

  node = NULL;

  element_container = NULL;
  element_properties = NULL;

  MassMatrix_ij = NULL;

  SolRest = NULL;

}

CGradientSmoothingSolver::CGradientSmoothingSolver(CGeometry *geometry, CConfig *config) : CSolver() {

  unsigned long iPoint;
  unsigned short iVar, jVar, iDim, jDim;
  unsigned short iTerm, iKind;

  nElement      = geometry->GetnElem();
  nDim          = geometry->GetnDim();
  nMarker       = geometry->GetnMarker();

  nPoint        = geometry->GetnPoint();
  nPointDomain  = geometry->GetnPointDomain();


  /*--- Here is where we assign the kind of each element ---*/

  /*--- First level: different possible terms of the equations ---*/
  element_container = new CElement** [MAX_TERMS];
  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++)
    element_container[iTerm] = new CElement* [MAX_FE_KINDS];

  for (iTerm = 0; iTerm < MAX_TERMS; iTerm++) {
    for (iKind = 0; iKind < MAX_FE_KINDS; iKind++) {
      element_container[iTerm][iKind] = NULL;
    }
  }

  if (nDim == 2) {

      /*--- Basic terms ---*/
      element_container[FEA_TERM][EL_TRIA] = new CTRIA1(nDim, config);
      element_container[FEA_TERM][EL_QUAD] = new CQUAD4(nDim, config);

      if (de_effects){
        element_container[DE_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[DE_TERM][EL_QUAD] = new CQUAD4(nDim, config);
      }

      if (incompressible){
        element_container[INC_TERM][EL_TRIA] = new CTRIA1(nDim, config);
        element_container[INC_TERM][EL_QUAD] = new CQUAD1(nDim, config);
      }

  }
  else if (nDim == 3) {

      element_container[FEA_TERM][EL_TETRA] = new CTETRA1(nDim, config);
      element_container[FEA_TERM][EL_HEXA] = new CHEXA8(nDim, config);

      if (de_effects){
        element_container[DE_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[DE_TERM][EL_HEXA] = new CHEXA8(nDim, config);
      }

      if (incompressible) {
        element_container[INC_TERM][EL_TETRA] = new CTETRA1(nDim, config);
        element_container[INC_TERM][EL_HEXA] = new CHEXA1(nDim, config);
      }


  }

  node = new CVariable*[nPoint];

  /*--- Set element properties ---*/
  elProperties = new unsigned long[4];
  for (iVar = 0; iVar < 4; iVar++)
    elProperties[iVar] = 0;
  Set_ElementProperties(geometry, config);

  GradN_X = new su2double [nDim];
  GradN_x = new su2double [nDim];

  nVar = nDim;

  /*--- The length of the solution vector depends on whether the problem is static or dynamic ---*/

  unsigned short nSolVar;

  nSolVar = nVar;

  su2double* SolRest = new su2double[nSolVar];

  /*--- Initialize from zero everywhere. ---*/

  for (iVar = 0; iVar < nSolVar; iVar++) SolRest[iVar] = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new CGradientSmoothingVariable(SolRest, nDim, nVar, config);
  }


  /*--- Initialization of matrix structures ---*/
  if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (Non-Linear Elasticity)." << endl;

  StiffnessMatrix.Initialize(nPoint, nPointDomain, nVar, nVar, false, geometry, config);

  /*--- Initialization of linear solver structures ---*/
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);


  /*--- Perform the MPI communication of the solution ---*/

  Set_MPI_Solution(geometry, config);

}

CGradientSmoothingSolver::~CGradientSmoothingSolver(void) {

  unsigned short iVar, jVar;
  unsigned long iElem;

  if (element_container != NULL) {
    for (iVar = 0; iVar < MAX_TERMS; iVar++) {
      for (jVar = 0; jVar < MAX_FE_KINDS; jVar++) {
        if (element_container[iVar][jVar] != NULL) delete element_container[iVar][jVar];
      }
      delete [] element_container[iVar];
    }
    delete [] element_container;
  }

  if (element_properties != NULL){
    for (iElem = 0; iElem < nElement; iElem++)
      if (element_properties[iElem] != NULL) delete element_properties[iElem];
    delete [] element_properties;
  }

  delete [] GradN_X;
  delete [] GradN_x;

}

void CGradientSmoothingSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {


  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

  bool dynamic = (config->GetDynamic_Analysis() == DYNAMIC);              // Dynamic simulations.

  unsigned short nSolVar;

  if (dynamic) nSolVar = 3 * nVar;
  else nSolVar = nVar;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nSolVar;     nBufferR_Vector = nVertexR*nSolVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution that should be sent ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
        if (dynamic) {
          for (iVar = 0; iVar < nVar; iVar++) {
            Buffer_Send_U[(iVar+nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Vel(iVar);
            Buffer_Send_U[(iVar+2*nVar)*nVertexS+iVertex] = node[iPoint]->GetSolution_Accel(iVar);
          }
        }
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                        Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);

#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        if (dynamic) {
          for (iVar = nVar; iVar < 3*nVar; iVar++)
            Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
        }
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy solution variables. ---*/
        for (iVar = 0; iVar < nSolVar; iVar++)
          SolRest[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];

        /*--- Store received values back into the variable. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, SolRest[iVar]);

        if (dynamic) {

          for (iVar = 0; iVar < nVar; iVar++) {
            node[iPoint]->SetSolution_Vel(iVar, SolRest[iVar+nVar]);
            node[iPoint]->SetSolution_Accel(iVar, SolRest[iVar+2*nVar]);
          }

        }

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }

}

void CGradientSmoothingSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, CNumerics **numerics, unsigned short iMesh, unsigned long Iteration, unsigned short RunTime_EqSystem, bool Output) {


  unsigned long iPoint;

  /*--- Set vector entries to zero ---*/
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
    LinSysSol.SetBlock_Zero(iPoint);
  }

  /*--- Set matrix entries to zero ---*/

  StiffnessMatrix.SetValZero();

}

void CGradientSmoothingSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {

  unsigned long iPoint, nPoint;
  bool incremental_load = config->GetIncrementalLoad();              // If an incremental load is applied

  nPoint = geometry[MESH_0]->GetnPoint();

  /*--- We store the current solution as "Solution Old", for the case that we need to retrieve it ---*/

  if (incremental_load) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) node[iPoint]->Set_OldSolution();
  }


}

void CGradientSmoothingSolver::Compute_StiffMatrix(CGeometry *geometry, CSolver **solver_container, CNumerics **numerics, CConfig *config) {

  unsigned long iElem, iVar, jVar;
  unsigned short iNode, iDim, nNodes = 0;
  unsigned long indexNode[8]={0,0,0,0,0,0,0,0};
  su2double val_Coord, val_Sol;
  int EL_KIND = 0;

  su2double *Kab = NULL, *Ta  = NULL;
  unsigned short NelNodes, jNode;

  /*--- Loops over all the elements ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {

    if (geometry->elem[iElem]->GetVTK_Type() == TRIANGLE)      {nNodes = 3; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == QUADRILATERAL) {nNodes = 4; EL_KIND = EL_QUAD;}
    if (geometry->elem[iElem]->GetVTK_Type() == TETRAHEDRON)   {nNodes = 4; EL_KIND = EL_TETRA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PYRAMID)       {nNodes = 5; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == PRISM)         {nNodes = 6; EL_KIND = EL_TRIA;}
    if (geometry->elem[iElem]->GetVTK_Type() == HEXAHEDRON)    {nNodes = 8; EL_KIND = EL_HEXA;}

    /*--- For the number of nodes, we get the coordinates from the connectivity matrix ---*/

    for (iNode = 0; iNode < nNodes; iNode++) {

      indexNode[iNode] = geometry->elem[iElem]->GetNode(iNode);

      for (iDim = 0; iDim < nDim; iDim++) {
        val_Coord = geometry->node[indexNode[iNode]]->GetCoord(iDim);
        val_Sol = node[indexNode[iNode]]->GetSolution(iDim) + val_Coord;
        element_container[FEA_TERM][EL_KIND]->SetRef_Coord(val_Coord, iNode, iDim);
        element_container[FEA_TERM][EL_KIND]->SetCurr_Coord(val_Sol, iNode, iDim);
      }
    }

    /*--- Set the properties of the element ---*/
    element_container[FEA_TERM][EL_KIND]->Set_ElProperties(element_properties[iElem]);

    /*--- Compute the components of the jacobian and the stress term ---*/
    numerics[FEA_TERM]->Compute_Tangent_Matrix(element_container[FEA_TERM][EL_KIND], config);


    NelNodes = element_container[FEA_TERM][EL_KIND]->GetnNodes();

    for (iNode = 0; iNode < NelNodes; iNode++) {

      Ta = element_container[FEA_TERM][EL_KIND]->Get_Kt_a(iNode);
      for (iVar = 0; iVar < nVar; iVar++) Res_Stress_i[iVar] = Ta[iVar];

      LinSysRes.SubtractBlock(indexNode[iNode], Res_Stress_i);

      for (jNode = 0; jNode < NelNodes; jNode++) {

        Kab = element_container[FEA_TERM][EL_KIND]->Get_Kab(iNode, jNode);

        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            Jacobian_ij[iVar][jVar] = Kab[iVar*nVar+jVar];
          }
        }

        StiffnessMatrix.AddBlock(indexNode[iNode], indexNode[jNode], Jacobian_ij);
      }

    }

  }


}

void CGradientSmoothingSolver::Initialize_SystemMatrix(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

}

void CGradientSmoothingSolver::Compute_IntegrationConstants(CConfig *config) {

  su2double Delta_t= config->GetDelta_DynTime();

  su2double gamma = config->GetNewmark_gamma(), beta = config->GetNewmark_beta();

  switch (config->GetKind_TimeIntScheme_FEA()) {
    case (CD_EXPLICIT):
      cout << "NOT IMPLEMENTED YET" << endl;
      break;
    case (NEWMARK_IMPLICIT):

      /*--- Integration constants for Newmark scheme ---*/

      a_dt[0]= 1 / (beta*pow(Delta_t,2.0));
      a_dt[1]= gamma / (beta*Delta_t);
      a_dt[2]= 1 / (beta*Delta_t);
      a_dt[3]= 1 /(2*beta) - 1;
      a_dt[4]= gamma/beta - 1;
      a_dt[5]= (Delta_t/2) * (gamma/beta - 2);
      a_dt[6]= Delta_t * (1-gamma);
      a_dt[7]= gamma * Delta_t;
      a_dt[8]= 0.0;

      break;

    case (GENERALIZED_ALPHA):

      /*--- Integration constants for Generalized Alpha ---*/
      /*--- Needs to be updated if accounting for structural damping ---*/

      //      su2double beta = config->Get_Int_Coeffs(0);
      //      //  su2double gamma =  config->Get_Int_Coeffs(1);
      //      su2double alpha_f = config->Get_Int_Coeffs(2), alpha_m =  config->Get_Int_Coeffs(3);
      //
      //      a_dt[0]= (1 / (beta*pow(Delta_t,2.0))) * ((1 - alpha_m) / (1 - alpha_f)) ;
      //      a_dt[1]= 0.0 ;
      //      a_dt[2]= (1 - alpha_m) / (beta*Delta_t);
      //      a_dt[3]= ((1 - 2*beta)*(1-alpha_m) / (2*beta)) - alpha_m;
      //      a_dt[4]= 0.0;
      //      a_dt[5]= 0.0;
      //      a_dt[6]= Delta_t * (1-delta);
      //      a_dt[7]= delta * Delta_t;
      //      a_dt[8]= (1 - alpha_m) / (beta*pow(Delta_t,2.0));

      break;
  }


}

void CGradientSmoothingSolver::BC_Impose(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {}

void CGradientSmoothingSolver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,  CNumerics **numerics,
                                           unsigned short iMesh) {

  unsigned short iVar;
  unsigned long iPoint, total_index;

  bool first_iter = (config->GetIntIter() == 0);
  bool nonlinear_analysis = (config->GetGeometricConditions() == LARGE_DEFORMATIONS);    // Nonlinear analysis.
  bool disc_adj_fem = (config->GetKind_Solver() == DISC_ADJ_FEM);

  su2double solNorm = 0.0, solNorm_recv = 0.0;

  if (disc_adj_fem) {

    if (nonlinear_analysis) {

      /*--- For nonlinear discrete adjoint, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)): ABSOLUTE, norm of the incremental displacements ------------*/
      /*--- RTOL = norm(Residual(k): ABSOLUTE, norm of the residual (T-F) -----------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k): ABSOLUTE, norm of the product disp * res -----------*/

      Conv_Check[0] = LinSysSol.norm();               // Norm of the delta-solution vector
      Conv_Check[1] = LinSysRes.norm();               // Norm of the residual
      Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);
    }
    else {
      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/
      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);

    }

  }
  else {
    if (nonlinear_analysis){

      /*--- If the problem is nonlinear, we have 3 convergence criteria ---*/

      /*--- UTOL = norm(Delta_U(k)) / norm(U(k)) --------------------------*/
      /*--- RTOL = norm(Residual(k)) / norm(Residual(0)) ------------------*/
      /*--- ETOL = Delta_U(k) * Residual(k) / Delta_U(0) * Residual(0) ----*/

      if (first_iter){
        Conv_Ref[0] = 1.0;                                        // Position for the norm of the solution
        Conv_Ref[1] = max(LinSysRes.norm(), EPS);                 // Position for the norm of the residual
        Conv_Ref[2] = max(dotProd(LinSysSol, LinSysRes), EPS);    // Position for the energy tolerance

        /*--- Make sure the computation runs at least 2 iterations ---*/
        Conv_Check[0] = 1.0;
        Conv_Check[1] = 1.0;
        Conv_Check[2] = 1.0;

        /*--- If absolute, we check the norms ---*/
        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
        }
      }
      else {
        /*--- Compute the norm of the solution vector Uk ---*/
        for (iPoint = 0; iPoint < nPointDomain; iPoint++){
          for (iVar = 0; iVar < nVar; iVar++){
            solNorm += node[iPoint]->GetSolution(iVar) * node[iPoint]->GetSolution(iVar);
          }
        }

        // We need to communicate the norm of the solution and compute the RMS throughout the different processors

#ifdef HAVE_MPI
        /*--- We sum the squares of the norms across the different processors ---*/
        SU2_MPI::Allreduce(&solNorm, &solNorm_recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        solNorm_recv         = solNorm;
#endif

        Conv_Ref[0] = max(sqrt(solNorm_recv), EPS);           // Norm of the solution vector

        switch (config->GetResidual_Criteria_FEM()) {
          case RESFEM_RELATIVE:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
          case RESFEM_ABSOLUTE:
            Conv_Check[0] = LinSysSol.norm();         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm();         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes);  // Position for the energy tolerance
            break;
          default:
            Conv_Check[0] = LinSysSol.norm() / Conv_Ref[0];         // Norm of the delta-solution vector
            Conv_Check[1] = LinSysRes.norm() / Conv_Ref[1];         // Norm of the residual
            Conv_Check[2] = dotProd(LinSysSol, LinSysRes) / Conv_Ref[2];  // Position for the energy tolerance
            break;
        }

      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

    } else {

      /*--- If the problem is linear, the only check we do is the RMS of the displacements ---*/

      /*---  Compute the residual Ax-f ---*/

      Jacobian.ComputeResidual(LinSysSol, LinSysRes, LinSysAux);

      /*--- Set maximum residual to zero ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        SetRes_RMS(iVar, 0.0);
        SetRes_Max(iVar, 0.0, 0);
      }

      /*--- Compute the residual ---*/

      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          AddRes_RMS(iVar, LinSysAux[total_index]*LinSysAux[total_index]);
          AddRes_Max(iVar, fabs(LinSysAux[total_index]), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
        }
      }

      /*--- MPI solution ---*/

      Set_MPI_Solution(geometry, config);

      /*--- Compute the root mean square residual ---*/

      SetResidual_RMS(geometry, config);
    }

  }

}

void CGradientSmoothingSolver::Solve_System(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned long IterLinSol = 0, iPoint, total_index;
  unsigned short iVar;

  /*--- Initialize residual and solution at the ghost points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar + iVar;
      LinSysRes[total_index] = 0.0;
      LinSysSol[total_index] = 0.0;
    }

  }

  CSysSolve LinSystem;
  IterLinSol = LinSystem.Solve(StiffnessMatrix, LinSysRes, LinSysSol, geometry, config);

  /*--- The the number of iterations of the linear solver ---*/

  SetIterLinSolver(IterLinSol);

}

su2double CGradientSmoothingSolver::GetRes(unsigned short val_var) {}

void CGradientSmoothingSolver::Set_ElementProperties(CGeometry *geometry, CConfig *config) {

  unsigned long iElem;
  unsigned long index;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = geometry->GetnZone();

  bool topology_mode = config->GetTopology_Optimization();

  string filename;
  ifstream properties_file;

  element_properties = new CElementProperty*[nElement];

  /*--- Restart the solution from file information ---*/

  filename = config->GetFEA_FileName();

  /*--- If multizone, append zone name ---*/
  if (nZone > 1)
    filename = config->GetMultizone_FileName(filename, iZone);

  if (rank == MASTER_NODE) cout << "Filename: " << filename << "." << endl;

  properties_file.open(filename.data(), ios::in);

  /*--- In case there is no file, all elements get the same property (0) ---*/

  if (properties_file.fail()) {
    if (rank == MASTER_NODE){
      cout << "There is no element-based properties file." << endl;
      cout << "The structural domain has uniform properties." << endl;

      if (topology_mode)
        SU2_MPI::Error("Topology mode requires an element-based properties file.",CURRENT_FUNCTION);
    }

    for (iElem = 0; iElem < nElement; iElem++){
      element_properties[iElem] = new CElementProperty(0, 0, 0, 0);
    }

    element_based = false;

  }
  else{

    element_based = true;

    /*--- In case this is a parallel simulation, we need to perform the
       Global2Local index transformation first. ---*/

    long *Global2Local = new long[geometry->GetGlobal_nElemDomain()];

    /*--- First, set all indices to a negative value by default ---*/

    for (iElem = 0; iElem < geometry->GetGlobal_nElemDomain(); iElem++)
      Global2Local[iElem] = -1;

    /*--- Now fill array with the transform values only for the points in the rank (including halos) ---*/

    for (iElem = 0; iElem < nElement; iElem++)
      Global2Local[geometry->elem[iElem]->GetGlobalIndex()] = iElem;

    /*--- Read all lines in the restart file ---*/

    long iElem_Local;
    unsigned long iElem_Global_Local = 0, iElem_Global = 0; string text_line;
    unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

    /*--- The first line is the header ---*/

    getline (properties_file, text_line);

    for (iElem_Global = 0; iElem_Global < geometry->GetGlobal_nElemDomain(); iElem_Global++ ) {

      getline (properties_file, text_line);

      istringstream point_line(text_line);

      /*--- Retrieve local index. If this element from the restart file lives
         only on a different processor, the value of iPoint_Local will be -1.
         Otherwise, the local index for this node on the current processor
         will be returned and used to instantiate the vars. ---*/

      iElem_Local = Global2Local[iElem_Global];

      if (iElem_Local >= 0) {

        if (config->GetDE_Effects())
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];
        else
          point_line >> index >> elProperties[0] >> elProperties[1] >> elProperties[2] >> elProperties[3];

        element_properties[iElem_Local] = new CElementProperty(elProperties[0],
                                                         elProperties[1],
                                                         elProperties[2],
                                                         elProperties[3]);

        /*--- For backwards compatibility we only read a fifth column in topology mode ---*/
        if (topology_mode) {
          su2double elDensity;
          point_line >> elDensity;
          element_properties[iElem_Local]->SetDesignDensity(elDensity);
        }

        iElem_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    if (iElem_Global_Local < nElement) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The properties file ") + filename + string(" doesn't match with the mesh file!\n")  +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    properties_file.close();

    /*--- Free memory needed for the transformation ---*/

    delete [] Global2Local;


  }


}

void CGradientSmoothingSolver::ExtractAdjoint_Variables(CGeometry *geometry, CSolver **solver_container, CConfig *config)
{

    /* possible MPI reduce?? */

    /* get the adjoint solution from solver_container */
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++){
        node[iPoint]->SetSensitivity(iDim,solver_container[ADJFLOW_SOL]->node[iPoint]->GetSensitivity(iDim);
    }
}


