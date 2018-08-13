/*!
 * \file solution_direct_radiation.cpp
 * \brief Main subrotuines for solving radiation problems (P1, M1, Discrete Ordinates).
 * \author R. Sanchez
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

#include "../include/solver_structure.hpp"

CRadSolver::CRadSolver(void) : CSolver() {

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

}

CRadSolver::CRadSolver(CGeometry* geometry, CConfig *config) : CSolver() {

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

  Absorption_Coeff = config->GetAbsorption_Coeff();
  Scattering_Coeff = config->GetScattering_Coeff();

  Absorption_Coeff = max(Absorption_Coeff,0.01);

}

CRadSolver::~CRadSolver(void) {

  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;

}

void CRadSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector, nBufferS_Scalar, nBufferR_Scalar;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL, *Buffer_Receive_muT = NULL, *Buffer_Send_muT = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      nBufferS_Scalar = nVertexS;             nBufferR_Scalar = nVertexR;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      Buffer_Receive_muT = new su2double [nBufferR_Scalar];
      Buffer_Send_muT = new su2double[nBufferS_Scalar];

      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_U, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_U, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      SU2_MPI::Sendrecv(Buffer_Send_muT, nBufferS_Scalar, MPI_DOUBLE, send_to, 1,
                   Buffer_Receive_muT, nBufferR_Scalar, MPI_DOUBLE, receive_from, 1, MPI_COMM_WORLD, &status);
#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        Buffer_Receive_muT[iVertex] = node[iPoint]->GetmuT();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_U[iVar*nVertexR+iVertex] = Buffer_Send_U[iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;
      delete [] Buffer_Send_muT;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy conservative variables. ---*/
        node[iPoint]->SetmuT(Buffer_Receive_muT[iVertex]);
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_muT;
      delete [] Buffer_Receive_U;

    }

  }

}

void CRadSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution_Old(iVar);
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
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_U;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();

        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Buffer_Receive_U[iVar*nVertexR+iVertex]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;

    }

  }
}

void CRadSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;

  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];

#ifdef HAVE_MPI
  int send_to, receive_from;
  SU2_MPI::Status status;
#endif

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
        (config->GetMarker_All_SendRecv(iMarker) > 0)) {

      MarkerS = iMarker;  MarkerR = iMarker+1;

#ifdef HAVE_MPI
      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
#endif

      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
      nBufferS_Vector = nVertexS*nVar*nDim;        nBufferR_Vector = nVertexR*nVar*nDim;

      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
      Buffer_Send_Gradient = new su2double[nBufferS_Vector];

      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Gradient[iDim*nVar*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient(iVar, iDim);
      }

#ifdef HAVE_MPI

      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else

      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];
      }

#endif

      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Gradient;

      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {

        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();

        /*--- Retrieve the supplied periodic information. ---*/
        angles = config->GetPeriodicRotation(iPeriodic_Index);

        /*--- Store angles separately for clarity. ---*/
        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis,
         then z-axis. Note that this is the transpose of the matrix
         used during the preprocessing stage. ---*/
        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Copy conserved variables before performing transformation. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nVar*nVertexR+iVar*nVertexR+iVertex];

        /*--- Need to rotate the gradients for all conserved variables. ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          if (nDim == 2) {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
          else {
            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nVar*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nVar*nVertexR+iVar*nVertexR+iVertex];
          }
        }

        /*--- Store the received information ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            node[iPoint]->SetGradient(iVar, iDim, Gradient[iVar][iDim]);

      }

      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Gradient;

    }

  }

  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Gradient[iVar];
  delete [] Gradient;

}

void CRadSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, Vol, *local_Res_TruncError;
  bool flow = ((config->GetKind_Solver() == NAVIER_STOKES)
               || (config->GetKind_Solver() == RANS)
               || (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)
               || (config->GetKind_Solver() == DISC_ADJ_RANS));


  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }

  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the residual ---*/

    local_Res_TruncError = node[iPoint]->GetResTruncError();

    /*--- Read the volume ---*/

    Vol = geometry->node[iPoint]->GetVolume();

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    if (node[iPoint]->GetDelta_Time() != 0.0) {

      if(flow) {
        Delta = Vol / node[iPoint]->GetDelta_Time();
        Jacobian.AddVal2Diag(iPoint, Delta);
      }
      else {
        Delta = Vol / node[iPoint]->GetDelta_Time();
        Jacobian.AddVal2Diag(iPoint, Delta);
      }

    } else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
        local_Res_TruncError[iVar] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index] + local_Res_TruncError[iVar]);
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

  CSysSolve system;
  system.Solve(Jacobian, LinSysRes, LinSysSol, geometry, config);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddSolution(iVar, LinSysSol[iPoint*nVar+iVar]);
    }
  }

  /*--- MPI solution ---*/

  Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);
}

void CRadSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  unsigned short iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Solution_Fine;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string UnstExt, text_line;
  ifstream restart_file;
  string restart_filename = config->GetSolution_FlowFileName();

  /*--- Modify file name for multizone problems ---*/
  if (nZone >1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time|| time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

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
      node[iPoint_Local]->SetSolution(Solution);
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
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- MPI solution and compute the eddy viscosity ---*/

//TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart.
  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);

  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
    solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
  }

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

CRadP1Solver::CRadP1Solver(void) : CRadSolver() {

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

}

CRadP1Solver::CRadP1Solver(CGeometry* geometry, CConfig *config) : CRadSolver() {

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

  Absorption_Coeff = config->GetAbsorption_Coeff();
  Scattering_Coeff = config->GetScattering_Coeff();

}

CRadP1Solver::~CRadP1Solver(void) {

  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;

}

void CRadP1Solver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;

  /*--- Initialize the residual vector ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    LinSysRes.SetBlock_Zero(iPoint);
  }

  /*--- Initialize the Jacobian matrix ---*/
  Jacobian.SetValZero();

  /*--- Compute the Solution gradients ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);

}

void CRadP1Solver::Postprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh) {

}

void CRadP1Solver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iEdge, iPoint, jPoint;

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    /*--- Points in edge ---*/

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Points coordinates, and normal vector ---*/

    numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
                       geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

    /*--- Radiation variables w/o reconstruction, and its gradients ---*/

    numerics->SetRadVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    numerics->SetRadVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());

    /*--- Compute residual, and Jacobians ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

    /*--- Add and subtract residual, and update Jacobians ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);
    LinSysRes.AddBlock(jPoint, Residual);

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
    Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_j);
    Jacobian.AddBlock(jPoint, iPoint, Jacobian_i);
    Jacobian.AddBlock(jPoint, jPoint, Jacobian_j);

  }

}

void CRadP1Solver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                    CConfig *config, unsigned short iMesh) {
  unsigned long iPoint;

  bool harmonic_balance = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  bool transition    = (config->GetKind_Trans_Model() == LM);
  bool transition_BC = (config->GetKind_Trans_Model() == BC);

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Conservative variables w/o reconstruction ---*/

    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);

    /*--- Radiation variables w/o reconstruction ---*/

    numerics->SetRadVar(node[iPoint]->GetSolution(), NULL);

    /*--- Set volume ---*/

    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Compute the source term ---*/

    numerics->ComputeResidual(Residual, Jacobian_i, config);

    /*--- Subtract residual and the Jacobian ---*/

    LinSysRes.SubtractBlock(iPoint, Residual);

    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }

}

void CRadP1Solver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//  unsigned long iPoint, iVertex;
//  unsigned short iVar;
//
//  /*--- The dirichlet condition is used only without wall function, otherwise the
//   convergence is compromised as we are providing nu tilde values for the
//   first point of the wall  ---*/
//
//  if (!config->GetWall_Functions()) {
//
//    for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//      iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//      /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//
//      if (geometry->node[iPoint]->GetDomain()) {
//
//        /*--- Get the velocity vector ---*/
//
//        for (iVar = 0; iVar < nVar; iVar++)
//          Solution[iVar] = 0.0;
//
//        node[iPoint]->SetSolution_Old(Solution);
//        LinSysRes.SetBlock_Zero(iPoint);
//
//        /*--- Includes 1 in the diagonal ---*/
//
//        Jacobian.DeleteValsRowi(iPoint);
//      }
//    }
//  }
//  else {
//
//    /*--- Evaluate nu tilde at the closest point to the surface using the wall functions ---*/
//
//    SetNuTilde_WF(geometry, solver_container, conv_numerics, visc_numerics, config, val_marker);
//
//  }

}

void CRadP1Solver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {
//  unsigned long iPoint, iVertex;
//  unsigned short iVar;
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Get the velocity vector ---*/
//      for (iVar = 0; iVar < nVar; iVar++)
//        Solution[iVar] = 0.0;
//
//      node[iPoint]->SetSolution_Old(Solution);
//      LinSysRes.SetBlock_Zero(iPoint);
//
//      /*--- Includes 1 in the diagonal ---*/
//
//      Jacobian.DeleteValsRowi(iPoint);
//    }
//  }

}

void CRadP1Solver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//
//  unsigned long iPoint, iVertex;
//  unsigned short iVar, iDim;
//  su2double *Normal, *V_infty, *V_domain;
//
//  bool grid_movement  = config->GetGrid_Movement();
//
//  Normal = new su2double[nDim];
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Allocate the value at the infinity ---*/
//
//      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
//
//      /*--- Retrieve solution at the farfield boundary node ---*/
//
//      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
//
//      /*--- Grid Movement ---*/
//
//      if (grid_movement)
//        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
//
//      conv_numerics->SetPrimitive(V_domain, V_infty);
//
//      /*--- Set turbulent variable at the wall, and at infinity ---*/
//
//      for (iVar = 0; iVar < nVar; iVar++)
//        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//      Solution_j[0] = nu_tilde_Inf;
//      conv_numerics->SetTurbVar(Solution_i, Solution_j);
//
//      /*--- Set Normal (it is necessary to change the sign) ---*/
//
//      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//      for (iDim = 0; iDim < nDim; iDim++)
//        Normal[iDim] = -Normal[iDim];
//      conv_numerics->SetNormal(Normal);
//
//      /*--- Compute residuals and Jacobians ---*/
//
//      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Add residuals and Jacobians ---*/
//
//      LinSysRes.AddBlock(iPoint, Residual);
//      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
//    }
//  }
//
//  delete [] Normal;

}

void CRadP1Solver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
//
//  unsigned short iDim;
//  unsigned long iVertex, iPoint;
//  su2double *V_inlet, *V_domain, *Normal;
//
//  Normal = new su2double[nDim];
//
//  bool grid_movement  = config->GetGrid_Movement();
//  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);
//
//  /*--- Loop over all the vertices on this boundary marker ---*/
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Normal vector for this vertex (negate for outward convention) ---*/
//
//      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
//
//      /*--- Allocate the value at the inlet ---*/
//
//      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
//
//      /*--- Retrieve solution at the farfield boundary node ---*/
//
//      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
//
//      /*--- Set various quantities in the solver class ---*/
//
//      conv_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Set the turbulent variable states (prescribed for an inflow) ---*/
//
//      Solution_i[0] = node[iPoint]->GetSolution(0);
//
//      /*--- Load the inlet turbulence variable (uniform by default). ---*/
//
//      Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
//
//      conv_numerics->SetTurbVar(Solution_i, Solution_j);
//
//      /*--- Set various other quantities in the conv_numerics class ---*/
//
//      conv_numerics->SetNormal(Normal);
//
//      if (grid_movement)
//        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
//                                  geometry->node[iPoint]->GetGridVel());
//
//      /*--- Compute the residual using an upwind scheme ---*/
//
//      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//      LinSysRes.AddBlock(iPoint, Residual);
//
//      /*--- Jacobian contribution for implicit integration ---*/
//
//      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
////      /*--- Viscous contribution, commented out because serious convergence problems ---*/
////
////      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
////      visc_numerics->SetNormal(Normal);
////
////      /*--- Conservative variables w/o reconstruction ---*/
////
////      visc_numerics->SetPrimitive(V_domain, V_inlet);
////
////      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
////
////      visc_numerics->SetTurbVar(Solution_i, Solution_j);
////      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
////
////      /*--- Compute residual, and Jacobians ---*/
////
////      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
////
////      /*--- Subtract residual, and update Jacobians ---*/
////
////      LinSysRes.SubtractBlock(iPoint, Residual);
////      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//    }
//  }
//
//  /*--- Free locally allocated memory ---*/
//  delete[] Normal;
//
}

void CRadP1Solver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {
//  unsigned long iPoint, iVertex;
//  unsigned short iVar, iDim;
//  su2double *V_outlet, *V_domain, *Normal;
//
//  bool grid_movement  = config->GetGrid_Movement();
//
//  Normal = new su2double[nDim];
//
//  /*--- Loop over all the vertices on this boundary marker ---*/
//
//  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
//    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
//
//    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
//
//    if (geometry->node[iPoint]->GetDomain()) {
//
//      /*--- Allocate the value at the outlet ---*/
//
//      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
//
//      /*--- Retrieve solution at the farfield boundary node ---*/
//
//      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
//
//      /*--- Set various quantities in the solver class ---*/
//
//      conv_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Set the turbulent variables. Here we use a Neumann BC such
//       that the turbulent variable is copied from the interior of the
//       domain to the outlet before computing the residual.
//       Solution_i --> TurbVar_internal,
//       Solution_j --> TurbVar_outlet ---*/
//
//      for (iVar = 0; iVar < nVar; iVar++) {
//        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
//        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
//      }
//      conv_numerics->SetTurbVar(Solution_i, Solution_j);
//
//      /*--- Set Normal (negate for outward convention) ---*/
//
//      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
//      for (iDim = 0; iDim < nDim; iDim++)
//        Normal[iDim] = -Normal[iDim];
//      conv_numerics->SetNormal(Normal);
//
//      if (grid_movement)
//        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
//                                  geometry->node[iPoint]->GetGridVel());
//
//      /*--- Compute the residual using an upwind scheme ---*/
//
//      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//      LinSysRes.AddBlock(iPoint, Residual);
//
//      /*--- Jacobian contribution for implicit integration ---*/
//
//      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
//
////      /*--- Viscous contribution, commented out because serious convergence problems ---*/
////
////      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
////      visc_numerics->SetNormal(Normal);
////
////      /*--- Conservative variables w/o reconstruction ---*/
////
////      visc_numerics->SetPrimitive(V_domain, V_outlet);
////
////      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
////
////      visc_numerics->SetTurbVar(Solution_i, Solution_j);
////      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
////
////      /*--- Compute residual, and Jacobians ---*/
////
////      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
////
////      /*--- Subtract residual, and update Jacobians ---*/
////
////      LinSysRes.SubtractBlock(iPoint, Residual);
////      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//
//    }
//  }
//
//  /*--- Free locally allocated memory ---*/
//
//  delete[] Normal;
//
}

void CRadP1Solver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}
