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
  Refractive_Index = config->GetRefractive_Index();

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

void CRadP1Solver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {

  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Vol;
  su2double Delta_time = 0.01, Delta;

  /*--- Set maximum residual to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  /*--- Build implicit system ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Read the volume ---*/

    Vol = geometry->node[iPoint]->GetVolume();

    /*--- Modify matrix diagonal to assure diagonal dominance ---*/

    if (node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
      Jacobian.SetVal2Diag(iPoint, 1.0);
      for (iVar = 0; iVar < nVar; iVar++) {
        total_index = iPoint*nVar + iVar;
        LinSysRes[total_index] = 0.0;
      }
    }

    /*--- Right hand side of the system (-Residual) and initial guess (x = 0) ---*/

    for (iVar = 0; iVar < nVar; iVar++) {
      total_index = iPoint*nVar+iVar;
      LinSysRes[total_index] = - (LinSysRes[total_index]);
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

  //Set_MPI_Solution(geometry, config);

  /*--- Compute the root mean square residual ---*/

  SetResidual_RMS(geometry, config);

}

void CRadSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

}


//void CRadSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {
//
//  /*--- Restart the solution from file information ---*/
//
//  unsigned short iVar, iMesh;
//  unsigned long iPoint, index, iChildren, Point_Fine;
//  su2double Area_Children, Area_Parent, *Solution_Fine;
//  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
//                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
//  bool time_stepping = (config->GetUnsteady_Simulation() == TIME_STEPPING);
//  unsigned short iZone = config->GetiZone();
//  unsigned short nZone = config->GetnZone();
//
//  string UnstExt, text_line;
//  ifstream restart_file;
//  string restart_filename = config->GetSolution_FlowFileName();
//
//  /*--- Modify file name for multizone problems ---*/
//  if (nZone >1)
//    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);
//
//  /*--- Modify file name for an unsteady restart ---*/
//
//  if (dual_time|| time_stepping)
//    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);
//
//  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/
//
//  if (config->GetRead_Binary_Restart()) {
//    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
//  } else {
//    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
//  }
//
//  int counter = 0;
//  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
//  unsigned long iPoint_Global_Local = 0;
//  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
//
//  /*--- Skip flow variables ---*/
//
//  unsigned short skipVars = 0;
//
//  if (nDim == 2) skipVars += 6;
//  if (nDim == 3) skipVars += 8;
//
//  /*--- Load data from the restart into correct containers. ---*/
//
//  counter = 0;
//  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++ ) {
//
//
//    /*--- Retrieve local index. If this node from the restart file lives
//     on the current processor, we will load and instantiate the vars. ---*/
//
//    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);
//
//    if (iPoint_Local > -1) {
//
//      /*--- We need to store this point's data, so jump to the correct
//       offset in the buffer of data from the restart file and load it. ---*/
//
//      index = counter*Restart_Vars[1] + skipVars;
//      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
//      node[iPoint_Local]->SetSolution(Solution);
//      iPoint_Global_Local++;
//
//      /*--- Increment the overall counter for how many points have been loaded. ---*/
//      counter++;
//    }
//
//  }
//
//  /*--- Detect a wrong solution file ---*/
//
//  if (iPoint_Global_Local < nPointDomain) { sbuf_NotMatching = 1; }
//
//#ifndef HAVE_MPI
//  rbuf_NotMatching = sbuf_NotMatching;
//#else
//  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
//#endif
//  if (rbuf_NotMatching != 0) {
//    SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
//                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
//  }
//
//  /*--- MPI solution and compute the eddy viscosity ---*/
//
////TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart.
//  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
//  solver[MESH_0][TURB_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
//
//  solver[MESH_0][FLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
//  solver[MESH_0][TURB_SOL]->Postprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0);
//
//  /*--- Interpolate the solution down to the coarse multigrid levels ---*/
//
//  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
//    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
//      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
//      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
//      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
//        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
//        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
//        Solution_Fine = solver[iMesh-1][TURB_SOL]->node[Point_Fine]->GetSolution();
//        for (iVar = 0; iVar < nVar; iVar++) {
//          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
//        }
//      }
//      solver[iMesh][TURB_SOL]->node[iPoint]->SetSolution(Solution);
//    }
//    solver[iMesh][TURB_SOL]->Set_MPI_Solution(geometry[iMesh], config);
//    solver[iMesh][FLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
//    solver[iMesh][TURB_SOL]->Postprocessing(geometry[iMesh], solver[iMesh], config, iMesh);
//  }
//
//  /*--- Delete the class memory that is used to load the restart. ---*/
//
//  if (Restart_Vars != NULL) delete [] Restart_Vars;
//  if (Restart_Data != NULL) delete [] Restart_Data;
//  Restart_Vars = NULL; Restart_Data = NULL;
//
//}

CRadP1Solver::CRadP1Solver(void) : CRadSolver() {

  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

}

CRadP1Solver::CRadP1Solver(CGeometry* geometry, CConfig *config) : CRadSolver(geometry, config) {

  unsigned long iPoint;
  unsigned short iVar, iDim;

  nDim =          geometry->GetnDim();
  nPoint =        geometry->GetnPoint();
  nPointDomain =  geometry->GetnPointDomain();
  nVar =          1;
  node =          new CVariable*[nPoint];

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  Residual = new su2double[nVar]; Residual_RMS = new su2double[nVar];
  Solution = new su2double[nVar]; Residual_Max = new su2double[nVar];

  Res_Visc = new su2double[nVar];

  /*--- Define some structures for locating max residuals ---*/

  Point_Max = new unsigned long[nVar];
  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Jacobians and vector structures for implicit computations ---*/

  if (config->GetKind_TimeIntScheme_Radiation() == EULER_IMPLICIT) {

    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    if (rank == MASTER_NODE) cout << "Initialize Jacobian structure (P1 radiation equation)." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);

  }

  /*--- Solution and residual vectors ---*/

  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysAux.Initialize(nPoint, nPointDomain, nVar, 0.0);

  /*--- Define some auxiliary vectors for computing flow variable
   gradients by least squares, S matrix := inv(R)*traspose(inv(R)),
   c vector := transpose(WA)*(Wb) ---*/

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];

  }

  /*--- Always instantiate and initialize the variable to a zero value. ---*/
  su2double init_val;
  switch(config->GetKind_P1_Init()){
    case P1_INIT_ZERO: init_val = 0.0; break;
    case P1_INIT_TEMP: init_val = 4.0*pow(Refractive_Index,2.0)*STEFAN_BOLTZMANN*pow(config->GetInc_Temperature_Init(),4.0); break;
    default: init_val = 0.0; break;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CRadP1Variable(init_val, nDim, nVar, config);

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

  unsigned long iPoint;
  su2double Energy, Temperature;
  su2double SourceTerm, SourceTerm_Derivative;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Retrieve the radiative energy ---*/
    Energy = node[iPoint]->GetSolution(0);

    /*--- Retrieve temperature from the flow solver ---*/
    Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive()[nDim+1];

    /*--- Compute the divergence of the radiative flux ---*/
    SourceTerm = Absorption_Coeff*(Energy - 4.0*pow(Refractive_Index,2.0)*STEFAN_BOLTZMANN*pow(Temperature,4.0));

    /*--- Compute the derivative of the source term with respect to the temperature ---*/
    SourceTerm_Derivative =  - 16.0*Absorption_Coeff*pow(Refractive_Index,2.0)*STEFAN_BOLTZMANN*pow(Temperature,3.0);

    /*--- Store the source term and its derivative ---*/
    node[iPoint]->SetRadiative_SourceTerm(0, SourceTerm);
    node[iPoint]->SetRadiative_SourceTerm(1, SourceTerm_Derivative);

  }

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

  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, total_index;

  su2double Theta, Ib_w, Temperature, Radiative_Energy;
  su2double *Normal, Area, Wall_Emissivity;
  su2double Radiative_Heat_Flux;
  su2double *Gradient, *Unit_Normal;

  Unit_Normal = new su2double[nDim];

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool energy        = config->GetEnergy_Equation();

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall emissivity from config ---*/
  Wall_Emissivity = config->GetWall_Emissivity(Marker_Tag);

  /*--- Compute the constant for the wall theta ---*/
  Theta = Wall_Emissivity / (2.0*(2.0 - Wall_Emissivity));

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      // Weak application of the boundary condition

      /*--- Initialize the viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Apply a weak boundary condition for the radiative transfer equation. ---*/

      /*--- Retrieve temperature from the flow solver ---*/
      Temperature = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive()[nDim+1];

      /*--- Compute the blackbody intensity at the wall. ---*/
      Ib_w = 4.0*pow(Refractive_Index,2.0)*STEFAN_BOLTZMANN*pow(Temperature,4.0);

      /*--- Compute the radiative heat flux. ---*/
      Radiative_Energy = node[iPoint]->GetSolution(0);
      Radiative_Heat_Flux = Theta*(Ib_w - Radiative_Energy);

      /*--- Compute the Viscous contribution to the residual ---*/
      Res_Visc[0] = Radiative_Heat_Flux*Area;

      /*--- Apply to the residual vector ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Compute the Jacobian contribution. ---*/
      if (implicit) {
        Jacobian_i[0][0] = - Theta;
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }

    }
  }

}

void CRadP1Solver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config,
                                       unsigned short val_marker) {

  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, total_index;

  su2double Theta, Ib_w, Temperature, Radiative_Energy;
  su2double *Normal, *Unit_Normal, Area, Wall_Emissivity;
  su2double Radiative_Heat_Flux;
  su2double Twall;
  su2double *Gradient;

  Unit_Normal = new su2double[nDim];

  bool implicit      = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  bool energy        = config->GetEnergy_Equation();

  /*--- Identify the boundary by string name ---*/
  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Get the specified wall emissivity from config ---*/
  Wall_Emissivity = config->GetWall_Emissivity(Marker_Tag);

  /*--- Compute the constant for the wall theta ---*/
  Theta = Wall_Emissivity / (2.0*(2.0 - Wall_Emissivity));

    /*--- Retrieve the specified wall temperature ---*/
  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  /*--- Loop over all of the vertices on this boundary marker ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute dual-grid area and boundary normal ---*/
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);

      // Weak application of the boundary condition

      /*--- Initialize the viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Visc[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar++)
            Jacobian_i[iVar][jVar] = 0.0;
        }
      }

      /*--- Apply a weak boundary condition for the radiative transfer equation. ---*/

      /*--- Compute the blackbody intensity at the wall. ---*/
      Ib_w = 4.0*pow(Refractive_Index,2.0)*STEFAN_BOLTZMANN*pow(Twall,4.0);

      /*--- Compute the radiative heat flux. ---*/
      Radiative_Energy = node[iPoint]->GetSolution(0);
      Radiative_Heat_Flux = 1.0*Theta*(Ib_w - Radiative_Energy);

      /*--- Compute the Viscous contribution to the residual ---*/
      Res_Visc[0] = Radiative_Heat_Flux*Area;

      /*--- Apply to the residual vector ---*/
      LinSysRes.SubtractBlock(iPoint, Res_Visc);

      /*--- Compute the Jacobian contribution. ---*/
      if (implicit) {
        Jacobian_i[0][0] = - Theta;
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
      }
    }
  }

}

void CRadP1Solver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

}

void CRadP1Solver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {

}

void CRadP1Solver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                              CConfig *config, unsigned short val_marker) {

}

void CRadP1Solver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container,
                                CNumerics *numerics, CConfig *config, unsigned short val_marker) {

  /*--- Convective fluxes across euler wall are equal to zero. ---*/

}

void CRadP1Solver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                               unsigned short iMesh, unsigned long Iteration) {

  unsigned short iDim, iMarker;
  unsigned long iEdge, iVertex, iPoint = 0, jPoint = 0;
  su2double *Normal, Area, Vol, Lambda;
  su2double Global_Delta_Time = 1E6, Global_Delta_UnstTimeND = 0.0, Local_Delta_Time = 0.0, Local_Delta_Time_Inv, Local_Delta_Time_Visc, CFL_Reduction, K_v = 0.25;
  su2double CFL = config->GetCFL_Rad();
  su2double GammaP1 = 1.0 / (3.0*(Absorption_Coeff + Scattering_Coeff));

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  /*--- Compute spectral radius based on thermal conductivity ---*/

  Min_Delta_Time = 1.E6; Max_Delta_Time = 0.0;

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    node[iPoint]->SetMax_Lambda_Visc(0.0);
  }

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Get the edge's normal vector to compute the edge's area ---*/
    Normal = geometry->edge[iEdge]->GetNormal();
    Area = 0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

    /*--- Viscous contribution ---*/

    Lambda = GammaP1*Area*Area;
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddMax_Lambda_Visc(Lambda);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {

      /*--- Point identification, Normal vector and area ---*/

      iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim]; Area = sqrt(Area);

      /*--- Viscous contribution ---*/

      Lambda = GammaP1*Area*Area;
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->AddMax_Lambda_Visc(Lambda);

    }
  }

  /*--- Each element uses their own speed, steady state simulation ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    Vol = geometry->node[iPoint]->GetVolume();

    if (Vol != 0.0) {

      /*--- Time step setting method ---*/

       Local_Delta_Time = CFL*K_v*Vol*Vol/ node[iPoint]->GetMax_Lambda_Visc();

      /*--- Min-Max-Logic ---*/

      Global_Delta_Time = min(Global_Delta_Time, Local_Delta_Time);
      Min_Delta_Time = min(Min_Delta_Time, Local_Delta_Time);
      Max_Delta_Time = max(Max_Delta_Time, Local_Delta_Time);
      if (Local_Delta_Time > config->GetMax_DeltaTime())
        Local_Delta_Time = config->GetMax_DeltaTime();

      node[iPoint]->SetDelta_Time(Local_Delta_Time);
    }
    else {
      node[iPoint]->SetDelta_Time(0.0);
    }
  }

  /*--- Compute the max and the min dt (in parallel) ---*/
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Min_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Min_Delta_Time = rbuf_time;

    sbuf_time = Max_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Max_Delta_Time = rbuf_time;
#endif
  }

  /*--- For exact time solution use the minimum delta time of the whole mesh ---*/
  if (config->GetUnsteady_Simulation() == TIME_STEPPING) {
#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_Time;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_Time = rbuf_time;
#endif
    for (iPoint = 0; iPoint < nPointDomain; iPoint++)
      node[iPoint]->SetDelta_Time(Global_Delta_Time);
  }

  /*--- Recompute the unsteady time step for the dual time strategy
   if the unsteady CFL is diferent from 0 ---*/
  if ((dual_time) && (Iteration == 0) && (config->GetUnst_CFL() != 0.0) && (iMesh == MESH_0)) {
    Global_Delta_UnstTimeND = config->GetUnst_CFL()*Global_Delta_Time/config->GetCFL(iMesh);

#ifdef HAVE_MPI
    su2double rbuf_time, sbuf_time;
    sbuf_time = Global_Delta_UnstTimeND;
    SU2_MPI::Reduce(&sbuf_time, &rbuf_time, 1, MPI_DOUBLE, MPI_MIN, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(&rbuf_time, 1, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    Global_Delta_UnstTimeND = rbuf_time;
#endif
    config->SetDelta_UnstTimeND(Global_Delta_UnstTimeND);
  }

  /*--- The pseudo local time (explicit integration) cannot be greater than the physical time ---*/
  if (dual_time)
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      if (!implicit) {
        cout << "Using unsteady time: " << config->GetDelta_UnstTimeND() << endl;
        Local_Delta_Time = min((2.0/3.0)*config->GetDelta_UnstTimeND(), node[iPoint]->GetDelta_Time());
        node[iPoint]->SetDelta_Time(Local_Delta_Time);
      }
  }
}
