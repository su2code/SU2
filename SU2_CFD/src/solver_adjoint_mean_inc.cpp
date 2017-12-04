/*!
 * \file solution_adjoint_mean_inc.cpp
 * \brief Main subrotuines for solving adjoint incompressible flow (Euler, Navier-Stokes, etc.).
 * \author F. Palacios, T. Economon
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

CAdjIncEulerSolver::CAdjIncEulerSolver(void) : CSolver() {
  
  /*--- Array initialization ---*/
  Phi_Inf = NULL;
  Sens_Mach = NULL;
  Sens_AoA = NULL;
  Sens_Geo = NULL;
  Sens_Press = NULL;
  Sens_Temp = NULL;
  Sens_BPress = NULL;
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Jacobian_Axisymmetric = NULL;
  CSensitivity = NULL;
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;
  
}

CAdjIncEulerSolver::CAdjIncEulerSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CSolver() {
  unsigned long iPoint, iVertex;
  string text_line, mesh_filename;
  unsigned short iDim, iVar, iMarker, nLineLets;
  ifstream restart_file;
  string filename, AdjExt;
  su2double myArea_Monitored, Area, *Normal;
  bool axisymmetric = config->GetAxisymmetric();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Array initialization ---*/
  Phi_Inf = NULL;
  Sens_Mach = NULL;
  Sens_AoA = NULL;
  Sens_Geo = NULL;
  Sens_Press = NULL;
  Sens_Temp = NULL;
  Sens_BPress = NULL;
  iPoint_UndLapl = NULL;
  jPoint_UndLapl = NULL;
  Jacobian_Axisymmetric = NULL;
  CSensitivity = NULL;
  FlowPrimVar_i = NULL;
  FlowPrimVar_j = NULL;

  /*--- Set the gamma value ---*/
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constans in the solver structure ---*/
  nDim = geometry->GetnDim();
  nMarker = config->GetnMarker_All();
  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  nVar = nDim + 1;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliary vectors related to the residual ---*/
  Residual = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 0.0;
  Point_Max = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  Residual_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]    = 0.0;
  Residual_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]    = 0.0;
  Res_Conv_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]    = 0.0;
  Res_Visc_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]    = 0.0;
  Res_Conv_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]    = 0.0;
  Res_Visc_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]    = 0.0;
  
  /*--- Define some auxiliary vectors related to the solution ---*/
  Solution   = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar]   = 0.0;
  Solution_j = new su2double[nVar];  for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary arrays related to the flow solution ---*/
  FlowPrimVar_i = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_i[iVar] = 0.0;
  FlowPrimVar_j = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Define some auxiliary vectors related to the undivided lapalacian ---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) {
    iPoint_UndLapl = new su2double [nPoint];
    jPoint_UndLapl = new su2double [nPoint];
  }
  
  /*--- Define some auxiliary vectors related to the geometry ---*/
  Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
    Jacobian_ii = new su2double* [nVar];
    Jacobian_ij = new su2double* [nVar];
    Jacobian_ji = new su2double* [nVar];
    Jacobian_jj = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_ii[iVar] = new su2double [nVar];
      Jacobian_ij[iVar] = new su2double [nVar];
      Jacobian_ji[iVar] = new su2double [nVar];
      Jacobian_jj[iVar] = new su2double [nVar];
    }
    
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
    if (axisymmetric) {
      Jacobian_Axisymmetric = new su2double* [nVar];
      for (iVar = 0; iVar < nVar; iVar++)
        Jacobian_Axisymmetric[iVar] = new su2double [nVar];
    }
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint Euler). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Sensitivity definition and coefficient in all the markers ---*/
  CSensitivity = new su2double* [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
  }
  Sens_Geo  = new su2double[nMarker];
  Sens_Mach = new su2double[nMarker];
  Sens_AoA  = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];
  Sens_BPress = new su2double[nMarker];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker]  = 0.0;
    Sens_Mach[iMarker] = 0.0;
    Sens_AoA[iMarker]  = 0.0;
    Sens_Press[iMarker] = 0.0;
    Sens_Temp[iMarker]  = 0.0;
    Sens_BPress[iMarker] = 0.0;
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
      CSensitivity[iMarker][iVertex] = 0.0;
  }
  
  /*--- Adjoint flow at the inifinity, initialization stuff ---*/
  PsiRho_Inf = 0.0; PsiE_Inf   = 0.0;
  Phi_Inf    = new su2double [nDim];
  Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
  if (nDim == 3) Phi_Inf[2] = 0.0;
  
  /*--- If outflow objective, nonzero initialization ---*/
  if ((config->GetKind_ObjFunc() == SURFACE_TOTAL_PRESSURE)){
    su2double SoundSpeed,*vel_inf,R,vel2,vel;
    R = config->GetGas_ConstantND();
    vel_inf = config->GetVelocity_FreeStreamND();
    vel2=0;
    for (iDim=0; iDim<nDim; iDim++)
      vel2 +=vel_inf[iDim]*vel_inf[iDim];
    vel = pow(vel2,0.5);
    SoundSpeed= pow(Gamma*config->GetTemperature_FreeStreamND()*R, 0.5);
    PsiE_Inf = Gamma_Minus_One*vel2/(vel2-pow(SoundSpeed,2.0))*0.5/vel;
    PsiRho_Inf += PsiE_Inf*(2*SoundSpeed*SoundSpeed+vel2*Gamma_Minus_One)/(2.0*Gamma_Minus_One);
    // Assumes +x flow direction
    // Assume v.n = |v|, n = -v/|v|

    for (iDim=0; iDim<nDim; iDim++){
      Phi_Inf[iDim] +=PsiE_Inf*(SoundSpeed*SoundSpeed/Gamma_Minus_One/vel2-1)*vel_inf[iDim];
      // Assumes n in direction of v
      Phi_Inf[iDim]+=vel_inf[iDim]/vel*(0.5);
    }

  }

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CAdjIncEulerVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
  
  /*--- Define solver parameters needed for execution of destructor ---*/
  if (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED) space_centered = true;
  else space_centered = false;
  
  /*--- Calculate area monitored for area-averaged-outflow-quantity-based objectives ---*/
  myArea_Monitored = 0.0;
  if (config->GetKind_ObjFunc()==SURFACE_TOTAL_PRESSURE || config->GetKind_ObjFunc()==SURFACE_STATIC_PRESSURE){
    for (iMarker =0; iMarker < config->GetnMarker_All();  iMarker++){
      if (config->GetMarker_All_Monitoring(iMarker) == YES){
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            myArea_Monitored += sqrt (Area);
          }
        }
      }
    }
  }
#ifdef HAVE_MPI
  Area_Monitored = 0.0;
  SU2_MPI::Allreduce(&myArea_Monitored, &Area_Monitored, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  Area_Monitored = myArea_Monitored;
#endif

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjIncEulerSolver::~CAdjIncEulerSolver(void) {
  unsigned short iVar, iMarker;
  
  if (Phi_Inf != NULL) delete [] Phi_Inf;
  if (Sens_Mach != NULL) delete [] Sens_Mach;
  if (Sens_AoA != NULL) delete [] Sens_AoA;
  if (Sens_Geo != NULL) delete [] Sens_Geo;
  if (Sens_Press != NULL) delete [] Sens_Press;
  if (Sens_Temp != NULL) delete [] Sens_Temp;
  if (Sens_BPress != NULL) delete [] Sens_BPress;
  if (iPoint_UndLapl != NULL) delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL) delete [] jPoint_UndLapl;
  if (FlowPrimVar_i != NULL) delete [] FlowPrimVar_i;
  if (FlowPrimVar_j != NULL) delete [] FlowPrimVar_j;
  
  if (Jacobian_Axisymmetric != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete Jacobian_Axisymmetric[iVar];
    delete [] Jacobian_Axisymmetric;
  }
  
  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      delete [] CSensitivity[iMarker];
    delete [] CSensitivity;
  }
  
}

void CAdjIncEulerSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration) {

  /*--- Use the flow solution to update the time step
   *    The time step depends on the characteristic velocity, which is the same
   *    for the adjoint and flow solutions, albeit in the opposite direction. ---*/
  solver_container[FLOW_SOL]->SetTime_Step(geometry, solver_container, config, iMesh, Iteration);
}


void CAdjIncEulerSolver::Set_MPI_Solution(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi, *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_U = new su2double [nBufferR_Vector];
      Buffer_Send_U = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_U[iVar*nVertexS+iVertex] = node[iPoint]->GetSolution(iVar);
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
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
  
}

void CAdjIncEulerSolver::Set_MPI_Solution_Old(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_U = NULL, *Buffer_Send_U = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
          Solution[iVar] = Buffer_Receive_U[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_U[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_U[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_U[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_U[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetSolution_Old(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_U;
      
    }
    
  }
}

void CAdjIncEulerSolver::Set_MPI_Solution_Limiter(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Limit = NULL, *Buffer_Send_Limit = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Limit = new su2double [nBufferR_Vector];
      Buffer_Send_Limit = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Limit[iVar*nVertexS+iVertex] = node[iPoint]->GetLimiter(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Limit, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Limit, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
      
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Limit[iVar*nVertexR+iVertex] = Buffer_Send_Limit[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Limit;
      
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
          Solution[iVar] = Buffer_Receive_Limit[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Limit[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Limit[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Limit[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetLimiter(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Limit;
      
    }
    
  }
}

void CAdjIncEulerSolver::Set_MPI_Solution_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
#endif
  
  su2double **Gradient = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Gradient[iVar] = new su2double[nDim];
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
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


void CAdjIncEulerSolver::Set_MPI_Undivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned short iVar, iMarker, iPeriodic_Index, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
  *Buffer_Receive_Undivided_Laplacian = NULL, *Buffer_Send_Undivided_Laplacian = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
      nBufferS_Vector = nVertexS*nVar;        nBufferR_Vector = nVertexR*nVar;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Undivided_Laplacian = new su2double [nBufferR_Vector];
      Buffer_Send_Undivided_Laplacian = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Send_Undivided_Laplacian[iVar*nVertexS+iVertex] = node[iPoint]->GetUndivided_Laplacian(iVar);
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Undivided_Laplacian, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Undivided_Laplacian, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        for (iVar = 0; iVar < nVar; iVar++)
          Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex] = Buffer_Send_Undivided_Laplacian[iVar*nVertexR+iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Undivided_Laplacian;
      
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
          Solution[iVar] = Buffer_Receive_Undivided_Laplacian[iVar*nVertexR+iVertex];
        
        /*--- Rotate the momentum components. ---*/
        if (nDim == 2) {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex];
        }
        else {
          Solution[1] = rotMatrix[0][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[0][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[0][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[2] = rotMatrix[1][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[1][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[1][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
          Solution[3] = rotMatrix[2][0]*Buffer_Receive_Undivided_Laplacian[1*nVertexR+iVertex] +
          rotMatrix[2][1]*Buffer_Receive_Undivided_Laplacian[2*nVertexR+iVertex] +
          rotMatrix[2][2]*Buffer_Receive_Undivided_Laplacian[3*nVertexR+iVertex];
        }
        
        /*--- Copy transformed conserved variables back into buffer. ---*/
        for (iVar = 0; iVar < nVar; iVar++)
          node[iPoint]->SetUndivided_Laplacian(iVar, Solution[iVar]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Undivided_Laplacian;
      
    }
    
  }
  
}

void CAdjIncEulerSolver::Set_MPI_Dissipation_Switch(CGeometry *geometry, CConfig *config) {
  unsigned short iMarker, MarkerS, MarkerR;
  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
  su2double *Buffer_Receive_Lambda = NULL, *Buffer_Send_Lambda = NULL;
  
#ifdef HAVE_MPI
  int send_to, receive_from;
  MPI_Status status;
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
      nBufferS_Vector = nVertexS;        nBufferR_Vector = nVertexR;
      
      /*--- Allocate Receive and send buffers  ---*/
      Buffer_Receive_Lambda = new su2double [nBufferR_Vector];
      Buffer_Send_Lambda = new su2double[nBufferS_Vector];
      
      /*--- Copy the solution old that should be sended ---*/
      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
        Buffer_Send_Lambda[iVertex] = node[iPoint]->GetSensor();
      }
      
#ifdef HAVE_MPI
      
      /*--- Send/Receive information using Sendrecv ---*/
      SU2_MPI::Sendrecv(Buffer_Send_Lambda, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
                   Buffer_Receive_Lambda, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
#else
      
      /*--- Receive information without MPI ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        Buffer_Receive_Lambda[iVertex] = Buffer_Send_Lambda[iVertex];
      }
      
#endif
      
      /*--- Deallocate send buffer ---*/
      delete [] Buffer_Send_Lambda;
      
      /*--- Do the coordinate transformation ---*/
      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
        
        /*--- Find point and its type of transformation ---*/
        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
        node[iPoint]->SetSensor(Buffer_Receive_Lambda[iVertex]);
        
      }
      
      /*--- Deallocate receive buffer ---*/
      delete [] Buffer_Receive_Lambda;
      
    }
    
  }
}

void CAdjIncEulerSolver::SetForceProj_Vector(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  su2double *ForceProj_Vector, x = 0.0, y = 0.0, z = 0.0, *Normal, CD, CL, Cp, CpTarget,
  CT, CQ, x_origin, y_origin, z_origin, WDrag, Area, invCD, CLCD2, invCQ, CTRCQ2;
  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  
  int rank = MASTER_NODE;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  su2double Alpha            = (config->GetAoA()*PI_NUMBER)/180.0;
  su2double Beta             = (config->GetAoS()*PI_NUMBER)/180.0;
  su2double RefLength  = config->GetRefLength();
  su2double *RefOriginMoment = config->GetRefOriginMoment(0);

  ForceProj_Vector = new su2double[nDim];
  
  /*--- Compute coefficients needed for objective function evaluation. ---*/
  
  CD = solver_container[FLOW_SOL]->GetTotal_CD();
  CL = solver_container[FLOW_SOL]->GetTotal_CL();
  CT = solver_container[FLOW_SOL]->GetTotal_CT();
  CQ = solver_container[FLOW_SOL]->GetTotal_CQ();
  invCD  = 1.0/CD; CLCD2  = CL/(CD*CD);
  invCQ  = 1.0/CQ; CTRCQ2 = CT/(RefLength*CQ*CQ);
  
  x_origin = RefOriginMoment[0]; y_origin = RefOriginMoment[1]; z_origin = RefOriginMoment[2];
  
  /*--- Evaluate the boundary condition coefficients. ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    
    if ((config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE) &&
        (config->GetMarker_All_Monitoring(iMarker) == YES))
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        x = geometry->node[iPoint]->GetCoord(0);
        y = geometry->node[iPoint]->GetCoord(1);
        if (nDim == 3) z = geometry->node[iPoint]->GetCoord(2);
        
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        switch (config->GetKind_ObjFunc()) {
          case DRAG_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha); ForceProj_Vector[1] = sin(Alpha); }
            if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta); ForceProj_Vector[1] = sin(Beta); ForceProj_Vector[2] = sin(Alpha)*cos(Beta); }
            break;
          case LIFT_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = cos(Alpha); }
            if (nDim == 3) { ForceProj_Vector[0] = -sin(Alpha); ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = cos(Alpha); }
            break;
          case SIDEFORCE_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl;
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) { ForceProj_Vector[0] = -sin(Beta) * cos(Alpha); ForceProj_Vector[1] = cos(Beta); ForceProj_Vector[2] = -sin(Beta) * sin(Alpha); }
            break;
          case INVERSE_DESIGN_PRESSURE :
            Cp = solver_container[FLOW_SOL]->GetCPressure(iMarker, iVertex);
            CpTarget = solver_container[FLOW_SOL]->GetCPressureTarget(iMarker, iVertex);
            Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1]);
            if (nDim == 3) Area = sqrt(Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);
            ForceProj_Vector[0] = -2.0*(Cp-CpTarget)*Normal[0]/Area; ForceProj_Vector[1] = -2.0*(Cp-CpTarget)*Normal[1]/Area;
            if (nDim == 3) ForceProj_Vector[2] = -2.0*(Cp-CpTarget)*Normal[2]/Area;
            break;
          case MOMENT_X_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(EXIT_FAILURE); }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = -(z - z_origin)/RefLength; ForceProj_Vector[2] = (y - y_origin)/RefLength; }
            break;
          case MOMENT_Y_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) { cout << "This functional is not possible in 2D!!" << endl; exit(EXIT_FAILURE); }
            if (nDim == 3) { ForceProj_Vector[0] = (z - z_origin)/RefLength; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = -(x - x_origin)/RefLength; }
            break;
          case MOMENT_Z_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = -(y - y_origin)/RefLength; ForceProj_Vector[1] = (x - x_origin)/RefLength; }
            if (nDim == 3) { ForceProj_Vector[0] = -(y - y_origin)/RefLength; ForceProj_Vector[1] = (x - x_origin)/RefLength; ForceProj_Vector[2] = 0; }
            break;
          case EFFICIENCY :
            if (nDim == 2) { ForceProj_Vector[0] = -(invCD*sin(Alpha)+CLCD2*cos(Alpha)); ForceProj_Vector[1] = (invCD*cos(Alpha)-CLCD2*sin(Alpha)); }
            if (nDim == 3) { ForceProj_Vector[0] = -(invCD*sin(Alpha)+CLCD2*cos(Alpha)*cos(Beta)); ForceProj_Vector[1] = -CLCD2*sin(Beta); ForceProj_Vector[2] = (invCD*cos(Alpha)-CLCD2*sin(Alpha)*cos(Beta)); }
            break;
          case EQUIVALENT_AREA :
            WDrag = config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha)*WDrag; ForceProj_Vector[1] = sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = sin(Beta)*WDrag; ForceProj_Vector[2] = sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case NEARFIELD_PRESSURE :
            WDrag = config->GetWeightCd();
            if (nDim == 2) { ForceProj_Vector[0] = cos(Alpha)*WDrag; ForceProj_Vector[1] = sin(Alpha)*WDrag; }
            if (nDim == 3) { ForceProj_Vector[0] = cos(Alpha)*cos(Beta)*WDrag; ForceProj_Vector[1] = sin(Beta)*WDrag; ForceProj_Vector[2] = sin(Alpha)*cos(Beta)*WDrag; }
            break;
          case FORCE_X_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 1.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
            break;
          case FORCE_Y_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 1.0; ForceProj_Vector[2] = 0.0; }
            break;
          case FORCE_Z_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 1.0; }
            break;
          case THRUST_COEFFICIENT :
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 1.0; }
            break;
          case TORQUE_COEFFICIENT :
            if (nDim == 2) { ForceProj_Vector[0] = (y - y_origin)/RefLength; ForceProj_Vector[1] = -(x - x_origin)/RefLength; }
            if (nDim == 3) { ForceProj_Vector[0] = (y - y_origin)/RefLength; ForceProj_Vector[1] = -(x - x_origin)/RefLength; ForceProj_Vector[2] = 0; }
            break;
          case FIGURE_OF_MERIT :
            if ((nDim == 2) && (rank == MASTER_NODE)) {cout << "This functional is not possible in 2D!!" << endl;
              exit(EXIT_FAILURE);
            }
            if (nDim == 3) {
              ForceProj_Vector[0] = -invCQ;
              ForceProj_Vector[1] = -CTRCQ2*(z - z_origin);
              ForceProj_Vector[2] =  CTRCQ2*(y - y_origin);
            }
            break;
          default :
            if (nDim == 2) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; }
            if (nDim == 3) { ForceProj_Vector[0] = 0.0; ForceProj_Vector[1] = 0.0; ForceProj_Vector[2] = 0.0; }
            break;
        }
        
        /*--- Store the force projection vector at this node ---*/
        
        node[iPoint]->SetForceProj_Vector(ForceProj_Vector);
        
      }
  
  delete [] ForceProj_Vector;
  
}

void CAdjIncEulerSolver::SetIntBoundary_Jump(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  unsigned short iMarker, iVar, jVar, kVar, iDim, jDim, iIndex;
  unsigned long iVertex, iPoint, iPointNearField, nPointNearField = 0;
  su2double factor = 1.0, AngleDouble, data, aux, *IntBound_Vector, *coord, *FlowSolution, WeightSB, MinDist = 1E6, Dist, DerivativeOF = 0.0, *Normal, Area, UnitNormal[3], velocity[3], Energy, Rho, sqvel, proj_vel, phi, a1, a2;
  su2double **A, **M, **AM, *b;
  short AngleInt = 0, IndexNF_inv[180], iColumn;
  ifstream index_file;
  string text_line;
  vector<vector<su2double> > NearFieldWeight;
  vector<su2double> CoordNF;
  vector<short> IndexNF;
  
  IntBound_Vector = new su2double [nVar];
  
  /*--- Allocate vectors and matrices ---*/
  
  b = new su2double [nVar];
  A = new su2double* [nVar];
  M = new su2double* [nVar];
  AM = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    A[iVar] = new su2double [nVar];
    M[iVar] = new su2double [nVar];
    AM[iVar] = new su2double [nVar];
  }
  
  /*--- If equivalent area objective function, read the value of
   the derivative from a file, this is a preprocess of the direct solution ---*/
  
  if (config->GetKind_ObjFunc() == EQUIVALENT_AREA) {
    
    /*--- Read derivative of the objective function at the NearField from file ---*/
    index_file.open("WeightNF.dat", ios::in);
    if (index_file.fail()) {
      cout << "There is no Weight Nearfield Pressure file (WeightNF.dat)." << endl;
      exit(EXIT_FAILURE);
    }
    
    nPointNearField = 0;
    
    while (index_file) {
      string line;
      getline(index_file, line);
      istringstream is(line);
      
      /*--- The first row provides the azimuthal angle ---*/
      
      if (nPointNearField == 0) {
        is >> data; // The first column is related with the coordinate
        while (is.good()) { is >> data; IndexNF.push_back(SU2_TYPE::Int(data)); }
      }
      else {
        is >> data; CoordNF.push_back(data); // The first column is the point coordinate
        vector<su2double> row;
        while (is.good()) { is >> data; row.push_back(data); }
        NearFieldWeight.push_back(row);
      }
      nPointNearField++;
    }
    
    /*--- Note tha the first row is the azimuthal angle ---*/
    
    nPointNearField = nPointNearField - 1;
    
    for (AngleInt = 0; AngleInt < 180; AngleInt++)
      IndexNF_inv[AngleInt] = -1;
    
  if (IndexNF.size() <= 180) {
    for (iIndex = 0; iIndex < IndexNF.size(); iIndex++)
      IndexNF_inv[IndexNF[iIndex]] = iIndex;
  }
  else {
    #ifndef HAVE_MPI
        exit(EXIT_FAILURE);
    #else
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, 1);
        MPI_Finalize();
    #endif
  }
    
  }
  
  /*--- Compute the jump on the adjoint variables for the upper and the lower side ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt (Area);
        
        for (iDim = 0; iDim < nDim; iDim++)
          UnitNormal[iDim] = Normal[iDim]/Area;
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          coord = geometry->node[iPoint]->GetCoord();
          DerivativeOF = 0.0;
          
          /*--- Just in case the functional depend also on the surface pressure ---*/
          
          WeightSB = 1.0-config->GetWeightCd();
          
          su2double AoA, XcoordRot = 0.0, YcoordRot = 0.0, ZcoordRot = 0.0;
          
          if (nDim == 2) XcoordRot = coord[0];
          if (nDim == 3) {
            
            /*--- Rotate the nearfield cylinder  ---*/
            
            AoA = -(config->GetAoA()*PI_NUMBER/180.0);
            XcoordRot = coord[0]*cos(AoA) - coord[2]*sin(AoA);
            YcoordRot = coord[1];
            ZcoordRot = coord[0]*sin(AoA) + coord[2]*cos(AoA);
          }
          
          switch (config->GetKind_ObjFunc()) {
            case EQUIVALENT_AREA :
              
              if (nDim == 2) AngleInt = 0;
              
              if (nDim == 3) {
                
                /*--- Compute the azimuthal angle of the iPoint ---*/
                
                AngleDouble = fabs(atan(-YcoordRot/ZcoordRot)*180.0/PI_NUMBER);
                
                /*--- Fix an azimuthal line due to misalignments of the near-field ---*/
                
                su2double FixAzimuthalLine = config->GetFixAzimuthalLine();
                
                if ((AngleDouble >= FixAzimuthalLine - 0.1) && (AngleDouble <= FixAzimuthalLine + 0.1)) AngleDouble = FixAzimuthalLine - 0.1;
                
                AngleInt = SU2_TYPE::Short(floor(AngleDouble + 0.5));
                if (AngleInt < 0) AngleInt = 180 + AngleInt;
                
              }
              
              if (AngleInt <= 60) {
                iColumn = IndexNF_inv[AngleInt];
                
                /*--- An azimuthal angle is not defined... this happens with MG levels ---*/
                
                if (iColumn < 0.0) {
                  if (IndexNF_inv[AngleInt+1] > 0) { iColumn = IndexNF_inv[AngleInt+1]; goto end; }
                  if (IndexNF_inv[AngleInt-1] > 0) { iColumn = IndexNF_inv[AngleInt-1]; goto end; }
                  if (IndexNF_inv[AngleInt+2] > 0) { iColumn = IndexNF_inv[AngleInt+2]; goto end; }
                  if (IndexNF_inv[AngleInt-2] > 0) { iColumn = IndexNF_inv[AngleInt-2]; goto end; }
                  if (IndexNF_inv[AngleInt+3] > 0) { iColumn = IndexNF_inv[AngleInt+3]; goto end; }
                  if (IndexNF_inv[AngleInt-3] > 0) { iColumn = IndexNF_inv[AngleInt-3]; goto end; }
                  if (IndexNF_inv[AngleInt+4] > 0) { iColumn = IndexNF_inv[AngleInt+4]; goto end; }
                  if (IndexNF_inv[AngleInt-4] > 0) { iColumn = IndexNF_inv[AngleInt-4]; goto end; }
                }
                
              end:
                
                if (iColumn < 0.0) { cout <<" An azimuthal angle is not defined..." << endl; }
                
                /*--- Find the value of the weight in the table, using the azimuthal angle  ---*/
                
                MinDist = 1E6;
                for (iPointNearField = 0; iPointNearField < nPointNearField; iPointNearField++) {
                  Dist = fabs(CoordNF[iPointNearField] - XcoordRot);
                  if (Dist <= MinDist) {
                    MinDist = Dist;
                    DerivativeOF = factor*WeightSB*NearFieldWeight[iPointNearField][iColumn];
                  }
                }
              }
              else DerivativeOF = 0.0;
              
              if ((MinDist > 1E-6) || (coord[nDim-1] > 0.0)) DerivativeOF = 0.0;
              
              break;
              
            case NEARFIELD_PRESSURE :
              
              DerivativeOF = factor*WeightSB*(solver_container[FLOW_SOL]->node[iPoint]->GetPressure()
                                              - solver_container[FLOW_SOL]->GetPressure_Inf());
              
              break;
              
          }
          
          /*--- Compute the jump of the adjoint variables (2D, and 3D problems) --*/
          
          FlowSolution = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
          
          Rho = FlowSolution[0];
          Energy = FlowSolution[nVar-1]/FlowSolution[0];
          
          sqvel = 0.0; proj_vel = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            velocity[iDim] = FlowSolution[iDim+1]/FlowSolution[0];
            sqvel    += velocity[iDim]*velocity[iDim];
            proj_vel += velocity[iDim]*UnitNormal[iDim];
          }
          
          if (nDim == 2) {
            
            /*--- Compute the projected Jacobian ---*/
            
            A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 1.0; A[0][3] = 0.0;
            A[1][0] = -velocity[0]*velocity[1]; A[1][1] = velocity[1]; A[1][2] = velocity[0]; A[1][3] = 0.0;
            A[2][0] = 0.5*(Gamma-3.0)*velocity[1]*velocity[1]+0.5*Gamma_Minus_One*velocity[0]*velocity[0]; A[2][1] = -Gamma_Minus_One*velocity[0];
            A[2][2] = (3.0-Gamma)*velocity[1]; A[2][3] = Gamma_Minus_One; A[3][0] = -Gamma*velocity[1]*Energy+Gamma_Minus_One*velocity[1]*sqvel;
            A[3][1] = -Gamma_Minus_One*velocity[0]*velocity[1]; A[3][2] = Gamma*Energy-0.5*Gamma_Minus_One*(velocity[0]*velocity[0]+3.0*velocity[1]*velocity[1]);  A[3][3] = Gamma*velocity[1];
            
            /*--- Compute the transformation matrix ---*/
            
            M[0][0] = 1.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = 0.0;
            M[1][0] = velocity[0]; M[1][1] = Rho; M[1][2] = 0.0; M[1][3] = 0.0;
            M[2][0] = velocity[1]; M[2][1] = 0.0; M[2][2] = Rho; M[2][3] = 0.0;
            M[3][0] = 0.5*sqvel;  M[3][1] = Rho*velocity[0]; M[3][2] = Rho*velocity[1]; M[3][3] = 1.0/Gamma_Minus_One;
            
            /*--- Create the soruce term (AM)^T X = b ---*/
            
            b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = DerivativeOF;
            
          }
          
          if (nDim == 3) {
            
            
            /*--- Compute the projected Jacobian ---*/
            
            phi = 0.5*Gamma_Minus_One*sqvel;
            a1 = Gamma*Energy-phi; a2 = Gamma-1.0;
            
            A[0][0] = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) A[0][iDim+1] = UnitNormal[iDim];
            A[0][nDim+1] = 0.0;
            
            for (iDim = 0; iDim < nDim; iDim++) {
              A[iDim+1][0] = (UnitNormal[iDim]*phi - velocity[iDim]*proj_vel);
              for (jDim = 0; jDim < nDim; jDim++)
                A[iDim+1][jDim+1] = (UnitNormal[jDim]*velocity[iDim]-a2*UnitNormal[iDim]*velocity[jDim]);
              A[iDim+1][iDim+1] += proj_vel;
              A[iDim+1][nDim+1] = a2*UnitNormal[iDim];
            }
            
            A[nDim+1][0] = proj_vel*(phi-a1);
            for (iDim = 0; iDim < nDim; iDim++)
              A[nDim+1][iDim+1] = (UnitNormal[iDim]*a1-a2*velocity[iDim]*proj_vel);
            A[nDim+1][nDim+1] = Gamma*proj_vel;
            
            /*--- Compute the transformation matrix ---*/
            
            M[0][0] = 1.0; M[0][1] = 0.0; M[0][2] = 0.0; M[0][3] = 0.0; M[0][4] = 0.0;
            M[1][0] = velocity[0]; M[1][1] = Rho; M[1][2] = 0.0; M[1][3] = 0.0; M[1][4] = 0.0;
            M[2][0] = velocity[1]; M[2][1] = 0.0; M[2][2] = Rho; M[2][3] = 0.0; M[2][4] = 0.0;
            M[3][0] = velocity[2]; M[3][1] = 0.0; M[3][2] = 0.0; M[3][3] = Rho; M[3][4] = 0.0;
            M[4][0] = 0.5*sqvel; M[4][1] = Rho*velocity[0]; M[4][2] = Rho*velocity[1];
            M[4][3] = Rho*velocity[2]; M[4][4] = 1.0/Gamma_Minus_One;
            
            /*--- Create the soruce term (AM)^T X = b ---*/
            
            b[0] = 0.0; b[1] = 0.0; b[2] = 0.0; b[3] = 0.0; b[4] = DerivativeOF;
            
          }
          
          /*--- Compute A times M ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++) {
              aux = 0.0;
              for (kVar = 0; kVar < nVar; kVar++)
                aux += A[iVar][kVar]*M[kVar][jVar];
              AM[iVar][jVar] = aux;
            }
          
          /*--- Compute the transpose matrix ---*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            for (jVar = 0; jVar < nVar; jVar++)
              A[iVar][jVar] = AM[jVar][iVar];
          
          /*--- Solve the linear system using a LU descomposition --*/
          
          Gauss_Elimination(A, b, nVar);
          
          /*--- Update the internal boundary jump --*/
          
          for (iVar = 0; iVar < nVar; iVar++)
            IntBound_Vector[iVar] = b[iVar];
          
          node[iPoint]->SetIntBoundary_Jump(IntBound_Vector);
          
        }
      }
  
  delete [] IntBound_Vector;
  
  /*--- Deallocate the linear system ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    delete [] A[iVar];
    delete [] M[iVar];
    delete [] AM[iVar];
  }
  delete [] A;
  delete [] M;
  delete [] AM;
  delete [] b;
  
}

void CAdjIncEulerSolver::SetInitialCondition(CGeometry **geometry, CSolver ***solver_container, CConfig *config, unsigned long ExtIter) {
  unsigned long iPoint, Point_Fine;
  unsigned short iMesh, iChildren, iVar;
  su2double Area_Children, Area_Parent, *Solution, *Solution_Fine;
  
  bool restart = config->GetRestart();
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  
  /*--- If restart solution, then interpolate the flow solution to
   all the multigrid levels, this is important with the dual time strategy ---*/
  if (restart) {
    Solution = new su2double[nVar];
    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
        Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
        for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
        for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
          Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
          Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
          Solution_Fine = solver_container[iMesh-1][ADJFLOW_SOL]->node[Point_Fine]->GetSolution();
          for (iVar = 0; iVar < nVar; iVar++) {
            Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
          }
        }
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution(Solution);
        
      }
      solver_container[iMesh][ADJFLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    }
    delete [] Solution;
  }
  
  /*--- The value of the solution for the first iteration of the dual time ---*/
  for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      if ((ExtIter == 0) && (dual_time)) {
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->Set_Solution_time_n();
        solver_container[iMesh][ADJFLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
      }
    }
  }
  
}

void CAdjIncEulerSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double SharpEdge_Distance;
  bool RightSol = true;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool muscl          = config->GetMUSCL_AdjFlow();
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool center         = (config->GetKind_ConvNumScheme_AdjFlow() == SPACE_CENTERED);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);

  /*--- Residual initialization ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Initialize the non-physical points vector ---*/

    node[iPoint]->SetNon_Physical(false);
    
    /*--- Set the primitive variables incompressible adjoint variables ---*/
    
    RightSol = node[iPoint]->SetPrimVar(SharpEdge_Distance, false, config);
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  
  if ((muscl) && (iMesh == MESH_0)) {
    
    /*--- Compute gradients for upwind second-order reconstruction ---*/

    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    
    /*--- Limiter computation ---*/
    
    if (limiter && !Output) SetSolution_Limiter(geometry, config);
    
  }
  
  /*--- Artificial dissipation for centered schemes ---*/
  
  if (center) {
    if ((center_jst) && (iMesh == MESH_0) && !Output) {
      SetDissipation_Switch(geometry, config);
      SetUndivided_Laplacian(geometry, config);
      if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
      if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
    }
  }
  
  /*--- Initialize the Jacobian for implicit integration ---*/
  
  if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CAdjIncEulerSolver::Centered_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                        CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  
  unsigned long iEdge, iPoint, jPoint;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool jst_scheme = ((config->GetKind_Centered_AdjFlow() == JST) && (iMesh == MESH_0));
  bool grid_movement  = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, normal, and neighbors---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    numerics->SetNeighbor(geometry->node[iPoint]->GetnNeighbor(), geometry->node[jPoint]->GetnNeighbor());
    
    /*--- Adjoint variables w/o reconstruction ---*/
    
    numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
    /*--- Conservative variables w/o reconstruction ---*/
    
    numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                              solver_container[FLOW_SOL]->node[jPoint]->GetSolution());

    numerics->SetDensity(solver_container[FLOW_SOL]->node[iPoint]->GetDensity(), solver_container[FLOW_SOL]->node[jPoint]->GetDensity());
    numerics->SetBetaInc2(solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2(), solver_container[FLOW_SOL]->node[jPoint]->GetBetaInc2());
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());

    
    numerics->SetLambda(solver_container[FLOW_SOL]->node[iPoint]->GetLambda(),
                        solver_container[FLOW_SOL]->node[jPoint]->GetLambda());
    
    if (jst_scheme) {
      numerics->SetUndivided_Laplacian(node[iPoint]->GetUndivided_Laplacian(), node[jPoint]->GetUndivided_Laplacian());
      numerics->SetSensor(node[iPoint]->GetSensor(), node[jPoint]->GetSensor());
    }
    
    /*--- Mesh motion ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- Compute residuals ---*/
    
    numerics->ComputeResidual(Res_Conv_i, Res_Visc_i, Res_Conv_j, Res_Visc_j,
                              Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Update convective and artificial dissipation residuals ---*/
    
    LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
    LinSysRes.SubtractBlock(jPoint, Res_Conv_j);
    LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
    LinSysRes.SubtractBlock(jPoint, Res_Visc_j);
    
    /*--- Implicit contribution to the residual ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
}


void CAdjIncEulerSolver::Upwind_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short iMesh) {
  
  su2double **Gradient_i, **Gradient_j, Project_Grad_i, Project_Grad_j, *Limiter_i = NULL,
  *Limiter_j = NULL, *Psi_i = NULL, *Psi_j = NULL, *V_i, *V_j, Non_Physical = 1.0;
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iDim, iVar;
  
  bool implicit         = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool muscl            = (config->GetMUSCL_AdjFlow() && (iMesh == MESH_0));
  bool limiter          = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool grid_movement    = config->GetGrid_Movement();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge and normal vectors ---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Adjoint variables w/o reconstruction ---*/
    
    Psi_i = node[iPoint]->GetSolution();
    Psi_j = node[jPoint]->GetSolution();
    numerics->SetAdjointVar(Psi_i, Psi_j);
    
    /*--- Primitive variables w/o reconstruction ---*/
    
    V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    V_j = solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();
    numerics->SetPrimitive(V_i, V_j);
    
    /*--- Grid velocities for dynamic meshes ---*/
    
    if (grid_movement) {
      numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[jPoint]->GetGridVel());
    }
    
    /*--- High order reconstruction using MUSCL strategy ---*/
    
    if (muscl) {
      
      for (iDim = 0; iDim < nDim; iDim++) {
        Vector_i[iDim] = 0.5*(geometry->node[jPoint]->GetCoord(iDim) - geometry->node[iPoint]->GetCoord(iDim));
        Vector_j[iDim] = 0.5*(geometry->node[iPoint]->GetCoord(iDim) - geometry->node[jPoint]->GetCoord(iDim));
      }
      
      /*--- Adjoint variables using gradient reconstruction and limiters ---*/

      Gradient_i = node[iPoint]->GetGradient(); Gradient_j = node[jPoint]->GetGradient();
      if (limiter) { Limiter_i = node[iPoint]->GetLimiter(); Limiter_j = node[jPoint]->GetLimiter(); }
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Project_Grad_i = 0; Project_Grad_j = 0;
        Non_Physical = node[iPoint]->GetNon_Physical()*node[jPoint]->GetNon_Physical();
        for (iDim = 0; iDim < nDim; iDim++) {
          Project_Grad_i += Vector_i[iDim]*Gradient_i[iVar][iDim]*Non_Physical;
          Project_Grad_j += Vector_j[iDim]*Gradient_j[iVar][iDim]*Non_Physical;
        }
        if (limiter) {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i*Limiter_i[iDim];
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j*Limiter_j[iDim];
        }
        else {
          Solution_i[iVar] = Psi_i[iVar] + Project_Grad_i;
          Solution_j[iVar] = Psi_j[iVar] + Project_Grad_j;
        }
      }
      
      numerics->SetAdjointVar(Solution_i, Solution_j);
      
    }
    
    /*--- Compute the residual---*/
    
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
    
    /*--- Add and Subtract Residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.SubtractBlock(jPoint, Residual_j);
    
    /*--- Implicit contribution to the residual ---*/
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.SubtractBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.SubtractBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
  
}

void CAdjIncEulerSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                      CConfig *config, unsigned short iMesh) {
  
  unsigned short iVar;
  unsigned long iPoint;
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  bool axisymmetric   = config->GetAxisymmetric();
  //  bool gravity        = (config->GetGravityForce() == YES);
  bool harmonic_balance  = (config->GetUnsteady_Simulation() == HARMONIC_BALANCE);
  
  /*--- Initialize the source residual to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++) Residual[iVar] = 0.0;
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the adjoint variables ---*/
      numerics->SetAdjointVar(node[iPoint]->GetSolution(),
                              node[iPoint]->GetSolution());
      
      /*--- Load the volume of the dual mesh cell ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the adjoint rotating frame source residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
  if (harmonic_balance) {
    
    su2double Volume, Source;
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Get control volume ---*/
      Volume = geometry->node[iPoint]->GetVolume();
      
      /*--- Get stored time spectral source term ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Source = node[iPoint]->GetHarmonicBalance_Source(iVar);
        Residual[iVar] = Source*Volume;
      }
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
    }
  }
  
  if (axisymmetric) {
    
    /*--- Zero out Jacobian structure ---*/
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar ++)
        for (unsigned short jVar = 0; jVar < nVar; jVar ++)
          Jacobian_i[iVar][jVar] = 0.0;
    }
    
    /*--- loop over points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Set solution ---*/
      numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(), solver_container[FLOW_SOL]->node[iPoint]->GetSolution());
      
      /*--- Set adjoint variables ---*/
      numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[iPoint]->GetSolution());
      
      /*--- Set control volume ---*/
      numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Set coordinate ---*/
      numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[iPoint]->GetCoord());
      
      /*--- Compute Source term Residual ---*/
      numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add Residual ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Implicit part ---*/
      if (implicit)
        Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CAdjIncEulerSolver::Source_Template(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                      CConfig *config, unsigned short iMesh) {
}

void CAdjIncEulerSolver::SetUndivided_Laplacian(CGeometry *geometry, CConfig *config) {
  unsigned long iPoint, jPoint, iEdge;
  unsigned short iVar;
  su2double *Diff;
  
  Diff = new su2double[nVar];
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    node[iPoint]->SetUnd_LaplZero();
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    for (iVar = 0; iVar < nVar; iVar++)
      Diff[iVar] = node[iPoint]->GetSolution(iVar) - node[jPoint]->GetSolution(iVar);
    
#ifdef STRUCTURED_GRID
    
    if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#else
    
    bool boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
    bool boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
    
    /*--- Both points inside the domain, or both in the boundary ---*/
    if ((!boundary_i && !boundary_j) || (boundary_i && boundary_j)) {
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    }
    
    /*--- iPoint inside the domain, jPoint on the boundary ---*/
    if (!boundary_i && boundary_j)
      if (geometry->node[iPoint]->GetDomain()) node[iPoint]->SubtractUnd_Lapl(Diff);
    
    /*--- jPoint inside the domain, iPoint on the boundary ---*/
    if (boundary_i && !boundary_j)
      if (geometry->node[jPoint]->GetDomain()) node[jPoint]->AddUnd_Lapl(Diff);
    
#endif
    
  }
  
#ifdef STRUCTURED_GRID
  
  unsigned long Point_Normal = 0, iVertex;
  unsigned short iMarker;
  su2double *Psi_mirror;
  
  Psi_mirror = new su2double[nVar];
  
  /*--- Loop over all boundaries and include an extra contribution
   from a halo node. Find the nearest normal, interior point
   for a boundary node and make a linear approximation. ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          Point_Normal = geometry->vertex[iMarker][iVertex]->GetNormal_Neighbor();
          
          /*--- Interpolate & compute difference in the conserved variables ---*/
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_mirror[iVar] = 2.0*node[iPoint]->GetSolution(iVar) - node[Point_Normal]->GetSolution(iVar);
            Diff[iVar]   = node[iPoint]->GetSolution(iVar) - Psi_mirror[iVar];
          }
          
          /*--- Subtract contribution at the boundary node only ---*/
          node[iPoint]->SubtractUnd_Lapl(Diff);
        }
      }
    }
  }
  
  delete [] Psi_mirror;
  
#endif
  
  delete [] Diff;
  
  /*--- MPI parallelization ---*/
  Set_MPI_Undivided_Laplacian(geometry, config);
  
}

void CAdjIncEulerSolver::SetDissipation_Switch(CGeometry *geometry, CConfig *config) {
  
  unsigned long iPoint;
  su2double SharpEdge_Distance, eps, ds, scale, Sensor, Param_Kappa_2, Param_Kappa_4;
  
  eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
  Param_Kappa_2 = config->GetKappa_2nd_AdjFlow();
  Param_Kappa_4 = config->GetKappa_4th_AdjFlow();
  
  if (Param_Kappa_2 != 0.0) scale = 2.0 * Param_Kappa_4 / Param_Kappa_2;
  else scale = 0.0;
  
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    
    SharpEdge_Distance = (geometry->node[iPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps);
    
    ds = 0.0;
    if (SharpEdge_Distance < -eps) ds = 1.0;
    if (fabs(SharpEdge_Distance) <= eps) ds = 1.0 - (0.5*(1.0+(SharpEdge_Distance/eps)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps)));
    if (SharpEdge_Distance > eps) ds = 0.0;
    
    Sensor = scale * ds;
    
    node[iPoint]->SetSensor(Sensor);
    
  }
  
  /*--- MPI parallelization ---*/
  Set_MPI_Dissipation_Switch(geometry, config);
  
}

void CAdjIncEulerSolver::ExplicitRK_Iteration(CGeometry *geometry, CSolver **solver_container,
                                           CConfig *config, unsigned short iRKStep) {
  su2double *Residual, *Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  su2double RK_AlphaCoeff = config->Get_Alpha_RKStep(iRKStep);
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
    
    Res_TruncError = node[iPoint]->GetResTruncError();
    Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = Residual[iVar] + Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta*RK_AlphaCoeff);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
    
  }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjIncEulerSolver::ExplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  su2double *local_Residual, *local_Res_TruncError, Vol, Delta, Res;
  unsigned short iVar;
  unsigned long iPoint;
  
  for (iVar = 0; iVar < nVar; iVar++) {
    SetRes_RMS(iVar, 0.0);
    SetRes_Max(iVar, 0.0, 0);
  }
  
  /*--- Update the solution ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Vol = geometry->node[iPoint]->GetVolume();
    Delta = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() / Vol;
    
    local_Res_TruncError = node[iPoint]->GetResTruncError();
    local_Residual = LinSysRes.GetBlock(iPoint);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      Res = local_Residual[iVar] + local_Res_TruncError[iVar];
      node[iPoint]->AddSolution(iVar, -Res*Delta);
      AddRes_RMS(iVar, Res*Res);
      AddRes_Max(iVar, fabs(Res), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
    
  }
  
  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  SetResidual_RMS(geometry, config);
  
}

void CAdjIncEulerSolver::ImplicitEuler_Iteration(CGeometry *geometry, CSolver **solver_container, CConfig *config) {
  
  unsigned short iVar;
  unsigned long iPoint, total_index;
  su2double Delta, *local_Res_TruncError, Vol;
  
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
    
    if (solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time() != 0.0) {
      Delta = Vol / solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();
      Jacobian.AddVal2Diag(iPoint, Delta);
    }
    else {
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
      LinSysRes[total_index] = -(LinSysRes[total_index] + local_Res_TruncError[iVar]);
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
  
  /*--- Update solution (system written in terms of increments) ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++)
    for (iVar = 0; iVar < nVar; iVar++) {
      node[iPoint]->AddSolution(iVar, config->GetRelaxation_Factor_AdjFlow()*LinSysSol[iPoint*nVar+iVar]);
    }
  
  /*--- MPI solution ---*/
  
  Set_MPI_Solution(geometry, config);
  
  /*--- Compute the root mean square residual ---*/
  
  SetResidual_RMS(geometry, config);
  
}

void CAdjIncEulerSolver::Inviscid_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned long iVertex, iPoint, Neigh;
  unsigned short iDim, iMarker, iNeigh;
  su2double *d = NULL, *Normal = NULL, *Psi = NULL, *U = NULL,  conspsi = 0.0,
  Area, **PrimVar_Grad = NULL, **ConsVar_Grad = NULL, *ConsPsi_Grad = NULL,
  ConsPsi, d_press, grad_v, Beta2, v_gradconspsi, *GridVel = NULL,
  eps, *USens, scale = 1.0;
  su2double RefVel2, RefDensity, Mach2Vel, *Velocity_Inf, factor;
  su2double *Velocity;
  
  USens = new su2double[nVar];
  Velocity = new su2double[nDim];

  su2double Gas_Constant    = config->GetGas_ConstantND();
  bool grid_movement     = config->GetGrid_Movement();
  su2double RefArea    = config->GetRefArea();
  su2double Mach_Motion     = config->GetMach_Motion();
  unsigned short ObjFunc = config->GetKind_ObjFunc();

  if (config->GetSystemMeasurements() == US) scale = 1.0/12.0;
  else scale = 1.0;
  
  /*--- Compute non-dimensional factor. For dynamic meshes, use the motion Mach 
   number as a reference value for computing the force coefficients. 
   Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = config->GetDensity_FreeStreamND();

  factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);
  
  if ((ObjFunc == INVERSE_DESIGN_HEATFLUX) ||
      (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
      (ObjFunc == SURFACE_MASSFLOW) ) factor = 1.0;

 if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE))
   factor = 1.0/Area_Monitored;
  
  /*--- Initialize sensitivities to zero ---*/
  
  Total_Sens_Geo = 0.0;
  Total_Sens_Mach = 0.0;
  Total_Sens_AoA = 0.0;
  Total_Sens_Press = 0.0;
  Total_Sens_Temp = 0.0;
  Total_Sens_BPress = 0.0;
  
  /*--- Loop over boundary markers to select those for Euler walls ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++)

    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL)
      
      /*--- Loop over points on the surface to store the auxiliary variable ---*/
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          Psi = node[iPoint]->GetSolution();
          U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
          Beta2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();
          conspsi = Beta2*Psi[0];

          for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
          
          node[iPoint]->SetAuxVar(conspsi);
          
          /*--- Also load the auxiliary variable for first neighbors ---*/

          for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
            Neigh = geometry->node[iPoint]->GetPoint(iNeigh);
            Psi = node[Neigh]->GetSolution();
            U = solver_container[FLOW_SOL]->node[Neigh]->GetSolution();
            Beta2 = solver_container[FLOW_SOL]->node[Neigh]->GetBetaInc2();
            conspsi = Beta2*Psi[0];
            for (iDim = 0; iDim < nDim; iDim++) conspsi += U[iDim+1]*Psi[iDim+1];
            node[Neigh]->SetAuxVar(conspsi);
          }
        }
      }
  
  /*--- Compute surface gradients of the auxiliary variable ---*/
  
  SetAuxVar_Surface_Gradient(geometry, config);
  
  /*--- Evaluate the shape sensitivity ---*/
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker] = 0.0;
    
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (geometry->node[iPoint]->GetDomain()) {
          
          d = node[iPoint]->GetForceProj_Vector();
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0;
          for (iDim = 0; iDim < nDim; iDim++)
            Area += Normal[iDim]*Normal[iDim];
          Area = sqrt(Area);
          
          PrimVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
          ConsVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient();
          ConsPsi_Grad = node[iPoint]->GetAuxVarGradient();
          ConsPsi = node[iPoint]->GetAuxVar();
          
          d_press = 0.0; grad_v = 0.0; v_gradconspsi = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            
            /*-- Retrieve the value of the pressure gradient ---*/
            
            d_press += d[iDim]*ConsVar_Grad[0][iDim];
            
            /*-- Retrieve the value of the velocity gradient ---*/
            
            grad_v += PrimVar_Grad[iDim+1][iDim]*ConsPsi;
            
            /*-- Retrieve the value of the theta gradient ---*/
            
            v_gradconspsi += solver_container[FLOW_SOL]->node[iPoint]->GetVelocity(iDim) * ConsPsi_Grad[iDim];
            
            /*--- Additional sensitivity term for grid movement ---*/
            
            if (grid_movement) {
              GridVel = geometry->node[iPoint]->GetGridVel();
              v_gradconspsi -= GridVel[iDim] * ConsPsi_Grad[iDim];
            }
            
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          
          CSensitivity[iMarker][iVertex] = (d_press + grad_v + v_gradconspsi) * Area * scale * factor;
          
          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
          
          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
            if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
              CSensitivity[iMarker][iVertex] = 0.0;
          }
          
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex];
          
        }
      }
      
      Total_Sens_Geo += Sens_Geo[iMarker];
      
    }
  }

#ifdef HAVE_MPI
  
  su2double MyTotal_Sens_Geo   = Total_Sens_Geo;     Total_Sens_Geo = 0.0;
  su2double MyTotal_Sens_Mach  = Total_Sens_Mach;    Total_Sens_Mach = 0.0;
  su2double MyTotal_Sens_AoA   = Total_Sens_AoA;     Total_Sens_AoA = 0.0;
  su2double MyTotal_Sens_Press = Total_Sens_Press;   Total_Sens_Press = 0.0;
  su2double MyTotal_Sens_Temp  = Total_Sens_Temp;    Total_Sens_Temp = 0.0;
  su2double MyTotal_Sens_BPress  = Total_Sens_BPress;    Total_Sens_BPress = 0.0;
  
  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_BPress, &Total_Sens_BPress, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  delete [] USens;
  delete [] Velocity;
  
}

void CAdjIncEulerSolver::Smooth_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  unsigned short iMarker;
  unsigned long iVertex, jVertex, nVertex, iPoint;
  su2double **A, *b, Sens, *ArchLength, *Coord_begin, *Coord_end, dist;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == EULER_WALL) {
      nVertex = geometry->nVertex[iMarker];
      
      /*--- Allocate the linear system ---*/
      
      A = new su2double* [nVertex];
      b = new su2double [nVertex];
      ArchLength = new su2double [nVertex];
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        A[iVertex] = new su2double [nVertex];
      }
      
      /*--- Initialization ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        b[iVertex] = 0.0; ArchLength[iVertex] = 0.0;
        for (jVertex = 0; jVertex < nVertex; jVertex++)
          A[iVertex][jVertex] = 0.0;
      }
      
      /*--- Set the arch length ---*/
      
      ArchLength[0] = 0.0;
      for (iVertex = 1; iVertex < nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex-1]->GetNode();
        Coord_begin = geometry->node[iPoint]->GetCoord();
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Coord_end = geometry->node[iPoint]->GetCoord();
        dist = sqrt (pow( Coord_end[0]-Coord_begin[0], 2.0) + pow( Coord_end[1]-Coord_begin[1], 2.0));
        ArchLength[iVertex] = ArchLength[iVertex-1] + dist;
      }
      
      /*--- Remove the trailing edge effect ---*/
      
      su2double MinPosSens = 0.0; su2double MinNegSens = 0.0;
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        Sens = CSensitivity[iMarker][iVertex];
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.01) { MinNegSens = Sens; break; }
      }
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        Sens = CSensitivity[iMarker][iVertex];
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99) { MinPosSens = Sens; break; }
      }
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        if (ArchLength[iVertex] < ArchLength[nVertex-1]*0.01)
          CSensitivity[iMarker][iVertex] = MinNegSens;
        if (ArchLength[iVertex] > ArchLength[nVertex-1]*0.99)
          CSensitivity[iMarker][iVertex] = MinPosSens;
      }
      
      /*--- Set the right hand side of the system ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        b[iVertex] = CSensitivity[iMarker][iVertex];
      }
      
      /*--- Set the mass matrix ---*/
      
      su2double Coeff = 0.0, BackDiff = 0.0, ForwDiff = 0.0, CentDiff = 0.0;
      su2double epsilon = 5E-5;
      for (iVertex = 0; iVertex < nVertex; iVertex++) {
        
        if ((iVertex != nVertex-1) && (iVertex != 0)) {
          BackDiff = (ArchLength[iVertex]-ArchLength[iVertex-1]);
          ForwDiff = (ArchLength[iVertex+1]-ArchLength[iVertex]);
          CentDiff = (ArchLength[iVertex+1]-ArchLength[iVertex-1]);
        }
        if (iVertex == nVertex-1) {
          BackDiff = (ArchLength[nVertex-1]-ArchLength[nVertex-2]);
          ForwDiff = (ArchLength[0]-ArchLength[nVertex-1]);
          CentDiff = (ArchLength[0]-ArchLength[nVertex-2]);
        }
        if (iVertex == 0) {
          BackDiff = (ArchLength[0]-ArchLength[nVertex-1]);
          ForwDiff = (ArchLength[1]-ArchLength[0]);
          CentDiff = (ArchLength[1]-ArchLength[nVertex-1]);
        }
        
        Coeff = epsilon*2.0/(BackDiff*ForwDiff*CentDiff);
        
        A[iVertex][iVertex] = Coeff*CentDiff;
        
        if (iVertex != 0) A[iVertex][iVertex-1] = -Coeff*ForwDiff;
        else A[iVertex][nVertex-1] = -Coeff*ForwDiff;
        
        if (iVertex != nVertex-1) A[iVertex][iVertex+1] = -Coeff*BackDiff;
        else A[iVertex][0] = -Coeff*BackDiff;
        
      }
      
      /*--- Add the gradient value in the main diagonal ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        A[iVertex][iVertex] += 1.0;
      
      /*--- Dirichlet boundary condition ---*/
      
      unsigned long iVertex = SU2_TYPE::Int(nVertex/2);
      A[iVertex][iVertex] = 1.0;
      A[iVertex][iVertex+1] = 0.0;
      A[iVertex][iVertex-1] = 0.0;
      
      Gauss_Elimination(A, b, (unsigned short)nVertex);
      
      /*--- Set the new value of the sensitiviy ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        CSensitivity[iMarker][iVertex] = b[iVertex];
      
      /*--- Deallocate the linear system ---*/
      
      for (iVertex = 0; iVertex < nVertex; iVertex++)
        delete [] A[iVertex];
      delete [] A;
      delete [] b;
      delete [] ArchLength;
      
    }
  }
  
  
}

void CAdjIncEulerSolver::BC_Euler_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config, unsigned short val_marker) {
  unsigned long iVertex, iPoint;
  su2double *d = NULL, *Normal, *U, *Psi_Aux, bcn, Area, *UnitNormal;
  su2double *Velocity, *Psi, phin, phis1, phis2, DensityInc = 0.0, BetaInc2 = 0.0;
  unsigned short iDim, iVar, jDim;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi      = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();
      
      /*--- Create a copy of the adjoint solution ---*/
      Psi_Aux = node[iPoint]->GetSolution();
      for (iVar = 0; iVar < nVar; iVar++) Psi[iVar] = Psi_Aux[iVar];
      
      /*--- Flow solution ---*/
      U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
      
      /*--- Read the value of the objective function ---*/
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Normal vector computation ---*/
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      for (iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = -Normal[iDim]/Area;
      
      /*--- Incompressible solver ---*/

      DensityInc = solver_container[FLOW_SOL]->node[iPoint]->GetDensity();
      BetaInc2 = solver_container[FLOW_SOL]->node[iPoint]->GetBetaInc2();

      for (iDim = 0; iDim < nDim; iDim++)
        Velocity[iDim] = U[iDim+1] / solver_container[FLOW_SOL]->node[iPoint]->GetDensity();

      /*--- Compute projections ---*/
      bcn = 0.0; phin = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        bcn += d[iDim]*UnitNormal[iDim];
        phin += Psi[iDim+1]*UnitNormal[iDim];
      }

      /*--- Introduce the boundary condition ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        Psi[iDim+1] -= ( phin - bcn ) * UnitNormal[iDim];

      /*--- Inner products after introducing BC (Psi has changed) ---*/
      phis1 = 0.0; phis2 = Psi[0] * (BetaInc2 / DensityInc);
      for (iDim = 0; iDim < nDim; iDim++) {
        phis1 -= Normal[iDim]*Psi[iDim+1];
        phis2 += Velocity[iDim]*Psi[iDim+1];
      }

      /*--- Flux of the Euler wall ---*/
      Residual[0] = phis1;
      for (iDim = 0; iDim < nDim; iDim++)
        Residual[iDim+1] = - phis2 * Normal[iDim];

      /*--- Update residual ---*/
      LinSysRes.SubtractBlock(iPoint, Residual);

      if (implicit) {

        /*--- Adjoint density ---*/
        Jacobian_ii[0][0] = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          Jacobian_ii[0][iDim+1] = - Normal[iDim];

        /*--- Adjoint velocities ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          Jacobian_ii[iDim+1][0] = -Normal[iDim] * (BetaInc2 / DensityInc) ;
          for (jDim = 0; jDim < nDim; jDim++)
            Jacobian_ii[iDim+1][jDim+1] = - Normal[iDim] * Velocity[jDim];
        }

        /*--- Update Jacobian ---*/
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      }


      
    }
  }
  
  delete [] Velocity;
  delete [] UnitNormal;
  delete [] Psi;
  
}

void CAdjIncEulerSolver::BC_Sym_Plane(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics,
                                   CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint;
  su2double *Normal, Area, *UnitNormal,
  *V_domain, *V_sym, *Psi_domain, *Psi_sym, *Velocity, NormalAdjVel;
  unsigned short iDim, iVar;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);

  Normal = new su2double[nDim];
  UnitNormal = new su2double[nDim];
  Velocity = new su2double[nDim];
  Psi_domain = new su2double[nVar];
  Psi_sym = new su2double[nVar];
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      
      Area = 0;
      for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
      Area = sqrt(Area);
      
      for (iDim = 0; iDim < nDim; iDim++)
        UnitNormal[iDim]   = -Normal[iDim]/Area;
      

      /*--- Set the normal vector ---*/

      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Retrieve solution at boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      V_sym = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      conv_numerics->SetPrimitive(V_domain, V_sym);

      /*--- Adjoint flow solution at the wall ---*/

      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
        Psi_sym[iVar] = node[iPoint]->GetSolution(iVar);
      }

      /*--- Compute normal component of the adjoint velocity ---*/

      NormalAdjVel = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        NormalAdjVel += Psi_domain[iDim+1]*UnitNormal[iDim];

      /*--- Remove the normal component ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        Psi_sym[iDim+1] -= NormalAdjVel*UnitNormal[iDim];

      /*--- Set the value of the adjoint variables ---*/

      conv_numerics->SetAdjointVar(Psi_domain, Psi_sym);

      /*--- Compute the upwind flux ---*/

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Update residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit stuff ---*/
      
      if (implicit) {
        
         /*--- Update jacobian ---*/
        
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        
      }

    }
  }
  
  delete [] Velocity;
  delete [] Psi_domain;
  delete [] Psi_sym;
  delete [] Normal;
  delete [] UnitNormal;
  
}

void CAdjIncEulerSolver::BC_Interface_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config) {
  
  unsigned long iVertex, iPoint, jPoint;
  unsigned short iDim, iVar, iMarker;
  su2double *V_i, *V_j;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *Psi_i = new su2double[nVar];
  su2double *Psi_j = new su2double[nVar];
  
#ifndef HAVE_MPI
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Adjoint variables w/o reconstruction ---*/
          
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_i[iVar] = node[iPoint]->GetSolution(iVar);
            Psi_j[iVar] = node[jPoint]->GetSolution(iVar);
          }
          numerics->SetAdjointVar(Psi_i, Psi_j);
          
          /*--- Conservative variables w/o reconstruction ---*/
          
          V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
          V_j = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
          numerics->SetPrimitive(V_i, V_j);
          
          /*--- Set face vector, and area ---*/
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
          for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
          numerics->SetNormal(Normal);
          
          /*--- Compute residual ---*/
          
          numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
          
          /*--- Add Residuals and Jacobians ---*/
          
          LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
          if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
          
        }
      }
    }
  }
  
#else
  
  int rank, jProcessor;
  MPI_Status send_stat[1], recv_stat[1];
  MPI_Request send_req[1], recv_req[1];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  bool compute;
  su2double *Buffer_Send_Psi = new su2double[nVar];
  su2double *Buffer_Receive_Psi = new su2double[nVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == INTERFACE_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          /*--- We only send the information that belong to other boundary ---*/
          
          if (compute) {
            
            if (jProcessor != rank) {
              
              /*--- Copy the adjoint variable ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                Buffer_Send_Psi[iVar] = node[iPoint]->GetSolution(iVar);
              
              SU2_MPI::Isend(Buffer_Send_Psi, nVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD, &send_req[0]);
              
              /*--- Wait for this set of non-blocking comm. to complete ---*/
              
              SU2_MPI::Waitall(1, send_req, send_stat);
              
            }
            
          }
          
        }
      }
      
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          if (compute) {
            
            /*--- We only receive the information that belong to other boundary ---*/
            
            if (jProcessor != rank) {
              
              SU2_MPI::Irecv(Buffer_Receive_Psi, nVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req[0]);
              
              /*--- Wait for the this set of non-blocking recv's to complete ---*/
              
              SU2_MPI::Waitall(1, recv_req, recv_stat);
              
            } else {
              for (iVar = 0; iVar < nVar; iVar++)
                Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar);
            }
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Psi_i[iVar] = node[iPoint]->GetSolution(iVar);
              Psi_j[iVar] = Buffer_Receive_Psi[iVar];
            }
            
            /*--- Set adjoint Variables ---*/
            
            numerics->SetAdjointVar(Psi_i, Psi_j);
            
            /*--- Conservative variables w/o reconstruction (the same at both points) ---*/
            
            V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
            V_j = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
            numerics->SetPrimitive(V_i, V_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute the convective residual using an upwind scheme ---*/
            numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
            if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
            
          }
          
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  delete[] Buffer_Send_Psi;
  delete[] Buffer_Receive_Psi;
  
#endif
  
  delete[] Normal;
  delete[] Psi_i;
  delete[] Psi_j;
  
}

void CAdjIncEulerSolver::BC_NearField_Boundary(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                            CConfig *config) {
  
  unsigned long iVertex, iPoint, jPoint, Pin, Pout;
  unsigned short iDim, iVar, iMarker;
  su2double *V_i, *V_j, *IntBoundary_Jump;
  
  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  
  su2double *Normal = new su2double[nDim];
  su2double *Psi_i = new su2double[nVar];
  su2double *Psi_j = new su2double[nVar];
  su2double *Psi_out = new su2double[nVar];
  su2double *Psi_in = new su2double[nVar];
  su2double *MeanPsi = new su2double[nVar];
  su2double *Psi_out_ghost = new su2double[nVar];
  su2double *Psi_in_ghost = new su2double[nVar];
  
  
#ifndef HAVE_MPI
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        jPoint = geometry->vertex[iMarker][iVertex]->GetDonorPoint();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Adjoint variables w/o reconstruction ---*/
          
          for (iVar = 0; iVar < nVar; iVar++) {
            Psi_i[iVar] = node[iPoint]->GetSolution(iVar);
            Psi_j[iVar] = node[jPoint]->GetSolution(iVar);
          }
          
          /*--- If equivalent area or nearfield pressure condition ---*/
          
          if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
              (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {
            
            /*--- Identify the inner and the outer point (based on the normal direction) ---*/
            
            if (Normal[nDim-1] < 0.0) { Pin = iPoint; Pout = jPoint; }
            else { Pout = iPoint; Pin = jPoint; }
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Psi_out[iVar] = node[Pout]->GetSolution(iVar);
              Psi_in[iVar] = node[Pin]->GetSolution(iVar);
              MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
            }
            
            IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();
            
            /*--- Inner point ---*/
            
            if (iPoint == Pin) {
              for (iVar = 0; iVar < nVar; iVar++)
                Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
              numerics->SetAdjointVar(Psi_in, Psi_in_ghost);
            }
            
            /*--- Outer point ---*/
            
            if (iPoint == Pout) {
              for (iVar = 0; iVar < nVar; iVar++)
                Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
              numerics->SetAdjointVar(Psi_out, Psi_out_ghost);
            }
            
          }
          else {
            
            /*--- Just do a periodic BC ---*/
            
            numerics->SetAdjointVar(Psi_i, Psi_j);
            
          }
          
          /*--- Conservative variables w/o reconstruction ---*/
          
          V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
          V_j = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
          numerics->SetPrimitive(V_i, V_j);
          
          /*--- Set Normal ---*/
          
          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
          for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
          numerics->SetNormal(Normal);
          
          
          /*--- Compute residual ---*/
          
          numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
          
          /*--- Add Residuals and Jacobians ---*/
          
          LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
          if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
          
        }
      }
    }
  }
  
#else
  
  int rank, jProcessor;
  MPI_Status status;
  //MPI_Status send_stat[1], recv_stat[1];
  //MPI_Request send_req[1], recv_req[1];
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  bool compute;
  su2double *Buffer_Send_Psi = new su2double[nVar];
  su2double *Buffer_Receive_Psi = new su2double[nVar];
  
  /*--- Do the send process, by the moment we are sending each
   node individually, this must be changed ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    
    if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          /*--- We only send the information that belong to other boundary ---*/
          if (compute) {
            
            if (jProcessor != rank) {
              
              /*--- Copy the adjoint variable ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                Buffer_Send_Psi[iVar] = node[iPoint]->GetSolution(iVar);
              
              SU2_MPI::Bsend(Buffer_Send_Psi, nVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD);
              
              //          SU2_MPI::Isend(Buffer_Send_Psi, nVar, MPI_DOUBLE, jProcessor, iPoint, MPI_COMM_WORLD, &send_req[0]);
              
              /*--- Wait for this set of non-blocking comm. to complete ---*/
              
              //          SU2_MPI::Waitall(1, send_req, send_stat);
              
            }
            
          }
          
        }
      }
      
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          /*--- Find the associate pair to the original node ---*/
          
          jPoint = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[0];
          jProcessor = geometry->vertex[iMarker][iVertex]->GetPeriodicPointDomain()[1];
          
          if ((iPoint == jPoint) && (jProcessor == rank)) compute = false;
          else compute = true;
          
          if (compute) {
            
            /*--- We only receive the information that belong to other boundary ---*/
            
            if (jProcessor != rank) {
              
              SU2_MPI::Recv(Buffer_Receive_Psi, nVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &status);
              
              //          SU2_MPI::Irecv(Buffer_Receive_Psi, nVar, MPI_DOUBLE, jProcessor, jPoint, MPI_COMM_WORLD, &recv_req[0]);
              
              /*--- Wait for the this set of non-blocking recv's to complete ---*/
              
              //          SU2_MPI::Waitall(1, recv_req, recv_stat);
              
            }
            else {
              for (iVar = 0; iVar < nVar; iVar++)
                Buffer_Receive_Psi[iVar] = node[jPoint]->GetSolution(iVar);
            }
            
            /*--- Store the solution for both points ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Psi_i[iVar] = node[iPoint]->GetSolution(iVar);
              Psi_j[iVar] = Buffer_Receive_Psi[iVar];
            }
            
            /*--- If equivalent area or nearfield pressure condition ---*/
            
            if ((config->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
                (config->GetKind_ObjFunc() == NEARFIELD_PRESSURE)) {
              
              /*--- Identify the inner and the outer point (based on the normal direction) ---*/
              
              if (Normal[nDim-1] < 0.0) { Pin = iPoint; Pout = jPoint; }
              else { Pout = iPoint; Pin = jPoint; }
              
              IntBoundary_Jump = node[iPoint]->GetIntBoundary_Jump();
              
              /*--- Inner point ---*/
              
              if (iPoint == Pin) {
                for (iVar = 0; iVar < nVar; iVar++) {
                  Psi_in[iVar] = Psi_i[iVar]; Psi_out[iVar] = Psi_j[iVar];
                  MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
                  Psi_in_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_in[iVar] - IntBoundary_Jump[iVar];
                }
                numerics->SetAdjointVar(Psi_in, Psi_in_ghost);
              }
              
              /*--- Outer point ---*/
              
              if (iPoint == Pout) {
                for (iVar = 0; iVar < nVar; iVar++) {
                  Psi_in[iVar] = Psi_j[iVar]; Psi_out[iVar] = Psi_i[iVar];
                  MeanPsi[iVar] = 0.5*(Psi_out[iVar] + Psi_in[iVar]);
                  Psi_out_ghost[iVar] = 2.0*MeanPsi[iVar] - Psi_out[iVar] + IntBoundary_Jump[iVar];
                }
                numerics->SetAdjointVar(Psi_out, Psi_out_ghost);
              }
            }
            else {
              
              /*--- Just do a periodic BC ---*/
              
              numerics->SetAdjointVar(Psi_i, Psi_j);
              
            }
            
            /*--- Conservative variables w/o reconstruction (the same at both points) ---*/
            
            V_i = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
            V_j = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
            numerics->SetPrimitive(V_i, V_j);
            
            /*--- Set Normal ---*/
            
            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
            numerics->SetNormal(Normal);
            
            /*--- Compute residual ---*/
            
            numerics->ComputeResidual(Res_Conv_i, Res_Conv_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
            
            /*--- Add Residuals and Jacobians ---*/
            
            LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
            if (implicit) Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
            
          }
        }
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);

  delete[] Buffer_Send_Psi;
  delete[] Buffer_Receive_Psi;
  
#endif
  
  delete[] Normal;
  delete[] Psi_i;
  delete[] Psi_j;
  delete[] Psi_out;
  delete[] Psi_in;
  delete[] MeanPsi;
  delete[] Psi_out_ghost;
  delete[] Psi_in_ghost;
  
  
}

void CAdjIncEulerSolver::BC_Far_Field(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, Point_Normal;
  unsigned short iVar, iDim;
  su2double *Normal, *V_domain, *V_infty, *Psi_domain, *Psi_infty;
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_infty = new su2double[nVar];
  
  /*--- Loop over all the vertices ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- If the node belongs to the domain ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Index of the closest interior node ---*/
      
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Set the normal vector ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      /*--- Allocate the value at the infinity ---*/
      
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the farfield boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      conv_numerics->SetPrimitive(V_domain, V_infty);
      
      /*--- Adjoint flow solution at the wall ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
        Psi_infty[iVar] = 0.0;
      }
      conv_numerics->SetAdjointVar(Psi_domain, Psi_infty);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(), geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the upwind flux ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
      /*--- Viscous residual contribution, it doesn't work ---*/
      
      if (config->GetViscous()) {
        
        /*--- Points in edge, coordinates and normal vector---*/
        
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);
        
        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/
        
        visc_numerics->SetPrimitive(V_domain, V_infty);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_infty);
        
        /*--- Gradient and limiter of Adjoint Variables ---*/
        
        visc_numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
        
        /*--- Compute residual ---*/
        
        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);
        
        /*--- Update adjoint viscous residual ---*/
        
        LinSysRes.SubtractBlock(iPoint, Residual_i);
        
        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        
      }
      
    }
  }
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_infty;
}

void CAdjIncEulerSolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double Area;
  su2double *V_inlet, *V_domain, *Normal, *Psi_domain, *Psi_inlet;

  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement = config->GetGrid_Movement();
  
  Normal = new su2double[nDim];
  Psi_domain = new su2double[nVar]; Psi_inlet = new su2double[nVar];
  
  /*--- Loop over all the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check that the node belongs to the domain (i.e., not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);
      
      Area = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        Area += Normal[iDim]*Normal[iDim];
      Area = sqrt (Area);
      
      /*--- Allocate the value at the inlet ---*/
      
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);
      
      /*--- Retrieve solution at the boundary node ---*/
      
      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
      
      /*--- Adjoint flow solution at the boundary ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);
        
      /*--- Adjoint solution at the inlet ---*/

      Psi_inlet[0] = node[iPoint]->GetSolution(0);
      for (iDim = 0; iDim < nDim; iDim++)
        Psi_inlet[iDim+1] = 0.0;


      /*--- Set the flow and adjoint states in the solver ---*/
      
      conv_numerics->SetPrimitive(V_domain, V_inlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_inlet);
      
      /*--- Grid Movement ---*/
      
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());
      
      /*--- Compute the residual using an upwind scheme ---*/
      
      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
                                     Jacobian_ji, Jacobian_jj, config);
      
      /*--- Add and Subtract Residual ---*/
      
      LinSysRes.SubtractBlock(iPoint, Residual_i);
      
      /*--- Implicit contribution to the residual ---*/
      
      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      
      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {
        /*--- Index of the closest interior node ---*/

        Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

        /*--- Points in edge, coordinates and normal vector---*/

        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());
        visc_numerics->SetNormal(Normal);

        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_inlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_inlet);

        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_inlet;
  
}

void CAdjIncEulerSolver::BC_Outlet(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_outlet, *V_domain, *Psi_domain, *Psi_outlet, *Normal;

  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  Psi_domain = new su2double [nVar]; Psi_outlet = new su2double [nVar];
  Normal = new su2double[nDim];

  /*--- Loop over all the vertices ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- If the node belong to the domain ---*/

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Normal vector for this vertex (negate for outward convention) ---*/

      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Set the normal point ---*/

      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the outlet ---*/

      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the boundary node ---*/

      V_domain = solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();

      /*--- Adjoint flow solution at the boundary ---*/

      for (iVar = 0; iVar < nVar; iVar++)
        Psi_domain[iVar] = node[iPoint]->GetSolution(iVar);

      /*--- Imposed pressure and density ---*/

      V_outlet[0] = solver_container[FLOW_SOL]->GetPressure_Inf();


      /*--- Neumann condition for the velocity ---*/

      for (iDim = 0; iDim < nDim; iDim++)
        V_outlet[iDim+1] = node[Point_Normal]->GetSolution(iDim+1);

      /*--- Adjoint flow solution at the outlet (hard-coded for size[3] again?) ---*/

      Psi_outlet[2] = 0.0;
      su2double coeff = (2.0*V_domain[1])/ solver_container[FLOW_SOL]->node[Point_Normal]->GetBetaInc2();
      Psi_outlet[1] = node[Point_Normal]->GetSolution(1);
      Psi_outlet[0] = -coeff*Psi_outlet[1];

      /*--- Set the flow and adjoint states in the solver ---*/

      conv_numerics->SetPrimitive(V_domain, V_outlet);
      conv_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

      /*--- Grid Movement ---*/

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
            geometry->node[iPoint]->GetGridVel());

      conv_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij,
          Jacobian_ji, Jacobian_jj, config);

      /*--- Add and Subtract Residual ---*/

      LinSysRes.SubtractBlock(iPoint, Residual_i);

      /*--- Implicit contribution to the residual ---*/

      if (implicit)
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      /*--- Viscous residual contribution, it doesn't work ---*/

      if (config->GetViscous()) {

        /*--- Set laminar and eddy viscosity at the infinity ---*/
        V_outlet[nDim+3] = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
        V_outlet[nDim+4] = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity();

        /*--- Points in edge, coordinates and normal vector---*/
        visc_numerics->SetNormal(Normal);
        visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[Point_Normal]->GetCoord());


        /*--- Conservative variables w/o reconstruction and adjoint variables w/o reconstruction---*/

        visc_numerics->SetPrimitive(V_domain, V_outlet);
        visc_numerics->SetAdjointVar(Psi_domain, Psi_outlet);

        /*--- Turbulent kinetic energy ---*/
        if (config->GetKind_Turb_Model() == SST)
          visc_numerics->SetTurbKineticEnergy(solver_container[TURB_SOL]->node[iPoint]->GetSolution(0), solver_container[TURB_SOL]->node[iPoint]->GetSolution(0));


        /*--- Gradient and limiter of Adjoint Variables ---*/

        visc_numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());

        /*--- Compute residual ---*/

        visc_numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

        /*--- Update adjoint viscous residual ---*/

        LinSysRes.SubtractBlock(iPoint, Residual_i);

        if (implicit)
          Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);

      }
    }
  }
  
  /*--- Free locally allocated memory ---*/
  
  delete [] Normal;
  delete [] Psi_domain; delete [] Psi_outlet;
  
}

void CAdjIncEulerSolver::SetResidual_DualTime(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iRKStep,
                                           unsigned short iMesh, unsigned short RunTime_EqSystem) {
  unsigned short iVar, jVar;
  unsigned long iPoint;
  su2double *U_time_nM1, *U_time_n, *U_time_nP1, Volume_nM1, Volume_n, Volume_nP1, TimeStep;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool FlowEq = (RunTime_EqSystem == RUNTIME_FLOW_SYS);
  bool AdjEq = (RunTime_EqSystem == RUNTIME_ADJFLOW_SYS);
  bool Grid_Movement = config->GetGrid_Movement();
  
  /*--- loop over points ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Solution at time n-1, n and n+1 ---*/
    U_time_nM1 = node[iPoint]->GetSolution_time_n1();
    U_time_n   = node[iPoint]->GetSolution_time_n();
    U_time_nP1 = node[iPoint]->GetSolution();
    
    /*--- Volume at time n-1 and n ---*/
    if (Grid_Movement) {
      Volume_nM1 = geometry->node[iPoint]->GetVolume_nM1();
      Volume_n = geometry->node[iPoint]->GetVolume_n();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }
    else {
      Volume_nM1 = geometry->node[iPoint]->GetVolume();
      Volume_n = geometry->node[iPoint]->GetVolume();
      Volume_nP1 = geometry->node[iPoint]->GetVolume();
    }
    
    /*--- Time Step ---*/
    TimeStep = config->GetDelta_UnstTimeND();
    
    /*--- Compute Residual ---*/
    for (iVar = 0; iVar < nVar; iVar++) {
      if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
        Residual[iVar] = ( U_time_nP1[iVar]*Volume_nP1 - U_time_n[iVar]*Volume_n ) / TimeStep;
      if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
        Residual[iVar] = ( 3.0*U_time_nP1[iVar]*Volume_nP1 - 4.0*U_time_n[iVar]*Volume_n
                          +  1.0*U_time_nM1[iVar]*Volume_nM1 ) / (2.0*TimeStep);
    }
    
    if (FlowEq || AdjEq) Residual[0] = 0.0;
    
    /*--- Add Residual ---*/
    LinSysRes.AddBlock(iPoint, Residual);
    
    if (implicit) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++)
          Jacobian_i[iVar][jVar] = 0.0;
        
        if (config->GetUnsteady_Simulation() == DT_STEPPING_1ST)
          Jacobian_i[iVar][iVar] = Volume_nP1 / TimeStep;
        if (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)
          Jacobian_i[iVar][iVar] = (Volume_nP1*3.0)/(2.0*TimeStep);
      }
      if (FlowEq || AdjEq) Jacobian_i[0][0] = 0.0;
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
    }
  }
  
}

void CAdjIncEulerSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iDim, iVar, iMesh;
  unsigned long iPoint, index, iChildren, Point_Fine;
  su2double Area_Children, Area_Parent, *Coord, *Solution_Fine;
  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));
  bool time_stepping = config->GetUnsteady_Simulation() == TIME_STEPPING;
  string UnstExt, text_line, filename, restart_filename;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  /*--- Restart the solution from file information ---*/

  filename         = config->GetSolution_AdjFileName();
  restart_filename = config->GetObjFunc_Extension(filename);

  Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[MESH_0]->GetnDim();

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    restart_filename = config->GetMultizone_FileName(restart_filename, iZone);

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    restart_filename = config->GetUnsteady_FileName(restart_filename, val_iter);

  /*--- Read all lines in the restart file ---*/

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned long iPoint_Global_Local = 0;
  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

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
    if (rank == MASTER_NODE) {
      cout << endl << "The solution file " << restart_filename.data() << " doesn't match with the mesh file!" << endl;
      cout << "It could be empty lines at the end of the file." << endl << endl;
    }
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /*--- Communicate the loaded solution on the fine grid before we transfer
   it down to the coarse levels. We also call the preprocessing routine
   on the fine level in order to have all necessary quantities updated. ---*/

  solver[MESH_0][ADJFLOW_SOL]->Set_MPI_Solution(geometry[MESH_0], config);
  solver[MESH_0][ADJFLOW_SOL]->Preprocessing(geometry[MESH_0], solver[MESH_0], config, MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

  /*--- Interpolate the solution down to the coarse multigrid levels ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iPoint = 0; iPoint < geometry[iMesh]->GetnPoint(); iPoint++) {
      Area_Parent = geometry[iMesh]->node[iPoint]->GetVolume();
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;
      for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
        Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
        Area_Children = geometry[iMesh-1]->node[Point_Fine]->GetVolume();
        Solution_Fine = solver[iMesh-1][ADJFLOW_SOL]->node[Point_Fine]->GetSolution();
        for (iVar = 0; iVar < nVar; iVar++) {
          Solution[iVar] += Solution_Fine[iVar]*Area_Children/Area_Parent;
        }
      }
      solver[iMesh][ADJFLOW_SOL]->node[iPoint]->SetSolution(Solution);
    }
    solver[iMesh][ADJFLOW_SOL]->Set_MPI_Solution(geometry[iMesh], config);
    solver[iMesh][ADJFLOW_SOL]->Preprocessing(geometry[iMesh], solver[iMesh], config, iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
  }

  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

CAdjIncNSSolver::CAdjIncNSSolver(void) : CAdjIncEulerSolver() { }

CAdjIncNSSolver::CAdjIncNSSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh) : CAdjIncEulerSolver() {
  unsigned long iPoint, iVertex;
  string text_line, mesh_filename;
  unsigned short iDim, iVar, iMarker, nLineLets;
  ifstream restart_file;
  string filename, AdjExt;
  su2double Area=0.0, *Normal = NULL, myArea_Monitored;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Norm heat flux objective test ---*/
  pnorm = 1.0;
  if (config->GetKind_ObjFunc()==MAXIMUM_HEATFLUX)
    pnorm = 8.0; // Matches MaxNorm defined in solver_direct_mean.
  
  /*--- Set the gamma value ---*/
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Define geometry constants in the solver structure ---*/
  
  nDim         = geometry->GetnDim();
  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();
  
  nVar = nDim + 1;
  
  /*--- Initialize nVarGrad for deallocation ---*/
  
  nVarGrad = nVar;
  
  node = new CVariable*[nPoint];
  
  /*--- Define some auxiliary arrays related to the residual ---*/
  
  Point_Max    = new unsigned long[nVar]; for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]  = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }
  Residual     = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]     = 0.0;
  Residual_RMS = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar] = 0.0;
  Residual_Max = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar] = 0.0;
  Residual_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_i[iVar]   = 0.0;
  Residual_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Residual_j[iVar]   = 0.0;
  Res_Conv_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_i[iVar]   = 0.0;
  Res_Visc_i   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_i[iVar]   = 0.0;
  Res_Conv_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Conv_j[iVar]   = 0.0;
  Res_Visc_j   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Res_Visc_j[iVar]   = 0.0;
  
  /*--- Define some auxiliary arrays related to the solution ---*/
  
  Solution   = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 0.0;
  Solution_i = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_i[iVar] = 0.0;
  Solution_j = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) Solution_j[iVar] = 0.0;

  /*--- Define some auxiliary arrays related to the flow solution ---*/
  
  FlowPrimVar_i = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_i[iVar] = 0.0;
  FlowPrimVar_j = new su2double[nDim+7]; for (iVar = 0; iVar < nDim+7; iVar++) FlowPrimVar_j[iVar] = 0.0;

  /*--- Define some auxiliary vectors related to the geometry ---*/
  
  Vector   = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector[iDim]   = 0.0;
  Vector_i = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_i[iDim] = 0.0;
  Vector_j = new su2double[nDim]; for (iDim = 0; iDim < nDim; iDim++) Vector_j[iDim] = 0.0;
  
  /*--- Point to point Jacobians. These are always defined because
   they are also used for sensitivity calculations. ---*/
  
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar];
    Jacobian_j[iVar] = new su2double [nVar];
  }
  
  /*--- Solution and residual vectors ---*/
  
  LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
  LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  
  /*--- Jacobians and vector structures for implicit computations ---*/
  
  if (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT) {
    Jacobian_ii = new su2double*[nVar];
    Jacobian_ij = new su2double*[nVar];
    Jacobian_ji = new su2double*[nVar];
    Jacobian_jj = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_ii[iVar] = new su2double[nVar];
      Jacobian_ij[iVar] = new su2double[nVar];
      Jacobian_ji[iVar] = new su2double[nVar];
      Jacobian_jj[iVar] = new su2double[nVar];
    }
    if (rank == MASTER_NODE)
      cout << "Initialize Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar, true, geometry, config);
    
    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) cout << "Compute linelet structure. " << nLineLets << " elements in each line (average)." << endl;
    }
    
  } else {
    if (rank == MASTER_NODE)
      cout << "Explicit scheme. No Jacobian structure (Adjoint N-S). MG level: " << iMesh <<"." << endl;
  }
  
  /*--- Array structures for computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Smatrix[iDim] = new su2double [nDim];
    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Cvector[iVar] = new su2double [nDim];
  }
  
  /*--- Sensitivity definition and coefficient on all markers ---*/
  CSensitivity = new su2double* [nMarker];
  for (iMarker=0; iMarker<nMarker; iMarker++) {
    CSensitivity[iMarker] = new su2double [geometry->nVertex[iMarker]];
  }
  Sens_Geo   = new su2double[nMarker];
  Sens_Mach  = new su2double[nMarker];
  Sens_AoA   = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];
  Sens_BPress = new su2double[nMarker];
  
  /*--- Initialize sensitivities to zero ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Sens_Geo[iMarker]   = 0.0;
    Sens_Mach[iMarker]  = 0.0;
    Sens_AoA[iMarker]   = 0.0;
    Sens_Press[iMarker] = 0.0;
    Sens_Temp[iMarker]  = 0.0;
    Sens_BPress[iMarker] = 0.0;
    for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++)
      CSensitivity[iMarker][iVertex] = 0.0;
  }
  
  /*--- Initialize the adjoint variables to zero (infinity state) ---*/
  PsiRho_Inf = 0.0;
  if ((config->GetKind_ObjFunc() == TOTAL_HEATFLUX) ||
      (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX) ||
      (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX))
    PsiE_Inf = 1.0;
  else
    PsiE_Inf = 0.0;
  Phi_Inf = new su2double [nDim];
  Phi_Inf[0] = 0.0; Phi_Inf[1] = 0.0;
  if (nDim == 3) Phi_Inf[2] = 0.0;

  /*--- Initialize the solution to the far-field state everywhere. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CAdjIncNSVariable(PsiRho_Inf, Phi_Inf, PsiE_Inf, nDim, nVar, config);
  
  /*--- Calculate area monitored for area-averaged-outflow-quantity-based objectives ---*/
  myArea_Monitored = 0.0;
  if (config->GetKind_ObjFunc()==SURFACE_TOTAL_PRESSURE ||  config->GetKind_ObjFunc()==SURFACE_STATIC_PRESSURE){
    for (iMarker =0; iMarker < config->GetnMarker_All();  iMarker++){
      if (config->GetMarker_All_Monitoring(iMarker) == YES){
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            myArea_Monitored += sqrt (Area);
          }
        }
      }
    }
  }
#ifdef HAVE_MPI
  Area_Monitored = 0.0;
  SU2_MPI::Allreduce(&myArea_Monitored, &Area_Monitored, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  Area_Monitored = myArea_Monitored;
#endif

  /*--- MPI solution ---*/
  Set_MPI_Solution(geometry, config);
  
}

CAdjIncNSSolver::~CAdjIncNSSolver(void) {
  
}


void CAdjIncNSSolver::SetTime_Step(CGeometry *geometry, CSolver **solver_container, CConfig *config,
                            unsigned short iMesh, unsigned long Iteration) {

  /*--- Use the flow solution to update the time step
   *    The time step depends on the characteristic velocity, which is the same
   *    for the adjoint and flow solutions, albeit in the opposite direction. ---*/
  solver_container[FLOW_SOL]->SetTime_Step(geometry, solver_container, config, iMesh, Iteration);
}

void CAdjIncNSSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  
  unsigned long iPoint, ErrorCounter = 0;
  su2double SharpEdge_Distance;
  bool RightSol = true;
  
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Retrieve information about the spatial and temporal integration for the
   adjoint equations (note that the flow problem may use different methods). ---*/
  
  bool implicit       = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool limiter        = (config->GetKind_SlopeLimit_AdjFlow() != NO_LIMITER);
  bool center_jst     = (config->GetKind_Centered_AdjFlow() == JST);

  /*--- Residual initialization ---*/
  
  for (iPoint = 0; iPoint < nPoint; iPoint ++) {
    
    /*--- Get the distance form a sharp edge ---*/
    
    SharpEdge_Distance = geometry->node[iPoint]->GetSharpEdge_Distance();
    
    /*--- Initialize the non-physical points vector ---*/
    
    node[iPoint]->SetNon_Physical(false);
    
    /*--- Set the primitive variables incompressible and compressible
     adjoint variables ---*/
    
    RightSol = node[iPoint]->SetPrimVar(SharpEdge_Distance, false, config);
    if (!RightSol) { node[iPoint]->SetNon_Physical(true); ErrorCounter++; }
    
    /*--- Initialize the convective residual vector ---*/
    
    if (!Output) LinSysRes.SetBlock_Zero(iPoint);
    
  }
  
  /*--- Compute gradients adj for solution reconstruction and viscous term ---*/
  
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
  
  /*--- Limiter computation (upwind reconstruction) ---*/
  
  if (limiter && !Output) SetSolution_Limiter(geometry, config);

  /*--- Compute gradients adj for viscous term coupling ---*/

  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
    if (config->GetKind_Gradient_Method() == GREEN_GAUSS) solver_container[ADJTURB_SOL]->SetSolution_Gradient_GG(geometry, config);
    if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) solver_container[ADJTURB_SOL]->SetSolution_Gradient_LS(geometry, config);
  }
  
  /*--- Artificial dissipation for centered schemes ---*/
  
  if (center_jst && (iMesh == MESH_0) && !Output) {
    SetDissipation_Switch(geometry, config);
    SetUndivided_Laplacian(geometry, config);
  }
  
  /*--- Initialize the Jacobian for implicit integration ---*/
  
  if (implicit) Jacobian.SetValZero();
  
  /*--- Error message ---*/
  
  if (config->GetConsole_Output_Verb() == VERB_HIGH) {
#ifdef HAVE_MPI
    unsigned long MyErrorCounter = ErrorCounter; ErrorCounter = 0;
    SU2_MPI::Allreduce(&MyErrorCounter, &ErrorCounter, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (iMesh == MESH_0) config->SetNonphysical_Points(ErrorCounter);
  }
  
}

void CAdjIncNSSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                    CConfig *config, unsigned short iMesh, unsigned short iRKStep) {
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    
    /*--- Points in edge, coordinates and normal vector---*/
    
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), geometry->node[jPoint]->GetCoord());
    numerics->SetNormal(geometry->edge[iEdge]->GetNormal());
    
    /*--- Primitive variables w/o reconstruction and adjoint variables w/o reconstruction---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(),
                           solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive());
    
    numerics->SetAdjointVar(node[iPoint]->GetSolution(), node[jPoint]->GetSolution());
    
    /*--- Gradient and limiter of Adjoint Variables ---*/
    
    numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), node[jPoint]->GetGradient());
    
    /*--- Compute residual ---*/
    
    numerics->ComputeResidual(Residual_i, Residual_j, Jacobian_ii, Jacobian_ij, Jacobian_ji, Jacobian_jj, config);

    /*--- Update adjoint viscous residual ---*/
    
    LinSysRes.SubtractBlock(iPoint, Residual_i);
    LinSysRes.AddBlock(jPoint, Residual_j);
    
    if (implicit) {
      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      Jacobian.SubtractBlock(iPoint, jPoint, Jacobian_ij);
      Jacobian.AddBlock(jPoint, iPoint, Jacobian_ji);
      Jacobian.AddBlock(jPoint, jPoint, Jacobian_jj);
    }
    
  }
  
}

void CAdjIncNSSolver::Source_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CNumerics *second_numerics,
                                   CConfig *config, unsigned short iMesh) {
  
  unsigned long iPoint, jPoint, iEdge;
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool rotating_frame = config->GetRotating_Frame();
  
  /*--- Loop over all the points, note that we are supposing that primitive and
   adjoint gradients have been computed previously ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Primitive variables w/o reconstruction, and its gradient ---*/
    
    numerics->SetPrimitive(solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive(), NULL);
    
    numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(), NULL);

    /*--- Gradient of adjoint variables ---*/
    
    numerics->SetAdjointVarGradient(node[iPoint]->GetGradient(), NULL);

    /*--- Set volume ---*/
    
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());
    
    /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
    
    if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {
      
      /*--- Turbulent variables w/o reconstruction and its gradient ---*/
      
      numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(), NULL);
      
      numerics->SetTurbVarGradient(solver_container[TURB_SOL]->node[iPoint]->GetGradient(), NULL);
      
      /*--- Turbulent adjoint variables w/o reconstruction and its gradient ---*/
      
      numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->node[iPoint]->GetSolution(), NULL);
      
      numerics->SetTurbAdjointGradient(solver_container[ADJTURB_SOL]->node[iPoint]->GetGradient(), NULL);
      
      /*--- Set distance to the surface ---*/
      
      numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), 0.0);
      
    }
    
    /*--- Compute residual ---*/
    
    numerics->ComputeResidual(Residual, config);
    
    /*--- Add to the residual ---*/
    
    LinSysRes.AddBlock(iPoint, Residual);
    
  }
  
  /*--- If turbulence computation we must add some coupling terms to the NS adjoint eq. ---*/
  
  if ((config->GetKind_Solver() == ADJ_RANS) && (!config->GetFrozen_Visc_Cont())) {

    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

      /*--- Points in edge, and normal vector ---*/

      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      second_numerics->SetNormal(geometry->edge[iEdge]->GetNormal());

      /*--- Conservative variables w/o reconstruction ---*/

      second_numerics->SetConservative(solver_container[FLOW_SOL]->node[iPoint]->GetSolution(),
                                     solver_container[FLOW_SOL]->node[jPoint]->GetSolution());

      /*--- Gradient of primitive variables w/o reconstruction ---*/
      
      second_numerics->SetPrimVarGradient(solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive(),
                                        solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive());

      /*--- Viscosity ---*/

      second_numerics->SetLaminarViscosity(solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity(),
                                         solver_container[FLOW_SOL]->node[jPoint]->GetLaminarViscosity());

      /*--- Turbulent variables w/o reconstruction ---*/

      second_numerics->SetTurbVar(solver_container[TURB_SOL]->node[iPoint]->GetSolution(),
                                solver_container[TURB_SOL]->node[jPoint]->GetSolution());

      /*--- Turbulent adjoint variables w/o reconstruction ---*/

      second_numerics->SetTurbAdjointVar(solver_container[ADJTURB_SOL]->node[iPoint]->GetSolution(),
                                       solver_container[ADJTURB_SOL]->node[jPoint]->GetSolution());

      /*--- Set distance to the surface ---*/

      second_numerics->SetDistance(geometry->node[iPoint]->GetWall_Distance(), geometry->node[jPoint]->GetWall_Distance());
      
      /*--- Update adjoint viscous residual ---*/
      
      second_numerics->ComputeResidual(Residual, config);
      
      LinSysRes.AddBlock(iPoint, Residual);
      LinSysRes.SubtractBlock(jPoint, Residual);
    }

  }
  
  // WARNING: The rotating frame source term has been placed in the second
  // source term container since the section below is commented. This needs a
  // permanent fix asap!
  
  if (rotating_frame) {
    
    /*--- Loop over all points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      
      /*--- Load the adjoint variables ---*/
      second_numerics->SetAdjointVar(node[iPoint]->GetSolution(),
                                     node[iPoint]->GetSolution());
      
      /*--- Load the volume of the dual mesh cell ---*/
      second_numerics->SetVolume(geometry->node[iPoint]->GetVolume());
      
      /*--- Compute the adjoint rotating frame source residual ---*/
      second_numerics->ComputeResidual(Residual, Jacobian_i, config);
      
      /*--- Add the source residual to the total ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      
      /*--- Add the implicit Jacobian contribution ---*/
      if (implicit) Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      
    }
  }
  
}

void CAdjIncNSSolver::Viscous_Sensitivity(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics, CConfig *config) {
  
  unsigned long iVertex, iPoint;
  unsigned short iDim, jDim, iMarker;
  su2double **PsiVar_Grad = NULL, **PrimVar_Grad = NULL, *Normal = NULL, Area,
  sigma_partial, Laminar_Viscosity = 0.0,
  temp_sens = 0.0, *Psi = NULL, *U = NULL, Enthalpy, **GridVel_Grad, gradPsi5_v,
  psi5_tau_partial, psi5_tau_grad_vel, source_v_1, Density, Pressure = 0.0, div_vel, val_turb_ke,
  vartheta, vartheta_partial, psi5_p_div_vel, Omega[3], rho_v[3] = {0.0,0.0,0.0},
  CrossProduct[3], delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}, eps, scale = 1.0;
  su2double RefVel2, RefDensity, Mach2Vel, *Velocity_Inf, factor;

  su2double *USens = new su2double[nVar];
  su2double *UnitNormal = new su2double[nDim];
  su2double *normal_grad_vel = new su2double[nDim];
  su2double *tang_deriv_psi5 = new su2double[nDim];
  su2double *tang_deriv_T = new su2double[nDim];
  su2double **Sigma = new su2double* [nDim];
  
  for (iDim = 0; iDim < nDim; iDim++)
    Sigma[iDim] = new su2double [nDim];
  
  su2double *normal_grad_gridvel = new su2double[nDim];
  su2double *normal_grad_v_ux =new su2double[nDim];
  su2double **Sigma_Psi5v = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Sigma_Psi5v[iDim] = new su2double [nDim];
  su2double **tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  
  bool rotating_frame    = config->GetRotating_Frame();
  bool grid_movement     = config->GetGrid_Movement();
  su2double RefArea    = config->GetRefArea();
  su2double Mach_Motion     = config->GetMach_Motion();
  unsigned short ObjFunc = config->GetKind_ObjFunc();
  su2double Gas_Constant    = config->GetGas_ConstantND();
  
  if (config->GetSystemMeasurements() == US) scale = 1.0/12.0;
  else scale = 1.0;
  
  /*--- Compute non-dimensional factor. For dynamic meshes, use the motion Mach
   number as a reference value for computing the force coefficients.
   Otherwise, use the freestream values,
   which is the standard convention. ---*/
  
  if (grid_movement) {
    Mach2Vel = sqrt(Gamma*Gas_Constant*config->GetTemperature_FreeStreamND());
    RefVel2 = (Mach_Motion*Mach2Vel)*(Mach_Motion*Mach2Vel);
  }
  else {
    Velocity_Inf = config->GetVelocity_FreeStreamND();
    RefVel2 = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      RefVel2  += Velocity_Inf[iDim]*Velocity_Inf[iDim];
  }
  
  RefDensity  = config->GetDensity_FreeStreamND();
  
  factor = 1.0/(0.5*RefDensity*RefArea*RefVel2);
  
  if ((ObjFunc == INVERSE_DESIGN_HEATFLUX) ||
      (ObjFunc == TOTAL_HEATFLUX) || (ObjFunc == MAXIMUM_HEATFLUX) ||
      (ObjFunc == SURFACE_MASSFLOW) ) factor = 1.0;

 if ((ObjFunc == SURFACE_TOTAL_PRESSURE) || (ObjFunc == SURFACE_STATIC_PRESSURE)) factor = 1.0/Area_Monitored;


  /*--- Compute gradient of the grid velocity, if applicable ---*/
  
  if (grid_movement)
    SetGridVel_Gradient(geometry, config);
  
  Total_Sens_Geo = 0.0;
  Total_Sens_Mach = 0.0;
  Total_Sens_AoA = 0.0;
  Total_Sens_Press = 0.0;
  Total_Sens_Temp = 0.0;
  
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    
    Sens_Geo[iMarker] = 0.0;
    
    if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
        (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)) {
      
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
        
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        
        if (geometry->node[iPoint]->GetDomain()) {
          
          PsiVar_Grad = node[iPoint]->GetGradient();
          PrimVar_Grad = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
          
          Laminar_Viscosity = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();
          
          /*--- Compute face area and the unit normal to the surface ---*/
          
          Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
          Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) { Area += Normal[iDim]*Normal[iDim]; } Area = sqrt(Area);
          for (iDim = 0; iDim < nDim; iDim++) { UnitNormal[iDim] = Normal[iDim] / Area; }

          /*--- Incompressible case ---*/

          temp_sens = 0.0;
          
          /*--- Term: sigma_partial = \Sigma_{ji} n_i \partial_n v_j ---*/
          
          for (iDim = 0; iDim < nDim; iDim++) {
            for (jDim = 0; jDim < nDim; jDim++)
              Sigma[iDim][jDim] = Laminar_Viscosity * PsiVar_Grad[jDim+1][iDim];
          }

          
          for (iDim = 0; iDim < nDim; iDim++) {
            normal_grad_vel[iDim] = 0.0;
            for (jDim = 0; jDim < nDim; jDim++)
              normal_grad_vel[iDim] += PrimVar_Grad[iDim+1][jDim]*UnitNormal[jDim];
          }
          
          sigma_partial = 0.0;
          for (iDim = 0; iDim < nDim; iDim++)
            for (jDim = 0; jDim < nDim; jDim++)
              sigma_partial += UnitNormal[iDim]*Sigma[iDim][jDim]*normal_grad_vel[jDim];
          
          /*--- Compute additional terms in the surface sensitivity for
           moving walls in a rotating frame or dynamic mesh problem. ---*/
          
          if (grid_movement) {
            
            Psi = node[iPoint]->GetSolution();
            U = solver_container[FLOW_SOL]->node[iPoint]->GetSolution();
            Density = U[0];
            Pressure = solver_container[FLOW_SOL]->node[iPoint]->GetPressure();
            Enthalpy = solver_container[FLOW_SOL]->node[iPoint]->GetEnthalpy();
            
            /*--- Turbulent kinetic energy ---*/
            
            if (config->GetKind_Turb_Model() == SST)
              val_turb_ke = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
            else
              val_turb_ke = 0.0;
            
            div_vel = 0.0;
            for (iDim = 0 ; iDim < nDim; iDim++) {
              Velocity[iDim] = U[iDim+1]/Density;
              div_vel += PrimVar_Grad[iDim+1][iDim];
            }
            
            for (iDim = 0 ; iDim < nDim; iDim++)
              for (jDim = 0 ; jDim < nDim; jDim++)
                tau[iDim][jDim] = Laminar_Viscosity*(PrimVar_Grad[jDim+1][iDim] + PrimVar_Grad[iDim+1][jDim])
                - TWO3*Laminar_Viscosity*div_vel*delta[iDim][jDim]
                - TWO3*Density*val_turb_ke*delta[iDim][jDim];
            
            /*--- Form normal_grad_gridvel = \partial_n (u_omega) ---*/
            
            GridVel_Grad = geometry->node[iPoint]->GetGridVel_Grad();
            for (iDim = 0; iDim < nDim; iDim++) {
              normal_grad_gridvel[iDim] = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                normal_grad_gridvel[iDim] += GridVel_Grad[iDim][jDim]*UnitNormal[jDim];
            }
            
            /*--- Form normal_grad_v_ux = \partial_n (v - u_omega) ---*/
            
            for (iDim = 0; iDim < nDim; iDim++)
              normal_grad_v_ux[iDim] = normal_grad_vel[iDim] - normal_grad_gridvel[iDim];
            
            /*--- Form Sigma_Psi5v ---*/
            
            gradPsi5_v = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) {
              gradPsi5_v += PsiVar_Grad[nDim+1][iDim]*Velocity[iDim];
              for (jDim = 0; jDim < nDim; jDim++)
                Sigma_Psi5v[iDim][jDim] = Laminar_Viscosity * (PsiVar_Grad[nDim+1][iDim]*Velocity[jDim]+PsiVar_Grad[nDim+1][jDim]*Velocity[iDim]);
            }
            for (iDim = 0; iDim < nDim; iDim++)
              Sigma_Psi5v[iDim][iDim] -= TWO3*Laminar_Viscosity * gradPsi5_v;
            
            
            /*--- Now compute terms of the surface sensitivity ---*/
            
            /*--- Form vartheta_partial = \vartheta * \partial_n (v - u_x) . n ---*/
            vartheta = Density*Psi[0] + Density*Enthalpy*Psi[nDim+1];
            for (iDim = 0; iDim < nDim; iDim++) {
              vartheta += U[iDim+1]*Psi[iDim+1];
            }
            vartheta_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              vartheta_partial += vartheta * normal_grad_v_ux[iDim] * UnitNormal[iDim];
            
            /*--- Form sigma_partial = n_i ( \Sigma_Phi_{ij} + \Sigma_Psi5v_{ij} ) \partial_n (v - u_x)_j ---*/
            
            sigma_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                sigma_partial += UnitNormal[iDim]*(Sigma[iDim][jDim]+Sigma_Psi5v[iDim][jDim])*normal_grad_v_ux[jDim];
            
            /*--- Form psi5_tau_partial = \Psi_5 * \partial_n (v - u_x)_i * tau_{ij} * n_j ---*/
            
            psi5_tau_partial = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                psi5_tau_partial -= Psi[nDim+1]*normal_grad_v_ux[iDim]*tau[iDim][jDim]*UnitNormal[jDim];
            
            /*--- Form psi5_p_div_vel = ---*/
            
            psi5_p_div_vel = -Psi[nDim+1]*Pressure*div_vel;
            
            /*--- Form psi5_tau_grad_vel = \Psi_5 * tau_{ij} : \nabla( v ) ---*/
            
            psi5_tau_grad_vel = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              for (jDim = 0; jDim < nDim; jDim++)
                psi5_tau_grad_vel += Psi[nDim+1]*tau[iDim][jDim]*PrimVar_Grad[iDim+1][jDim];
            
            /*--- Retrieve the angular velocity vector ---*/
            
            source_v_1 = 0.0;
            if (rotating_frame) {
              
              Omega[0]  = (config->GetRotation_Rate_X(ZONE_0)/config->GetOmega_Ref());
              Omega[1]  = (config->GetRotation_Rate_Y(ZONE_0)/config->GetOmega_Ref());
              Omega[2]  = (config->GetRotation_Rate_Z(ZONE_0)/config->GetOmega_Ref());
              
              /*--- Calculate momentum source terms as: rho * ( Omega X V ) ---*/
              
              for (iDim = 0; iDim < nDim; iDim++)
                rho_v[iDim] = U[iDim+1];
              if (nDim == 2) rho_v[2] = 0.0;
              
              CrossProduct[0] = Omega[1]*rho_v[2] - Omega[2]*rho_v[1];
              CrossProduct[1] = Omega[2]*rho_v[0] - Omega[0]*rho_v[2];
              CrossProduct[2] = Omega[0]*rho_v[1] - Omega[1]*rho_v[0];
              
              
              for (iDim = 0; iDim < nDim; iDim++) {
                source_v_1 += Psi[iDim+1]*CrossProduct[iDim];
              }
            }
            
            /*--- For simplicity, store all additional terms within sigma_partial ---*/
            
            sigma_partial = sigma_partial + vartheta_partial + psi5_tau_partial + psi5_p_div_vel + psi5_tau_grad_vel + source_v_1;
            
          }
          
          /*--- Compute sensitivity for each surface point ---*/
          
          CSensitivity[iMarker][iVertex] = (sigma_partial - temp_sens) * Area * scale * factor;
            
          /*--- If sharp edge, set the sensitivity to 0 on that region ---*/
          
          if (config->GetSens_Remove_Sharp()) {
            eps = config->GetVenkat_LimiterCoeff()*config->GetRefElemLength();
            if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetAdjSharp_LimiterCoeff()*eps )
              CSensitivity[iMarker][iVertex] = 0.0;
          }
          
          Sens_Geo[iMarker] -= CSensitivity[iMarker][iVertex];
          
        }
      }
      
      Total_Sens_Geo += Sens_Geo[iMarker];
      
    }
  }

#ifdef HAVE_MPI
  
  su2double MyTotal_Sens_Geo   = Total_Sens_Geo;     Total_Sens_Geo = 0.0;
  su2double MyTotal_Sens_Mach  = Total_Sens_Mach;    Total_Sens_Mach = 0.0;
  su2double MyTotal_Sens_AoA   = Total_Sens_AoA;     Total_Sens_AoA = 0.0;
  su2double MyTotal_Sens_Press = Total_Sens_Press;   Total_Sens_Press = 0.0;
  su2double MyTotal_Sens_Temp  = Total_Sens_Temp;    Total_Sens_Temp = 0.0;
  
  SU2_MPI::Allreduce(&MyTotal_Sens_Geo, &Total_Sens_Geo, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Mach, &Total_Sens_Mach, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_AoA, &Total_Sens_AoA, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MyTotal_Sens_Temp, &Total_Sens_Temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
#endif
  
  delete [] USens;
  delete [] UnitNormal;
  delete [] normal_grad_vel;
  delete [] tang_deriv_psi5;
  delete [] tang_deriv_T;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Sigma[iDim];
  delete [] Sigma;
  delete [] normal_grad_gridvel;
  delete [] normal_grad_v_ux;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Sigma_Psi5v[iDim];
  delete [] Sigma_Psi5v;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] tau[iDim];
  delete [] tau;
  delete [] Velocity;
  
}

void CAdjIncNSSolver::BC_HeatFlux_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned short iDim, iVar, jVar;
  unsigned long iVertex, iPoint, total_index;
  
  su2double *d, l1psi, phi[] = {0.0,0.0,0.0}, *GridVel;
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  
  su2double *Psi = new su2double[nVar];
  su2double **Tau = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  su2double *Normal = new su2double[nDim];
  su2double *Edge_Vector = new su2double[nDim];
  su2double **GradPhi = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  su2double *GradPsiE = new su2double [nDim];
  
  /*--- Loop over all of the vertices on this boundary marker ---*/
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0; Res_Visc_i[iVar] = 0.0;
        if (implicit) { for (jVar = 0; jVar < nVar; jVar ++) Jacobian_ii[iVar][jVar] = 0.0; }
      }
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
      
      for (iVar = 0; iVar < nVar; iVar++)
        Psi[iVar] = node[iPoint]->GetSolution(iVar);
      
      /*--- Get the force projection vector (based on the objective function) ---*/
      
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Set the adjoint velocity BC ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) { phi[iDim] = d[iDim]; }
      
      /*--- Correct the adjoint velocity BC for dynamic meshes ---*/
      
      if (grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++)
          phi[iDim] -= Psi[nDim+1]*GridVel[iDim];
      }
      
      /*--- Impose the value of the adjoint velocity as a strong boundary
       condition (Dirichlet). Fix the adjoint velocity and remove any addtional
       contribution to the residual at this node. ---*/
      
      for (iDim = 0; iDim < nDim; iDim++)
        node[iPoint]->SetSolution_Old(iDim+1, phi[iDim]);
      
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();

      /*--- Pressure residual due to the convective term ---*/

      l1psi = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        l1psi += Normal[iDim]*d[iDim];
      Res_Conv_i[0] = l1psi;

      /*--- Convective contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Conv_i);
      
      /*--- Viscous contribution to the residual at the wall ---*/
      
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      
      /*--- Enforce the no-slip boundary condition in a strong way by
       modifying the velocity-rows of the Jacobian (1 on the diagonal). ---*/
      
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
    }
    
  }
  
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Tau[iDim];
  delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] Edge_Vector;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  
}


void CAdjIncNSSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics, CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  
  unsigned long iVertex, iPoint, total_index;
  unsigned short iDim, iVar, jVar;
  su2double *d, q, l1psi, phi[3] = {0.0,0.0,0.0};
  su2double *GridVel;
  su2double Laminar_Viscosity;
  su2double kGTdotn=0.0, Area=0.0, Xi=0.0;
  
  su2double *Psi = new su2double[nVar];
  su2double **Tau = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Tau[iDim] = new su2double [nDim];
  su2double *Velocity = new su2double[nDim];
  su2double *Normal = new su2double[nDim];
  
  su2double **GradPhi = new su2double* [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    GradPhi[iDim] = new su2double [nDim];
  su2double *GradPsiE = new su2double [nDim];
  su2double *GradT;// = new su2double[nDim];
  su2double *dPoRho2 = new su2double[nDim];
  
  bool implicit = (config->GetKind_TimeIntScheme_AdjFlow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();
  bool heat_flux_obj  = ((config->GetKind_ObjFunc() == TOTAL_HEATFLUX) ||
                         (config->GetKind_ObjFunc() == MAXIMUM_HEATFLUX) ||
                         (config->GetKind_ObjFunc() == INVERSE_DESIGN_HEATFLUX));
  
  su2double Prandtl_Lam  = config->GetPrandtl_Lam();
  su2double Gas_Constant = config->GetGas_ConstantND();
  su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();
    
    if (geometry->node[iPoint]->GetDomain()) {
      
      /*--- Initialize the convective & viscous residuals to zero ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Res_Conv_i[iVar] = 0.0;
        Res_Visc_i[iVar] = 0.0;
        if (implicit) {
          for (jVar = 0; jVar < nVar; jVar ++)
            Jacobian_ii[iVar][jVar] = 0.0;
        }
      }
      
      /*--- Retrieve adjoint solution at the wall boundary node ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Psi[iVar] = node[iPoint]->GetSolution(iVar);
      
      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];
      
      /*--- Get the force projection vector (based on the objective function) ---*/
      d = node[iPoint]->GetForceProj_Vector();
      
      /*--- Adjustments to strong boundary condition for dynamic meshes ---*/
      if ( grid_movement) {
        GridVel = geometry->node[iPoint]->GetGridVel();
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim] - Psi[nVar-1]*GridVel[iDim];
        }
      } else {
        for (iDim = 0; iDim < nDim; iDim++) {
          phi[iDim] = d[iDim];
        }
      }
      
      /*--- Strong BC imposition for the adjoint velocity equations ---*/
      for (iDim = 0; iDim < nDim; iDim++)
        LinSysRes.SetBlock_Zero(iPoint, iDim+1);
      node[iPoint]->SetVel_ResTruncError_Zero();
      for (iDim = 0; iDim < nDim; iDim++)
        node[iPoint]->SetSolution_Old(iDim+1, phi[iDim]);
      if (implicit) {
        for (iVar = 1; iVar <= nDim; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }
      }
      
      /*--- Get transport coefficient information ---*/
      Laminar_Viscosity    = solver_container[FLOW_SOL]->node[iPoint]->GetLaminarViscosity();

//      GradV = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      
      /*--- Calculate Dirichlet condition for energy equation ---*/
      if (!heat_flux_obj) {
        q = 0.0;
      }
      else {

        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) Area += Normal[iDim]*Normal[iDim];
        Area = sqrt(Area);

        /* --- Temperature gradient term ---*/
        GradT = solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive()[0];
        kGTdotn = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          kGTdotn += Cp * Laminar_Viscosity/Prandtl_Lam*GradT[iDim]*Normal[iDim]/Area;
        // Cp * Viscosity/Prandtl_Lam matches term used in solver_direct_mean
        /*--- constant term to multiply max heat flux objective ---*/
        Xi = solver_container[FLOW_SOL]->GetTotal_HeatFlux(); // versions for max heat flux
        Xi = pow(Xi, 1.0/pnorm-1.0)/pnorm;

        /*--- Boundary condition value ---*/
        q = Xi * pnorm * pow(kGTdotn, pnorm-1.0)*Area;
      }
      
      /*--- Strong BC enforcement of the energy equation ---*/
      LinSysRes.SetBlock_Zero(iPoint, nVar-1);
      node[iPoint]->SetEnergy_ResTruncError_Zero();
      node[iPoint]->SetSolution_Old(nDim+1, q);
      if (implicit) {
        iVar = nDim+1;
        total_index = iPoint*nVar+iVar;
        Jacobian.DeleteValsRowi(total_index);
      }

      /*--- Pressure residual due to the convective term ---*/
      l1psi = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        l1psi += Normal[iDim]*d[iDim];
      Res_Conv_i[0] = l1psi;
      
      /*--- Update convective and viscous residuals ---*/
      LinSysRes.AddBlock(iPoint, Res_Conv_i);
      LinSysRes.SubtractBlock(iPoint, Res_Visc_i);
      if (implicit) {
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_ii);
      }
      
    }
    
  }
  
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Tau[iDim];
  delete [] Tau;
  delete [] Psi;
  delete [] Velocity;
  delete [] Normal;
  delete [] GradPsiE;
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] GradPhi[iDim];
  delete [] GradPhi;
  delete [] dPoRho2;
}

