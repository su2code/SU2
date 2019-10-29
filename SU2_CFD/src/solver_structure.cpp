/*!
 * \file solver_structure.cpp
 * \brief Main subroutines for solving primal and adjoint problems.
 * \author F. Palacios, T. Economon
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
#include "../include/variables/CBaselineVariable.hpp"
#include "../../Common/include/toolboxes/MMS/CIncTGVSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CInviscidVortexSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSIncEulerSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSIncNSSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSNSUnitQuadSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.hpp"
#include "../../Common/include/toolboxes/MMS/CNSUnitQuadSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CRinglebSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CTGVSolution.hpp"
#include "../../Common/include/toolboxes/MMS/CUserDefinedSolution.hpp"
#include "../../Common/include/toolboxes/printing_toolbox.hpp"


CSolver::CSolver(bool mesh_deform_mode) : System(mesh_deform_mode) {

  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
  
  adjoint = false;

  /*--- Set the multigrid level to the finest grid. This can be
        overwritten in the constructors of the derived classes. ---*/
  MGLevel = MESH_0;
  
  /*--- Array initialization ---*/
  
  OutputHeadingNames = NULL;
  Residual_RMS       = NULL;
  Residual_Max       = NULL;
  Residual_BGS       = NULL;
  Residual_Max_BGS   = NULL;
  Residual           = NULL;
  Residual_i         = NULL;
  Residual_j         = NULL;
  Point_Max          = NULL;
  Point_Max_Coord    = NULL;
  Point_Max_BGS      = NULL;
  Point_Max_Coord_BGS = NULL;
  Solution           = NULL;
  Solution_i         = NULL;
  Solution_j         = NULL;
  Vector             = NULL;
  Vector_i           = NULL;
  Vector_j           = NULL;
  Res_Conv           = NULL;
  Res_Visc           = NULL;
  Res_Sour           = NULL;
  Res_Conv_i         = NULL;
  Res_Visc_i         = NULL;
  Res_Conv_j         = NULL;
  Res_Visc_j         = NULL;
  Jacobian_i         = NULL;
  Jacobian_j         = NULL;
  Jacobian_ii        = NULL;
  Jacobian_ij        = NULL;
  Jacobian_ji        = NULL;
  Jacobian_jj        = NULL;
  iPoint_UndLapl     = NULL;
  jPoint_UndLapl     = NULL;
  Smatrix            = NULL;
  Cvector            = NULL;
  Restart_Vars       = NULL;
  Restart_Data       = NULL;
  base_nodes         = nullptr;
  nOutputVariables   = 0;
  valResidual        = 0.0;
  

  /*--- Inlet profile data structures. ---*/

  nRowCum_InletFile = NULL;
  nRow_InletFile    = NULL;
  nCol_InletFile    = NULL;
  Inlet_Data        = NULL;

  /*--- Variable initialization to avoid valgrid warnings when not used. ---*/
  
  IterLinSolver = 0;

  /*--- Initialize pointer for any verification solution. ---*/
  VerificationSolution  = NULL;
  
  /*--- Flags for the periodic BC communications. ---*/
  
  rotate_periodic   = false;
  implicit_periodic = false;

  /*--- Containers to store the markers. ---*/
  nMarker = 0;
  nVertex = nullptr;

  /*--- Flags for the dynamic grid (rigid movement or unsteady deformation). ---*/
  dynamic_grid = false;

  /*--- Container to store the vertex tractions. ---*/
  VertexTraction = NULL;
  VertexTractionAdjoint = NULL;

  
  nPrimVarGrad = 0;
  nPrimVar     = 0;
  
}

CSolver::~CSolver(void) {

  unsigned short iVar, iDim;
  unsigned long iMarker, iVertex;

  /*--- Public variables, may be accessible outside ---*/

  if ( OutputHeadingNames != NULL) {
    delete [] OutputHeadingNames;
  }

  /*--- Private ---*/

  if (Residual_RMS != NULL) delete [] Residual_RMS;
  if (Residual_Max != NULL) delete [] Residual_Max;
  if (Residual != NULL) delete [] Residual;
  if (Residual_i != NULL) delete [] Residual_i;
  if (Residual_j != NULL) delete [] Residual_j;
  if (Point_Max != NULL) delete [] Point_Max;

  if (Residual_BGS != NULL) delete [] Residual_BGS;
  if (Residual_Max_BGS != NULL) delete [] Residual_Max_BGS;
  if (Point_Max_BGS != NULL) delete [] Point_Max_BGS;

  if (Point_Max_Coord != NULL) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Point_Max_Coord[iVar];
    }
    delete [] Point_Max_Coord;
  }

  if (Point_Max_Coord_BGS != NULL) {
    for (iVar = 0; iVar < nVar; iVar++) {
      delete [] Point_Max_Coord_BGS[iVar];
    }
    delete [] Point_Max_Coord_BGS;
  }

  if (Solution != NULL) delete [] Solution;
  if (Solution_i != NULL) delete [] Solution_i;
  if (Solution_j != NULL) delete [] Solution_j;
  if (Vector != NULL) delete [] Vector;
  if (Vector_i != NULL) delete [] Vector_i;
  if (Vector_j != NULL) delete [] Vector_j;
  if (Res_Conv != NULL) delete [] Res_Conv;
  if (Res_Visc != NULL) delete [] Res_Visc;
  if (Res_Sour != NULL) delete [] Res_Sour;
  if (Res_Conv_i != NULL) delete [] Res_Conv_i;
  if (Res_Visc_i != NULL) delete [] Res_Visc_i;
  if (Res_Visc_j != NULL) delete [] Res_Visc_j;

  if (iPoint_UndLapl != NULL) delete [] iPoint_UndLapl;
  if (jPoint_UndLapl != NULL) delete [] jPoint_UndLapl;

  if (Jacobian_i != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_i[iVar];
    delete [] Jacobian_i;
  }

  if (Jacobian_j != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_j[iVar];
    delete [] Jacobian_j;
  }

  if (Jacobian_ii != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ii[iVar];
    delete [] Jacobian_ii;
  }

  if (Jacobian_ij != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ij[iVar];
    delete [] Jacobian_ij;
  }

  if (Jacobian_ji != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_ji[iVar];
    delete [] Jacobian_ji;
  }

  if (Jacobian_jj != NULL) {
    for (iVar = 0; iVar < nVar; iVar++)
      delete [] Jacobian_jj[iVar];
    delete [] Jacobian_jj;
  }

  if (Smatrix != NULL) {
    for (iDim = 0; iDim < nDim; iDim++)
      delete [] Smatrix[iDim];
    delete [] Smatrix;
  }

  if (Cvector != NULL) {
    for (iVar = 0; iVar < nVarGrad; iVar++)
      delete [] Cvector[iVar];
    delete [] Cvector;
  }

  if (VertexTraction != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] VertexTraction[iMarker][iVertex];
      delete [] VertexTraction[iMarker];
    }
    delete [] VertexTraction;
  }

  if (VertexTractionAdjoint != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        delete [] VertexTractionAdjoint[iMarker][iVertex];
      delete [] VertexTractionAdjoint[iMarker];
    }
    delete [] VertexTractionAdjoint;
  }

  if (nVertex != nullptr) delete [] nVertex;

  if (Restart_Vars != NULL) {delete [] Restart_Vars; Restart_Vars = NULL;}
  if (Restart_Data != NULL) {delete [] Restart_Data; Restart_Data = NULL;}

  if (nRowCum_InletFile != NULL) {delete [] nRowCum_InletFile; nRowCum_InletFile = NULL;}
  if (nRow_InletFile    != NULL) {delete [] nRow_InletFile;    nRow_InletFile    = NULL;}
  if (nCol_InletFile    != NULL) {delete [] nCol_InletFile;    nCol_InletFile    = NULL;}
  if (Inlet_Data        != NULL) {delete [] Inlet_Data;        Inlet_Data        = NULL;}

  if (VerificationSolution != NULL) {delete VerificationSolution; VerificationSolution = NULL;}
  
}

void CSolver::InitiatePeriodicComms(CGeometry *geometry,
                                    CConfig *config,
                                    unsigned short val_periodic_index,
                                    unsigned short commType) {
  
  /*--- Local variables ---*/
  
  bool boundary_i, boundary_j;
  
  unsigned short iVar, jVar, iDim;
  unsigned short iNeighbor, nNeighbor = 0;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;
  unsigned short ICOUNT          = nVar;
  unsigned short JCOUNT          = nVar;
  
  int iMessage, iSend, nSend;

  unsigned long iPoint, jPoint, msg_offset, buf_offset, iPeriodic, Neighbor_Point;
  
  su2double *Diff      = new su2double[nVar];
  su2double *Und_Lapl  = new su2double[nVar];
  su2double *Sol_Min   = new su2double[nPrimVarGrad];
  su2double *Sol_Max   = new su2double[nPrimVarGrad];
  su2double *rotPrim_i = new su2double[nPrimVar];
  su2double *rotPrim_j = new su2double[nPrimVar];
  
  su2double Sensor_i = 0.0, Sensor_j = 0.0, Pressure_i, Pressure_j;
  su2double *Coord_i, *Coord_j, r11, r12, r13, r22, r23_a, r23_b, r33, weight;
  su2double *center, *angles, translation[3]={0.0,0.0,0.0}, *trans, dx, dy, dz;
  su2double rotMatrix[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double Theta, Phi, Psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  su2double rotCoord_i[3] = {0.0, 0.0, 0.0}, rotCoord_j[3] = {0.0, 0.0, 0.0};
  
  string Marker_Tag;
  
  /*--- Set the size of the data packet and type depending on quantity. ---*/
  
  switch (commType) {
    case PERIODIC_VOLUME:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_NEIGHBORS:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_UNSIGNED_SHORT;
      break;
    case PERIODIC_RESIDUAL:
      COUNT_PER_POINT  = nVar + nVar*nVar + 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_IMPLICIT:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LAPLACIAN:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_MAX_EIG:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_SENSOR:
      COUNT_PER_POINT  = 2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_SOL_GG:
      COUNT_PER_POINT  = nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nVar;
      JCOUNT           = nDim;
      break;
    case PERIODIC_PRIM_GG:
      COUNT_PER_POINT  = nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      ICOUNT           = nPrimVarGrad;
      JCOUNT           = nDim;
      break;
    case PERIODIC_SOL_LS:
      COUNT_PER_POINT  = nDim*nDim + nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_PRIM_LS:
      COUNT_PER_POINT  = nDim*nDim + nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LIM_PRIM_1:
      COUNT_PER_POINT  = nPrimVarGrad*2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LIM_PRIM_2:
      COUNT_PER_POINT  = nPrimVarGrad;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LIM_SOL_1:
      COUNT_PER_POINT  = nVar*2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_LIM_SOL_2:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                     CURRENT_FUNCTION);
      break;
  }
  
  su2double **jacBlock = new su2double*[ICOUNT];
  su2double **rotBlock = new su2double*[ICOUNT];
  for (iVar = 0; iVar < ICOUNT; iVar++) {
    jacBlock[iVar] = new su2double[JCOUNT];
    rotBlock[iVar] = new su2double[JCOUNT];
  }
  
  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. It will be reallocated whenever
   we find a larger count per point than currently exists. After the
   first cycle of comms, this should be inactive. ---*/
  
  if (COUNT_PER_POINT > geometry->countPerPeriodicPoint) {
    geometry->AllocatePeriodicComms(COUNT_PER_POINT);
  }
  
  /*--- Set some local pointers to make access simpler. ---*/
  
  su2double *bufDSend = geometry->bufD_PeriodicSend;
  
  unsigned short *bufSSend = geometry->bufS_PeriodicSend;
  
  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/
  
  if (geometry->nPeriodicSend > 0) {
    
    /*--- Post all non-blocking recvs first before sends. ---*/
    
    geometry->PostPeriodicRecvs(geometry, config, MPI_TYPE);
    
    for (iMessage = 0; iMessage < geometry->nPeriodicSend; iMessage++) {
      
      /*--- Get the offset in the buffer for the start of this message. ---*/
      
      msg_offset = geometry->nPoint_PeriodicSend[iMessage];
      
      /*--- Get the number of periodic points we need to
       communicate on the current periodic marker. ---*/
      
      nSend = (geometry->nPoint_PeriodicSend[iMessage+1] -
               geometry->nPoint_PeriodicSend[iMessage]);
      
      for (iSend = 0; iSend < nSend; iSend++) {
        
        /*--- Get the local index for this communicated data. We need
         both the node and periodic face index (for rotations). ---*/
        
        iPoint    = geometry->Local_Point_PeriodicSend[msg_offset  + iSend];
        iPeriodic = geometry->Local_Marker_PeriodicSend[msg_offset + iSend];
        
        /*--- Retrieve the supplied periodic information. ---*/
        
        Marker_Tag = config->GetMarker_All_TagBound(iPeriodic);
        center     = config->GetPeriodicRotCenter(Marker_Tag);
        angles     = config->GetPeriodicRotAngles(Marker_Tag);
        trans      = config->GetPeriodicTranslation(Marker_Tag);
        
        /*--- Store (center+trans) as it is constant and will be added. ---*/
        
        translation[0] = center[0] + trans[0];
        translation[1] = center[1] + trans[1];
        translation[2] = center[2] + trans[2];
        
        /*--- Store angles separately for clarity. Compute sines/cosines. ---*/
        
        Theta    = angles[0];      Phi = angles[1];     Psi = angles[2];
        cosTheta = cos(Theta);  cosPhi = cos(Phi);   cosPsi = cos(Psi);
        sinTheta = sin(Theta);  sinPhi = sin(Phi);   sinPsi = sin(Psi);
        
        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. ---*/
        
        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = cosPhi*sinPsi;
        rotMatrix[2][0] = -sinPhi;
        
        rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = sinTheta*cosPhi;
        
        rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[2][2] = cosTheta*cosPhi;
        
        /*--- Compute the offset in the recv buffer for this point. ---*/
        
        buf_offset = (msg_offset + iSend)*geometry->countPerPeriodicPoint;
        
        /*--- Load the send buffers depending on the particular value
         that has been requested for communication. ---*/
        
        switch (commType) {
            
          case PERIODIC_VOLUME:
            
            /*--- Load the volume of the current periodic CV so that
             we can accumulate the total control volume size on all
             periodic faces. ---*/
            
            bufDSend[buf_offset] = geometry->node[iPoint]->GetVolume() +
            geometry->node[iPoint]->GetPeriodicVolume();
            
            break;
            
          case PERIODIC_NEIGHBORS:
            
            nNeighbor = 0;
            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              Neighbor_Point = geometry->node[iPoint]->GetPoint(iNeighbor);
              
              /*--- Check if this neighbor lies on the periodic face so
               that we avoid double counting neighbors on both sides. If
               not, increment the count of neighbors for the donor. ---*/
              
              if (!geometry->node[Neighbor_Point]->GetPeriodicBoundary())
              nNeighbor++;
              
            }
            
            /*--- Store the number of neighbors in bufffer. ---*/
            
            bufSSend[buf_offset] = nNeighbor;
            
            break;
            
          case PERIODIC_RESIDUAL:
            
            /*--- Communicate the residual from our partial control
             volume to the other side of the periodic face. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = LinSysRes.GetBlock(iPoint, iVar);
            }
            
            /*--- Rotate the momentum components of the residual array. ---*/
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*LinSysRes.GetBlock(iPoint, 1) +
                                          rotMatrix[0][1]*LinSysRes.GetBlock(iPoint, 2));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*LinSysRes.GetBlock(iPoint, 1) +
                                          rotMatrix[1][1]*LinSysRes.GetBlock(iPoint, 2));
              } else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*LinSysRes.GetBlock(iPoint, 1) +
                                          rotMatrix[0][1]*LinSysRes.GetBlock(iPoint, 2) +
                                          rotMatrix[0][2]*LinSysRes.GetBlock(iPoint, 3));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*LinSysRes.GetBlock(iPoint, 1) +
                                          rotMatrix[1][1]*LinSysRes.GetBlock(iPoint, 2) +
                                          rotMatrix[1][2]*LinSysRes.GetBlock(iPoint, 3));
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*LinSysRes.GetBlock(iPoint, 1) +
                                          rotMatrix[2][1]*LinSysRes.GetBlock(iPoint, 2) +
                                          rotMatrix[2][2]*LinSysRes.GetBlock(iPoint, 3));
              }
            }
            buf_offset += nVar;
            
            /*--- Load the time step for the current point. ---*/
            
            bufDSend[buf_offset] = base_nodes->GetDelta_Time(iPoint);
            buf_offset++;
            
            /*--- For implicit calculations, we will communicate the
             contributions to the Jacobian block diagonal, i.e., the
             impact of the point upon itself, J_ii. ---*/
            
            if (implicit_periodic) {
              
              for (iVar = 0; iVar < nVar; iVar++) {
                for (jVar = 0; jVar < nVar; jVar++) {
                  jacBlock[iVar][jVar] = Jacobian.GetBlock(iPoint, iPoint, iVar, jVar);
                }
              }
              
              /*--- Rotate the momentum columns of the Jacobian. ---*/
              
              if (rotate_periodic) {
                for (iVar = 0; iVar < nVar; iVar++) {
                  if (nDim == 2) {
                    jacBlock[1][iVar] = (rotMatrix[0][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix[0][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar));
                    jacBlock[2][iVar] = (rotMatrix[1][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix[1][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar));
                  } else {
                    
                    jacBlock[1][iVar] = (rotMatrix[0][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix[0][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix[0][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                    jacBlock[2][iVar] = (rotMatrix[1][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix[1][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix[1][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                    jacBlock[3][iVar] = (rotMatrix[2][0]*Jacobian.GetBlock(iPoint, iPoint, 1, iVar) +
                                         rotMatrix[2][1]*Jacobian.GetBlock(iPoint, iPoint, 2, iVar) +
                                         rotMatrix[2][2]*Jacobian.GetBlock(iPoint, iPoint, 3, iVar));
                  }
                }
              }
              
              /*--- Load the Jacobian terms into the buffer for sending. ---*/
              
              for (iVar = 0; iVar < nVar; iVar++) {
                for (jVar = 0; jVar < nVar; jVar++) {
                  bufDSend[buf_offset] = jacBlock[iVar][jVar];
                  buf_offset++;
                }
              }
            }
            
            break;
            
          case PERIODIC_IMPLICIT:
            
            /*--- Communicate the solution from our master set of periodic
             nodes (from the linear solver perspective) to the passive
             periodic nodes on the matching face. This is done at the
             end of the iteration to synchronize the solution after the
             linear solve. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            }
            
            /*--- Rotate the momentum components of the solution array. ---*/
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetSolution(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetSolution(iPoint,2));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetSolution(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetSolution(iPoint,2));
              } else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetSolution(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetSolution(iPoint,2) +
                                          rotMatrix[0][2]*base_nodes->GetSolution(iPoint,3));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetSolution(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetSolution(iPoint,2) +
                                          rotMatrix[1][2]*base_nodes->GetSolution(iPoint,3));
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*base_nodes->GetSolution(iPoint,1) +
                                          rotMatrix[2][1]*base_nodes->GetSolution(iPoint,2) +
                                          rotMatrix[2][2]*base_nodes->GetSolution(iPoint,3));
              }
            }
            
            break;
            
          case PERIODIC_LAPLACIAN:
            
            /*--- For JST, the undivided Laplacian must be computed
             consistently by using the complete control volume info
             from both sides of the periodic face. ---*/
            
            for (iVar = 0; iVar< nVar; iVar++)
            Und_Lapl[iVar] = 0.0;
            
            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              
              /*--- Avoid periodic boundary points so that we do not
               duplicate edges on both sides of the periodic BC. ---*/
              
              if (!geometry->node[jPoint]->GetPeriodicBoundary()) {
                
                /*--- Solution differences ---*/
                
                for (iVar = 0; iVar < nVar; iVar++)
                Diff[iVar] = (base_nodes->GetSolution(iPoint, iVar) -
                              base_nodes->GetSolution(jPoint,iVar));
                
                /*--- Correction for compressible flows (use enthalpy) ---*/
                
                if (!(config->GetKind_Regime() == INCOMPRESSIBLE)) {
                  Pressure_i   = base_nodes->GetPressure(iPoint);
                  Pressure_j   = base_nodes->GetPressure(jPoint);
                  Diff[nVar-1] = ((base_nodes->GetSolution(iPoint,nVar-1) + Pressure_i) -
                                  (base_nodes->GetSolution(jPoint,nVar-1) + Pressure_j));
                }
                
                boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
                boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
                
                /*--- Both points inside the domain, or both in the boundary ---*/
                
                if ((!boundary_i && !boundary_j) ||
                    ( boundary_i &&  boundary_j)) {
                  if (geometry->node[iPoint]->GetDomain()) {
                    for (iVar = 0; iVar< nVar; iVar++)
                    Und_Lapl[iVar] -= Diff[iVar];
                  }
                }
                
                /*--- iPoint inside the domain, jPoint on the boundary ---*/
                
                if (!boundary_i && boundary_j)
                if (geometry->node[iPoint]->GetDomain()){
                  for (iVar = 0; iVar< nVar; iVar++)
                  Und_Lapl[iVar] -= Diff[iVar];
                }
                
              }
            }
            
            /*--- Store the components to be communicated in the buffer. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
            bufDSend[buf_offset+iVar] = Und_Lapl[iVar];
            
            /*--- Rotate the momentum components of the Laplacian. ---*/
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Und_Lapl[1] +
                                          rotMatrix[0][1]*Und_Lapl[2]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Und_Lapl[1] +
                                          rotMatrix[1][1]*Und_Lapl[2]);
              }
              else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Und_Lapl[1] +
                                          rotMatrix[0][1]*Und_Lapl[2] +
                                          rotMatrix[0][2]*Und_Lapl[3]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Und_Lapl[1] +
                                          rotMatrix[1][1]*Und_Lapl[2] +
                                          rotMatrix[1][2]*Und_Lapl[3]);
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*Und_Lapl[1] +
                                          rotMatrix[2][1]*Und_Lapl[2] +
                                          rotMatrix[2][2]*Und_Lapl[3]);
              }
            }
            
            break;
            
          case PERIODIC_MAX_EIG:
            
            /*--- Simple summation of eig calc on both periodic faces. ---*/
            
            bufDSend[buf_offset] = base_nodes->GetLambda(iPoint);
            
            break;
            
          case PERIODIC_SENSOR:
            
            /*--- For the centered schemes, the sensor must be computed
             consistently using info from the entire control volume
             on both sides of the periodic face. ---*/
            
            Sensor_i = 0.0; Sensor_j = 0.0;
            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              
              /*--- Avoid halos and boundary points so that we don't
               duplicate edges on both sides of the periodic BC. ---*/
              
              if (!geometry->node[jPoint]->GetPeriodicBoundary()) {
                
                /*--- Use density instead of pressure for incomp. flows. ---*/
                
                if ((config->GetKind_Regime() == INCOMPRESSIBLE)) {
                  Pressure_i = base_nodes->GetDensity(iPoint);
                  Pressure_j = base_nodes->GetDensity(jPoint);
                } else {
                  Pressure_i = base_nodes->GetPressure(iPoint);
                  Pressure_j = base_nodes->GetPressure(jPoint);
                }
                
                boundary_i = geometry->node[iPoint]->GetPhysicalBoundary();
                boundary_j = geometry->node[jPoint]->GetPhysicalBoundary();
                
                /*--- Both points inside domain, or both on boundary ---*/
                
                if ((!boundary_i && !boundary_j) ||
                    (boundary_i && boundary_j)) {
                  if (geometry->node[iPoint]->GetDomain()) {
                    Sensor_i += Pressure_j - Pressure_i;
                    Sensor_j += Pressure_i + Pressure_j;
                  }
                }
                
                /*--- iPoint inside the domain, jPoint on the boundary ---*/
                
                if (!boundary_i && boundary_j) {
                  if (geometry->node[iPoint]->GetDomain()) {
                    Sensor_i += (Pressure_j - Pressure_i);
                    Sensor_j += (Pressure_i + Pressure_j);
                    
                  }
                }
                
              }
            }
            
            /*--- Store the sensor increments to buffer. After summing
             all contributions, these will be divided. ---*/
            
            bufDSend[buf_offset] = Sensor_i;
            buf_offset++;
            bufDSend[buf_offset] = Sensor_j;
            
            break;
            
          case PERIODIC_SOL_GG:
            
            /*--- Access and rotate the partial G-G gradient. These will be
             summed on both sides of the periodic faces before dividing
             by the volume to complete the Green-Gauss gradient calc. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                jacBlock[iVar][iDim] = base_nodes->GetGradient(iPoint, iVar, iDim);
                rotBlock[iVar][iDim] = base_nodes->GetGradient(iPoint, iVar, iDim);
              }
            }
            
            /*--- Rotate the gradients in x,y,z space for all variables. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              if (nDim == 2) {
                rotBlock[iVar][0] = (rotMatrix[0][0]*jacBlock[iVar][0] +
                                     rotMatrix[0][1]*jacBlock[iVar][1]);
                rotBlock[iVar][1] = (rotMatrix[1][0]*jacBlock[iVar][0] +
                                     rotMatrix[1][1]*jacBlock[iVar][1]);
              } else {
                
                rotBlock[iVar][0] = (rotMatrix[0][0]*jacBlock[iVar][0] +
                                     rotMatrix[0][1]*jacBlock[iVar][1] +
                                     rotMatrix[0][2]*jacBlock[iVar][2]);
                rotBlock[iVar][1] = (rotMatrix[1][0]*jacBlock[iVar][0] +
                                     rotMatrix[1][1]*jacBlock[iVar][1] +
                                     rotMatrix[1][2]*jacBlock[iVar][2]);
                rotBlock[iVar][2] = (rotMatrix[2][0]*jacBlock[iVar][0] +
                                     rotMatrix[2][1]*jacBlock[iVar][1] +
                                     rotMatrix[2][2]*jacBlock[iVar][2]);
              }
            }
            
            /*--- Store the partial gradient in the buffer. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+iVar*nDim+iDim] = rotBlock[iVar][iDim];
              }
            }
            
            break;
            
          case PERIODIC_PRIM_GG:
            
            /*--- Access and rotate the partial G-G gradient. These will be
             summed on both sides of the periodic faces before dividing
             by the volume to complete the Green-Gauss gradient calc. ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++){
                jacBlock[iVar][iDim] = base_nodes->GetGradient_Primitive(iPoint, iVar, iDim);
                rotBlock[iVar][iDim] = base_nodes->GetGradient_Primitive(iPoint, iVar, iDim);
              }
            }
            
            /*--- Rotate the partial gradients in space for all variables. ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              if (nDim == 2) {
                rotBlock[iVar][0] = (rotMatrix[0][0]*jacBlock[iVar][0] +
                                     rotMatrix[0][1]*jacBlock[iVar][1]);
                rotBlock[iVar][1] = (rotMatrix[1][0]*jacBlock[iVar][0] +
                                     rotMatrix[1][1]*jacBlock[iVar][1]);
              } else {
                rotBlock[iVar][0] = (rotMatrix[0][0]*jacBlock[iVar][0] +
                                     rotMatrix[0][1]*jacBlock[iVar][1] +
                                     rotMatrix[0][2]*jacBlock[iVar][2]);
                rotBlock[iVar][1] = (rotMatrix[1][0]*jacBlock[iVar][0] +
                                     rotMatrix[1][1]*jacBlock[iVar][1] +
                                     rotMatrix[1][2]*jacBlock[iVar][2]);
                rotBlock[iVar][2] = (rotMatrix[2][0]*jacBlock[iVar][0] +
                                     rotMatrix[2][1]*jacBlock[iVar][1] +
                                     rotMatrix[2][2]*jacBlock[iVar][2]);
              }
            }
            
            /*--- Store the partial gradient in the buffer. ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+iVar*nDim+iDim] = rotBlock[iVar][iDim];
              }
            }
            
            break;
            
          case PERIODIC_SOL_LS:
            
            /*--- For L-S gradient calculations with rotational periodicity,
             we will need to rotate the x,y,z components. To make the process
             easier, we choose to rotate the initial periodic point and their
             neighbor points into their location on the donor marker before
             computing the terms that we need to communicate. ---*/
            
            /*--- Get coordinates for the current point. ---*/
            
            Coord_i = geometry->node[iPoint]->GetCoord();
            
            /*--- Get the position vector from rotation center to point. ---*/
            
            dx = Coord_i[0] - center[0];
            dy = Coord_i[1] - center[1];
            if (nDim == 3) dz = Coord_i[2] - center[2];
            else           dz = 0.0;
            
            /*--- Compute transformed point coordinates. ---*/
            
            rotCoord_i[0] = (rotMatrix[0][0]*dx +
                             rotMatrix[0][1]*dy +
                             rotMatrix[0][2]*dz + translation[0]);
            
            rotCoord_i[1] = (rotMatrix[1][0]*dx +
                             rotMatrix[1][1]*dy +
                             rotMatrix[1][2]*dz + translation[1]);
            
            rotCoord_i[2] = (rotMatrix[2][0]*dx +
                             rotMatrix[2][1]*dy +
                             rotMatrix[2][2]*dz + translation[2]);
            
            /*--- Get conservative solution and rotate if necessary. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
            rotPrim_i[iVar] = base_nodes->GetSolution(iPoint, iVar);
            
            if (rotate_periodic) {
              if (nDim == 2) {
                rotPrim_i[1] = (rotMatrix[0][0]*base_nodes->GetSolution(iPoint,1) +
                                rotMatrix[0][1]*base_nodes->GetSolution(iPoint,2));
                rotPrim_i[2] = (rotMatrix[1][0]*base_nodes->GetSolution(iPoint,1) +
                                rotMatrix[1][1]*base_nodes->GetSolution(iPoint,2));
              }
              else {
                rotPrim_i[1] = (rotMatrix[0][0]*base_nodes->GetSolution(iPoint,1) +
                                rotMatrix[0][1]*base_nodes->GetSolution(iPoint,2) +
                                rotMatrix[0][2]*base_nodes->GetSolution(iPoint,3));
                rotPrim_i[2] = (rotMatrix[1][0]*base_nodes->GetSolution(iPoint,1) +
                                rotMatrix[1][1]*base_nodes->GetSolution(iPoint,2) +
                                rotMatrix[1][2]*base_nodes->GetSolution(iPoint,3));
                rotPrim_i[3] = (rotMatrix[2][0]*base_nodes->GetSolution(iPoint,1) +
                                rotMatrix[2][1]*base_nodes->GetSolution(iPoint,2) +
                                rotMatrix[2][2]*base_nodes->GetSolution(iPoint,3));
              }
            }
            
            /*--- Inizialization of variables ---*/
            
            for (iVar = 0; iVar < nVar; iVar++)
            for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] = 0.0;
            
            r11 = 0.0;   r12 = 0.0;   r22 = 0.0;
            r13 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
            
            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              
              /*--- Avoid periodic boundary points so that we do not
               duplicate edges on both sides of the periodic BC. ---*/
              
              if (!geometry->node[jPoint]->GetPeriodicBoundary()) {
                
                /*--- Get coordinates for the neighbor point. ---*/
                
                Coord_j = geometry->node[jPoint]->GetCoord();
                
                /*--- Get the position vector from rotation center. ---*/
                
                dx = Coord_j[0] - center[0];
                dy = Coord_j[1] - center[1];
                if (nDim == 3) dz = Coord_j[2] - center[2];
                else           dz = 0.0;
                
                /*--- Compute transformed point coordinates. ---*/
                
                rotCoord_j[0] = (rotMatrix[0][0]*dx +
                                 rotMatrix[0][1]*dy +
                                 rotMatrix[0][2]*dz + translation[0]);
                
                rotCoord_j[1] = (rotMatrix[1][0]*dx +
                                 rotMatrix[1][1]*dy +
                                 rotMatrix[1][2]*dz + translation[1]);
                
                rotCoord_j[2] = (rotMatrix[2][0]*dx +
                                 rotMatrix[2][1]*dy +
                                 rotMatrix[2][2]*dz + translation[2]);
                
                /*--- Get conservative solution and rotte if necessary. ---*/
                
                for (iVar = 0; iVar < nVar; iVar++)
                rotPrim_j[iVar] = base_nodes->GetSolution(jPoint,iVar);
                
                if (rotate_periodic) {
                  if (nDim == 2) {
                    rotPrim_j[1] = (rotMatrix[0][0]*base_nodes->GetSolution(jPoint,1) +
                                    rotMatrix[0][1]*base_nodes->GetSolution(jPoint,2));
                    rotPrim_j[2] = (rotMatrix[1][0]*base_nodes->GetSolution(jPoint,1) +
                                    rotMatrix[1][1]*base_nodes->GetSolution(jPoint,2));
                  }
                  else {
                    rotPrim_j[1] = (rotMatrix[0][0]*base_nodes->GetSolution(jPoint,1) +
                                    rotMatrix[0][1]*base_nodes->GetSolution(jPoint,2) +
                                    rotMatrix[0][2]*base_nodes->GetSolution(jPoint,3));
                    rotPrim_j[2] = (rotMatrix[1][0]*base_nodes->GetSolution(jPoint,1) +
                                    rotMatrix[1][1]*base_nodes->GetSolution(jPoint,2) +
                                    rotMatrix[1][2]*base_nodes->GetSolution(jPoint,3));
                    rotPrim_j[3] = (rotMatrix[2][0]*base_nodes->GetSolution(jPoint,1) +
                                    rotMatrix[2][1]*base_nodes->GetSolution(jPoint,2) +
                                    rotMatrix[2][2]*base_nodes->GetSolution(jPoint,3));
                  }
                }
                
                weight = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                  weight += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                             (rotCoord_j[iDim]-rotCoord_i[iDim]));
                }
                
                /*--- Sumations for entries of upper triangular matrix R ---*/
                
                if (weight != 0.0) {
                  
                  r11 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[0]-rotCoord_i[0])/weight);
                  r12 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);
                  r22 += ((rotCoord_j[1]-rotCoord_i[1])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);
                  
                  if (nDim == 3) {
                    r13   += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_a += ((rotCoord_j[1]-rotCoord_i[1])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_b += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r33   += ((rotCoord_j[2]-rotCoord_i[2])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                  }
                  
                  /*--- Entries of c:= transpose(A)*b ---*/
                  
                  for (iVar = 0; iVar < nVar; iVar++)
                  for (iDim = 0; iDim < nDim; iDim++)
                  Cvector[iVar][iDim] += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                                          (rotPrim_j[iVar]-rotPrim_i[iVar])/weight);
                  
                }
              }
            }
            
            /*--- We store and communicate the increments for the matching
             upper triangular matrix (weights) and the r.h.s. vector.
             These will be accumulated before completing the L-S gradient
             calculation for each periodic point. ---*/
            
            if (nDim == 2) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
            }
            if (nDim == 3) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = r13;   buf_offset++;
              
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
              bufDSend[buf_offset] = r23_a; buf_offset++;
              
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r23_b; buf_offset++;
              bufDSend[buf_offset] = r33;   buf_offset++;
            }
            
            for (iVar = 0; iVar < nVar; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset] = Cvector[iVar][iDim];
                buf_offset++;
              }
            }
            
            break;
            
          case PERIODIC_PRIM_LS:
            
            /*--- For L-S gradient calculations with rotational periodicity,
             we will need to rotate the x,y,z components. To make the process
             easier, we choose to rotate the initial periodic point and their
             neighbor points into their location on the donor marker before
             computing the terms that we need to communicate. ---*/
            
            /*--- Get coordinates ---*/
            
            Coord_i = geometry->node[iPoint]->GetCoord();
            
            /*--- Get the position vector from rot center to point. ---*/
            
            dx = Coord_i[0] - center[0];
            dy = Coord_i[1] - center[1];
            if (nDim == 3) dz = Coord_i[2] - center[2];
            else           dz = 0.0;
            
            /*--- Compute transformed point coordinates. ---*/
            
            rotCoord_i[0] = (rotMatrix[0][0]*dx +
                             rotMatrix[0][1]*dy +
                             rotMatrix[0][2]*dz + translation[0]);
            
            rotCoord_i[1] = (rotMatrix[1][0]*dx +
                             rotMatrix[1][1]*dy +
                             rotMatrix[1][2]*dz + translation[1]);
            
            rotCoord_i[2] = (rotMatrix[2][0]*dx +
                             rotMatrix[2][1]*dy +
                             rotMatrix[2][2]*dz + translation[2]);
            
            /*--- Get primitives and rotate if necessary. ---*/
            
            for (iVar = 0; iVar < nPrimVar; iVar++)
            rotPrim_i[iVar] = base_nodes->GetPrimitive(iPoint, iVar);
            
            if (rotate_periodic) {
              if (nDim == 2) {
                rotPrim_i[1] = (rotMatrix[0][0]*base_nodes->GetPrimitive(iPoint,1) +
                                rotMatrix[0][1]*base_nodes->GetPrimitive(iPoint,2));
                rotPrim_i[2] = (rotMatrix[1][0]*base_nodes->GetPrimitive(iPoint,1) +
                                rotMatrix[1][1]*base_nodes->GetPrimitive(iPoint,2));
              }
              else {
                rotPrim_i[1] = (rotMatrix[0][0]*base_nodes->GetPrimitive(iPoint,1) +
                                rotMatrix[0][1]*base_nodes->GetPrimitive(iPoint,2) +
                                rotMatrix[0][2]*base_nodes->GetPrimitive(iPoint,3));
                rotPrim_i[2] = (rotMatrix[1][0]*base_nodes->GetPrimitive(iPoint,1) +
                                rotMatrix[1][1]*base_nodes->GetPrimitive(iPoint,2) +
                                rotMatrix[1][2]*base_nodes->GetPrimitive(iPoint,3));
                rotPrim_i[3] = (rotMatrix[2][0]*base_nodes->GetPrimitive(iPoint,1) +
                                rotMatrix[2][1]*base_nodes->GetPrimitive(iPoint,2) +
                                rotMatrix[2][2]*base_nodes->GetPrimitive(iPoint,3));
              }
            }
            
            /*--- Inizialization of variables ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
            for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] = 0.0;
            
            r11 = 0.0;   r12 = 0.0;   r22 = 0.0;
            r13 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
            
            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              
              /*--- Avoid periodic boundary points so that we do not
               duplicate edges on both sides of the periodic BC. ---*/
              
              if (!geometry->node[jPoint]->GetPeriodicBoundary()) {
                
                /*--- Get coordinates for the neighbor point. ---*/
                
                Coord_j = geometry->node[jPoint]->GetCoord();
                
                /*--- Get the position vector from rotation center. ---*/
                
                dx = Coord_j[0] - center[0];
                dy = Coord_j[1] - center[1];
                if (nDim == 3) dz = Coord_j[2] - center[2];
                else           dz = 0.0;
                
                /*--- Compute transformed point coordinates. ---*/
                
                rotCoord_j[0] = (rotMatrix[0][0]*dx +
                                 rotMatrix[0][1]*dy +
                                 rotMatrix[0][2]*dz + translation[0]);
                
                rotCoord_j[1] = (rotMatrix[1][0]*dx +
                                 rotMatrix[1][1]*dy +
                                 rotMatrix[1][2]*dz + translation[1]);
                
                rotCoord_j[2] = (rotMatrix[2][0]*dx +
                                 rotMatrix[2][1]*dy +
                                 rotMatrix[2][2]*dz + translation[2]);
                
                /*--- Get primitives from CVariable ---*/
                
                for (iVar = 0; iVar < nPrimVar; iVar++)
                rotPrim_j[iVar] = base_nodes->GetPrimitive(jPoint,iVar);
                
                if (rotate_periodic) {
                  if (nDim == 2) {
                    rotPrim_j[1] = (rotMatrix[0][0]*base_nodes->GetPrimitive(jPoint,1) +
                                    rotMatrix[0][1]*base_nodes->GetPrimitive(jPoint,2));
                    rotPrim_j[2] = (rotMatrix[1][0]*base_nodes->GetPrimitive(jPoint,1) +
                                    rotMatrix[1][1]*base_nodes->GetPrimitive(jPoint,2));
                  }
                  else {
                    rotPrim_j[1] = (rotMatrix[0][0]*base_nodes->GetPrimitive(jPoint,1) +
                                    rotMatrix[0][1]*base_nodes->GetPrimitive(jPoint,2) +
                                    rotMatrix[0][2]*base_nodes->GetPrimitive(jPoint,3));
                    rotPrim_j[2] = (rotMatrix[1][0]*base_nodes->GetPrimitive(jPoint,1) +
                                    rotMatrix[1][1]*base_nodes->GetPrimitive(jPoint,2) +
                                    rotMatrix[1][2]*base_nodes->GetPrimitive(jPoint,3));
                    rotPrim_j[3] = (rotMatrix[2][0]*base_nodes->GetPrimitive(jPoint,1) +
                                    rotMatrix[2][1]*base_nodes->GetPrimitive(jPoint,2) +
                                    rotMatrix[2][2]*base_nodes->GetPrimitive(jPoint,3));
                  }
                }
                
                weight = 0.0;
                for (iDim = 0; iDim < nDim; iDim++)
                weight += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                           (rotCoord_j[iDim]-rotCoord_i[iDim]));
                
                /*--- Sumations for entries of upper triangular matrix R ---*/
                
                if (weight != 0.0) {
                  
                  r11 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[0]-rotCoord_i[0])/weight);
                  r12 += ((rotCoord_j[0]-rotCoord_i[0])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);
                  r22 += ((rotCoord_j[1]-rotCoord_i[1])*
                          (rotCoord_j[1]-rotCoord_i[1])/weight);
                  
                  if (nDim == 3) {
                    r13   += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_a += ((rotCoord_j[1]-rotCoord_i[1])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r23_b += ((rotCoord_j[0]-rotCoord_i[0])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                    r33   += ((rotCoord_j[2]-rotCoord_i[2])*
                              (rotCoord_j[2]-rotCoord_i[2])/weight);
                  }
                  
                  /*--- Entries of c:= transpose(A)*b ---*/
                  
                  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
                  for (iDim = 0; iDim < nDim; iDim++)
                  Cvector[iVar][iDim] += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                                          (rotPrim_j[iVar]-rotPrim_i[iVar])/weight);
                  
                }
              }
            }
            
            /*--- We store and communicate the increments for the matching
             upper triangular matrix (weights) and the r.h.s. vector.
             These will be accumulated before completing the L-S gradient
             calculation for each periodic point. ---*/
            
            if (nDim == 2) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
            }
            if (nDim == 3) {
              bufDSend[buf_offset] = r11;   buf_offset++;
              bufDSend[buf_offset] = r12;   buf_offset++;
              bufDSend[buf_offset] = r13;   buf_offset++;
              
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r22;   buf_offset++;
              bufDSend[buf_offset] = r23_a; buf_offset++;
              
              bufDSend[buf_offset] = 0.0;   buf_offset++;
              bufDSend[buf_offset] = r23_b; buf_offset++;
              bufDSend[buf_offset] = r33;   buf_offset++;
            }
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset] = Cvector[iVar][iDim];
                buf_offset++;
              }
            }
            
            break;
            
          case PERIODIC_LIM_PRIM_1:
            
            /*--- The first phase of the periodic limiter calculation
             ensures that the proper min and max of the solution are found
             among all nodes adjacent to periodic faces. ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              Sol_Min[iVar] = base_nodes->GetSolution_Min(iPoint, iVar);
              Sol_Max[iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
              
              bufDSend[buf_offset+iVar]              = base_nodes->GetSolution_Min(iPoint, iVar);
              bufDSend[buf_offset+nPrimVarGrad+iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
            }
            
            /*--- Rotate the momentum components of the min/max. ---*/
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Sol_Min[1] +
                                          rotMatrix[0][1]*Sol_Min[2]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Sol_Min[1] +
                                          rotMatrix[1][1]*Sol_Min[2]);
                
                bufDSend[buf_offset+nPrimVarGrad+1] = (rotMatrix[0][0]*Sol_Max[1] +
                                                       rotMatrix[0][1]*Sol_Max[2]);
                bufDSend[buf_offset+nPrimVarGrad+2] = (rotMatrix[1][0]*Sol_Max[1] +
                                                       rotMatrix[1][1]*Sol_Max[2]);
                
              } else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Sol_Min[1] +
                                          rotMatrix[0][1]*Sol_Min[2] +
                                          rotMatrix[0][2]*Sol_Min[3]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Sol_Min[1] +
                                          rotMatrix[1][1]*Sol_Min[2] +
                                          rotMatrix[1][2]*Sol_Min[3]);
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*Sol_Min[1] +
                                          rotMatrix[2][1]*Sol_Min[2] +
                                          rotMatrix[2][2]*Sol_Min[3]);
                
                bufDSend[buf_offset+nPrimVarGrad+1] = (rotMatrix[0][0]*Sol_Max[1] +
                                                       rotMatrix[0][1]*Sol_Max[2] +
                                                       rotMatrix[0][2]*Sol_Max[3]);
                bufDSend[buf_offset+nPrimVarGrad+2] = (rotMatrix[1][0]*Sol_Max[1] +
                                                       rotMatrix[1][1]*Sol_Max[2] +
                                                       rotMatrix[1][2]*Sol_Max[3]);
                bufDSend[buf_offset+nPrimVarGrad+3] = (rotMatrix[2][0]*Sol_Max[1] +
                                                       rotMatrix[2][1]*Sol_Max[2] +
                                                       rotMatrix[2][2]*Sol_Max[3]);
              }
            }
            
            break;
            
          case PERIODIC_LIM_PRIM_2:
            
            /*--- The second phase of the periodic limiter calculation
             ensures that the correct minimum value of the limiter is
             found for a node on a periodic face and stores it. ---*/
            
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetLimiter_Primitive(iPoint, iVar);
            }
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetLimiter_Primitive(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetLimiter_Primitive(iPoint,2));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetLimiter_Primitive(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetLimiter_Primitive(iPoint,2));
                
              }
              else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetLimiter_Primitive(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetLimiter_Primitive(iPoint,2) +
                                          rotMatrix[0][2]*base_nodes->GetLimiter_Primitive(iPoint,3));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetLimiter_Primitive(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetLimiter_Primitive(iPoint,2) +
                                          rotMatrix[1][2]*base_nodes->GetLimiter_Primitive(iPoint,3));
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*base_nodes->GetLimiter_Primitive(iPoint,1) +
                                          rotMatrix[2][1]*base_nodes->GetLimiter_Primitive(iPoint,2) +
                                          rotMatrix[2][2]*base_nodes->GetLimiter_Primitive(iPoint,3));
              }
            }
            
            break;
            
          case PERIODIC_LIM_SOL_1:
            
            /*--- The first phase of the periodic limiter calculation
             ensures that the proper min and max of the solution are found
             among all nodes adjacent to periodic faces. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              Sol_Min[iVar] = base_nodes->GetSolution_Min(iPoint, iVar);
              Sol_Max[iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
              
              bufDSend[buf_offset+iVar]      = base_nodes->GetSolution_Min(iPoint, iVar);
              bufDSend[buf_offset+nVar+iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
            }
            
            /*--- Rotate the momentum components of the min/max. ---*/
            
            if (rotate_periodic) {
              
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Sol_Min[1] +
                                          rotMatrix[0][1]*Sol_Min[2]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Sol_Min[1] +
                                          rotMatrix[1][1]*Sol_Min[2]);
                
                bufDSend[buf_offset+nVar+1] = (rotMatrix[0][0]*Sol_Max[1] +
                                               rotMatrix[0][1]*Sol_Max[2]);
                bufDSend[buf_offset+nVar+2] = (rotMatrix[1][0]*Sol_Max[1] +
                                               rotMatrix[1][1]*Sol_Max[2]);
                
              }
              else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*Sol_Min[1] +
                                          rotMatrix[0][1]*Sol_Min[2] +
                                          rotMatrix[0][2]*Sol_Min[3]);
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*Sol_Min[1] +
                                          rotMatrix[1][1]*Sol_Min[2] +
                                          rotMatrix[1][2]*Sol_Min[3]);
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*Sol_Min[1] +
                                          rotMatrix[2][1]*Sol_Min[2] +
                                          rotMatrix[2][2]*Sol_Min[3]);
                
                bufDSend[buf_offset+nVar+1] = (rotMatrix[0][0]*Sol_Max[1] +
                                               rotMatrix[0][1]*Sol_Max[2] +
                                               rotMatrix[0][2]*Sol_Max[3]);
                bufDSend[buf_offset+nVar+2] = (rotMatrix[1][0]*Sol_Max[1] +
                                               rotMatrix[1][1]*Sol_Max[2] +
                                               rotMatrix[1][2]*Sol_Max[3]);
                bufDSend[buf_offset+nVar+3] = (rotMatrix[2][0]*Sol_Max[1] +
                                               rotMatrix[2][1]*Sol_Max[2] +
                                               rotMatrix[2][2]*Sol_Max[3]);
                
              }
            }
            
            break;
            
          case PERIODIC_LIM_SOL_2:
            
            /*--- The second phase of the periodic limiter calculation
             ensures that the correct minimum value of the limiter is
             found for a node on a periodic face and stores it. ---*/
            
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetLimiter(iPoint, iVar);
            }
            
            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetLimiter(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetLimiter(iPoint,2));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetLimiter(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetLimiter(iPoint,2));
                
              }
              else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*base_nodes->GetLimiter(iPoint,1) +
                                          rotMatrix[0][1]*base_nodes->GetLimiter(iPoint,2) +
                                          rotMatrix[0][2]*base_nodes->GetLimiter(iPoint,3));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*base_nodes->GetLimiter(iPoint,1) +
                                          rotMatrix[1][1]*base_nodes->GetLimiter(iPoint,2) +
                                          rotMatrix[1][2]*base_nodes->GetLimiter(iPoint,3));
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*base_nodes->GetLimiter(iPoint,1) +
                                          rotMatrix[2][1]*base_nodes->GetLimiter(iPoint,2) +
                                          rotMatrix[2][2]*base_nodes->GetLimiter(iPoint,3));
              }
            }
            
            break;
            
          default:
            SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                           CURRENT_FUNCTION);
            break;
        }
      }
      
      /*--- Launch the point-to-point MPI send for this message. ---*/
      
      geometry->PostPeriodicSends(geometry, config, MPI_TYPE, iMessage);
      
    }
  }
  
  delete [] Diff;
  delete [] Und_Lapl;
  delete [] Sol_Min;
  delete [] Sol_Max;
  delete [] rotPrim_i;
  delete [] rotPrim_j;
  
  for (iVar = 0; iVar < ICOUNT; iVar++) {
    delete [] jacBlock[iVar];
    delete [] rotBlock[iVar];
  }
  delete [] jacBlock;
  delete [] rotBlock;
  
}

void CSolver::CompletePeriodicComms(CGeometry *geometry,
                                    CConfig *config,
                                    unsigned short val_periodic_index,
                                    unsigned short commType) {
  
  /*--- Local variables ---*/
  
  unsigned short nPeriodic = config->GetnMarker_Periodic();
  unsigned short iDim, jDim, iVar, jVar, iPeriodic, nNeighbor;
  
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset, total_index;
  
  int source, iMessage, jRecv;
  
  SU2_MPI::Status status;
  
  su2double *Diff = new su2double[nVar];
  
  su2double Time_Step, Volume, Solution_Min, Solution_Max, Limiter_Min;
  
  /*--- Set some local pointers to make access simpler. ---*/
  
  su2double *bufDRecv = geometry->bufD_PeriodicRecv;
  
  unsigned short *bufSRecv = geometry->bufS_PeriodicRecv;
  
  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/
  
  if (geometry->nPeriodicRecv > 0) {
    
    for (iMessage = 0; iMessage < geometry->nPeriodicRecv; iMessage++) {
      
      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/
      
#ifdef HAVE_MPI
      /*--- Once we have recv'd a message, get the source rank. ---*/
      int ind;
      SU2_MPI::Waitany(geometry->nPeriodicRecv,
                       geometry->req_PeriodicRecv,
                       &ind, &status);
      source = status.MPI_SOURCE;
#else
      /*--- For serial calculations, we know the rank. ---*/
      source = rank;
#endif
      
      /*--- We know the offsets based on the source rank. ---*/
      
      jRecv = geometry->PeriodicRecv2Neighbor[source];
      
      /*--- Get the offset in the buffer for the start of this message. ---*/
      
      msg_offset = geometry->nPoint_PeriodicRecv[jRecv];
      
      /*--- Get the number of packets to be received in this message. ---*/
      
      nRecv = (geometry->nPoint_PeriodicRecv[jRecv+1] -
               geometry->nPoint_PeriodicRecv[jRecv]);
      
      for (iRecv = 0; iRecv < nRecv; iRecv++) {
        
        /*--- Get the local index for this communicated data. ---*/
        
        iPoint    = geometry->Local_Point_PeriodicRecv[msg_offset  + iRecv];
        iPeriodic = geometry->Local_Marker_PeriodicRecv[msg_offset + iRecv];
        
        /*--- While all periodic face data was accumulated, we only store
         the values for the current pair of periodic faces. This is slightly
         inefficient when we have multiple pairs of periodic faces, but
         it simplifies the communications. ---*/
        
        if ((iPeriodic == val_periodic_index) ||
            (iPeriodic == val_periodic_index + nPeriodic/2)) {
          
          /*--- Compute the offset in the recv buffer for this point. ---*/
          
          buf_offset = (msg_offset + iRecv)*geometry->countPerPeriodicPoint;
          
          /*--- Store the data correctly depending on the quantity. ---*/
          
          switch (commType) {
              
            case PERIODIC_VOLUME:
              
              /*--- The periodic points need to keep track of their
               total volume spread across the periodic faces. ---*/
              
              Volume = (bufDRecv[buf_offset] +
                        geometry->node[iPoint]->GetPeriodicVolume());
              geometry->node[iPoint]->SetPeriodicVolume(Volume);
              
              break;
              
            case PERIODIC_NEIGHBORS:
              
              /*--- Store the extra neighbors on the periodic face. ---*/
              
              nNeighbor = (geometry->node[iPoint]->GetnNeighbor() +
                           bufSRecv[buf_offset]);
              geometry->node[iPoint]->SetnNeighbor(nNeighbor);
              
              break;
              
            case PERIODIC_RESIDUAL:
              
              /*--- Access the residual from the donor. ---*/
              
              for (iVar = 0; iVar < nVar; iVar++) {
                Residual[iVar] = bufDRecv[buf_offset];
                buf_offset++;
              }
              
              /*--- Check the computed time step against the donor
               value and keep the minimum in order to be conservative. ---*/
              
              Time_Step = base_nodes->GetDelta_Time(iPoint);
              if (bufDRecv[buf_offset] < Time_Step)
                base_nodes->SetDelta_Time(iPoint,bufDRecv[buf_offset]);
              buf_offset++;
              
              /*--- Access the Jacobian from the donor if implicit. ---*/
              
              if (implicit_periodic) {
                for (iVar = 0; iVar < nVar; iVar++) {
                  for (jVar = 0; jVar < nVar; jVar++) {
                    Jacobian_i[iVar][jVar] = bufDRecv[buf_offset];
                    buf_offset++;
                  }
                }
              }
              
              /*--- Add contributions to total residual. ---*/
              
              LinSysRes.AddBlock(iPoint, Residual);
              
              /*--- For implicit integration, we choose the first
               periodic face of each pair to be the master/owner of
               the solution for the linear system while fixing the
               solution at the matching face during the solve. Here,
               we remove the Jacobian and residual contributions from
               the passive face such that it does not participate in
               the linear solve. ---*/
              
              if (implicit_periodic) {
                
                Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
                
                if (iPeriodic == val_periodic_index + nPeriodic/2) {
                  for (iVar = 0; iVar < nVar; iVar++) {
                    LinSysRes.SetBlock_Zero(iPoint, iVar);
                    total_index = iPoint*nVar+iVar;
                    Jacobian.DeleteValsRowi(total_index);
                  }
                }
                
              }
              
              break;
              
            case PERIODIC_IMPLICIT:
              
              /*--- For implicit integration, we choose the first
               periodic face of each pair to be the master/owner of
               the solution for the linear system while fixing the
               solution at the matching face during the solve. Here,
               we are updating the solution at the passive nodes
               using the new solution from the master. ---*/
              
              if ((implicit_periodic) &&
                  (iPeriodic == val_periodic_index + nPeriodic/2)) {
                
                /*--- Access the solution from the donor. ---*/
                
                for (iVar = 0; iVar < nVar; iVar++) {
                  Solution[iVar] = bufDRecv[buf_offset];
                  buf_offset++;
                }
                
                /*--- Directly set the solution on the passive periodic
                 face that is provided from the master. ---*/
                
                for (iVar = 0; iVar < nVar; iVar++) {
                  base_nodes->SetSolution(iPoint, iVar, Solution[iVar]);
                  base_nodes->SetSolution_Old(iPoint, iVar, Solution[iVar]);
                }
                
              }
              
              break;
              
            case PERIODIC_LAPLACIAN:
              
              /*--- Adjust the undivided Laplacian. The accumulation was
               with a subtraction before communicating, so now just add. ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                Diff[iVar] = bufDRecv[buf_offset+iVar];
              
              base_nodes->AddUnd_Lapl(iPoint,Diff);
              
              break;
              
            case PERIODIC_MAX_EIG:
              
              /*--- Simple accumulation of the max eig on periodic faces. ---*/
              
              base_nodes->AddLambda(iPoint,bufDRecv[buf_offset]);
              
              break;
              
            case PERIODIC_SENSOR:
              
              /*--- Simple accumulation of the sensors on periodic faces. ---*/
              
              iPoint_UndLapl[iPoint] += bufDRecv[buf_offset]; buf_offset++;
              jPoint_UndLapl[iPoint] += bufDRecv[buf_offset];
              
              break;
              
            case PERIODIC_SOL_GG:
              
              /*--- For G-G, we accumulate partial gradients then compute
               the final value using the entire volume of the periodic cell. ---*/
              
              for (iVar = 0; iVar < nVar; iVar++)
                for (iDim = 0; iDim < nDim; iDim++)
                  base_nodes->SetGradient(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim] + base_nodes->GetGradient(iPoint, iVar, iDim));
              
              break;
              
            case PERIODIC_PRIM_GG:
              
              /*--- For G-G, we accumulate partial gradients then compute
               the final value using the entire volume of the periodic cell. ---*/
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++)
                for (iDim = 0; iDim < nDim; iDim++)
                  base_nodes->SetGradient_Primitive(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim] + base_nodes->GetGradient_Primitive(iPoint, iVar, iDim));
              break;
              
            case PERIODIC_SOL_LS:
              
              /*--- For L-S, we build the upper triangular matrix and the
               r.h.s. vector by accumulating from all periodic partial
               control volumes. ---*/
              
              for (iDim = 0; iDim < nDim; iDim++) {
                for (jDim = 0; jDim < nDim; jDim++) {
                  base_nodes->AddRmatrix(iPoint, iDim,jDim,bufDRecv[buf_offset]);
                  buf_offset++;
                }
              }
              for (iVar = 0; iVar < nVar; iVar++) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  base_nodes->AddGradient(iPoint, iVar, iDim, bufDRecv[buf_offset]);
                  buf_offset++;
                }
              }
              
              break;
              
            case PERIODIC_PRIM_LS:
              
              /*--- For L-S, we build the upper triangular matrix and the
               r.h.s. vector by accumulating from all periodic partial
               control volumes. ---*/
              
              for (iDim = 0; iDim < nDim; iDim++) {
                for (jDim = 0; jDim < nDim; jDim++) {
                  base_nodes->AddRmatrix(iPoint, iDim,jDim,bufDRecv[buf_offset]);
                  buf_offset++;
                }
              }
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                for (iDim = 0; iDim < nDim; iDim++) {
                  base_nodes->AddGradient_Primitive(iPoint, iVar, iDim, bufDRecv[buf_offset]);
                  buf_offset++;
                }
              }
              
              break;
              
            case PERIODIC_LIM_PRIM_1:
              
              /*--- Check the min and max values found on the matching
               perioic faces for the solution, and store the proper min
               and max for this point.  ---*/
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                base_nodes->SetSolution_Min(iPoint, iVar, min(base_nodes->GetSolution_Min(iPoint, iVar), bufDRecv[buf_offset+iVar]));
                base_nodes->SetSolution_Max(iPoint, iVar, max(base_nodes->GetSolution_Max(iPoint, iVar), bufDRecv[buf_offset+nPrimVarGrad+iVar]));
              }
              
              break;
              
            case PERIODIC_LIM_PRIM_2:
              
              /*--- Check the min values found on the matching periodic
               faces for the limiter, and store the proper min value. ---*/
              
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                base_nodes->SetLimiter_Primitive(iPoint, iVar, min(base_nodes->GetLimiter_Primitive(iPoint, iVar), bufDRecv[buf_offset+iVar]));
              }
              
              break;
              
            case PERIODIC_LIM_SOL_1:
              
              /*--- Check the min and max values found on the matching
               perioic faces for the solution, and store the proper min
               and max for this point.  ---*/
              
              for (iVar = 0; iVar < nVar; iVar++) {
                
                /*--- Solution minimum. ---*/
                
                Solution_Min = min(base_nodes->GetSolution_Min(iPoint, iVar),
                                   bufDRecv[buf_offset+iVar]);
                base_nodes->SetSolution_Min(iPoint, iVar, Solution_Min);
                
                /*--- Solution maximum. ---*/
                
                Solution_Max = max(base_nodes->GetSolution_Max(iPoint, iVar),
                                   bufDRecv[buf_offset+nVar+iVar]);
                base_nodes->SetSolution_Max(iPoint, iVar, Solution_Max);
                
              }
              
              break;
              
            case PERIODIC_LIM_SOL_2:
              
              /*--- Check the min values found on the matching periodic
               faces for the limiter, and store the proper min value. ---*/
              
              for (iVar = 0; iVar < nVar; iVar++) {
                Limiter_Min = min(base_nodes->GetLimiter_Primitive(iPoint, iVar),
                                  bufDRecv[buf_offset+iVar]);
                base_nodes->SetLimiter_Primitive(iPoint, iVar, Limiter_Min);
              }
              
              break;
              
            default:
              
              SU2_MPI::Error("Unrecognized quantity for periodic communication.",
                             CURRENT_FUNCTION);
              break;
              
          }
        }
      }
    }
    
    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Waitall(geometry->nPeriodicSend,
                     geometry->req_PeriodicSend,
                     MPI_STATUS_IGNORE);
#endif
    
  }
  
  delete [] Diff;
  
}

void CSolver::InitiateComms(CGeometry *geometry,
                            CConfig *config,
                            unsigned short commType) {
  
  /*--- Local variables ---*/
  
  unsigned short iVar, iDim;
  unsigned short COUNT_PER_POINT = 0;
  unsigned short MPI_TYPE        = 0;
  
  unsigned long iPoint, msg_offset, buf_offset;
  
  int iMessage, iSend, nSend;
  
  /*--- Set the size of the data packet and type depending on quantity. ---*/
  
  switch (commType) {
    case SOLUTION:
    case SOLUTION_OLD:
    case UNDIVIDED_LAPLACIAN:
    case SOLUTION_LIMITER:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case MAX_EIGENVALUE:
    case SENSOR:
      COUNT_PER_POINT  = 1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_GRADIENT:
      COUNT_PER_POINT  = nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PRIMITIVE_GRADIENT:
      COUNT_PER_POINT  = nPrimVarGrad*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PRIMITIVE_LIMITER:
      COUNT_PER_POINT  = nPrimVarGrad;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_EDDY:
      COUNT_PER_POINT  = nVar+1;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_FEA:
      if (config->GetTime_Domain())
        COUNT_PER_POINT  = nVar*3;
      else
        COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_FEA_OLD:
      COUNT_PER_POINT  = nVar*3;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_DISPONLY:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_PRED:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_PRED_OLD:
      COUNT_PER_POINT  = nVar*3;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case AUXVAR_GRADIENT:
      COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case MESH_DISPLACEMENTS:
      COUNT_PER_POINT  = nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_TIME_N:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case SOLUTION_TIME_N1:
      COUNT_PER_POINT  = nVar;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    default:
      SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                     CURRENT_FUNCTION);
      break;
  }
  
  /*--- Check to make sure we have created a large enough buffer
   for these comms during preprocessing. This is only for the su2double
   buffer. It will be reallocated whenever we find a larger count
   per point. After the first cycle of comms, this should be inactive. ---*/
  
  if (COUNT_PER_POINT > geometry->countPerPoint) {
    geometry->AllocateP2PComms(COUNT_PER_POINT);
  }
  
  /*--- Set some local pointers to make access simpler. ---*/
  
  su2double *bufDSend = geometry->bufD_P2PSend;
  
  /*--- Load the specified quantity from the solver into the generic
   communication buffer in the geometry class. ---*/
  
  if (geometry->nP2PSend > 0) {
    
    /*--- Post all non-blocking recvs first before sends. ---*/
    
    geometry->PostP2PRecvs(geometry, config, MPI_TYPE, false);
    
    for (iMessage = 0; iMessage < geometry->nP2PSend; iMessage++) {
      
      /*--- Get the offset in the buffer for the start of this message. ---*/
      
      msg_offset = geometry->nPoint_P2PSend[iMessage];
      
      /*--- Total count can include multiple pieces of data per element. ---*/
      
      nSend = (geometry->nPoint_P2PSend[iMessage+1] -
               geometry->nPoint_P2PSend[iMessage]);
      
      for (iSend = 0; iSend < nSend; iSend++) {
        
        /*--- Get the local index for this communicated data. ---*/
        
        iPoint = geometry->Local_Point_P2PSend[msg_offset + iSend];
        
        /*--- Compute the offset in the recv buffer for this point. ---*/
        
        buf_offset = (msg_offset + iSend)*geometry->countPerPoint;
        
        switch (commType) {
          case SOLUTION:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            break;
          case SOLUTION_OLD:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_Old(iPoint, iVar);
            break;
          case SOLUTION_EDDY:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            bufDSend[buf_offset+nVar]   = base_nodes->GetmuT(iPoint);
            break;
          case UNDIVIDED_LAPLACIAN:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetUndivided_Laplacian(iPoint, iVar);
            break;
          case SOLUTION_LIMITER:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetLimiter(iPoint, iVar);
            break;
          case MAX_EIGENVALUE:
            bufDSend[buf_offset] = base_nodes->GetLambda(iPoint);
            break;
          case SENSOR:
            bufDSend[buf_offset] = base_nodes->GetSensor(iPoint);
            break;
          case SOLUTION_GRADIENT:
            for (iVar = 0; iVar < nVar; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                bufDSend[buf_offset+iVar*nDim+iDim] = base_nodes->GetGradient(iPoint, iVar, iDim);
            break;
          case PRIMITIVE_GRADIENT:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                bufDSend[buf_offset+iVar*nDim+iDim] = base_nodes->GetGradient_Primitive(iPoint, iVar, iDim);
            break;
          case PRIMITIVE_LIMITER:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetLimiter_Primitive(iPoint, iVar);
            break;
          case AUXVAR_GRADIENT:
            for (iDim = 0; iDim < nDim; iDim++)
              bufDSend[buf_offset+iDim] = base_nodes->GetAuxVarGradient(iPoint, iDim);
            break;
          case SOLUTION_FEA:
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
              if (config->GetTime_Domain()) {
                bufDSend[buf_offset+nVar+iVar]   = base_nodes->GetSolution_Vel(iPoint, iVar);
                bufDSend[buf_offset+nVar*2+iVar] = base_nodes->GetSolution_Accel(iPoint, iVar);
              }
            }
            break;
          case SOLUTION_FEA_OLD:
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar]        = base_nodes->GetSolution_time_n(iPoint, iVar);
              bufDSend[buf_offset+nVar+iVar]   = base_nodes->GetSolution_Vel_time_n(iPoint, iVar);
              bufDSend[buf_offset+nVar*2+iVar] = base_nodes->GetSolution_Accel_time_n(iPoint, iVar);
            }
            break;
          case SOLUTION_DISPONLY:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution(iPoint, iVar);
            break;
          case SOLUTION_PRED:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_Pred(iPoint, iVar);
            break;
          case SOLUTION_PRED_OLD:
            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar]        = base_nodes->GetSolution_Old(iPoint, iVar);
              bufDSend[buf_offset+nVar+iVar]   = base_nodes->GetSolution_Pred(iPoint, iVar);
              bufDSend[buf_offset+nVar*2+iVar] = base_nodes->GetSolution_Pred_Old(iPoint, iVar);
            }
            break;
          case MESH_DISPLACEMENTS:
            for (iDim = 0; iDim < nDim; iDim++)
              bufDSend[buf_offset+iDim] = base_nodes->GetBound_Disp(iPoint, iDim);
            break;
          case SOLUTION_TIME_N:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_time_n(iPoint, iVar);
            break;
          case SOLUTION_TIME_N1:
            for (iVar = 0; iVar < nVar; iVar++)
              bufDSend[buf_offset+iVar] = base_nodes->GetSolution_time_n1(iPoint, iVar);
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }
      
      /*--- Launch the point-to-point MPI send for this message. ---*/
      
      geometry->PostP2PSends(geometry, config, MPI_TYPE, iMessage, false);
      
    }
  }
  
}
void CSolver::CompleteComms(CGeometry *geometry,
                            CConfig *config,
                            unsigned short commType) {
  
  /*--- Local variables ---*/
  
  unsigned short iDim, iVar;
  unsigned long iPoint, iRecv, nRecv, msg_offset, buf_offset;
  
  int ind, source, iMessage, jRecv;
  SU2_MPI::Status status;
  
  /*--- Set some local pointers to make access simpler. ---*/
  
  su2double *bufDRecv = geometry->bufD_P2PRecv;
  
  /*--- Store the data that was communicated into the appropriate
   location within the local class data structures. ---*/
  
  if (geometry->nP2PRecv > 0) {
    
    for (iMessage = 0; iMessage < geometry->nP2PRecv; iMessage++) {
      
      /*--- For efficiency, recv the messages dynamically based on
       the order they arrive. ---*/
      
      SU2_MPI::Waitany(geometry->nP2PRecv, geometry->req_P2PRecv,
                       &ind, &status);
      
      /*--- Once we have recv'd a message, get the source rank. ---*/
      
      source = status.MPI_SOURCE;
      
      /*--- We know the offsets based on the source rank. ---*/
      
      jRecv = geometry->P2PRecv2Neighbor[source];
      
      /*--- Get the offset in the buffer for the start of this message. ---*/
      
      msg_offset = geometry->nPoint_P2PRecv[jRecv];
      
      /*--- Get the number of packets to be received in this message. ---*/
      
      nRecv = (geometry->nPoint_P2PRecv[jRecv+1] -
               geometry->nPoint_P2PRecv[jRecv]);
      
      for (iRecv = 0; iRecv < nRecv; iRecv++) {
        
        /*--- Get the local index for this communicated data. ---*/
        
        iPoint = geometry->Local_Point_P2PRecv[msg_offset + iRecv];
        
        /*--- Compute the offset in the recv buffer for this point. ---*/
        
        buf_offset = (msg_offset + iRecv)*geometry->countPerPoint;
        
        /*--- Store the data correctly depending on the quantity. ---*/
        
        switch (commType) {
          case SOLUTION:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_OLD:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution_Old(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_EDDY:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            base_nodes->SetmuT(iPoint,bufDRecv[buf_offset+nVar]);
            break;
          case UNDIVIDED_LAPLACIAN:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetUndivided_Laplacian(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_LIMITER:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetLimiter(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case MAX_EIGENVALUE:
            base_nodes->SetLambda(iPoint,bufDRecv[buf_offset]);
            break;
          case SENSOR:
            base_nodes->SetSensor(iPoint,bufDRecv[buf_offset]);
            break;
          case SOLUTION_GRADIENT:
            for (iVar = 0; iVar < nVar; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                base_nodes->SetGradient(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim]);
            break;
          case PRIMITIVE_GRADIENT:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
              for (iDim = 0; iDim < nDim; iDim++)
                base_nodes->SetGradient_Primitive(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim]);
            break;
          case PRIMITIVE_LIMITER:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++)
              base_nodes->SetLimiter_Primitive(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case AUXVAR_GRADIENT:
            for (iDim = 0; iDim < nDim; iDim++)
              base_nodes->SetAuxVarGradient(iPoint, iDim, bufDRecv[buf_offset+iDim]);
            break;
          case SOLUTION_FEA:
            for (iVar = 0; iVar < nVar; iVar++) {
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
              if (config->GetTime_Domain()) {
                base_nodes->SetSolution_Vel(iPoint, iVar, bufDRecv[buf_offset+nVar+iVar]);
                base_nodes->SetSolution_Accel(iPoint, iVar, bufDRecv[buf_offset+nVar*2+iVar]);
              }
            }
            break;
          case SOLUTION_FEA_OLD:
            for (iVar = 0; iVar < nVar; iVar++) {
              base_nodes->Set_Solution_time_n(iPoint, iVar, bufDRecv[buf_offset+iVar]);
              base_nodes->SetSolution_Vel_time_n(iPoint, iVar, bufDRecv[buf_offset+nVar+iVar]);
              base_nodes->SetSolution_Accel_time_n(iPoint, iVar, bufDRecv[buf_offset+nVar*2+iVar]);
            }
            break;
          case SOLUTION_DISPONLY:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_PRED:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->SetSolution_Pred(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_PRED_OLD:
            for (iVar = 0; iVar < nVar; iVar++) {
              base_nodes->SetSolution_Old(iPoint, iVar, bufDRecv[buf_offset+iVar]);
              base_nodes->SetSolution_Pred(iPoint, iVar, bufDRecv[buf_offset+nVar+iVar]);
              base_nodes->SetSolution_Pred_Old(iPoint, iVar, bufDRecv[buf_offset+nVar*2+iVar]);
            }
            break;
          case MESH_DISPLACEMENTS:
            for (iDim = 0; iDim < nDim; iDim++)
              base_nodes->SetBound_Disp(iPoint, iDim, bufDRecv[buf_offset+iDim]);
            break;
          case SOLUTION_TIME_N:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->Set_Solution_time_n(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          case SOLUTION_TIME_N1:
            for (iVar = 0; iVar < nVar; iVar++)
              base_nodes->Set_Solution_time_n1(iPoint, iVar, bufDRecv[buf_offset+iVar]);
            break;
          default:
            SU2_MPI::Error("Unrecognized quantity for point-to-point MPI comms.",
                           CURRENT_FUNCTION);
            break;
        }
      }
    }
    
    /*--- Verify that all non-blocking point-to-point sends have finished.
     Note that this should be satisfied, as we have received all of the
     data in the loop above at this point. ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Waitall(geometry->nP2PSend, geometry->req_P2PSend, MPI_STATUS_IGNORE);
#endif
    
  }
  
}

void CSolver::SetResidual_RMS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;
  
#ifndef HAVE_MPI
  
  for (iVar = 0; iVar < nVar; iVar++) {
    
    if (GetRes_RMS(iVar) != GetRes_RMS(iVar)) {
        SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);
    }
    if (log10(sqrt(GetRes_RMS(iVar)/geometry->GetnPoint())) > 20 ){
      SU2_MPI::Error("SU2 has diverged. (Residual > 10^20 detected)", CURRENT_FUNCTION);
    }

    SetRes_RMS(iVar, max(EPS*EPS, sqrt(GetRes_RMS(iVar)/geometry->GetnPoint())));
    
  }
  
#else
  
  int nProcessor = size, iProcessor;

  su2double *sbuf_residual, *rbuf_residual, *sbuf_coord, *rbuf_coord, *Coord;
  unsigned long *sbuf_point, *rbuf_point, Global_nPointDomain;
  unsigned short iDim;
  
  /*--- Set the L2 Norm residual in all the processors ---*/
  
  sbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;
  
  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_RMS(iVar);
  
  if (config->GetComm_Level() == COMM_FULL) {
    
    unsigned long Local_nPointDomain = geometry->GetnPointDomain();
    SU2_MPI::Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    
  } else {
    
    /*--- Reduced MPI comms have been requested. Use a local residual only. ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = sbuf_residual[iVar];
    Global_nPointDomain = geometry->GetnPointDomain();
    
  }
  
  
  for (iVar = 0; iVar < nVar; iVar++) {
    
    if (rbuf_residual[iVar] != rbuf_residual[iVar]) {
      SU2_MPI::Error("SU2 has diverged. (NaN detected)", CURRENT_FUNCTION);
    }
    
    SetRes_RMS(iVar, max(EPS*EPS, sqrt(rbuf_residual[iVar]/Global_nPointDomain)));
    
  }

  delete [] sbuf_residual;
  delete [] rbuf_residual;
  
  /*--- Set the Maximum residual in all the processors ---*/
  
  if (config->GetComm_Level() == COMM_FULL) {
    
    sbuf_residual = new su2double [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
    sbuf_point = new unsigned long [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_point[iVar] = 0;
    sbuf_coord = new su2double[nVar*nDim]; for (iVar = 0; iVar < nVar*nDim; iVar++) sbuf_coord[iVar] = 0.0;
    
    rbuf_residual = new su2double [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_residual[iVar] = 0.0;
    rbuf_point = new unsigned long [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_point[iVar] = 0;
    rbuf_coord = new su2double[nProcessor*nVar*nDim]; for (iVar = 0; iVar < nProcessor*nVar*nDim; iVar++) rbuf_coord[iVar] = 0.0;
    
    for (iVar = 0; iVar < nVar; iVar++) {
      sbuf_residual[iVar] = GetRes_Max(iVar);
      sbuf_point[iVar] = GetPoint_Max(iVar);
      Coord = GetPoint_Max_Coord(iVar);
      for (iDim = 0; iDim < nDim; iDim++)
        sbuf_coord[iVar*nDim+iDim] = Coord[iDim];
    }
    
    SU2_MPI::Allgather(sbuf_residual, nVar, MPI_DOUBLE, rbuf_residual, nVar, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(sbuf_point, nVar, MPI_UNSIGNED_LONG, rbuf_point, nVar, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(sbuf_coord, nVar*nDim, MPI_DOUBLE, rbuf_coord, nVar*nDim, MPI_DOUBLE, MPI_COMM_WORLD);
    
    for (iVar = 0; iVar < nVar; iVar++) {
      for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
        AddRes_Max(iVar, rbuf_residual[iProcessor*nVar+iVar], rbuf_point[iProcessor*nVar+iVar], &rbuf_coord[iProcessor*nVar*nDim+iVar*nDim]);
      }
    }
    
    delete [] sbuf_residual;
    delete [] rbuf_residual;
    
    delete [] sbuf_point;
    delete [] rbuf_point;
    
    delete [] sbuf_coord;
    delete [] rbuf_coord;
    
  }
  
#endif
  
}

void CSolver::SetResidual_BGS(CGeometry *geometry, CConfig *config) {
  unsigned short iVar;

#ifndef HAVE_MPI

  for (iVar = 0; iVar < nVar; iVar++) {

//    if (GetRes_BGS(iVar) != GetRes_BGS(iVar)) {
//      SU2_MPI::Error("SU2 has diverged.", CURRENT_FUNCTION);
//    }

    SetRes_BGS(iVar, max(EPS*EPS, sqrt(GetRes_BGS(iVar)/geometry->GetnPoint())));

  }

#else

  int nProcessor = size, iProcessor;

  su2double *sbuf_residual, *rbuf_residual, *sbuf_coord, *rbuf_coord, *Coord;
  unsigned long *sbuf_point, *rbuf_point, Local_nPointDomain, Global_nPointDomain;
  unsigned short iDim;

  /*--- Set the L2 Norm residual in all the processors ---*/

  sbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  rbuf_residual  = new su2double[nVar]; for (iVar = 0; iVar < nVar; iVar++) rbuf_residual[iVar] = 0.0;

  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_BGS(iVar);
  Local_nPointDomain = geometry->GetnPointDomain();


  SU2_MPI::Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);


  for (iVar = 0; iVar < nVar; iVar++) {

//    if (rbuf_residual[iVar] != rbuf_residual[iVar]) {

//      SU2_MPI::Error("SU2 has diverged (NaN detected)", CURRENT_FUNCTION);

//    }

    SetRes_BGS(iVar, max(EPS*EPS, sqrt(rbuf_residual[iVar]/Global_nPointDomain)));

  }

  delete [] sbuf_residual;
  delete [] rbuf_residual;

  /*--- Set the Maximum residual in all the processors ---*/
  sbuf_residual = new su2double [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = 0.0;
  sbuf_point = new unsigned long [nVar]; for (iVar = 0; iVar < nVar; iVar++) sbuf_point[iVar] = 0;
  sbuf_coord = new su2double[nVar*nDim]; for (iVar = 0; iVar < nVar*nDim; iVar++) sbuf_coord[iVar] = 0.0;

  rbuf_residual = new su2double [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_residual[iVar] = 0.0;
  rbuf_point = new unsigned long [nProcessor*nVar]; for (iVar = 0; iVar < nProcessor*nVar; iVar++) rbuf_point[iVar] = 0;
  rbuf_coord = new su2double[nProcessor*nVar*nDim]; for (iVar = 0; iVar < nProcessor*nVar*nDim; iVar++) rbuf_coord[iVar] = 0.0;

  for (iVar = 0; iVar < nVar; iVar++) {
    sbuf_residual[iVar] = GetRes_Max_BGS(iVar);
    sbuf_point[iVar] = GetPoint_Max_BGS(iVar);
    Coord = GetPoint_Max_Coord_BGS(iVar);
    for (iDim = 0; iDim < nDim; iDim++)
      sbuf_coord[iVar*nDim+iDim] = Coord[iDim];
  }

  SU2_MPI::Allgather(sbuf_residual, nVar, MPI_DOUBLE, rbuf_residual, nVar, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_point, nVar, MPI_UNSIGNED_LONG, rbuf_point, nVar, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(sbuf_coord, nVar*nDim, MPI_DOUBLE, rbuf_coord, nVar*nDim, MPI_DOUBLE, MPI_COMM_WORLD);

  for (iVar = 0; iVar < nVar; iVar++) {
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
      AddRes_Max_BGS(iVar, rbuf_residual[iProcessor*nVar+iVar], rbuf_point[iProcessor*nVar+iVar], &rbuf_coord[iProcessor*nVar*nDim+iVar*nDim]);
    }
  }

  delete [] sbuf_residual;
  delete [] rbuf_residual;

  delete [] sbuf_point;
  delete [] rbuf_point;

  delete [] sbuf_coord;
  delete [] rbuf_coord;

#endif

}

void CSolver::SetRotatingFrame_GCL(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, nDim = geometry->GetnDim(), iVar, nVar = GetnVar(), iMarker;
  unsigned long iVertex, iEdge;
  su2double ProjGridVel, *Normal;

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    const unsigned long iPoint = geometry->edge[iEdge]->GetNode(0);
    const unsigned long jPoint = geometry->edge[iEdge]->GetNode(1);

    /*--- Solution at each edge point ---*/

    su2double *Solution_i = base_nodes->GetSolution(iPoint);
    su2double *Solution_j = base_nodes->GetSolution(jPoint);

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = 0.5* (Solution_i[iVar] + Solution_j[iVar]);

    /*--- Grid Velocity at each edge point ---*/

    su2double *GridVel_i = geometry->node[iPoint]->GetGridVel();
    su2double *GridVel_j = geometry->node[jPoint]->GetGridVel();
    for (iDim = 0; iDim < nDim; iDim++)
      Vector[iDim] = 0.5* (GridVel_i[iDim] + GridVel_j[iDim]);

    Normal = geometry->edge[iEdge]->GetNormal();

    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += Vector[iDim]*Normal[iDim];

    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = ProjGridVel*Solution_i[iVar];

    LinSysRes.AddBlock(iPoint, Residual);

    for (iVar = 0; iVar < nVar; iVar++)
      Residual[iVar] = ProjGridVel*Solution_j[iVar];

    LinSysRes.SubtractBlock(jPoint, Residual);

  }

  /*--- Loop boundary edges ---*/

  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY)  &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
        const unsigned long Point = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Solution at each edge point ---*/

        su2double *Solution = base_nodes->GetSolution(Point);

        /*--- Grid Velocity at each edge point ---*/

        su2double *GridVel = geometry->node[Point]->GetGridVel();

        /*--- Summed normal components ---*/

        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();

        ProjGridVel = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          ProjGridVel += GridVel[iDim]*Normal[iDim];

        for (iVar = 0; iVar < nVar; iVar++)
          Residual[iVar] = ProjGridVel*Solution[iVar];

        LinSysRes.SubtractBlock(Point, Residual);
        
      }
    }
  }
  
}

void CSolver::SetAuxVar_Gradient_GG(CGeometry *geometry, CConfig *config) {
  
  unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
  unsigned short nDim = geometry->GetnDim(), iDim, iMarker;
  
  su2double AuxVar_Vertex, AuxVar_i, AuxVar_j, AuxVar_Average;
  su2double *Gradient, DualArea, Partial_Res, Grad_Val, *Normal;
  
  /*--- Set Gradient to Zero ---*/

  base_nodes->SetAuxVarGradientZero();
  
  /*--- Loop interior edges ---*/
  
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    AuxVar_i = base_nodes->GetAuxVar(iPoint);
    AuxVar_j = base_nodes->GetAuxVar(jPoint);
    
    Normal = geometry->edge[iEdge]->GetNormal();
    AuxVar_Average =  0.5 * ( AuxVar_i + AuxVar_j);
    for (iDim = 0; iDim < nDim; iDim++) {
      Partial_Res = AuxVar_Average*Normal[iDim];
      base_nodes->AddAuxVarGradient(iPoint, iDim, Partial_Res);
      base_nodes->SubtractAuxVarGradient(jPoint,iDim, Partial_Res);
    }
  }
  
  /*--- Loop boundary edges ---*/
  
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++)
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      Point = geometry->vertex[iMarker][iVertex]->GetNode();
      AuxVar_Vertex = base_nodes->GetAuxVar(Point);
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = AuxVar_Vertex*Normal[iDim];
        base_nodes->SubtractAuxVarGradient(Point,iDim, Partial_Res);
      }
    }
    }
  
  for (iPoint=0; iPoint<geometry->GetnPoint(); iPoint++)
    for (iDim = 0; iDim < nDim; iDim++) {
      Gradient = base_nodes->GetAuxVarGradient(iPoint);
      DualArea = geometry->node[iPoint]->GetVolume();
      Grad_Val = Gradient[iDim]/(DualArea+EPS);
      base_nodes->SetAuxVarGradient(iPoint, iDim, Grad_Val);
    }
  
  /*--- Gradient MPI ---*/
  
  InitiateComms(geometry, config, AUXVAR_GRADIENT);
  CompleteComms(geometry, config, AUXVAR_GRADIENT);

}

void CSolver::SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iNeigh;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j, weight, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, z11, z12, z13, z22, z23, z33, detR2, product;
  bool singular = false;
  
  su2double *Cvector = new su2double [nDim];
  
  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    AuxVar_i = base_nodes->GetAuxVar(iPoint);
    
    /*--- Inizialization of variables ---*/
    for (iDim = 0; iDim < nDim; iDim++)
      Cvector[iDim] = 0.0;
    
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      AuxVar_j = base_nodes->GetAuxVar(jPoint);
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (fabs(weight) > EPS) {
        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iDim = 0; iDim < nDim; iDim++)
          Cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/(weight);
      }
      
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (fabs(r11) < EPS) r11 = EPS;
    r11 = sqrt(r11);
    r12 = r12/r11;
    r22 = sqrt(r22-r12*r12);
    if (fabs(r22) < EPS) r22 = EPS;
    if (nDim == 3) {
      r13 = r13/r11;
      r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
      r33 = sqrt(r33-r23*r23-r13*r13);
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (fabs(detR2) < EPS) singular = true;
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    
    for (iDim = 0; iDim < nDim; iDim++) {
      product = 0.0;
      for (jDim = 0; jDim < nDim; jDim++)
        product += Smatrix[iDim][jDim]*Cvector[jDim];
      if (geometry->node[iPoint]->GetDomain())
        base_nodes->SetAuxVarGradient(iPoint, iDim, product);
    }
  }
  
  delete [] Cvector;
  
  /*--- Gradient MPI ---*/
  
  InitiateComms(geometry, config, AUXVAR_GRADIENT);
  CompleteComms(geometry, config, AUXVAR_GRADIENT);
  
}

void CSolver::SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config) {
  unsigned long Point = 0, iPoint = 0, jPoint = 0, iEdge, iVertex;
  unsigned short iVar, iDim, iMarker;
  su2double *Solution_Vertex, *Solution_i, *Solution_j, Solution_Average, **Gradient,
  Partial_Res, Grad_Val, *Normal, Vol;
  
  /*--- Set Gradient to Zero ---*/
  base_nodes->SetGradientZero();
  
  /*--- Loop interior edges ---*/
  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
    iPoint = geometry->edge[iEdge]->GetNode(0);
    jPoint = geometry->edge[iEdge]->GetNode(1);
    
    Solution_i = base_nodes->GetSolution(iPoint);
    Solution_j = base_nodes->GetSolution(jPoint);
    Normal = geometry->edge[iEdge]->GetNormal();
    for (iVar = 0; iVar< nVar; iVar++) {
      Solution_Average =  0.5 * (Solution_i[iVar] + Solution_j[iVar]);
      for (iDim = 0; iDim < nDim; iDim++) {
        Partial_Res = Solution_Average*Normal[iDim];
        if (geometry->node[iPoint]->GetDomain())
          base_nodes->AddGradient(iPoint, iVar, iDim, Partial_Res);
        if (geometry->node[jPoint]->GetDomain())
          base_nodes->SubtractGradient(jPoint,iVar, iDim, Partial_Res);
      }
    }
  }
  
  /*--- Loop boundary edges ---*/
  for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)) {
    for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++) {
      Point = geometry->vertex[iMarker][iVertex]->GetNode();
      Solution_Vertex = base_nodes->GetSolution(Point);
      Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
      for (iVar = 0; iVar < nVar; iVar++)
        for (iDim = 0; iDim < nDim; iDim++) {
          Partial_Res = Solution_Vertex[iVar]*Normal[iDim];
          if (geometry->node[Point]->GetDomain())
            base_nodes->SubtractGradient(Point,iVar, iDim, Partial_Res);
        }
    }
  }
  }
  
  /*--- Correct the gradient values for any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SOL_GG);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SOL_GG);
  }
  
  /*--- Compute gradient ---*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    /*--- Get the volume, which may include periodic components. ---*/
    
    Vol = (geometry->node[iPoint]->GetVolume() +
           geometry->node[iPoint]->GetPeriodicVolume());
    
    for (iVar = 0; iVar < nVar; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Gradient = base_nodes->GetGradient(iPoint);
        Grad_Val = Gradient[iVar][iDim] / (Vol+EPS);
        base_nodes->SetGradient(iPoint, iVar, iDim, Grad_Val);
      }
    }
    
  }
  
  /*--- Gradient MPI ---*/
  
  InitiateComms(geometry, config, SOLUTION_GRADIENT);
  CompleteComms(geometry, config, SOLUTION_GRADIENT);
  
}

void CSolver::SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iVar, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, *Solution_i, *Solution_j;
  su2double r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight;
  su2double detR2, z11, z12, z13, z22, z23, z33;
  bool singular = false;
  
  su2double **Cvector = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  
  /*--- Clear Rmatrix, which could eventually be computed once
     and stored for static meshes, as well as the gradient. ---*/

  base_nodes->SetRmatrixZero();
  base_nodes->SetGradientZero();

  /*--- Loop over points of the grid ---*/
  
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    /*--- Set the value of the singular ---*/
    singular = false;
    
    /*--- Get coordinates ---*/
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    
    /*--- Get consevative solution ---*/
    
    Solution_i = base_nodes->GetSolution(iPoint);
    
    /*--- Inizialization of variables ---*/
    
    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      
      Solution_j = base_nodes->GetSolution(jPoint);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      
      if (weight != 0.0) {
        
        base_nodes->AddRmatrix(iPoint,0, 0, (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight);
        base_nodes->AddRmatrix(iPoint,0, 1, (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight);
        base_nodes->AddRmatrix(iPoint,1, 1, (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight);
        
        if (nDim == 3) {
          base_nodes->AddRmatrix(iPoint,0, 2, (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight);
          base_nodes->AddRmatrix(iPoint,1, 2, (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight);
          base_nodes->AddRmatrix(iPoint,2, 1, (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight);
          base_nodes->AddRmatrix(iPoint,2, 2, (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight);
        }
        
        /*--- Entries of c:= transpose(A)*b ---*/
        
        for (iVar = 0; iVar < nVar; iVar++) {
          for (iDim = 0; iDim < nDim; iDim++) {
            base_nodes->AddGradient(iPoint, iVar,iDim, (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/weight);
          }
        }
        
      }
    }
  }
  
  /*--- Correct the gradient values for any periodic boundaries. ---*/
  
  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_SOL_LS);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_SOL_LS);
  }
  
  /*--- Second loop over points of the grid to compute final gradient ---*/
  
  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    
    /*--- Set the value of the singular ---*/
    
    singular = false;
    
    /*--- Entries of upper triangular matrix R ---*/
    
    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;
    
    r11 = base_nodes->GetRmatrix(iPoint,0,0);
    r12 = base_nodes->GetRmatrix(iPoint,0,1);
    r22 = base_nodes->GetRmatrix(iPoint,1,1);
    
    /*--- Entries of upper triangular matrix R ---*/
    
    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;
    
    if (nDim == 3) {
      r13   = base_nodes->GetRmatrix(iPoint,0,2);
      r23_a = base_nodes->GetRmatrix(iPoint,1,2);
      r23_b = base_nodes->GetRmatrix(iPoint,2,1);
      r33   = base_nodes->GetRmatrix(iPoint,2,2);
      
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }
    
    /*--- Compute determinant ---*/
    
    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);
    
    /*--- Detect singular matrices ---*/
    
    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }
    
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    
    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }
    
    /*--- Computation of the gradient: S*c ---*/
    
    for (iVar = 0; iVar < nVar; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        Cvector[iVar][iDim] = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          Cvector[iVar][iDim] += Smatrix[iDim][jDim]*base_nodes->GetGradient(iPoint, iVar, jDim);
        }
      }
    }
    
    for (iVar = 0; iVar < nVar; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        base_nodes->SetGradient(iPoint, iVar, iDim, Cvector[iVar][iDim]);
      }
    }
    
  }
  
  /*--- Deallocate memory ---*/
  
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Cvector[iVar];
  delete [] Cvector;
  
  /*--- Gradient MPI ---*/
  
  InitiateComms(geometry, config, SOLUTION_GRADIENT);
  CompleteComms(geometry, config, SOLUTION_GRADIENT);
  
}

void CSolver::Add_ExternalOld_To_Solution(CGeometry *geometry) {
  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    base_nodes->AddSolution(iPoint, base_nodes->Get_ExternalOld(iPoint));
  }
}

void CSolver::Add_Solution_To_External(CGeometry *geometry) {
  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    base_nodes->Add_External(iPoint, base_nodes->GetSolution(iPoint));
  }
}

void CSolver::Add_Solution_To_ExternalOld(CGeometry *geometry) {
  for (unsigned long iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    base_nodes->Add_ExternalOld(iPoint, base_nodes->GetSolution(iPoint));
  }
}

void CSolver::SetGridVel_Gradient(CGeometry *geometry, CConfig *config) {
  unsigned short iDim, jDim, iVar, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, *Solution_i, *Solution_j, Smatrix[3][3],
  r11, r12, r13, r22, r23, r23_a, r23_b, r33, weight, detR2, z11, z12, z13,
  z22, z23, z33, product;
  su2double **Cvector;
  
  /*--- Note that all nVar entries in this routine have been changed to nDim ---*/
  Cvector = new su2double* [nDim];
  for (iVar = 0; iVar < nDim; iVar++)
    Cvector[iVar] = new su2double [nDim];
  
  /*--- Loop over points of the grid ---*/
  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    
    Coord_i = geometry->node[iPoint]->GetCoord();
    Solution_i = geometry->node[iPoint]->GetGridVel();
    
    /*--- Inizialization of variables ---*/
    for (iVar = 0; iVar < nDim; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;
    r11 = 0.0; r12 = 0.0; r13 = 0.0; r22 = 0.0; r23 = 0.0; r23_a = 0.0; r23_b = 0.0; r33 = 0.0;
    
    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();
      Solution_j = geometry->node[jPoint]->GetGridVel();
      
      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      
      /*--- Sumations for entries of upper triangular matrix R ---*/
      r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/(weight);
      r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/(weight);
      r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/(weight);
      if (nDim == 3) {
        r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
        r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/(weight);
        r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/(weight);
        r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/(weight);
      }
      
      /*--- Entries of c:= transpose(A)*b ---*/
      for (iVar = 0; iVar < nDim; iVar++)
        for (iDim = 0; iDim < nDim; iDim++)
          Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(Solution_j[iVar]-Solution_i[iVar])/(weight);
    }
    
    /*--- Entries of upper triangular matrix R ---*/
    r11 = sqrt(r11);
    r12 = r12/(r11);
    r22 = sqrt(r22-r12*r12);
    if (nDim == 3) {
      r13 = r13/(r11);
      r23 = r23_a/(r22) - r23_b*r12/(r11*r22);
      r33 = sqrt(r33-r23*r23-r13*r13);
    }
    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    if (nDim == 2) {
      detR2 = (r11*r22)*(r11*r22);
      Smatrix[0][0] = (r12*r12+r22*r22)/(detR2);
      Smatrix[0][1] = -r11*r12/(detR2);
      Smatrix[1][0] = Smatrix[0][1];
      Smatrix[1][1] = r11*r11/(detR2);
    }
    else {
      detR2 = (r11*r22*r33)*(r11*r22*r33);
      z11 = r22*r33;
      z12 = -r12*r33;
      z13 = r12*r23-r13*r22;
      z22 = r11*r33;
      z23 = -r11*r23;
      z33 = r11*r22;
      Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/(detR2);
      Smatrix[0][1] = (z12*z22+z13*z23)/(detR2);
      Smatrix[0][2] = (z13*z33)/(detR2);
      Smatrix[1][0] = Smatrix[0][1];
      Smatrix[1][1] = (z22*z22+z23*z23)/(detR2);
      Smatrix[1][2] = (z23*z33)/(detR2);
      Smatrix[2][0] = Smatrix[0][2];
      Smatrix[2][1] = Smatrix[1][2];
      Smatrix[2][2] = (z33*z33)/(detR2);
    }
    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nDim; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++)
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        geometry->node[iPoint]->SetGridVel_Grad(iVar, iDim, product);
      }
    }
  }
  
  /*--- Deallocate memory ---*/
  for (iVar = 0; iVar < nDim; iVar++)
    delete [] Cvector[iVar];
  delete [] Cvector;
  
}

void CSolver::SetAuxVar_Surface_Gradient(CGeometry *geometry, CConfig *config) {
  
  unsigned short iDim, jDim, iNeigh, iMarker, Boundary;
  unsigned short nDim = geometry->GetnDim();
  unsigned long iPoint, jPoint, iVertex;
  su2double *Coord_i, *Coord_j, AuxVar_i, AuxVar_j;
  su2double **Smatrix, *Cvector;
  
  Smatrix = new su2double* [nDim];
  Cvector = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];
  
  
  /*--- Loop over boundary markers to select those for Euler or NS walls ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    switch (Boundary) {
      case EULER_WALL:
      case HEAT_FLUX:
      case ISOTHERMAL:
      case CHT_WALL_INTERFACE:
        
        /*--- Loop over points on the surface (Least-Squares approximation) ---*/
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          if (geometry->node[iPoint]->GetDomain()) {
            Coord_i = geometry->node[iPoint]->GetCoord();
            AuxVar_i = base_nodes->GetAuxVar(iPoint);
            
            /*--- Inizialization of variables ---*/
            for (iDim = 0; iDim < nDim; iDim++)
              Cvector[iDim] = 0.0;
            su2double r11 = 0.0, r12 = 0.0, r13 = 0.0, r22 = 0.0, r23 = 0.0, r23_a = 0.0, r23_b = 0.0, r33 = 0.0;
            
            for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
              Coord_j = geometry->node[jPoint]->GetCoord();
              AuxVar_j = base_nodes->GetAuxVar(jPoint);
              
              su2double weight = 0;
              for (iDim = 0; iDim < nDim; iDim++)
                weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
              
              /*--- Sumations for entries of upper triangular matrix R ---*/
              r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
              r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
              r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;
              if (nDim == 3) {
                r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
                r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
                r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
                r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
              }
              
              /*--- Entries of c:= transpose(A)*b ---*/
              for (iDim = 0; iDim < nDim; iDim++)
                Cvector[iDim] += (Coord_j[iDim]-Coord_i[iDim])*(AuxVar_j-AuxVar_i)/weight;
            }
            
            /*--- Entries of upper triangular matrix R ---*/
            r11 = sqrt(r11);
            r12 = r12/r11;
            r22 = sqrt(r22-r12*r12);
            if (nDim == 3) {
              r13 = r13/r11;
              r23 = r23_a/r22 - r23_b*r12/(r11*r22);
              r33 = sqrt(r33-r23*r23-r13*r13);
            }
            /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
            if (nDim == 2) {
              su2double detR2 = (r11*r22)*(r11*r22);
              Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
              Smatrix[0][1] = -r11*r12/detR2;
              Smatrix[1][0] = Smatrix[0][1];
              Smatrix[1][1] = r11*r11/detR2;
            }
            else {
              su2double detR2 = (r11*r22*r33)*(r11*r22*r33);
              su2double z11, z12, z13, z22, z23, z33; // aux vars
              z11 = r22*r33;
              z12 = -r12*r33;
              z13 = r12*r23-r13*r22;
              z22 = r11*r33;
              z23 = -r11*r23;
              z33 = r11*r22;
              Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
              Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
              Smatrix[0][2] = (z13*z33)/detR2;
              Smatrix[1][0] = Smatrix[0][1];
              Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
              Smatrix[1][2] = (z23*z33)/detR2;
              Smatrix[2][0] = Smatrix[0][2];
              Smatrix[2][1] = Smatrix[1][2];
              Smatrix[2][2] = (z33*z33)/detR2;
            }
            /*--- Computation of the gradient: S*c ---*/
            su2double product;
            for (iDim = 0; iDim < nDim; iDim++) {
              product = 0.0;
              for (jDim = 0; jDim < nDim; jDim++)
                product += Smatrix[iDim][jDim]*Cvector[jDim];
              base_nodes->SetAuxVarGradient(iPoint, iDim, product);
            }
          }
        } /*--- End of loop over surface points ---*/
        break;
      default:
        break;
    }
  }
  
  /*--- Memory deallocation ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    delete [] Smatrix[iDim];
  delete [] Cvector;
  delete [] Smatrix;
}

void CSolver::SetSolution_Limiter(CGeometry *geometry, CConfig *config) {
  
  unsigned long iEdge, iPoint, jPoint;
  unsigned short iVar, iDim;
  su2double **Gradient_i, **Gradient_j, *Coord_i, *Coord_j,
  *Solution, *Solution_i, *Solution_j,
  *LocalMinSolution = NULL, *LocalMaxSolution = NULL,
  *GlobalMinSolution = NULL, *GlobalMaxSolution = NULL,
  dave, LimK, eps1, eps2, dm, dp, du, ds, y, limiter, SharpEdge_Distance;
  
#ifdef CODI_REVERSE_TYPE
  bool TapeActive = false;

  if (config->GetDiscrete_Adjoint() && config->GetFrozen_Limiter_Disc()) {
    /*--- If limiters are frozen do not record the computation ---*/
    TapeActive = AD::globalTape.isActive();
    AD::StopRecording();
  }
#endif
  
  dave = config->GetRefElemLength();
  LimK = config->GetVenkat_LimiterCoeff();
  
  if (config->GetKind_SlopeLimit() == NO_LIMITER) {
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        base_nodes->SetLimiter(iPoint, iVar, 1.0);
      }
    }
    
  }
  
  else {
    
    /*--- Initialize solution max and solution min and the limiter in the entire domain --*/
    
    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        base_nodes->SetSolution_Max(iPoint, iVar, -EPS);
        base_nodes->SetSolution_Min(iPoint, iVar, EPS);
        base_nodes->SetLimiter(iPoint, iVar, 2.0);
      }
    }
    
    /*--- Establish bounds for Spekreijse monotonicity by finding max & min values of neighbor variables --*/
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      /*--- Point identification, Normal vector and area ---*/
      
      iPoint = geometry->edge[iEdge]->GetNode(0);
      jPoint = geometry->edge[iEdge]->GetNode(1);
      
      /*--- Get the conserved variables ---*/
      
      Solution_i = base_nodes->GetSolution(iPoint);
      Solution_j = base_nodes->GetSolution(jPoint);
      
      /*--- Compute the maximum, and minimum values for nodes i & j ---*/
      
      for (iVar = 0; iVar < nVar; iVar++) {
        du = (Solution_j[iVar] - Solution_i[iVar]);
        base_nodes->SetSolution_Min(iPoint, iVar, min(base_nodes->GetSolution_Min(iPoint, iVar), du));
        base_nodes->SetSolution_Max(iPoint, iVar, max(base_nodes->GetSolution_Max(iPoint, iVar), du));
        base_nodes->SetSolution_Min(jPoint, iVar, min(base_nodes->GetSolution_Min(jPoint, iVar), -du));
        base_nodes->SetSolution_Max(jPoint, iVar, max(base_nodes->GetSolution_Max(jPoint, iVar), -du));
      }
      
    }
    
    /*--- Correct the limiter values across any periodic boundaries. ---*/
    
    for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
      InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_SOL_1);
      CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_SOL_1);
    }
    
  }
  
  /*--- Barth-Jespersen limiter with Venkatakrishnan modification ---*/
  
  if (config->GetKind_SlopeLimit_Flow() == BARTH_JESPERSEN) {
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = base_nodes->GetGradient(iPoint);
      Gradient_j = base_nodes->GetGradient(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nVar, nDim);
      AD::SetPreaccIn(Gradient_j, nVar, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nVar; iVar++) {
        
        AD::SetPreaccIn(base_nodes->GetSolution_Max(iPoint, iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Min(iPoint, iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Max(jPoint,iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Min(jPoint,iVar));

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(iPoint, iVar);
          else dp = base_nodes->GetSolution_Min(iPoint, iVar);
          limiter = dp/dm;
        }
        
        if (limiter < base_nodes->GetLimiter(iPoint, iVar)) {
          base_nodes->SetLimiter(iPoint, iVar, limiter);
          AD::SetPreaccOut(base_nodes->GetLimiter(iPoint)[iVar]);
        }
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if (dm == 0.0) { limiter = 2.0; }
        else {
          if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(jPoint,iVar);
          else dp = base_nodes->GetSolution_Min(jPoint,iVar);
          limiter = dp/dm;
        }
        
        if (limiter < base_nodes->GetLimiter(jPoint,iVar)) {
          base_nodes->SetLimiter(jPoint,iVar, limiter);
          AD::SetPreaccOut(base_nodes->GetLimiter(jPoint)[iVar]);
        }

      }
      
      AD::EndPreacc();
      
    }


    for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        y =  base_nodes->GetLimiter(iPoint, iVar);
        limiter = (y*y + 2.0*y) / (y*y + y + 2.0);
        base_nodes->SetLimiter(iPoint, iVar, limiter);
      }
    }
    
  }

  /*--- Venkatakrishnan limiter ---*/
  
  if ((config->GetKind_SlopeLimit() == VENKATAKRISHNAN) || (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG)) {
    
    if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {

      /*--- Allocate memory for the max and min solution value --*/
      
      LocalMinSolution = new su2double [nVar]; GlobalMinSolution = new su2double [nVar];
      LocalMaxSolution = new su2double [nVar]; GlobalMaxSolution = new su2double [nVar];
      
      /*--- Compute the max value and min value of the solution ---*/
      
      Solution = base_nodes->GetSolution(iPoint);
      for (iVar = 0; iVar < nVar; iVar++) {
        LocalMinSolution[iVar] = Solution[iVar];
        LocalMaxSolution[iVar] = Solution[iVar];
      }
      
      for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
        
        /*--- Get the solution variables ---*/
        
        Solution = base_nodes->GetSolution(iPoint);
        
        for (iVar = 0; iVar < nVar; iVar++) {
          LocalMinSolution[iVar] = min (LocalMinSolution[iVar], Solution[iVar]);
          LocalMaxSolution[iVar] = max (LocalMaxSolution[iVar], Solution[iVar]);
        }
        
      }
      
#ifdef HAVE_MPI
      SU2_MPI::Allreduce(LocalMinSolution, GlobalMinSolution, nVar, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(LocalMaxSolution, GlobalMaxSolution, nVar, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
      for (iVar = 0; iVar < nVar; iVar++) {
        GlobalMinSolution[iVar] = LocalMinSolution[iVar];
        GlobalMaxSolution[iVar] = LocalMaxSolution[iVar];
      }
#endif
    }
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = base_nodes->GetGradient(iPoint);
      Gradient_j = base_nodes->GetGradient(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      AD::StartPreacc();
      AD::SetPreaccIn(Gradient_i, nVar, nDim);
      AD::SetPreaccIn(Gradient_j, nVar, nDim);
      AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);

      for (iVar = 0; iVar < nVar; iVar++) {
          
        AD::StartPreacc();
        AD::SetPreaccIn(Gradient_i[iVar], nDim);
        AD::SetPreaccIn(Gradient_j[iVar], nDim);
        AD::SetPreaccIn(Coord_i, nDim);
        AD::SetPreaccIn(Coord_j, nDim);
        AD::SetPreaccIn(base_nodes->GetSolution_Max(iPoint, iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Min(iPoint, iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Max(jPoint,iVar));
        AD::SetPreaccIn(base_nodes->GetSolution_Min(jPoint,iVar));
        
        if (config->GetKind_SlopeLimit_Flow() == VENKATAKRISHNAN_WANG) {
          AD::SetPreaccIn(GlobalMaxSolution[iVar]);
          AD::SetPreaccIn(GlobalMinSolution[iVar]);
          eps1 = LimK * (GlobalMaxSolution[iVar] - GlobalMinSolution[iVar]);
          eps2 = eps1*eps1;
        }
        else {
          eps1 = LimK*dave;
          eps2 = eps1*eps1*eps1;
        }

        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(iPoint, iVar);
        else dp = base_nodes->GetSolution_Min(iPoint, iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(iPoint, iVar)) {
          base_nodes->SetLimiter(iPoint, iVar, limiter);
          AD::SetPreaccOut(base_nodes->GetLimiter(iPoint)[iVar]);
        }
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(jPoint,iVar);
        else dp = base_nodes->GetSolution_Min(jPoint,iVar);
        
        limiter = ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(jPoint,iVar)) {
          base_nodes->SetLimiter(jPoint,iVar, limiter);
          AD::SetPreaccOut(base_nodes->GetLimiter(jPoint)[iVar]);
        }
        
        AD::EndPreacc();
      }
    }
    
    if (LocalMinSolution  != NULL) delete [] LocalMinSolution;
    if (LocalMaxSolution  != NULL) delete [] LocalMaxSolution;
    if (GlobalMinSolution != NULL) delete [] GlobalMinSolution;
    if (GlobalMaxSolution != NULL) delete [] GlobalMaxSolution;

  }
  
  /*--- Sharp edges limiter ---*/
  
  if (config->GetKind_SlopeLimit() == SHARP_EDGES) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetVenkat_LimiterCoeff();
    eps1 = LimK*dave;
    eps2 = eps1*eps1*eps1;
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = base_nodes->GetGradient(iPoint);
      Gradient_j = base_nodes->GetGradient(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(iPoint, iVar);
        else dp = base_nodes->GetSolution_Min(iPoint, iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[iPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(iPoint, iVar))
          base_nodes->SetLimiter(iPoint, iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(jPoint,iVar);
        else dp = base_nodes->GetSolution_Min(jPoint,iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[jPoint]->GetSharpEdge_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(jPoint,iVar))
          base_nodes->SetLimiter(jPoint,iVar, limiter);
        
      }
    }
  }
  
  /*--- Sharp edges limiter ---*/
  
  if (config->GetKind_SlopeLimit() == WALL_DISTANCE) {
    
    /*-- Get limiter parameters from the configuration file ---*/
    
    dave = config->GetRefElemLength();
    LimK = config->GetVenkat_LimiterCoeff();
    eps1 = LimK*dave;
    eps2 = eps1*eps1*eps1;
    
    for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {
      
      iPoint     = geometry->edge[iEdge]->GetNode(0);
      jPoint     = geometry->edge[iEdge]->GetNode(1);
      Gradient_i = base_nodes->GetGradient(iPoint);
      Gradient_j = base_nodes->GetGradient(jPoint);
      Coord_i    = geometry->node[iPoint]->GetCoord();
      Coord_j    = geometry->node[jPoint]->GetCoord();
      
      for (iVar = 0; iVar < nVar; iVar++) {
        
        /*--- Calculate the interface left gradient, delta- (dm) ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_j[iDim]-Coord_i[iDim])*Gradient_i[iVar][iDim];
        
        /*--- Calculate the interface right gradient, delta+ (dp) ---*/
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(iPoint, iVar);
        else dp = base_nodes->GetSolution_Min(iPoint, iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[iPoint]->GetWall_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(iPoint, iVar))
          base_nodes->SetLimiter(iPoint, iVar, limiter);
        
        /*-- Repeat for point j on the edge ---*/
        
        dm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++)
          dm += 0.5*(Coord_i[iDim]-Coord_j[iDim])*Gradient_j[iVar][iDim];
        
        if ( dm > 0.0 ) dp = base_nodes->GetSolution_Max(jPoint,iVar);
        else dp = base_nodes->GetSolution_Min(jPoint,iVar);
        
        /*--- Compute the distance to a sharp edge ---*/
        
        SharpEdge_Distance = (geometry->node[jPoint]->GetWall_Distance() - config->GetAdjSharp_LimiterCoeff()*eps1);
        ds = 0.0;
        if (SharpEdge_Distance < -eps1) ds = 0.0;
        if (fabs(SharpEdge_Distance) <= eps1) ds = 0.5*(1.0+(SharpEdge_Distance/eps1)+(1.0/PI_NUMBER)*sin(PI_NUMBER*SharpEdge_Distance/eps1));
        if (SharpEdge_Distance > eps1) ds = 1.0;
        
        limiter = ds * ( dp*dp + 2.0*dp*dm + eps2 )/( dp*dp + dp*dm + 2.0*dm*dm + eps2);
        
        if (limiter < base_nodes->GetLimiter(jPoint,iVar))
          base_nodes->SetLimiter(jPoint,iVar, limiter);
        
      }
    }
  }

  /*--- Correct the limiter values across any periodic boundaries. ---*/

  for (unsigned short iPeriodic = 1; iPeriodic <= config->GetnMarker_Periodic()/2; iPeriodic++) {
    InitiatePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_SOL_2);
    CompletePeriodicComms(geometry, config, iPeriodic, PERIODIC_LIM_SOL_2);
  }
  
  /*--- Limiter MPI ---*/
  
  InitiateComms(geometry, config, SOLUTION_LIMITER);
  CompleteComms(geometry, config, SOLUTION_LIMITER);

#ifdef CODI_REVERSE_TYPE
  if (TapeActive) AD::StartRecording();
#endif
}

void CSolver::Gauss_Elimination(su2double** A, su2double* rhs, unsigned short nVar) {
  
  short iVar, jVar, kVar;
  su2double weight, aux;
  
  if (nVar == 1)
    rhs[0] /= A[0][0];
  else {
    
    /*--- Transform system in Upper Matrix ---*/
    
    for (iVar = 1; iVar < (short)nVar; iVar++) {
      for (jVar = 0; jVar < iVar; jVar++) {
        weight = A[iVar][jVar]/A[jVar][jVar];
        for (kVar = jVar; kVar < (short)nVar; kVar++)
          A[iVar][kVar] -= weight*A[jVar][kVar];
        rhs[iVar] -= weight*rhs[jVar];
      }
    }
    
    /*--- Backwards substitution ---*/
    
    rhs[nVar-1] = rhs[nVar-1]/A[nVar-1][nVar-1];
    for (iVar = (short)nVar-2; iVar >= 0; iVar--) {
      aux = 0;
      for (jVar = iVar+1; jVar < (short)nVar; jVar++)
        aux += A[iVar][jVar]*rhs[jVar];
      rhs[iVar] = (rhs[iVar]-aux)/A[iVar][iVar];
      if (iVar == 0) break;
    }
  }
  
}

void CSolver::Aeroelastic(CSurfaceMovement *surface_movement, CGeometry *geometry, CConfig *config, unsigned long TimeIter) {
  
  /*--- Variables used for Aeroelastic case ---*/
  
  su2double Cl, Cd, Cn, Ct, Cm, Cn_rot;
  su2double Alpha = config->GetAoA()*PI_NUMBER/180.0;
  vector<su2double> structural_solution(4,0.0); //contains solution(displacements and rates) of typical section wing model.
  
  unsigned short iMarker, iMarker_Monitoring, Monitoring;
  string Marker_Tag, Monitoring_Tag;
  
  /*--- Loop over markers and find the ones being monitored. ---*/
  
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    Monitoring = config->GetMarker_All_Monitoring(iMarker);
    if (Monitoring == YES) {
      
      /*--- Find the particular marker being monitored and get the forces acting on it. ---*/
      
      for (iMarker_Monitoring = 0; iMarker_Monitoring < config->GetnMarker_Monitoring(); iMarker_Monitoring++) {
        Monitoring_Tag = config->GetMarker_Monitoring_TagBound(iMarker_Monitoring);
        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        if (Marker_Tag == Monitoring_Tag) {
          
          Cl = GetSurface_CL(iMarker_Monitoring);
          Cd = GetSurface_CD(iMarker_Monitoring);
          
          /*--- For typical section wing model want the force normal to the airfoil (in the direction of the spring) ---*/
          Cn = Cl*cos(Alpha) + Cd*sin(Alpha);
          Ct = -Cl*sin(Alpha) + Cd*cos(Alpha);
          
          Cm = GetSurface_CMz(iMarker_Monitoring);
          
          /*--- Calculate forces for the Typical Section Wing Model taking into account rotation ---*/
          
          /*--- Note that the calculation of the forces and the subsequent displacements ...
           is only correct for the airfoil that starts at the 0 degree position ---*/
          
          if (config->GetKind_GridMovement() == AEROELASTIC_RIGID_MOTION) {
            su2double Omega, dt, psi;
            dt = config->GetDelta_UnstTimeND();
            Omega  = (config->GetRotation_Rate(2)/config->GetOmega_Ref());
            psi = Omega*(dt*TimeIter);
            
            /*--- Correct for the airfoil starting position (This is hardcoded in here) ---*/
            if (Monitoring_Tag == "Airfoil1") {
              psi = psi + 0.0;
            }
            else if (Monitoring_Tag == "Airfoil2") {
              psi = psi + 2.0/3.0*PI_NUMBER;
            }
            else if (Monitoring_Tag == "Airfoil3") {
              psi = psi + 4.0/3.0*PI_NUMBER;
            }
            else
              cout << "WARNING: There is a marker that we are monitoring that doesn't match the values hardcoded above!" << endl;
            
            cout << Monitoring_Tag << " position " << psi*180.0/PI_NUMBER << " degrees. " << endl;
            
            Cn_rot = Cn*cos(psi) - Ct*sin(psi); //Note the signs are different for accounting for the AOA.
            Cn = Cn_rot;
          }
          
          /*--- Solve the aeroelastic equations for the particular marker(surface) ---*/
          
          SolveTypicalSectionWingModel(geometry, Cn, Cm, config, iMarker_Monitoring, structural_solution);
          
          break;
        }
      }
      
      /*--- Compute the new surface node locations ---*/
      surface_movement->AeroelasticDeform(geometry, config, TimeIter, iMarker, iMarker_Monitoring, structural_solution);
      
    }
    
  }
  
}

void CSolver::SetUpTypicalSectionWingModel(vector<vector<su2double> >& Phi, vector<su2double>& omega, CConfig *config) {
  
  /*--- Retrieve values from the config file ---*/
  su2double w_h = config->GetAeroelastic_Frequency_Plunge();
  su2double w_a = config->GetAeroelastic_Frequency_Pitch();
  su2double x_a = config->GetAeroelastic_CG_Location();
  su2double r_a = sqrt(config->GetAeroelastic_Radius_Gyration_Squared());
  su2double w = w_h/w_a;
  
  // Mass Matrix
  vector<vector<su2double> > M(2,vector<su2double>(2,0.0));
  M[0][0] = 1;
  M[0][1] = x_a;
  M[1][0] = x_a;
  M[1][1] = r_a*r_a;
  
  // Stiffness Matrix
  //  vector<vector<su2double> > K(2,vector<su2double>(2,0.0));
  //  K[0][0] = (w_h/w_a)*(w_h/w_a);
  //  K[0][1] = 0.0;
  //  K[1][0] = 0.0;
  //  K[1][1] = r_a*r_a;
  
  /* Eigenvector and Eigenvalue Matrices of the Generalized EigenValue Problem. */
  
  vector<vector<su2double> > Omega2(2,vector<su2double>(2,0.0));
  su2double aux; // auxiliary variable
  aux = sqrt(pow(r_a,2)*pow(w,4) - 2*pow(r_a,2)*pow(w,2) + pow(r_a,2) + 4*pow(x_a,2)*pow(w,2));
  Phi[0][0] = (r_a * (r_a - r_a*pow(w,2) + aux)) / (2*x_a*pow(w, 2));
  Phi[0][1] = (r_a * (r_a - r_a*pow(w,2) - aux)) / (2*x_a*pow(w, 2));
  Phi[1][0] = 1.0;
  Phi[1][1] = 1.0;
  
  Omega2[0][0] = (r_a * (r_a + r_a*pow(w,2) - aux)) / (2*(pow(r_a, 2) - pow(x_a, 2)));
  Omega2[0][1] = 0;
  Omega2[1][0] = 0;
  Omega2[1][1] = (r_a * (r_a + r_a*pow(w,2) + aux)) / (2*(pow(r_a, 2) - pow(x_a, 2)));
  
  /* Nondimesionalize the Eigenvectors such that Phi'*M*Phi = I and PHI'*K*PHI = Omega */
  // Phi'*M*Phi = D
  // D^(-1/2)*Phi'*M*Phi*D^(-1/2) = D^(-1/2)*D^(1/2)*D^(1/2)*D^(-1/2) = I,  D^(-1/2) = inv(sqrt(D))
  // Phi = Phi*D^(-1/2)
  
  vector<vector<su2double> > Aux(2,vector<su2double>(2,0.0));
  vector<vector<su2double> > D(2,vector<su2double>(2,0.0));
  // Aux = M*Phi
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      Aux[i][j] = 0;
      for (int k=0; k<2; k++) {
        Aux[i][j] += M[i][k]*Phi[k][j];
      }
    }
  }
  
  // D = Phi'*Aux
  for (int i=0; i<2; i++) {
    for (int j=0; j<2; j++) {
      D[i][j] = 0;
      for (int k=0; k<2; k++) {
        D[i][j] += Phi[k][i]*Aux[k][j]; //PHI transpose
      }
    }
  }
  
  //Modify the first column
  Phi[0][0] = Phi[0][0] * 1/sqrt(D[0][0]);
  Phi[1][0] = Phi[1][0] * 1/sqrt(D[0][0]);
  //Modify the second column
  Phi[0][1] = Phi[0][1] * 1/sqrt(D[1][1]);
  Phi[1][1] = Phi[1][1] * 1/sqrt(D[1][1]);
  
  // Sqrt of the eigenvalues (frequency of vibration of the modes)
  omega[0] = sqrt(Omega2[0][0]);
  omega[1] = sqrt(Omega2[1][1]);
  
}

void CSolver::SolveTypicalSectionWingModel(CGeometry *geometry, su2double Cl, su2double Cm, CConfig *config, unsigned short iMarker, vector<su2double>& displacements) {
  
  /*--- The aeroelastic model solved in this routine is the typical section wing model
   The details of the implementation are similar to those found in J.J. Alonso 
   "Fully-Implicit Time-Marching Aeroelastic Solutions" 1994. ---*/
  
  /*--- Retrieve values from the config file ---*/
  su2double w_alpha = config->GetAeroelastic_Frequency_Pitch();
  su2double vf      = config->GetAeroelastic_Flutter_Speed_Index();
  su2double b       = config->GetLength_Reynolds()/2.0; // airfoil semichord, Reynolds length is by defaul 1.0
  su2double dt      = config->GetDelta_UnstTimeND();
  dt = dt*w_alpha; //Non-dimensionalize the structural time.
  
  /*--- Structural Equation damping ---*/
  vector<su2double> xi(2,0.0);
  
  /*--- Eigenvectors and Eigenvalues of the Generalized EigenValue Problem. ---*/
  vector<vector<su2double> > Phi(2,vector<su2double>(2,0.0));   // generalized eigenvectors.
  vector<su2double> w(2,0.0);        // sqrt of the generalized eigenvalues (frequency of vibration of the modes).
  SetUpTypicalSectionWingModel(Phi, w, config);
  
  /*--- Solving the Decoupled Aeroelastic Problem with second order time discretization Eq (9) ---*/
  
  /*--- Solution variables description. //x[j][i], j-entry, i-equation. // Time (n+1)->np1, n->n, (n-1)->n1 ---*/
  vector<vector<su2double> > x_np1(2,vector<su2double>(2,0.0));
  
  /*--- Values from previous movement of spring at true time step n+1
   We use this values because we are solving for delta changes not absolute changes ---*/
  vector<vector<su2double> > x_np1_old = config->GetAeroelastic_np1(iMarker);
  
  /*--- Values at previous timesteps. ---*/
  vector<vector<su2double> > x_n = config->GetAeroelastic_n(iMarker);
  vector<vector<su2double> > x_n1 = config->GetAeroelastic_n1(iMarker);
  
  /*--- Set up of variables used to solve the structural problem. ---*/
  vector<su2double> f_tilde(2,0.0);
  vector<vector<su2double> > A_inv(2,vector<su2double>(2,0.0));
  su2double detA;
  su2double s1, s2;
  vector<su2double> rhs(2,0.0); //right hand side
  vector<su2double> eta(2,0.0);
  vector<su2double> eta_dot(2,0.0);
  
  /*--- Forcing Term ---*/
  su2double cons = vf*vf/PI_NUMBER;
  vector<su2double> f(2,0.0);
  f[0] = cons*(-Cl);
  f[1] = cons*(2*-Cm);
  
  //f_tilde = Phi'*f
  for (int i=0; i<2; i++) {
    f_tilde[i] = 0;
    for (int k=0; k<2; k++) {
      f_tilde[i] += Phi[k][i]*f[k]; //PHI transpose
    }
  }
  
  /*--- solve each decoupled equation (The inverse of the 2x2 matrix is provided) ---*/
  for (int i=0; i<2; i++) {
    /* Matrix Inverse */
    detA = 9.0/(4.0*dt*dt) + 3*w[i]*xi[i]/(dt) + w[i]*w[i];
    A_inv[0][0] = 1/detA * (3/(2.0*dt) + 2*xi[i]*w[i]);
    A_inv[0][1] = 1/detA * 1;
    A_inv[1][0] = 1/detA * -w[i]*w[i];
    A_inv[1][1] = 1/detA * 3/(2.0*dt);
    
    /* Source Terms from previous iterations */
    s1 = (-4*x_n[0][i] + x_n1[0][i])/(2.0*dt);
    s2 = (-4*x_n[1][i] + x_n1[1][i])/(2.0*dt);
    
    /* Problem Right Hand Side */
    rhs[0] = -s1;
    rhs[1] = f_tilde[i]-s2;
    
    /* Solve the equations */
    x_np1[0][i] = A_inv[0][0]*rhs[0] + A_inv[0][1]*rhs[1];
    x_np1[1][i] = A_inv[1][0]*rhs[0] + A_inv[1][1]*rhs[1];
    
    eta[i] = x_np1[0][i]-x_np1_old[0][i];  // For displacements, the change(deltas) is used.
    eta_dot[i] = x_np1[1][i]; // For velocities, absolute values are used.
  }
  
  /*--- Transform back from the generalized coordinates to get the actual displacements in plunge and pitch  q = Phi*eta ---*/
  vector<su2double> q(2,0.0);
  vector<su2double> q_dot(2,0.0);
  for (int i=0; i<2; i++) {
    q[i] = 0;
    q_dot[i] = 0;
    for (int k=0; k<2; k++) {
      q[i] += Phi[i][k]*eta[k];
      q_dot[i] += Phi[i][k]*eta_dot[k];
    }
  }
  
  su2double dh = b*q[0];
  su2double dalpha = q[1];
  
  su2double h_dot = w_alpha*b*q_dot[0];  //The w_a brings it back to actual time.
  su2double alpha_dot = w_alpha*q_dot[1];
  
  /*--- Set the solution of the structural equations ---*/
  displacements[0] = dh;
  displacements[1] = dalpha;
  displacements[2] = h_dot;
  displacements[3] = alpha_dot;
  
  /*--- Calculate the total plunge and total pitch displacements for the unsteady step by summing the displacement at each sudo time step ---*/
  su2double pitch, plunge;
  pitch = config->GetAeroelastic_pitch(iMarker);
  plunge = config->GetAeroelastic_plunge(iMarker);
  
  config->SetAeroelastic_pitch(iMarker , pitch+dalpha);
  config->SetAeroelastic_plunge(iMarker , plunge+dh/b);
  
  /*--- Set the Aeroelastic solution at time n+1. This gets update every sudo time step
   and after convering the sudo time step the solution at n+1 get moved to the solution at n
   in SetDualTime_Solver method ---*/
  
  config->SetAeroelastic_np1(iMarker, x_np1);
  
}

void CSolver::Restart_OldGeometry(CGeometry *geometry, CConfig *config) {

  /*--- This function is intended for dual time simulations ---*/

  int Unst_RestartIter;
  ifstream restart_file_n;

  string filename = config->GetSolution_FileName();
  string filename_n;

  /*--- Auxiliary vector for storing the coordinates ---*/
  su2double *Coord;
  Coord = new su2double[nDim];

  /*--- Variables for reading the restart files ---*/
  string text_line;
  long iPoint_Local;
  unsigned long iPoint_Global_Local = 0, iPoint_Global = 0;
  unsigned short rbuf_NotMatching, sbuf_NotMatching;

  /*--- First, we load the restart file for time n ---*/

  /*-------------------------------------------------------------------------------------------*/

  /*--- Modify file name for an unsteady restart ---*/
  Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
  filename_n = config->GetFilename(filename, ".csv", Unst_RestartIter);

  /*--- Open the restart file, throw an error if this fails. ---*/

  restart_file_n.open(filename_n.data(), ios::in);
  if (restart_file_n.fail()) {
    SU2_MPI::Error(string("There is no flow restart file ") + filename_n, CURRENT_FUNCTION);
  }

  /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
  iPoint_Global_Local = 0; iPoint_Global = 0;

  /*--- Read all lines in the restart file ---*/
  /*--- The first line is the header ---*/

  getline (restart_file_n, text_line);

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    
    getline (restart_file_n, text_line);
    
    vector<string> point_line = PrintingToolbox::split(text_line, ',');

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      
      Coord[0] = PrintingToolbox::stod(point_line[1]);
      Coord[1] = PrintingToolbox::stod(point_line[2]);
      if (nDim == 3){
        Coord[2] = PrintingToolbox::stod(point_line[3]);          
      } 
      geometry->node[iPoint_Local]->SetCoord_n(Coord);

      iPoint_Global_Local++;
    }
  }

  /*--- Detect a wrong solution file ---*/

  rbuf_NotMatching = 0; sbuf_NotMatching = 0;

  if (iPoint_Global_Local < geometry->GetnPointDomain()) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
    SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  /*--- Close the restart file ---*/

  restart_file_n.close();

  /*-------------------------------------------------------------------------------------------*/
  /*-------------------------------------------------------------------------------------------*/

  /*--- Now, we load the restart file for time n-1, if the simulation is 2nd Order ---*/

  if (config->GetTime_Marching() == DT_STEPPING_2ND) {

    ifstream restart_file_n1;
    string filename_n1;

    /*--- Modify file name for an unsteady restart ---*/
    Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    filename_n1 = config->GetFilename(filename, ".csv", Unst_RestartIter);

    /*--- Open the restart file, throw an error if this fails. ---*/

    restart_file_n1.open(filename_n1.data(), ios::in);
    if (restart_file_n1.fail()) {
        SU2_MPI::Error(string("There is no flow restart file ") + filename_n1, CURRENT_FUNCTION);

    }

    /*--- First, set all indices to a negative value by default, and Global n indices to 0 ---*/
    iPoint_Global_Local = 0; iPoint_Global = 0;

    /*--- Read all lines in the restart file ---*/
    /*--- The first line is the header ---*/

    getline (restart_file_n1, text_line);

    for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
      
      getline (restart_file_n1, text_line);

      vector<string> point_line = PrintingToolbox::split(text_line, ',');
      
      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {
        
        Coord[0] = PrintingToolbox::stod(point_line[1]);
        Coord[1] = PrintingToolbox::stod(point_line[2]);
        if (nDim == 3){
          Coord[2] = PrintingToolbox::stod(point_line[3]);          
        } 

        geometry->node[iPoint_Local]->SetCoord_n1(Coord);

        iPoint_Global_Local++;
      }

    }

    /*--- Detect a wrong solution file ---*/

    rbuf_NotMatching = 0; sbuf_NotMatching = 0;

    if (iPoint_Global_Local < geometry->GetnPointDomain()) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
    rbuf_NotMatching = sbuf_NotMatching;
#else
    SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
    if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
    }

    /*--- Close the restart file ---*/

    restart_file_n1.close();

  }

  /*--- It's necessary to communicate this information ---*/
  
  geometry->InitiateComms(geometry, config, COORDINATES_OLD);
  geometry->CompleteComms(geometry, config, COORDINATES_OLD);
  
  delete [] Coord;

}

void CSolver::Read_SU2_Restart_ASCII(CGeometry *geometry, CConfig *config, string val_filename) {

  ifstream restart_file;
  string text_line, Tag;
  unsigned short iVar;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  int counter = 0;
  fields.clear();

  Restart_Vars = new int[5];
  
  string error_string = "Note: ASCII restart files must be in CSV format since v7.0.\n"
                        "Check https://su2code.github.io/docs/Guide-to-v7 for more information.";

  /*--- First, check that this is not a binary restart file. ---*/

  char fname[100];
  val_filename += ".csv";
  strcpy(fname, val_filename.c_str());
  int magic_number;

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + fname, CURRENT_FUNCTION);
  }

  /*--- Attempt to read the first int, which should be our magic number. ---*/

  ret = fread(&magic_number, sizeof(int), 1, fhw);
  if (ret != 1) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  int ierr;

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("SU2 ASCII restart file ") + string(fname) + string(" not found.\n") + error_string, 
                   CURRENT_FUNCTION);
  }

  /*--- Have the master attempt to read the magic number. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  MPI_File_close(&fhw);

#endif

  /*--- Open the restart file ---*/

  restart_file.open(val_filename.data(), ios::in);

  /*--- In case there is no restart file ---*/

  if (restart_file.fail()) {
    SU2_MPI::Error(string("SU2 ASCII restart file ") + string(fname) + string(" not found.\n") + error_string, 
                   CURRENT_FUNCTION);
  }

  /*--- Identify the number of fields (and names) in the restart file ---*/

  getline (restart_file, text_line);
  
  char delimiter = ',';
  fields = PrintingToolbox::split(text_line, delimiter);
  
  if (fields.size() <= 1) {
    SU2_MPI::Error(string("Restart file does not seem to be a CSV file.\n") + error_string, CURRENT_FUNCTION);
  }
  
  for (unsigned short iField = 0; iField < fields.size(); iField++){
    PrintingToolbox::trim(fields[iField]);
  }

  /*--- Set the number of variables, one per field in the
   restart file (without including the PointID) ---*/

  Restart_Vars[1] = (int)fields.size() - 1;

  /*--- Allocate memory for the restart data. ---*/

  Restart_Data = new passivedouble[Restart_Vars[1]*geometry->GetnPointDomain()];

  /*--- Read all lines in the restart file and extract data. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    getline (restart_file, text_line);
    
    vector<string> point_line = PrintingToolbox::split(text_line, delimiter);

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      
      /*--- Store the solution (starting with node coordinates) --*/

      for (iVar = 0; iVar < Restart_Vars[1]; iVar++)
        Restart_Data[counter*Restart_Vars[1] + iVar] = SU2_TYPE::GetValue(PrintingToolbox::stod(point_line[iVar+1]));

      /*--- Increment our local point counter. ---*/

      counter++;

    }
  }

}

void CSolver::Read_SU2_Restart_Binary(CGeometry *geometry, CConfig *config, string val_filename) {

  char str_buf[CGNS_STRING_SIZE], fname[100];
  unsigned short iVar;
  val_filename += ".dat";  
  strcpy(fname, val_filename.c_str());
  int nRestart_Vars = 5, nFields;
  Restart_Vars = new int[5];
  fields.clear();

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, read the number of variables and points. ---*/

  ret = fread(Restart_Vars, sizeof(int), nRestart_Vars, fhw);
  if (ret != (unsigned long)nRestart_Vars) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  /*--- Store the number of fields to be read for clarity. ---*/

  nFields = Restart_Vars[1];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. We pad the beginning of the
   variable string vector with the Point_ID tag that wasn't written. ---*/

  fields.push_back("Point_ID");
  for (iVar = 0; iVar < nFields; iVar++) {
    ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, fhw);
    if (ret != (unsigned long)CGNS_STRING_SIZE) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }
    fields.push_back(str_buf);
  }

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  Restart_Data = new passivedouble[nFields*geometry->GetnPointDomain()];

  /*--- Read in the data for the restart at all local points. ---*/

  ret = fread(Restart_Data, sizeof(passivedouble), nFields*geometry->GetnPointDomain(), fhw);
  if (ret != (unsigned long)nFields*geometry->GetnPointDomain()) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Close the file. ---*/

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  SU2_MPI::Status status;
  MPI_Datatype etype, filetype;
  MPI_Offset disp;
  unsigned long iPoint_Global, index, iChar;
  string field_buf;

  int ierr;

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- First, read the number of variables and points (i.e., cols and rows),
   which we will need in order to read the file later. Also, read the
   variable string names here. Only the master rank reads the header. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, Restart_Vars, nRestart_Vars, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(Restart_Vars, nRestart_Vars, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Check that this is an SU2 binary file. SU2 binary files
   have the hex representation of "SU2" as the first int in the file. ---*/

  if (Restart_Vars[0] != 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                   string("SU2 reads/writes binary restart files by default.\n") +
                   string("Note that backward compatibility for ASCII restart files is\n") +
                   string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  /*--- Store the number of fields to be read for clarity. ---*/

  nFields = Restart_Vars[1];

  /*--- Read the variable names from the file. Note that we are adopting a
   fixed length of 33 for the string length to match with CGNS. This is
   needed for when we read the strings later. ---*/

  char *mpi_str_buf = new char[nFields*CGNS_STRING_SIZE];
  if (rank == MASTER_NODE) {
    disp = nRestart_Vars*sizeof(int);
    MPI_File_read_at(fhw, disp, mpi_str_buf, nFields*CGNS_STRING_SIZE,
                     MPI_CHAR, MPI_STATUS_IGNORE);
  }

  /*--- Broadcast the string names of the variables. ---*/

  SU2_MPI::Bcast(mpi_str_buf, nFields*CGNS_STRING_SIZE, MPI_CHAR,
                 MASTER_NODE, MPI_COMM_WORLD);

  /*--- Now parse the string names and load into the config class in case
   we need them for writing visualization files (SU2_SOL). ---*/

  fields.push_back("Point_ID");
  for (iVar = 0; iVar < nFields; iVar++) {
    index = iVar*CGNS_STRING_SIZE;
    field_buf.append("\"");
    for (iChar = 0; iChar < (unsigned long)CGNS_STRING_SIZE; iChar++) {
      str_buf[iChar] = mpi_str_buf[index + iChar];
    }
    field_buf.append(str_buf);
    field_buf.append("\"");
    fields.push_back(field_buf.c_str());
    field_buf.clear();
  }

  /*--- Free string buffer memory. ---*/

  delete [] mpi_str_buf;

  /*--- We're writing only su2doubles in the data portion of the file. ---*/

  etype = MPI_DOUBLE;

  /*--- We need to ignore the 4 ints describing the nVar_Restart and nPoints,
   along with the string names of the variables. ---*/

  disp = nRestart_Vars*sizeof(int) + CGNS_STRING_SIZE*nFields*sizeof(char);

  /*--- Define a derived datatype for this rank's set of non-contiguous data
   that will be placed in the restart. Here, we are collecting each one of the
   points which are distributed throughout the file in blocks of nVar_Restart data. ---*/

  int *blocklen = new int[geometry->GetnPointDomain()];
  int *displace = new int[geometry->GetnPointDomain()];
  int counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {
    if (geometry->GetGlobal_to_Local_Point(iPoint_Global) > -1) {
      blocklen[counter] = nFields;
      displace[counter] = iPoint_Global*nFields;
      counter++;
    }
  }
  MPI_Type_indexed(geometry->GetnPointDomain(), blocklen, displace, MPI_DOUBLE, &filetype);
  MPI_Type_commit(&filetype);

  /*--- Set the view for the MPI file write, i.e., describe the location in
   the file that this rank "sees" for writing its piece of the restart file. ---*/

  MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

  /*--- For now, create a temp 1D buffer to read the data from file. ---*/

  Restart_Data = new passivedouble[nFields*geometry->GetnPointDomain()];

  /*--- Collective call for all ranks to read from their view simultaneously. ---*/

  MPI_File_read_all(fhw, Restart_Data, nFields*geometry->GetnPointDomain(), MPI_DOUBLE, &status);

  /*--- All ranks close the file after writing. ---*/

  MPI_File_close(&fhw);

  /*--- Free the derived datatype and release temp memory. ---*/

  MPI_Type_free(&filetype);

  delete [] blocklen;
  delete [] displace;
  
#endif
  
}

void CSolver::Read_SU2_Restart_Metadata(CGeometry *geometry, CConfig *config, bool adjoint, string val_filename) {

  su2double AoA_ = config->GetAoA();
  su2double AoS_ = config->GetAoS();
  su2double BCThrust_ = config->GetInitial_BCThrust();
  su2double dCD_dCL_ = config->GetdCD_dCL();
  su2double dCMx_dCL_ = config->GetdCMx_dCL();
  su2double dCMy_dCL_ = config->GetdCMy_dCL();
  su2double dCMz_dCL_ = config->GetdCMz_dCL();
  string::size_type position;
  unsigned long InnerIter_ = 0;
  ifstream restart_file;
  
  /*--- Carry on with ASCII metadata reading. ---*/
  
  restart_file.open(val_filename.data(), ios::in);
  if (restart_file.fail()) {
    if (rank == MASTER_NODE) {
      cout << " Warning: There is no restart file (" << val_filename.data() << ")."<< endl;
      cout << " Computation will continue without updating metadata parameters." << endl;
    }
  } 
  else {
    
    string text_line;
    
    /*--- Space for extra info (if any) ---*/
    
    while (getline (restart_file, text_line)) {
      
      /*--- External iteration ---*/
      
      position = text_line.find ("ITER=",0);
      if (position != string::npos) {
        text_line.erase (0,9); InnerIter_ = atoi(text_line.c_str());
      }
      
      /*--- Angle of attack ---*/
      
      position = text_line.find ("AOA=",0);
      if (position != string::npos) {
        text_line.erase (0,4); AoA_ = atof(text_line.c_str());
      }
      
      /*--- Sideslip angle ---*/
      
      position = text_line.find ("SIDESLIP_ANGLE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); AoS_ = atof(text_line.c_str());
      }
      
      /*--- BCThrust angle ---*/
      
      position = text_line.find ("INITIAL_BCTHRUST=",0);
      if (position != string::npos) {
        text_line.erase (0,17); BCThrust_ = atof(text_line.c_str());
      }
             
      /*--- dCD_dCL coefficient ---*/
      
      position = text_line.find ("DCD_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,14); dCD_dCL_ = atof(text_line.c_str());
      }
      
      /*--- dCMx_dCL coefficient ---*/
      
      position = text_line.find ("DCMX_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMx_dCL_ = atof(text_line.c_str());
      }
      
      /*--- dCMy_dCL coefficient ---*/
      
      position = text_line.find ("DCMY_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMy_dCL_ = atof(text_line.c_str());
      }
      
      /*--- dCMz_dCL coefficient ---*/
      
      position = text_line.find ("DCMZ_DCL_VALUE=",0);
      if (position != string::npos) {
        text_line.erase (0,15); dCMz_dCL_ = atof(text_line.c_str());
      }     
      
    }
    
    /*--- Close the restart meta file. ---*/
    
    restart_file.close();
    
  }
  

  /*--- Load the metadata. ---*/
  
  /*--- Angle of attack ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
      cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
    }
    config->SetAoA(AoA_ + config->GetAoA_Offset());
  }

  else {
    if ((config->GetAoA() != AoA_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoA in the solution file." << endl;
  }

  /*--- Sideslip angle ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
      cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
    }
    config->SetAoS(AoS_ + config->GetAoS_Offset());
  }
  else {
    if ((config->GetAoS() != AoS_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoS in the solution file." << endl;
  }

  /*--- BCThrust ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
      config->SetInitial_BCThrust(BCThrust_);
    }
  else {
    if ((config->GetInitial_BCThrust() != BCThrust_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
  }


  if (config->GetDiscard_InFiles() == false) {

    if ((config->GetdCD_dCL() != dCD_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCD/dCL provided in the direct solution file: " << dCD_dCL_ << "." << endl;
    config->SetdCD_dCL(dCD_dCL_);

    if ((config->GetdCMx_dCL() != dCMx_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMx/dCL provided in the direct solution file: " << dCMx_dCL_ << "." << endl;
    config->SetdCMx_dCL(dCMx_dCL_);

    if ((config->GetdCMy_dCL() != dCMy_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMy/dCL provided in the direct solution file: " << dCMy_dCL_ << "." << endl;
    config->SetdCMy_dCL(dCMy_dCL_);

    if ((config->GetdCMz_dCL() != dCMz_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMz/dCL provided in the direct solution file: " << dCMz_dCL_ << "." << endl;
    config->SetdCMz_dCL(dCMz_dCL_);

  }
	
  else {

    if ((config->GetdCD_dCL() != dCD_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCD/dCL in the direct solution file." << endl;
    
    if ((config->GetdCMx_dCL() != dCMx_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMx/dCL in the direct solution file." << endl;
    
    if ((config->GetdCMy_dCL() != dCMy_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMy/dCL in the direct solution file." << endl;
    
    if ((config->GetdCMz_dCL() != dCMz_dCL_) &&  (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMz/dCL in the direct solution file." << endl;

  }
  
  /*--- External iteration ---*/

  if ((config->GetDiscard_InFiles() == false) && (!adjoint || (adjoint && config->GetRestart())))
    config->SetExtIter_OffSet(InnerIter_);

}

void CSolver::Read_InletFile_ASCII(CGeometry *geometry, CConfig *config, string val_filename) {

  ifstream inlet_file;
  string text_line;
  unsigned long iVar, iMarker, iChar, iRow;
  int counter = 0;
  string::size_type position;

  /*--- Open the inlet profile file (we have already error checked) ---*/

  inlet_file.open(val_filename.data(), ios::in);

  /*--- Identify the markers and data set in the inlet profile file ---*/

  while (getline (inlet_file, text_line)) {

    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {
      text_line.erase (0,6); nMarker_InletFile = atoi(text_line.c_str());

      nRow_InletFile    = new unsigned long[nMarker_InletFile];
      nRowCum_InletFile = new unsigned long[nMarker_InletFile+1];
      nCol_InletFile    = new unsigned long[nMarker_InletFile];

      for (iMarker = 0 ; iMarker < nMarker_InletFile; iMarker++) {

        getline (inlet_file, text_line);
        text_line.erase (0,11);
        for (iChar = 0; iChar < 20; iChar++) {
          position = text_line.find( " ", 0 );  if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\r", 0 ); if (position != string::npos) text_line.erase (position,1);
          position = text_line.find( "\n", 0 ); if (position != string::npos) text_line.erase (position,1);
        }
        Marker_Tags_InletFile.push_back(text_line.c_str());

        getline (inlet_file, text_line);
        text_line.erase (0,5); nRow_InletFile[iMarker] = atoi(text_line.c_str());

        getline (inlet_file, text_line);
        text_line.erase (0,5); nCol_InletFile[iMarker] = atoi(text_line.c_str());

        /*--- Skip the data. This is read in the next loop. ---*/

        for (iRow = 0; iRow < nRow_InletFile[iMarker]; iRow++) getline (inlet_file, text_line);

      }
    } else {
      SU2_MPI::Error("While opening inlet file, no \"NMARK=\" specification was found", CURRENT_FUNCTION);
    }
  }

  inlet_file.close();

  /*--- Compute array bounds and offsets. Allocate data structure. ---*/

  maxCol_InletFile = 0; nRowCum_InletFile[0] = 0;
  for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {
    if (nCol_InletFile[iMarker] > maxCol_InletFile)
      maxCol_InletFile = nCol_InletFile[iMarker];

    /*--- Put nRow into cumulative storage format. ---*/

    nRowCum_InletFile[iMarker+1] = nRowCum_InletFile[iMarker] + nRow_InletFile[iMarker];
    
  }

  Inlet_Data = new passivedouble[nRowCum_InletFile[nMarker_InletFile]*maxCol_InletFile];

  for (unsigned long iPoint = 0; iPoint < nRowCum_InletFile[nMarker_InletFile]*maxCol_InletFile; iPoint++)
    Inlet_Data[iPoint] = 0.0;

  /*--- Read all lines in the inlet profile file and extract data. ---*/

  inlet_file.open(val_filename.data(), ios::in);

  counter = 0;
  while (getline (inlet_file, text_line)) {

    position = text_line.find ("NMARK=",0);
    if (position != string::npos) {

      for (iMarker = 0; iMarker < nMarker_InletFile; iMarker++) {

        /*--- Skip the tag, nRow, and nCol lines. ---*/

        getline (inlet_file, text_line);
        getline (inlet_file, text_line);
        getline (inlet_file, text_line);

        /*--- Now read the data for each row and store. ---*/

        for (iRow = 0; iRow < nRow_InletFile[iMarker]; iRow++) {

          getline (inlet_file, text_line);

          istringstream point_line(text_line);

          /*--- Store the values (starting with node coordinates) --*/

          for (iVar = 0; iVar < nCol_InletFile[iMarker]; iVar++)
            point_line >> Inlet_Data[counter*maxCol_InletFile + iVar];

          /*--- Increment our local row counter. ---*/

          counter++;

        }
      }
    }
  }
  
  inlet_file.close();
  
}

void CSolver::LoadInletProfile(CGeometry **geometry,
                               CSolver ***solver,
                               CConfig *config,
                               int val_iter,
                               unsigned short val_kind_solver,
                               unsigned short val_kind_marker) {

  /*-- First, set the solver and marker kind for the particular problem at
   hand. Note that, in the future, these routines can be used for any solver
   and potentially any marker type (beyond inlets). ---*/

  unsigned short KIND_SOLVER = val_kind_solver;
  unsigned short KIND_MARKER = val_kind_marker;

  /*--- Local variables ---*/

  unsigned short iDim, iVar, iMesh, iMarker, jMarker;
  unsigned long iPoint, iVertex, index, iChildren, Point_Fine, iRow;
  su2double Area_Children, Area_Parent, *Coord, dist, min_dist;
  bool dual_time = ((config->GetTime_Marching() == DT_STEPPING_1ST) ||
                    (config->GetTime_Marching() == DT_STEPPING_2ND));
  bool time_stepping = config->GetTime_Marching() == TIME_STEPPING;

  string UnstExt, text_line;
  ifstream restart_file;

  unsigned short iZone = config->GetiZone();
  unsigned short nZone = config->GetnZone();

  string Marker_Tag;
  string profile_filename = config->GetInlet_FileName();
  ifstream inlet_file;

  su2double *Inlet_Values = NULL;
  su2double *Inlet_Fine   = NULL;
  su2double *Normal       = new su2double[nDim];

  unsigned long Marker_Counter = 0;

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    profile_filename = config->GetMultizone_FileName(profile_filename, iZone, ".dat");

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    profile_filename = config->GetUnsteady_FileName(profile_filename, val_iter, ".dat");

  /*--- Open the file and check for problems. If a file can not be found,
   then a warning will be printed, but the calculation will continue
   using the uniform inlet values. A template inlet file will be written
   at a later point using COutput. ---*/

  inlet_file.open(profile_filename.data(), ios::in);

  if (!inlet_file.fail()) {

    /*--- Close the file and start the loading. ---*/

    inlet_file.close();

    /*--- Read the profile data from an ASCII file. ---*/

    Read_InletFile_ASCII(geometry[MESH_0], config, profile_filename);

    /*--- Load data from the restart into correct containers. ---*/

    Marker_Counter = 0;

    Inlet_Values = new su2double[maxCol_InletFile];
    Inlet_Fine   = new su2double[maxCol_InletFile];

    unsigned short global_failure = 0, local_failure = 0;
    ostringstream error_msg;

    const su2double tolerance = config->GetInlet_Profile_Matching_Tolerance();

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == KIND_MARKER) {

        /*--- Get tag in order to identify the correct inlet data. ---*/

        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

        for (jMarker = 0; jMarker < nMarker_InletFile; jMarker++) {

          /*--- If we have found the matching marker string, continue. ---*/

          if (Marker_Tags_InletFile[jMarker] == Marker_Tag) {

            /*--- Increment our counter for marker matches. ---*/

            Marker_Counter++;

            /*--- Loop through the nodes on this marker. ---*/

            for (iVertex = 0; iVertex < geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

              iPoint   = geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
              Coord    = geometry[MESH_0]->node[iPoint]->GetCoord();
              min_dist = 1e16;

              /*--- Find the distance to the closest point in our inlet profile data. ---*/

              for (iRow = nRowCum_InletFile[jMarker]; iRow < nRowCum_InletFile[jMarker+1]; iRow++) {

                /*--- Get the coords for this data point. ---*/

                index = iRow*maxCol_InletFile;

                dist = 0.0;
                for (unsigned short iDim = 0; iDim < nDim; iDim++)
                  dist += pow(Inlet_Data[index+iDim] - Coord[iDim], 2);
                dist = sqrt(dist);

                /*--- Check is this is the closest point and store data if so. ---*/

                if (dist < min_dist) {
                  min_dist = dist;
                  for (iVar = 0; iVar < maxCol_InletFile; iVar++)
                    Inlet_Values[iVar] = Inlet_Data[index+iVar];
                }

              }

              /*--- If the diff is less than the tolerance, match the two.
               We could modify this to simply use the nearest neighbor, or
               eventually add something more elaborate here for interpolation. ---*/

              if (min_dist < tolerance) {

                solver[MESH_0][KIND_SOLVER]->SetInletAtVertex(Inlet_Values, iMarker, iVertex);

              } else {

                unsigned long GlobalIndex = geometry[MESH_0]->node[iPoint]->GetGlobalIndex();
                cout << "WARNING: Did not find a match between the points in the inlet file" << endl;
                cout << "and point " << GlobalIndex;
                cout << std::scientific;
                cout << " at location: [" << Coord[0] << ", " << Coord[1];
                if (nDim ==3) error_msg << ", " << Coord[2];
                cout << "]" << endl;
                cout << "Distance to closest point: " << min_dist << endl;
                cout << "Current tolerance:         " << tolerance << endl;
                cout << endl;
                cout << "You can widen the tolerance for point matching by changing the value" << endl;
                cout << "of the option INLET_MATCHING_TOLERANCE in your *.cfg file." << endl;
                local_failure++;
                break;

              }
            }
          }
        }
      }

      if (local_failure > 0) break;
    }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&local_failure, &global_failure, 1, MPI_UNSIGNED_SHORT,
                       MPI_SUM, MPI_COMM_WORLD);
#else
    global_failure = local_failure;
#endif

    if (global_failure > 0) {
      SU2_MPI::Error(string("Prescribed inlet data does not match markers within tolerance."), CURRENT_FUNCTION);
    }

    /*--- Copy the inlet data down to the coarse levels if multigrid is active.
     Here, we use a face area-averaging to restrict the values. ---*/

    for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == KIND_MARKER) {

          Marker_Tag = config->GetMarker_All_TagBound(iMarker);
          
          /*--- Loop through the nodes on this marker. ---*/

          for (iVertex = 0; iVertex < geometry[iMesh]->nVertex[iMarker]; iVertex++) {

            /*--- Get the coarse mesh point and compute the boundary area. ---*/

            iPoint = geometry[iMesh]->vertex[iMarker][iVertex]->GetNode();
            geometry[iMesh]->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area_Parent = 0.0;
            for (iDim = 0; iDim < nDim; iDim++) Area_Parent += Normal[iDim]*Normal[iDim];
            Area_Parent = sqrt(Area_Parent);

            /*--- Reset the values for the coarse point. ---*/

            for (iVar = 0; iVar < maxCol_InletFile; iVar++) Inlet_Values[iVar] = 0.0;

            /*-- Loop through the children and extract the inlet values
             from those nodes that lie on the boundary as well as their
             boundary area. We build a face area-averaged value for the
             coarse point values from the fine grid points. Note that
             children from the interior volume will not be included in
             the averaging. ---*/

            for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
              Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
              for (iVar = 0; iVar < maxCol_InletFile; iVar++) Inlet_Fine[iVar] = 0.0;
              Area_Children = solver[iMesh-1][KIND_SOLVER]->GetInletAtVertex(Inlet_Fine, Point_Fine, KIND_MARKER, Marker_Tag, geometry[iMesh-1], config);
              for (iVar = 0; iVar < maxCol_InletFile; iVar++) {
                Inlet_Values[iVar] += Inlet_Fine[iVar]*Area_Children/Area_Parent;
              }
            }

            /*--- Set the boundary area-averaged inlet values for the coarse point. ---*/

            solver[iMesh][KIND_SOLVER]->SetInletAtVertex(Inlet_Values, iMarker, iVertex);

          }
        }
      }
    }

    /*--- Delete the class memory that is used to load the inlets. ---*/

    Marker_Tags_InletFile.clear();

    if (nRowCum_InletFile != NULL) {delete [] nRowCum_InletFile; nRowCum_InletFile = NULL;}
    if (nRow_InletFile    != NULL) {delete [] nRow_InletFile;    nRow_InletFile    = NULL;}
    if (nCol_InletFile    != NULL) {delete [] nCol_InletFile;    nCol_InletFile    = NULL;}
    if (Inlet_Data        != NULL) {delete [] Inlet_Data;        Inlet_Data        = NULL;}

  } else {

    if (rank == MASTER_NODE) {
      cout << endl;
      cout << "WARNING: Could not find the input file for the inlet profile." << endl;
      cout << "Looked for: " << profile_filename << "." << endl;
      cout << "A template inlet profile file will be written, and the " << endl;
      cout << "calculation will continue with uniform inlets." << endl << endl;
    }

    /*--- Set the bit to write a template inlet profile file. ---*/

    config->SetWrt_InletFile(true);

    /*--- Set the mean flow inlets to uniform. ---*/

    for (iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          solver[iMesh][KIND_SOLVER]->SetUniformInlet(config, iMarker);
      }
    }

  }

  /*--- Deallocated local data. ---*/

  if (Inlet_Values != NULL) delete [] Inlet_Values;
  if (Inlet_Fine   != NULL) delete [] Inlet_Fine;
  delete [] Normal;
  
}

void CSolver::ComputeVertexTractions(CGeometry *geometry, CConfig *config){

  /*--- Compute the constant factor to dimensionalize pressure and shear stress. ---*/
  su2double *Velocity_ND, *Velocity_Real;
  su2double Density_ND,  Density_Real, Velocity2_Real, Velocity2_ND;
  su2double factor;

  unsigned short iDim, jDim;

  // Check whether the problem is viscous
  bool viscous_flow = ((config->GetKind_Solver() == NAVIER_STOKES) ||
                       (config->GetKind_Solver() == INC_NAVIER_STOKES) ||
                       (config->GetKind_Solver() == RANS) ||
                       (config->GetKind_Solver() == INC_RANS) ||
                       (config->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES) ||
                       (config->GetKind_Solver() == DISC_ADJ_INC_NAVIER_STOKES) ||
                       (config->GetKind_Solver() == DISC_ADJ_INC_RANS) ||
                       (config->GetKind_Solver() == DISC_ADJ_RANS));

  // Parameters for the calculations
  su2double Pn = 0.0, div_vel = 0.0;
  su2double Viscosity = 0.0;
  su2double Tau[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double Grad_Vel[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0, 1.0, 0.0},{0.0, 0.0, 1.0}};
  su2double auxForce[3] = {1.0, 0.0, 0.0};

  unsigned short iMarker;
  unsigned long iVertex, iPoint;
  su2double const *iNormal;

  su2double Pressure_Inf = config->GetPressure_FreeStreamND();

  Velocity_Real = config->GetVelocity_FreeStream();
  Density_Real  = config->GetDensity_FreeStream();

  Velocity_ND = config->GetVelocity_FreeStreamND();
  Density_ND  = config->GetDensity_FreeStreamND();

  Velocity2_Real = 0.0;
  Velocity2_ND   = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Velocity2_Real += Velocity_Real[iDim]*Velocity_Real[iDim];
    Velocity2_ND   += Velocity_ND[iDim]*Velocity_ND[iDim];
  }

  factor = Density_Real * Velocity2_Real / ( Density_ND * Velocity2_ND );

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as an interface marker ---*/
    if (config->GetMarker_All_Fluid_Load(iMarker) == YES) {

      // Loop over the vertices
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        // Recover the point index
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        // Get the normal at the vertex: this normal goes inside the fluid domain.
        iNormal = geometry->vertex[iMarker][iVertex]->GetNormal();

        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        if (geometry->node[iPoint]->GetDomain()) {

          // Retrieve the values of pressure
          Pn = base_nodes->GetPressure(iPoint);

          // Calculate tn in the fluid nodes for the inviscid term --> Units of force (non-dimensional).
          for (iDim = 0; iDim < nDim; iDim++)
            auxForce[iDim] = -(Pn-Pressure_Inf)*iNormal[iDim];

          // Calculate tn in the fluid nodes for the viscous term
          if (viscous_flow) {

            Viscosity = base_nodes->GetLaminarViscosity(iPoint);

            for (iDim = 0; iDim < nDim; iDim++) {
              for (jDim = 0 ; jDim < nDim; jDim++) {
                Grad_Vel[iDim][jDim] = base_nodes->GetGradient_Primitive(iPoint, iDim+1, jDim);
              }
            }

            // Divergence of the velocity
            div_vel = 0.0; for (iDim = 0; iDim < nDim; iDim++) div_vel += Grad_Vel[iDim][iDim];

            for (iDim = 0; iDim < nDim; iDim++) {
              for (jDim = 0 ; jDim < nDim; jDim++) {

                // Viscous stress
                Tau[iDim][jDim] = Viscosity*(Grad_Vel[jDim][iDim] + Grad_Vel[iDim][jDim])
                                 - TWO3*Viscosity*div_vel*delta[iDim][jDim];

                // Viscous component in the tn vector --> Units of force (non-dimensional).
                auxForce[iDim] += Tau[iDim][jDim]*iNormal[jDim];
              }
            }
          }

          // Redimensionalize the forces
          for (iDim = 0; iDim < nDim; iDim++) {
            VertexTraction[iMarker][iVertex][iDim] = factor * auxForce[iDim];
          }
        }
        else{
          for (iDim = 0; iDim < nDim; iDim++) {
            VertexTraction[iMarker][iVertex][iDim] = 0.0;
          }
        }
      }
    }
  }

}

void CSolver::RegisterVertexTractions(CGeometry *geometry, CConfig *config){

  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;

  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as an interface marker ---*/
    if (config->GetMarker_All_Fluid_Load(iMarker) == YES) {

      /*--- Loop over the vertices ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Recover the point index ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Register the vertex traction as output ---*/
          for (iDim = 0; iDim < nDim; iDim++) {
            AD::RegisterOutput(VertexTraction[iMarker][iVertex][iDim]);
          }

        }
      }
    }
  }

}

void CSolver::SetVertexTractionsAdjoint(CGeometry *geometry, CConfig *config){

  unsigned short iMarker, iDim;
  unsigned long iVertex, iPoint;

  /*--- Loop over all the markers ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- If this is defined as an interface marker ---*/
    if (config->GetMarker_All_Fluid_Load(iMarker) == YES) {

      /*--- Loop over the vertices ---*/
      for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

        /*--- Recover the point index ---*/
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

        /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
        if (geometry->node[iPoint]->GetDomain()) {

          /*--- Set the adjoint of the vertex traction from the value received ---*/
          for (iDim = 0; iDim < nDim; iDim++) {

            SU2_TYPE::SetDerivative(VertexTraction[iMarker][iVertex][iDim],
                                    SU2_TYPE::GetValue(VertexTractionAdjoint[iMarker][iVertex][iDim]));
          }

        }
      }
    }
  }

}


void CSolver::SetVerificationSolution(unsigned short nDim,
                                      unsigned short nVar,
                                      CConfig        *config) {

  /*--- Determine the verification solution to be set and
        allocate memory for the corresponding class. ---*/
  switch( config->GetVerification_Solution() ) {

    case NO_VERIFICATION_SOLUTION:
      VerificationSolution = NULL; break;
    case INVISCID_VORTEX:
      VerificationSolution = new CInviscidVortexSolution(nDim, nVar, MGLevel, config); break;
    case RINGLEB:
      VerificationSolution = new CRinglebSolution(nDim, nVar, MGLevel, config); break;
    case NS_UNIT_QUAD:
      VerificationSolution = new CNSUnitQuadSolution(nDim, nVar, MGLevel, config); break;
    case TAYLOR_GREEN_VORTEX:
      VerificationSolution = new CTGVSolution(nDim, nVar, MGLevel, config); break;
    case INC_TAYLOR_GREEN_VORTEX:
      VerificationSolution = new CIncTGVSolution(nDim, nVar, MGLevel, config); break;
    case MMS_NS_UNIT_QUAD:
      VerificationSolution = new CMMSNSUnitQuadSolution(nDim, nVar, MGLevel, config); break;
    case MMS_NS_UNIT_QUAD_WALL_BC:
      VerificationSolution = new CMMSNSUnitQuadSolutionWallBC(nDim, nVar, MGLevel, config); break;
    case MMS_NS_TWO_HALF_CIRCLES:
      VerificationSolution = new CMMSNSTwoHalfCirclesSolution(nDim, nVar, MGLevel, config); break;
    case MMS_NS_TWO_HALF_SPHERES:
      VerificationSolution = new CMMSNSTwoHalfSpheresSolution(nDim, nVar, MGLevel, config); break;
    case MMS_INC_EULER:
      VerificationSolution = new CMMSIncEulerSolution(nDim, nVar, MGLevel, config); break;
    case MMS_INC_NS:
      VerificationSolution = new CMMSIncNSSolution(nDim, nVar, MGLevel, config); break;
    case USER_DEFINED_SOLUTION:
      VerificationSolution = new CUserDefinedSolution(nDim, nVar, MGLevel, config); break;
  }
}

void CSolver::ComputeResidual_Multizone(CGeometry *geometry, CConfig *config){

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
    SetRes_BGS(iVar,0.0);
    SetRes_Max_BGS(iVar,0.0,0);
  }

  /*--- Set the residuals ---*/
  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      residual = base_nodes->GetSolution(iPoint,iVar) - base_nodes->Get_BGSSolution_k(iPoint,iVar);
      AddRes_BGS(iVar,residual*residual);
      AddRes_Max_BGS(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
    }
  }

  SetResidual_BGS(geometry, config);

}


void CSolver::UpdateSolution_BGS(CGeometry *geometry, CConfig *config){

  /*--- To nPoint: The solution must be communicated beforehand ---*/
  base_nodes->Set_BGSSolution_k();
}

CBaselineSolver::CBaselineSolver(void) : CSolver() { }

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config) {

  nPoint = geometry->GetnPoint();

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Routines to access the number of variables and string names. ---*/

  SetOutputVariables(geometry, config);
  
  /*--- Initialize a zero solution and instantiate the CVariable class. ---*/

  Solution = new su2double[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) Solution[iVar] = 0.0;

  nodes = new CBaselineVariable(nPoint, nVar, config);
  SetBaseClassPointerToNodes();

  dynamic_grid = config->GetDynamic_Grid();
}

CBaselineSolver::CBaselineSolver(CGeometry *geometry, CConfig *config, unsigned short val_nvar, vector<string> field_names) {

  /*--- Define geometry constants in the solver structure ---*/

  nPoint = geometry->GetnPoint();
  nDim = geometry->GetnDim();
  nVar = val_nvar;
  fields = field_names;

  /*--- Allocate the node variables ---*/

  nodes = new CBaselineVariable(nPoint, nVar, config);
  SetBaseClassPointerToNodes();

  dynamic_grid = config->GetDynamic_Grid();

}

void CBaselineSolver::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  /*--- Open the restart file and extract the nVar and field names. ---*/

  string Tag, text_line;

  ifstream restart_file;
  string filename;

  /*--- Retrieve filename from config ---*/

  if (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint()) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }


  /*--- Read only the number of variables in the restart file. ---*/

  if (config->GetRead_Binary_Restart()) {
    
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    filename = config->GetFilename(filename, ".dat", config->GetTimeIter());

    char fname[100];
    strcpy(fname, filename.c_str());
    int nVar_Buf = 5;
    int var_buf[5];

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }
    
    /*--- First, read the number of variables and points. ---*/

    ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }
    
    /*--- Close the file. ---*/

    fclose(fhw);
    
    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];
#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;
    MPI_Offset disp;
    unsigned short iVar;
    unsigned long index, iChar;
    string field_buf;
    char str_buf[CGNS_STRING_SIZE];
    
    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

    if (rank == MASTER_NODE) {
      MPI_File_read(fhw, var_buf, nVar_Buf, MPI_INT, MPI_STATUS_IGNORE);
    }

    /*--- Broadcast the number of variables to all procs and store more clearly. ---*/

    SU2_MPI::Bcast(var_buf, nVar_Buf, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }



    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];
    
    /*--- Read the variable names from the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/
  
    char *mpi_str_buf = new char[nVar*CGNS_STRING_SIZE];
    if (rank == MASTER_NODE) {
      disp = nVar_Buf*sizeof(int);
      MPI_File_read_at(fhw, disp, mpi_str_buf, nVar*CGNS_STRING_SIZE,
                       MPI_CHAR, MPI_STATUS_IGNORE);
    }
    
    /*--- Broadcast the string names of the variables. ---*/
  
    SU2_MPI::Bcast(mpi_str_buf, nVar*CGNS_STRING_SIZE, MPI_CHAR,
                   MASTER_NODE, MPI_COMM_WORLD);
    
    fields.push_back("Point_ID");
    
    for (iVar = 0; iVar < nVar; iVar++) {
      index = iVar*CGNS_STRING_SIZE;
      field_buf.append("\"");
      for (iChar = 0; iChar < (unsigned long)CGNS_STRING_SIZE; iChar++) {
        str_buf[iChar] = mpi_str_buf[index + iChar];
      }
      field_buf.append(str_buf);
      field_buf.append("\"");
      fields.push_back(field_buf.c_str());
      field_buf.clear();
    }
    
    /*--- All ranks close the file after writing. ---*/
    
    MPI_File_close(&fhw);
    
#endif
  } else {
    
    /*--- Multizone problems require the number of the zone to be appended. ---*/
    
    filename = config->GetFilename(filename, ".csv", config->GetTimeIter());

    /*--- First, check that this is not a binary restart file. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }
    
    MPI_File_close(&fhw);
    
#endif

    /*--- Open the restart file ---*/

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (restart_file.fail()) {
      SU2_MPI::Error(string("SU2 solution file ") + filename + string(" not found"), CURRENT_FUNCTION);
    }
    
    /*--- Identify the number of fields (and names) in the restart file ---*/

    getline (restart_file, text_line);

    fields = PrintingToolbox::split(text_line, ',');
    
    for (unsigned short iField = 0; iField < fields.size(); iField++){
      PrintingToolbox::trim(fields[iField]);
    }
    
    /*--- Close the file (the solution date is read later). ---*/
    
    restart_file.close();

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = fields.size() - 1;

  }

}

void CBaselineSolver::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/

  string filename;
  unsigned long index;
  unsigned short iDim, iVar;
  bool adjoint = ( config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint() ); 
  unsigned short iInst = config->GetiInst();
  bool steady_restart = config->GetSteadyRestart();
  unsigned short turb_model = config->GetKind_Turb_Model();

  su2double *Coord = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++)
    Coord[iDim] = 0.0;

  /*--- Skip coordinates ---*/

  unsigned short skipVars = geometry[iInst]->GetnDim();

  /*--- Retrieve filename from config ---*/

  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }

  filename = config->GetFilename(filename, "", val_iter);

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[iInst], config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[iInst], config, filename);
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/

  for (iPoint_Global = 0; iPoint_Global < geometry[iInst]->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[iInst]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      
      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);
     
      /*--- For dynamic meshes, read in and store the
       grid coordinates and grid velocities for each node. ---*/
      
      if (dynamic_grid && val_update_geo) {

        /*--- First, remove any variables for the turbulence model that
         appear in the restart file before the grid velocities. ---*/

        if (turb_model == SA || turb_model == SA_NEG) {
          index++;
        } else if (turb_model == SST) {
          index+=2;
        }
        
        /*--- Read in the next 2 or 3 variables which are the grid velocities ---*/
        /*--- If we are restarting the solution from a previously computed static calculation (no grid movement) ---*/
        /*--- the grid velocities are set to 0. This is useful for FSI computations ---*/
        
        su2double GridVel[3] = {0.0,0.0,0.0};
        if (!steady_restart) {

          /*--- Rewind the index to retrieve the Coords. ---*/
          index = counter*Restart_Vars[1];
          for (iDim = 0; iDim < nDim; iDim++) { Coord[iDim] = Restart_Data[index+iDim]; }

          /*--- Move the index forward to get the grid velocities. ---*/
          index = counter*Restart_Vars[1] + skipVars + nVar;
          for (iDim = 0; iDim < nDim; iDim++) { GridVel[iDim] = Restart_Data[index+iDim]; }
        }

        for (iDim = 0; iDim < nDim; iDim++) {
          geometry[iInst]->node[iPoint_Local]->SetCoord(iDim, Coord[iDim]);
          geometry[iInst]->node[iPoint_Local]->SetGridVel(iDim, GridVel[iDim]);
        }
      }

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }
    
  }

  /*--- MPI solution ---*/
  
  InitiateComms(geometry[iInst], config, SOLUTION);
  CompleteComms(geometry[iInst], config, SOLUTION);

  /*--- Update the geometry for flows on dynamic meshes ---*/
  
  if (dynamic_grid && val_update_geo) {
    
    /*--- Communicate the new coordinates and grid velocities at the halos ---*/
    
    geometry[iInst]->InitiateComms(geometry[iInst], config, COORDINATES);
    geometry[iInst]->CompleteComms(geometry[iInst], config, COORDINATES);
        
    geometry[iInst]->InitiateComms(geometry[iInst], config, GRID_VELOCITY);
    geometry[iInst]->CompleteComms(geometry[iInst], config, GRID_VELOCITY);

  }
  
  delete [] Coord;

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

void CBaselineSolver::LoadRestart_FSI(CGeometry *geometry, CConfig *config, int val_iter) {

  /*--- Restart the solution from file information ---*/
  string filename;
  unsigned long index;
  unsigned short iVar;
  bool adjoint = (config->GetContinuous_Adjoint() || config->GetDiscrete_Adjoint());

  /*--- Retrieve filename from config ---*/
  if (adjoint) {
    filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(filename);
  } else {
    filename = config->GetSolution_FileName();
  }

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  filename = config->GetFilename(filename, "", val_iter);

  /*--- Output the file name to the console. ---*/

  if (rank == MASTER_NODE)
    cout << "Reading and storing the solution from " << filename
    << "." << endl;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry, config, filename);
  } else {
    Read_SU2_Restart_ASCII(geometry, config, filename);
  }

  unsigned short nVar_Local = Restart_Vars[1];
  su2double *Solution_Local = new su2double[nVar_Local];

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;

  /*--- Load data from the restart into correct containers. ---*/
  
  for (iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar_Local; iVar++) Solution[iVar] = Restart_Data[index+iVar];
      nodes->SetSolution(iPoint_Local,Solution);

      /*--- Increment the overall counter for how many points have been loaded. ---*/

      counter++;

    }

  }

  delete [] Solution_Local;

}

CBaselineSolver::~CBaselineSolver(void) {
  if (nodes != nullptr) delete nodes;
}

CBaselineSolver_FEM::CBaselineSolver_FEM(void) : CSolver() { }

CBaselineSolver_FEM::CBaselineSolver_FEM(CGeometry *geometry, CConfig *config) {

  /*--- Define geometry constants in the solver structure ---*/

  nDim = geometry->GetnDim();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. If necessary, it is
   possible to increase nMatchingFacesWithHaloElem a bit, such that
   the computation of the external faces may be more efficient when
   using multiple threads. ---*/

  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  nVolElemTot   = DGGeometry->GetNVolElemTot();
  nVolElemOwned = DGGeometry->GetNVolElemOwned();
  volElem       = DGGeometry->GetVolElem();

  /*--- Routines to access the number of variables and string names. ---*/

  SetOutputVariables(geometry, config);

  /*--- Determine the total number of DOFs stored on this rank and allocate the memory
   to store the conservative variables. ---*/
  nDOFsLocOwned = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) nDOFsLocOwned += volElem[i].nDOFsSol;

  nDOFsLocTot = nDOFsLocOwned;
  for(unsigned long i=nVolElemOwned; i<nVolElemTot; ++i) nDOFsLocTot += volElem[i].nDOFsSol;

  VecSolDOFs.resize(nVar*nDOFsLocTot);

  /*--- Determine the global number of DOFs. ---*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nDOFsLocOwned, &nDOFsGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  nDOFsGlobal = nDOFsLocOwned;
#endif

  /*--- Store the number of DOFs in the geometry class in case of restart. ---*/
  geometry->SetnPointDomain(nDOFsLocOwned);
  geometry->SetGlobal_nPointDomain(nDOFsGlobal);

  /*--- Initialize the solution to zero. ---*/

  unsigned long ii = 0;
  for(unsigned long i=0; i<nDOFsLocTot; ++i) {
    for(unsigned short j=0; j<nVar; ++j, ++ii) {
      VecSolDOFs[ii] = 0.0;
    }
  }

}

void CBaselineSolver_FEM::SetOutputVariables(CGeometry *geometry, CConfig *config) {

  /*--- Open the restart file and extract the nVar and field names. ---*/

  string Tag, text_line, AdjExt, UnstExt;
  unsigned long TimeIter = config->GetTimeIter();

  ifstream restart_file;
  string filename;

  /*--- Retrieve filename from config ---*/

  filename = config->GetSolution_FileName();

  /*--- Unsteady problems require an iteration number to be appended. ---*/

  if (config->GetTime_Domain()) {
    filename = config->GetUnsteady_FileName(filename, SU2_TYPE::Int(TimeIter), ".dat");
  }

  /*--- Read only the number of variables in the restart file. ---*/

  if (config->GetRead_Binary_Restart()) {

    int nVar_Buf = 5;
    int var_buf[5];

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(filename.c_str(),"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- First, read the number of variables and points. ---*/

    ret = fread(var_buf, sizeof(int), nVar_Buf, fhw);
    if (ret != (unsigned long)nVar_Buf) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);

    /*--- Close the file. ---*/

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, var_buf, nVar_Buf, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store more clearly. ---*/

    SU2_MPI::Bcast(var_buf, nVar_Buf, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (var_buf[0] != 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);

    /*--- All ranks close the file after writing. ---*/

    MPI_File_close(&fhw);

#endif

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = var_buf[1];

  } else {

    /*--- First, check that this is not a binary restart file. ---*/

    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(filename.c_str(), "rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr)
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + filename,
                     CURRENT_FUNCTION);

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532)
      SU2_MPI::Error(string("File ") + filename + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);

    MPI_File_close(&fhw);
    
#endif

    /*--- Open the restart file ---*/

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no restart file ---*/

    if (restart_file.fail())
      SU2_MPI::Error(string("SU2 solution file ") + filename + string(" not found"), CURRENT_FUNCTION);

    /*--- Identify the number of fields (and names) in the restart file ---*/

    getline (restart_file, text_line);

    stringstream ss(text_line);
    while (ss >> Tag) {
      config->fields.push_back(Tag);
      if (ss.peek() == ',') ss.ignore();
    }

    /*--- Close the file (the solution date is read later). ---*/

    restart_file.close();

    /*--- Set the number of variables, one per field in the
     restart file (without including the PointID) ---*/

    nVar = config->fields.size() - 1;

    /*--- Clear the fields vector since we'll read it again. ---*/

    config->fields.clear();

  }

}

void CBaselineSolver_FEM::LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo) {

  /*--- Restart the solution from file information ---*/
  unsigned short iVar;
  unsigned long index;

  string restart_filename = config->GetSolution_FileName();

  if (config->GetTime_Domain()) {
    restart_filename = config->GetUnsteady_FileName(restart_filename, SU2_TYPE::Int(val_iter), "");
  }

  int counter = 0;
  long iPoint_Local = 0; unsigned long iPoint_Global = 0;
  unsigned short rbuf_NotMatching = 0;
  unsigned long nDOF_Read = 0;

  /*--- Read the restart data from either an ASCII or binary SU2 file. ---*/

  if (config->GetRead_Binary_Restart()) {
    Read_SU2_Restart_Binary(geometry[MESH_0], config, restart_filename);
  } else {
    Read_SU2_Restart_ASCII(geometry[MESH_0], config, restart_filename);
  }

  /*--- Load data from the restart into correct containers. ---*/

  counter = 0;
  for (iPoint_Global = 0; iPoint_Global < geometry[MESH_0]->GetGlobal_nPointDomain(); iPoint_Global++) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = geometry[MESH_0]->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      index = counter*Restart_Vars[1];
      for (iVar = 0; iVar < nVar; iVar++) {
        VecSolDOFs[nVar*iPoint_Local+iVar] = Restart_Data[index+iVar];
      }
      /*--- Update the local counter nDOF_Read. ---*/
      ++nDOF_Read;

      /*--- Increment the overall counter for how many points have been loaded. ---*/
      counter++;
    }

  }

  /*--- Detect a wrong solution file ---*/
  if(nDOF_Read < nDOFsLocOwned) rbuf_NotMatching = 1;

#ifdef HAVE_MPI
  unsigned short sbuf_NotMatching = rbuf_NotMatching;
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_MAX, MPI_COMM_WORLD);
#endif

  if (rbuf_NotMatching != 0)
    SU2_MPI::Error(string("The solution file ") + restart_filename +
                   string(" doesn't match with the mesh file!\n") +
                   string("It could be empty lines at the end of the file."),
                   CURRENT_FUNCTION);

  /*--- Delete the class memory that is used to load the restart. ---*/

  if (Restart_Vars != NULL) delete [] Restart_Vars;
  if (Restart_Data != NULL) delete [] Restart_Data;
  Restart_Vars = NULL; Restart_Data = NULL;

}

CBaselineSolver_FEM::~CBaselineSolver_FEM(void) { }
