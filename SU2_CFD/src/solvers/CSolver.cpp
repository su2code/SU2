/*!
 * \file CSolver.cpp
 * \brief Main subroutines for CSolver class.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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


#include "../../include/solvers/CSolver.hpp"
#include "../../include/gradients/computeGradientsGreenGauss.hpp"
#include "../../include/gradients/computeGradientsLeastSquares.hpp"
#include "../../include/limiters/computeLimiters.hpp"
#include "../../../Common/include/toolboxes/MMS/CIncTGVSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CInviscidVortexSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSIncEulerSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSIncNSSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSNSTwoHalfCirclesSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSNSTwoHalfSpheresSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSNSUnitQuadSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CMMSNSUnitQuadSolutionWallBC.hpp"
#include "../../../Common/include/toolboxes/MMS/CNSUnitQuadSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CRinglebSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CTGVSolution.hpp"
#include "../../../Common/include/toolboxes/MMS/CUserDefinedSolution.hpp"
#include "../../../Common/include/toolboxes/printing_toolbox.hpp"
#include "../../../Common/include/toolboxes/C1DInterpolation.hpp"
#include "../../include/CMarkerProfileReaderFVM.hpp"


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
  ResLinSolver       = 0.0;

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

  /*--- Auxiliary data needed for CFL adaption. ---*/

  NonLinRes_Value = 0;
  NonLinRes_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  NonLinRes_Counter = 0;

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

  if (VerificationSolution != NULL) {delete VerificationSolution; VerificationSolution = NULL;}

}

void CSolver::InitiatePeriodicComms(CGeometry *geometry,
                                    CConfig *config,
                                    unsigned short val_periodic_index,
                                    unsigned short commType) {

  /*--- Check for dummy communication. ---*/

  if (commType == PERIODIC_NONE) return;

  /*--- Local variables ---*/

  bool boundary_i, boundary_j;
  bool weighted = true;

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
    case PERIODIC_SOL_LS: case PERIODIC_SOL_ULS:
      COUNT_PER_POINT  = nDim*nDim + nVar*nDim;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PERIODIC_PRIM_LS: case PERIODIC_PRIM_ULS:
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
              bufDSend[buf_offset+iVar] = LinSysRes(iPoint, iVar);
            }

            /*--- Rotate the momentum components of the residual array. ---*/

            if (rotate_periodic) {
              if (nDim == 2) {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*LinSysRes(iPoint, 1) +
                                          rotMatrix[0][1]*LinSysRes(iPoint, 2));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*LinSysRes(iPoint, 1) +
                                          rotMatrix[1][1]*LinSysRes(iPoint, 2));
              } else {
                bufDSend[buf_offset+1] = (rotMatrix[0][0]*LinSysRes(iPoint, 1) +
                                          rotMatrix[0][1]*LinSysRes(iPoint, 2) +
                                          rotMatrix[0][2]*LinSysRes(iPoint, 3));
                bufDSend[buf_offset+2] = (rotMatrix[1][0]*LinSysRes(iPoint, 1) +
                                          rotMatrix[1][1]*LinSysRes(iPoint, 2) +
                                          rotMatrix[1][2]*LinSysRes(iPoint, 3));
                bufDSend[buf_offset+3] = (rotMatrix[2][0]*LinSysRes(iPoint, 1) +
                                          rotMatrix[2][1]*LinSysRes(iPoint, 2) +
                                          rotMatrix[2][2]*LinSysRes(iPoint, 3));
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

          case PERIODIC_SOL_LS: case PERIODIC_SOL_ULS:

            /*--- For L-S gradient calculations with rotational periodicity,
             we will need to rotate the x,y,z components. To make the process
             easier, we choose to rotate the initial periodic point and their
             neighbor points into their location on the donor marker before
             computing the terms that we need to communicate. ---*/

            /*--- Set a flag for unweighted or weighted least-squares. ---*/

            weighted = true;
            if (commType == PERIODIC_SOL_ULS) {
              weighted = false;
            }

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

                if (weighted) {
                  weight = 0.0;
                  for (iDim = 0; iDim < nDim; iDim++) {
                    weight += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                               (rotCoord_j[iDim]-rotCoord_i[iDim]));
                  }
                } else {
                  weight = 1.0;
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

          case PERIODIC_PRIM_LS: case PERIODIC_PRIM_ULS:

            /*--- For L-S gradient calculations with rotational periodicity,
             we will need to rotate the x,y,z components. To make the process
             easier, we choose to rotate the initial periodic point and their
             neighbor points into their location on the donor marker before
             computing the terms that we need to communicate. ---*/

            /*--- Set a flag for unweighted or weighted least-squares. ---*/

            weighted = true;
            if (commType == PERIODIC_PRIM_ULS) {
              weighted = false;
            }

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

                if (weighted) {
                  weight = 0.0;
                  for (iDim = 0; iDim < nDim; iDim++) {
                    weight += ((rotCoord_j[iDim]-rotCoord_i[iDim])*
                               (rotCoord_j[iDim]-rotCoord_i[iDim]));
                  }
                } else {
                  weight = 1.0;
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

            /*--- We send the min and max over "our" neighbours. ---*/

            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              Sol_Min[iVar] = base_nodes->GetSolution_Min(iPoint, iVar);
              Sol_Max[iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
            }

            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                Sol_Min[iVar] = min(Sol_Min[iVar], base_nodes->GetPrimitive(jPoint, iVar));
                Sol_Max[iVar] = max(Sol_Max[iVar], base_nodes->GetPrimitive(jPoint, iVar));
              }
            }

            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              bufDSend[buf_offset+iVar]              = Sol_Min[iVar];
              bufDSend[buf_offset+nPrimVarGrad+iVar] = Sol_Max[iVar];
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

            /*--- We send the min and max over "our" neighbours. ---*/

            for (iVar = 0; iVar < nVar; iVar++) {
              Sol_Min[iVar] = base_nodes->GetSolution_Min(iPoint, iVar);
              Sol_Max[iVar] = base_nodes->GetSolution_Max(iPoint, iVar);
            }

            for (iNeighbor = 0; iNeighbor < geometry->node[iPoint]->GetnPoint(); iNeighbor++) {
              jPoint = geometry->node[iPoint]->GetPoint(iNeighbor);
              for (iVar = 0; iVar < nVar; iVar++) {
                Sol_Min[iVar] = min(Sol_Min[iVar], base_nodes->GetSolution(jPoint, iVar));
                Sol_Max[iVar] = max(Sol_Max[iVar], base_nodes->GetSolution(jPoint, iVar));
              }
            }

            for (iVar = 0; iVar < nVar; iVar++) {
              bufDSend[buf_offset+iVar]      = Sol_Min[iVar];
              bufDSend[buf_offset+nVar+iVar] = Sol_Max[iVar];
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

  /*--- Check for dummy communication. ---*/

  if (commType == PERIODIC_NONE) return;

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

                Jacobian.AddBlock2Diag(iPoint, Jacobian_i);

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

            case PERIODIC_SOL_LS: case PERIODIC_SOL_ULS:

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

            case PERIODIC_PRIM_LS: case PERIODIC_PRIM_ULS:

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

              /*--- Update solution min/max with min/max between "us" and
               the periodic match plus its neighbors, computation will need to
               be concluded on "our" side to account for "our" neighbors. ---*/

              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {

                /*--- Solution minimum. ---*/

                Solution_Min = min(base_nodes->GetSolution_Min(iPoint, iVar),
                                   bufDRecv[buf_offset+iVar]);
                base_nodes->SetSolution_Min(iPoint, iVar, Solution_Min);

                /*--- Solution maximum. ---*/

                Solution_Max = max(base_nodes->GetSolution_Max(iPoint, iVar),
                                   bufDRecv[buf_offset+nPrimVarGrad+iVar]);
                base_nodes->SetSolution_Max(iPoint, iVar, Solution_Max);
              }

              break;

            case PERIODIC_LIM_PRIM_2:

              /*--- Check the min values found on the matching periodic
               faces for the limiter, and store the proper min value. ---*/

              for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
                Limiter_Min = min(base_nodes->GetLimiter_Primitive(iPoint, iVar),
                                  bufDRecv[buf_offset+iVar]);
                base_nodes->SetLimiter_Primitive(iPoint, iVar, Limiter_Min);
              }

              break;

            case PERIODIC_LIM_SOL_1:

              /*--- Update solution min/max with min/max between "us" and
               the periodic match plus its neighbors, computation will need to
               be concluded on "our" side to account for "our" neighbors. ---*/

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
      COUNT_PER_POINT  = nVar*nDim*2;
      MPI_TYPE         = COMM_TYPE_DOUBLE;
      break;
    case PRIMITIVE_GRADIENT:
      COUNT_PER_POINT  = nPrimVarGrad*nDim*2;
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
            for (iVar = 0; iVar < nVar; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+iVar*nDim+iDim] = base_nodes->GetGradient(iPoint, iVar, iDim);
                bufDSend[buf_offset+iVar*nDim+iDim+nDim*nVar] = base_nodes->GetGradient_Reconstruction(iPoint, iVar, iDim);
              }
            }
            break;
          case PRIMITIVE_GRADIENT:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                bufDSend[buf_offset+iVar*nDim+iDim] = base_nodes->GetGradient_Primitive(iPoint, iVar, iDim);
                bufDSend[buf_offset+iVar*nDim+iDim+nDim*nPrimVarGrad] = base_nodes->GetGradient_Reconstruction(iPoint, iVar, iDim);
              }
            }
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
              base_nodes->SetUnd_Lapl(iPoint, iVar, bufDRecv[buf_offset+iVar]);
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
            for (iVar = 0; iVar < nVar; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                base_nodes->SetGradient(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim]);
                base_nodes->SetGradient_Reconstruction(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim+nDim*nVar]);
              }
            }
            break;
          case PRIMITIVE_GRADIENT:
            for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
              for (iDim = 0; iDim < nDim; iDim++) {
                base_nodes->SetGradient_Primitive(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim]);
                base_nodes->SetGradient_Reconstruction(iPoint, iVar, iDim, bufDRecv[buf_offset+iVar*nDim+iDim+nDim*nPrimVarGrad]);
              }
            }
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

void CSolver::ResetCFLAdapt() {
  NonLinRes_Series.clear();
  NonLinRes_Value = 0;
  NonLinRes_Func = 0;
  Old_Func = 0;
  New_Func = 0;
  NonLinRes_Counter = 0;
}


void CSolver::AdaptCFLNumber(CGeometry **geometry,
                             CSolver   ***solver_container,
                             CConfig   *config) {

  /* Adapt the CFL number on all multigrid levels using an
   exponential progression with under-relaxation approach. */

  vector<su2double> MGFactor(config->GetnMGLevels()+1,1.0);
  const su2double CFLFactorDecrease = config->GetCFL_AdaptParam(0);
  const su2double CFLFactorIncrease = config->GetCFL_AdaptParam(1);
  const su2double CFLMin            = config->GetCFL_AdaptParam(2);
  const su2double CFLMax            = config->GetCFL_AdaptParam(3);
  const bool fullComms              = (config->GetComm_Level() == COMM_FULL);

  for (unsigned short iMesh = 0; iMesh <= config->GetnMGLevels(); iMesh++) {

    /* Store the mean flow, and turbulence solvers more clearly. */

    CSolver *solverFlow = solver_container[iMesh][FLOW_SOL];
    CSolver *solverTurb = solver_container[iMesh][TURB_SOL];

    /* Compute the reduction factor for CFLs on the coarse levels. */

    if (iMesh == MESH_0) {
      MGFactor[iMesh] = 1.0;
    } else {
      const su2double CFLRatio = config->GetCFL(iMesh)/config->GetCFL(iMesh-1);
      MGFactor[iMesh] = MGFactor[iMesh-1]*CFLRatio;
    }

    /* Check whether we achieved the requested reduction in the linear
     solver residual within the specified number of linear iterations. */

    bool reduceCFL = false;
    su2double linResFlow = solverFlow->GetResLinSolver();
    su2double linResTurb = -1.0;
    if ((iMesh == MESH_0) && (config->GetKind_Turb_Model() != NONE)) {
      linResTurb = solverTurb->GetResLinSolver();
    }

    su2double maxLinResid = max(linResFlow, linResTurb);
    if (maxLinResid > 0.5) {
      reduceCFL = true;
    }

    /* Check that we are meeting our nonlinear residual reduction target
     over time so that we do not get stuck in limit cycles. */

    SU2_OMP_MASTER
    { /* Only the master thread updates the shared variables. */

    Old_Func = New_Func;
    unsigned short Res_Count = 100;
    if (NonLinRes_Series.size() == 0) NonLinRes_Series.resize(Res_Count,0.0);

    /* Sum the RMS residuals for all equations. */

    New_Func = 0.0;
    for (unsigned short iVar = 0; iVar < solverFlow->GetnVar(); iVar++) {
      New_Func += solverFlow->GetRes_RMS(iVar);
    }
    if ((iMesh == MESH_0) && (config->GetKind_Turb_Model() != NONE)) {
      for (unsigned short iVar = 0; iVar < solverTurb->GetnVar(); iVar++) {
        New_Func += solverTurb->GetRes_RMS(iVar);
      }
    }

    /* Compute the difference in the nonlinear residuals between the
     current and previous iterations. */

    NonLinRes_Func = (New_Func - Old_Func);
    NonLinRes_Series[NonLinRes_Counter] = NonLinRes_Func;

    /* Increment the counter, if we hit the max size, then start over. */

    NonLinRes_Counter++;
    if (NonLinRes_Counter == Res_Count) NonLinRes_Counter = 0;

    /* Sum the total change in nonlinear residuals over the previous
     set of all stored iterations. */

    NonLinRes_Value = New_Func;
    if (config->GetTimeIter() >= Res_Count) {
      NonLinRes_Value = 0.0;
      for (unsigned short iCounter = 0; iCounter < Res_Count; iCounter++)
        NonLinRes_Value += NonLinRes_Series[iCounter];
    }

    /* If the sum is larger than a small fraction of the current nonlinear
     residual, then we are not decreasing the nonlinear residual at a high
     rate. In this situation, we force a reduction of the CFL in all cells.
     Reset the array so that we delay the next decrease for some iterations. */

    if (fabs(NonLinRes_Value) < 0.1*New_Func) {
      NonLinRes_Counter = 0;
      for (unsigned short iCounter = 0; iCounter < Res_Count; iCounter++)
        NonLinRes_Series[iCounter] = New_Func;
    }

    } /* End SU2_OMP_MASTER, now all threads update the CFL number. */
    SU2_OMP_BARRIER

    if (fabs(NonLinRes_Value) < 0.1*New_Func) {
      reduceCFL = true;
    }

    /* Loop over all points on this grid and apply CFL adaption. */

    su2double myCFLMin = 1e30, myCFLMax = 0.0, myCFLSum = 0.0;

    SU2_OMP_MASTER
    if ((iMesh == MESH_0) && fullComms) {
      Min_CFL_Local = 1e30;
      Max_CFL_Local = 0.0;
      Avg_CFL_Local = 0.0;
    }

    SU2_OMP_FOR_STAT(roundUpDiv(geometry[iMesh]->GetnPointDomain(),omp_get_max_threads()))
    for (unsigned long iPoint = 0; iPoint < geometry[iMesh]->GetnPointDomain(); iPoint++) {

      /* Get the current local flow CFL number at this point. */

      su2double CFL = solverFlow->GetNodes()->GetLocalCFL(iPoint);

      /* Get the current under-relaxation parameters that were computed
       during the previous nonlinear update. If we have a turbulence model,
       take the minimum under-relaxation parameter between the mean flow
       and turbulence systems. */

      su2double underRelaxationFlow = solverFlow->GetNodes()->GetUnderRelaxation(iPoint);
      su2double underRelaxationTurb = 1.0;
      if ((iMesh == MESH_0) && (config->GetKind_Turb_Model() != NONE))
        underRelaxationTurb = solverTurb->GetNodes()->GetUnderRelaxation(iPoint);
      const su2double underRelaxation = min(underRelaxationFlow,underRelaxationTurb);

      /* If we apply a small under-relaxation parameter for stability,
       then we should reduce the CFL before the next iteration. If we
       are able to add the entire nonlinear update (under-relaxation = 1)
       then we schedule an increase the CFL number for the next iteration. */

      su2double CFLFactor = 1.0;
      if ((underRelaxation < 0.1)) {
        CFLFactor = CFLFactorDecrease;
      } else if (underRelaxation >= 0.1 && underRelaxation < 1.0) {
        CFLFactor = 1.0;
      } else {
        CFLFactor = CFLFactorIncrease;
      }

      /* Check if we are hitting the min or max and adjust. */

      if (CFL*CFLFactor <= CFLMin) {
        CFL       = CFLMin;
        CFLFactor = MGFactor[iMesh];
      } else if (CFL*CFLFactor >= CFLMax) {
        CFL       = CFLMax;
        CFLFactor = MGFactor[iMesh];
      }

      /* If we detect a stalled nonlinear residual, then force the CFL
       for all points to the minimum temporarily to restart the ramp. */

      if (reduceCFL) {
        CFL       = CFLMin;
        CFLFactor = MGFactor[iMesh];
      }

      /* Apply the adjustment to the CFL and store local values. */

      CFL *= CFLFactor;
      solverFlow->GetNodes()->SetLocalCFL(iPoint, CFL);
      if ((iMesh == MESH_0) && (config->GetKind_Turb_Model() != NONE)) {
        solverTurb->GetNodes()->SetLocalCFL(iPoint, CFL);
      }

      /* Store min and max CFL for reporting on the fine grid. */

      if ((iMesh == MESH_0) && fullComms) {
        myCFLMin = min(CFL,myCFLMin);
        myCFLMax = max(CFL,myCFLMax);
        myCFLSum += CFL;
      }

    }

    /* Reduce the min/max/avg local CFL numbers. */

    if ((iMesh == MESH_0) && fullComms) {
      SU2_OMP_CRITICAL
      { /* OpenMP reduction. */
        Min_CFL_Local = min(Min_CFL_Local,myCFLMin);
        Max_CFL_Local = max(Max_CFL_Local,myCFLMax);
        Avg_CFL_Local += myCFLSum;
      }
      SU2_OMP_BARRIER

      SU2_OMP_MASTER
      { /* MPI reduction. */
        myCFLMin = Min_CFL_Local; myCFLMax = Max_CFL_Local; myCFLSum = Avg_CFL_Local;
        SU2_MPI::Allreduce(&myCFLMin, &Min_CFL_Local, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&myCFLMax, &Max_CFL_Local, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        SU2_MPI::Allreduce(&myCFLSum, &Avg_CFL_Local, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        Avg_CFL_Local /= su2double(geometry[iMesh]->GetGlobal_nPointDomain());
      }
      SU2_OMP_BARRIER
    }

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

  sbuf_residual = new su2double[nVar];
  rbuf_residual = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) sbuf_residual[iVar] = GetRes_RMS(iVar);

  if (config->GetComm_Level() == COMM_FULL) {

    unsigned long Local_nPointDomain = geometry->GetnPointDomain();
    SU2_MPI::Allreduce(sbuf_residual, rbuf_residual, nVar, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  }
  else {
    /*--- Reduced MPI comms have been requested. Use a local residual only. ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      rbuf_residual[iVar] = sbuf_residual[iVar];
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

    sbuf_residual = new su2double [nVar]();
    sbuf_point = new unsigned long [nVar]();
    sbuf_coord = new su2double[nVar*nDim]();

    rbuf_residual = new su2double [nProcessor*nVar]();
    rbuf_point = new unsigned long [nProcessor*nVar]();
    rbuf_coord = new su2double[nProcessor*nVar*nDim]();

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
  su2double ProjGridVel;
  const su2double* Normal;

  /*--- Loop interior edges ---*/

  for (iEdge = 0; iEdge < geometry->GetnEdge(); iEdge++) {

    const unsigned long iPoint = geometry->edges->GetNode(iEdge,0);
    const unsigned long jPoint = geometry->edges->GetNode(iEdge,1);

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

    Normal = geometry->edges->GetNormal(iEdge);

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

  const auto solution = base_nodes->GetAuxVar();
  auto gradient = base_nodes->GetAuxVarGradient();

  computeGradientsGreenGauss(this, AUXVAR_GRADIENT, PERIODIC_NONE, *geometry,
                             *config, solution, 0, 1, gradient);
}

void CSolver::SetAuxVar_Gradient_LS(CGeometry *geometry, CConfig *config) {

  bool weighted = true;
  const auto solution = base_nodes->GetAuxVar();
  auto gradient = base_nodes->GetAuxVarGradient();
  auto& rmatrix  = base_nodes->GetRmatrix();

  computeGradientsLeastSquares(this, AUXVAR_GRADIENT, PERIODIC_NONE, *geometry, *config,
                               weighted, solution, 0, 1, gradient, rmatrix);
}

void CSolver::SetSolution_Gradient_GG(CGeometry *geometry, CConfig *config, bool reconstruction) {

  const auto& solution = base_nodes->GetSolution();
  auto& gradient = reconstruction? base_nodes->GetGradient_Reconstruction() : base_nodes->GetGradient();

  computeGradientsGreenGauss(this, SOLUTION_GRADIENT, PERIODIC_SOL_GG, *geometry,
                             *config, solution, 0, nVar, gradient);
}

void CSolver::SetSolution_Gradient_LS(CGeometry *geometry, CConfig *config, bool reconstruction) {

  /*--- Set a flag for unweighted or weighted least-squares. ---*/
  bool weighted;

  if (reconstruction)
    weighted = (config->GetKind_Gradient_Method_Recon() == WEIGHTED_LEAST_SQUARES);
  else
    weighted = (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES);

  const auto& solution = base_nodes->GetSolution();
  auto& rmatrix = base_nodes->GetRmatrix();
  auto& gradient = reconstruction? base_nodes->GetGradient_Reconstruction() : base_nodes->GetGradient();
  PERIODIC_QUANTITIES kindPeriodicComm = weighted? PERIODIC_SOL_LS : PERIODIC_SOL_ULS;

  computeGradientsLeastSquares(this, SOLUTION_GRADIENT, kindPeriodicComm, *geometry, *config,
                               weighted, solution, 0, nVar, gradient, rmatrix);
}

void CSolver::Add_External_To_Solution() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    base_nodes->AddSolution(iPoint, base_nodes->Get_External(iPoint));
  }
}

void CSolver::Add_Solution_To_External() {
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    base_nodes->Add_External(iPoint, base_nodes->GetSolution(iPoint));
  }
}

void CSolver::Update_Cross_Term(CConfig *config, su2passivematrix &cross_term) {

  /*--- This method is for discrete adjoint solvers and it is used in multi-physics
   *    contexts, "cross_term" is the old value, the new one is in "Solution".
   *    We update "cross_term" and the sum of all cross terms (in "External")
   *    with a fraction of the difference between new and old.
   *    When "alpha" is 1, i.e. no relaxation, we effectively subtract the old
   *    value and add the new one to the total ("External"). ---*/

  passivedouble alpha = SU2_TYPE::GetValue(config->GetAitkenStatRelax());

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      passivedouble
      new_val = SU2_TYPE::GetValue(base_nodes->GetSolution(iPoint,iVar)),
      delta = alpha * (new_val - cross_term(iPoint,iVar));
      /*--- Update cross term. ---*/
      cross_term(iPoint,iVar) += delta;
      Solution[iVar] = delta;
    }
    /*--- Update the sum of all cross-terms. ---*/
    base_nodes->Add_External(iPoint, Solution);
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

  auto kindLimiter = static_cast<ENUM_LIMITER>(config->GetKind_SlopeLimit());
  const auto& solution = base_nodes->GetSolution();
  const auto& gradient = base_nodes->GetGradient_Reconstruction();
  auto& solMin = base_nodes->GetSolution_Min();
  auto& solMax = base_nodes->GetSolution_Max();
  auto& limiter = base_nodes->GetLimiter();

  computeLimiters(kindLimiter, this, SOLUTION_LIMITER, PERIODIC_LIM_SOL_1, PERIODIC_LIM_SOL_2,
                  *geometry, *config, 0, nVar, solution, gradient, solMin, solMax, limiter);
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

  /*--- First, we load the restart file for time n ---*/

  /*-------------------------------------------------------------------------------------------*/

  /*--- Modify file name for an unsteady restart ---*/
  if (config->GetRestart()) Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-1;
  else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-1;
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

  if (iPoint_Global_Local < geometry->GetnPointDomain()) {
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
    if (config->GetRestart()) Unst_RestartIter = SU2_TYPE::Int(config->GetRestart_Iter())-2;
    else Unst_RestartIter = SU2_TYPE::Int(config->GetUnst_AdjointIter())-2;
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

    if (iPoint_Global_Local < geometry->GetnPointDomain()) {
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
    if ((config->GetAoA() != AoA_) && (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoA in the solution file (" << AoA_ << " deg.) +" << endl;
      cout << "         AoA offset in mesh file (" << config->GetAoA_Offset() << " deg.) = " << AoA_ + config->GetAoA_Offset() << " deg." << endl;
    }
    config->SetAoA(AoA_ + config->GetAoA_Offset());
  }

  else {
    if ((config->GetAoA() != AoA_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoA in the solution file." << endl;
  }

  /*--- Sideslip angle ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetAoS() != AoS_) && (rank == MASTER_NODE)) {
      cout.precision(6);
      cout <<"WARNING: AoS in the solution file (" << AoS_ << " deg.) +" << endl;
      cout << "         AoS offset in mesh file (" << config->GetAoS_Offset() << " deg.) = " << AoS_ + config->GetAoS_Offset() << " deg." << endl;
    }
    config->SetAoS(AoS_ + config->GetAoS_Offset());
  }
  else {
    if ((config->GetAoS() != AoS_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the AoS in the solution file." << endl;
  }

  /*--- BCThrust ---*/

  if (config->GetDiscard_InFiles() == false) {
    if ((config->GetInitial_BCThrust() != BCThrust_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the initial BC Thrust provided in the solution file: " << BCThrust_ << " lbs." << endl;
    config->SetInitial_BCThrust(BCThrust_);
  }
  else {
    if ((config->GetInitial_BCThrust() != BCThrust_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the BC Thrust in the solution file." << endl;
  }


  if (config->GetDiscard_InFiles() == false) {

    if ((config->GetdCD_dCL() != dCD_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCD/dCL provided in the direct solution file: " << dCD_dCL_ << "." << endl;
    config->SetdCD_dCL(dCD_dCL_);

    if ((config->GetdCMx_dCL() != dCMx_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMx/dCL provided in the direct solution file: " << dCMx_dCL_ << "." << endl;
    config->SetdCMx_dCL(dCMx_dCL_);

    if ((config->GetdCMy_dCL() != dCMy_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMy/dCL provided in the direct solution file: " << dCMy_dCL_ << "." << endl;
    config->SetdCMy_dCL(dCMy_dCL_);

    if ((config->GetdCMz_dCL() != dCMz_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: SU2 will use the dCMz/dCL provided in the direct solution file: " << dCMz_dCL_ << "." << endl;
    config->SetdCMz_dCL(dCMz_dCL_);

  }

  else {

    if ((config->GetdCD_dCL() != dCD_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCD/dCL in the direct solution file." << endl;

    if ((config->GetdCMx_dCL() != dCMx_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMx/dCL in the direct solution file." << endl;

    if ((config->GetdCMy_dCL() != dCMy_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMy/dCL in the direct solution file." << endl;

    if ((config->GetdCMz_dCL() != dCMz_dCL_) && (rank == MASTER_NODE))
      cout <<"WARNING: Discarding the dCMz/dCL in the direct solution file." << endl;

  }

  /*--- External iteration ---*/

  if ((config->GetDiscard_InFiles() == false) && (!adjoint || (adjoint && config->GetRestart())))
    config->SetExtIter_OffSet(InnerIter_);

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
  su2double Area_Children, Area_Parent, dist, min_dist, Interp_Radius, Theta;
  const su2double *Coord = nullptr;
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
  string Interpolation_Function, Interpolation_Type;
  bool Interpolate = false;

  su2double *Normal = new su2double[nDim];

  unsigned long Marker_Counter = 0;

  bool turbulent = (config->GetKind_Solver() == RANS ||
                    config->GetKind_Solver() == INC_RANS ||
                    config->GetKind_Solver() == ADJ_RANS ||
                    config->GetKind_Solver() == DISC_ADJ_RANS ||
                    config->GetKind_Solver() == DISC_ADJ_INC_RANS);

  unsigned short nVar_Turb = 0;
  if (turbulent)
    switch (config->GetKind_Turb_Model()) {
      case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
        nVar_Turb = 1;
        break;
      case SST: case SST_SUST:
        nVar_Turb = 2;
        break;
      default:
        SU2_MPI::Error("Specified turbulence model unavailable or none selected", CURRENT_FUNCTION);
        break;
    }

  /*--- Count the number of columns that we have for this flow case,
   excluding the coordinates. Here, we have 2 entries for the total
   conditions or mass flow, another nDim for the direction vector, and
   finally entries for the number of turbulence variables. This is only
   necessary in case we are writing a template profile file or for Inlet
   Interpolation purposes. ---*/

  unsigned short nCol_InletFile = 2 + nDim + nVar_Turb;

  /*--- Multizone problems require the number of the zone to be appended. ---*/

  if (nZone > 1)
    profile_filename = config->GetMultizone_FileName(profile_filename, iZone, ".dat");

  /*--- Modify file name for an unsteady restart ---*/

  if (dual_time || time_stepping)
    profile_filename = config->GetUnsteady_FileName(profile_filename, val_iter, ".dat");

  /*--- Read the profile data from an ASCII file. ---*/

  CMarkerProfileReaderFVM profileReader(geometry[MESH_0], config, profile_filename, KIND_MARKER, nCol_InletFile);

  /*--- Load data from the restart into correct containers. ---*/

  Marker_Counter = 0;

  unsigned short global_failure = 0, local_failure = 0;
  ostringstream error_msg;

  const su2double tolerance = config->GetInlet_Profile_Matching_Tolerance();

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    /*--- Skip if this is the wrong type of marker. ---*/

    if (config->GetMarker_All_KindBC(iMarker) != KIND_MARKER) continue;

    /*--- Get tag in order to identify the correct inlet data. ---*/

    Marker_Tag = config->GetMarker_All_TagBound(iMarker);

    for (jMarker = 0; jMarker < profileReader.GetNumberOfProfiles(); jMarker++) {

      /*--- If we have not found the matching marker string, continue to next marker. ---*/

      if (profileReader.GetTagForProfile(jMarker) != Marker_Tag) continue;

      /*--- Increment our counter for marker matches. ---*/

      Marker_Counter++;

      /*--- Get data for this profile. ---*/

      vector<passivedouble> Inlet_Data = profileReader.GetDataForProfile(jMarker);
      unsigned short nColumns = profileReader.GetNumberOfColumnsInProfile(jMarker);
      vector<su2double> Inlet_Data_Interpolated ((nCol_InletFile+nDim)*geometry[MESH_0]->nVertex[iMarker]);

      /*--- Define Inlet Values vectors before and after interpolation (if needed) ---*/
      vector<su2double> Inlet_Values(nCol_InletFile+nDim);
      vector<su2double> Inlet_Interpolated(nColumns);

      unsigned long nRows = profileReader.GetNumberOfRowsInProfile(jMarker);

      /*--- Pointer to call Set and Evaluate functions. ---*/
      vector<C1DInterpolation*> interpolator (nColumns);
      string interpolation_function, interpolation_type;

      /*--- Define the reference for interpolation. ---*/
      unsigned short radius_index=0;
      vector<su2double> InletRadii = profileReader.GetColumnForProfile(jMarker, radius_index);
      vector<su2double> Interpolation_Column (nRows);

      switch(config->GetKindInletInterpolationFunction()){

        case (NONE):
          Interpolate = false;
          break;

        case (AKIMA_1D):
          for (unsigned short iCol=0; iCol < nColumns; iCol++){
            Interpolation_Column = profileReader.GetColumnForProfile(jMarker, iCol);
            interpolator[iCol] = new CAkimaInterpolation(InletRadii,Interpolation_Column);
            interpolation_function = "AKIMA";
            Interpolate = true;
          }
          break;

        case (LINEAR_1D):
          for (unsigned short iCol=0; iCol < nColumns; iCol++){
            Interpolation_Column = profileReader.GetColumnForProfile(jMarker, iCol);
            interpolator[iCol] = new CLinearInterpolation(InletRadii,Interpolation_Column);
            interpolation_function = "LINEAR";
            Interpolate = true;
          }
          break;

        default:
          SU2_MPI::Error("Error in the Kind_InletInterpolation Marker\n",CURRENT_FUNCTION);
          break;
      }

      if (Interpolate == true){
        switch(config->GetKindInletInterpolationType()){
          case(VR_VTHETA):
            interpolation_type="VR_VTHETA";
            break;
          case(ALPHA_PHI):
            interpolation_type="ALPHA_PHI";
            break;
        }
        cout<<"Inlet Interpolation being done using "<<interpolation_function
            <<" function and type "<<interpolation_type<<" for "<< Marker_Tag<<endl;
        if(nDim == 3)
          cout<<"Ensure the flow direction is in z direction"<<endl;
        else if (nDim == 2)
          cout<<"Ensure the flow direction is in x direction"<<endl;
      }
      else if(Interpolate == false) {
        cout<<"No Inlet Interpolation being used"<<endl;
      }

      /*--- Loop through the nodes on this marker. ---*/

      for (iVertex = 0; iVertex < geometry[MESH_0]->nVertex[iMarker]; iVertex++) {

        iPoint = geometry[MESH_0]->vertex[iMarker][iVertex]->GetNode();
        Coord = geometry[MESH_0]->node[iPoint]->GetCoord();

        if(Interpolate == false) {

          min_dist = 1e16;

          /*--- Find the distance to the closest point in our inlet profile data. ---*/

          for (iRow = 0; iRow < nRows; iRow++) {

            /*--- Get the coords for this data point. ---*/

            index = iRow*nColumns;

            dist = 0.0;
            for (unsigned short iDim = 0; iDim < nDim; iDim++)
            dist += pow(Inlet_Data[index+iDim] - Coord[iDim], 2);
            dist = sqrt(dist);

            /*--- Check is this is the closest point and store data if so. ---*/

            if (dist < min_dist) {
            min_dist = dist;
            for (iVar = 0; iVar < nColumns; iVar++)
              Inlet_Values[iVar] = Inlet_Data[index+iVar];
            }

          }

          /*--- If the diff is less than the tolerance, match the two.
          We could modify this to simply use the nearest neighbor, or
          eventually add something more elaborate here for interpolation. ---*/

          if (min_dist < tolerance) {

            solver[MESH_0][KIND_SOLVER]->SetInletAtVertex(Inlet_Values.data(), iMarker, iVertex);

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

        else if(Interpolate == true) {

          /* --- Calculating the radius and angle of the vertex ---*/
          /* --- Flow should be in z direction for 3D cases ---*/
          /* --- Or in x direction for 2D cases ---*/
          Interp_Radius = sqrt(pow(Coord[0],2)+ pow(Coord[1],2));
          Theta = atan2(Coord[1],Coord[0]);

          /* --- Evaluating and saving the final spline data ---*/
          for  (unsigned short iVar=0; iVar < nColumns; iVar++){

            /*---Evaluate spline will get the respective value of the Data set (column) specified
            for that interpolator[iVar], cycling through all columns to get all the
            data for that vertex ---*/
            Inlet_Interpolated[iVar]=interpolator[iVar]->EvaluateSpline(Interp_Radius);
            if (interpolator[iVar]->GetPointMatch() == false){
              cout << "WARNING: Did not find a match between the radius in the inlet file " ;
              cout << std::scientific;
              cout << "at location: [" << Coord[0] << ", " << Coord[1];
              if (nDim == 3) {cout << ", " << Coord[2];}
              cout << "]";
              cout << " with Radius: "<< Interp_Radius << endl;
              cout << "You can add a row for Radius: " << Interp_Radius <<" in the inlet file ";
              cout << "to eliminate this issue or give proper data" << endl;
              local_failure++;
              break;
            }
          }

          /* --- Correcting for Interpolation Type ---*/
          switch(config->GetKindInletInterpolationType()){
          case(VR_VTHETA):
            Inlet_Values = CorrectedInletValues(Inlet_Interpolated, Theta, nDim, Coord, nVar_Turb, VR_VTHETA);
          break;
          case(ALPHA_PHI):
            Inlet_Values = CorrectedInletValues(Inlet_Interpolated, Theta, nDim, Coord, nVar_Turb, ALPHA_PHI);
          break;
          }

          solver[MESH_0][KIND_SOLVER]->SetInletAtVertex(Inlet_Values.data(), iMarker, iVertex);

          for (unsigned short iVar=0; iVar < (nCol_InletFile+nDim); iVar++)
            Inlet_Data_Interpolated[iVertex*(nCol_InletFile+nDim)+iVar] = Inlet_Values[iVar];

        }

      } // end iVertex loop

      if(config->GetPrintInlet_InterpolatedData() == true) {
          PrintInletInterpolatedData(Inlet_Data_Interpolated, profileReader.GetTagForProfile(jMarker),
                                     geometry[MESH_0]->nVertex[iMarker], nDim, nCol_InletFile+nDim);
      }

      for (int i=0; i<nColumns;i++)
        delete interpolator[i];

    } // end jMarker loop

    if (local_failure > 0) break;

  } // end iMarker loop

  SU2_MPI::Allreduce(&local_failure, &global_failure, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);

  if (global_failure > 0) {
    SU2_MPI::Error("Prescribed inlet data does not match markers within tolerance.", CURRENT_FUNCTION);
  }

  /*--- Copy the inlet data down to the coarse levels if multigrid is active.
   Here, we use a face area-averaging to restrict the values. ---*/

  for (iMesh = 1; iMesh <= config->GetnMGLevels(); iMesh++) {
    for (iMarker=0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == KIND_MARKER) {

        Marker_Tag = config->GetMarker_All_TagBound(iMarker);

        /* Check the number of columns and allocate temp array. */

        unsigned short nColumns = 0;
        for (jMarker = 0; jMarker < profileReader.GetNumberOfProfiles(); jMarker++) {
          if (profileReader.GetTagForProfile(jMarker) == Marker_Tag) {
            nColumns = profileReader.GetNumberOfColumnsInProfile(jMarker);
          }
        }
        vector<su2double> Inlet_Values(nColumns);
        vector<su2double> Inlet_Fine(nColumns);

        /*--- Loop through the nodes on this marker. ---*/

        for (iVertex = 0; iVertex < geometry[iMesh]->nVertex[iMarker]; iVertex++) {

          /*--- Get the coarse mesh point and compute the boundary area. ---*/

          iPoint = geometry[iMesh]->vertex[iMarker][iVertex]->GetNode();
          geometry[iMesh]->vertex[iMarker][iVertex]->GetNormal(Normal);
          Area_Parent = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) Area_Parent += Normal[iDim]*Normal[iDim];
          Area_Parent = sqrt(Area_Parent);

          /*--- Reset the values for the coarse point. ---*/

          for (iVar = 0; iVar < nColumns; iVar++) Inlet_Values[iVar] = 0.0;

          /*-- Loop through the children and extract the inlet values
           from those nodes that lie on the boundary as well as their
           boundary area. We build a face area-averaged value for the
           coarse point values from the fine grid points. Note that
           children from the interior volume will not be included in
           the averaging. ---*/

          for (iChildren = 0; iChildren < geometry[iMesh]->node[iPoint]->GetnChildren_CV(); iChildren++) {
            Point_Fine = geometry[iMesh]->node[iPoint]->GetChildren_CV(iChildren);
            for (iVar = 0; iVar < nColumns; iVar++) Inlet_Fine[iVar] = 0.0;
            Area_Children = solver[iMesh-1][KIND_SOLVER]->GetInletAtVertex(Inlet_Fine.data(), Point_Fine, KIND_MARKER,
                                                                           Marker_Tag, geometry[iMesh-1], config);
            for (iVar = 0; iVar < nColumns; iVar++) {
              Inlet_Values[iVar] += Inlet_Fine[iVar]*Area_Children/Area_Parent;
            }
          }

          /*--- Set the boundary area-averaged inlet values for the coarse point. ---*/

          solver[iMesh][KIND_SOLVER]->SetInletAtVertex(Inlet_Values.data(), iMarker, iVertex);

        }
      }
    }
  }

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

  /*--- Set the residuals and BGSSolution_k to solution for next multizone outer iteration. ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    const su2double domain = (iPoint < nPointDomain);
    for (iVar = 0; iVar < nVar; iVar++) {
      residual = (base_nodes->Get_BGSSolution(iPoint,iVar) - base_nodes->Get_BGSSolution_k(iPoint,iVar))*domain;
      base_nodes->Set_BGSSolution_k(iPoint,iVar, base_nodes->Get_BGSSolution(iPoint,iVar));
      AddRes_BGS(iVar, residual*residual);
      AddRes_Max_BGS(iVar, fabs(residual), geometry->node[iPoint]->GetGlobalIndex(), geometry->node[iPoint]->GetCoord());
    }
  }

  SetResidual_BGS(geometry, config);

}
