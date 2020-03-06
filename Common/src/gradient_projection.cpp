/*!
 * \file gradient_projection.cpp
 * \brief Functions for projections between mesh coordinate sensitivities and design values
 *        similar to the functions in SU2_DOT
 * \author T.Dick
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../include/gradient_projection.hpp"

/* development comments:
 * Call the AD routine for each point to get the whole Jacobian
 * do we need to use iDV_Value=0.0 or do we use the current position???
 * 31.01.2020 T. Dick
 */
void GetParameterizationJacobianForward(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double **Jacobian) {

  su2double DV_Value = 0.0;

  unsigned short iDV = 0, iDV_Value = 0, jDV=0, jDV_Value=0;
  unsigned short iDV_Index = 0;
  unsigned long iPoint;
  su2double* VarCoord;
  unsigned short nDim = geometry->GetnDim();

  // loop over all design values and claculate the gradient for each of them
  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++) {

      // set the derivatives to the i-th unit vector
      for (jDV = 0; jDV < config->GetnDV(); jDV++) {
        for (jDV_Value = 0; jDV_Value < config->GetnDV_Value(jDV); jDV_Value++) {
          DV_Value = config->GetDV_Value(jDV, jDV_Value);
          DV_Value = 0.0;
          SU2_TYPE::SetDerivative(DV_Value, 0.0);
          config->SetDV_Value(jDV, jDV_Value, DV_Value);
        }
      }
      DV_Value = config->GetDV_Value(iDV, iDV_Value);
      DV_Value = 0.0;
      SU2_TYPE::SetDerivative(DV_Value, 1.0);
      config->SetDV_Value(iDV, iDV_Value, DV_Value);

      // run the function in AD forward mode
      surface_movement->SetSurface_Deformation(geometry, config);

      // extract derivatives and store them in the matrix
      for (auto iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (auto iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
          for (auto iDim = 0; iDim < nDim; iDim++) {
            auto total_index = iPoint*nDim + iDim;
            Jacobian[iDV_Index][total_index] = SU2_TYPE::GetDerivative(VarCoord[iDim]);
          }
        }
      }

      iDV_Index++;
    }
  }

}

/* development comments:
 * Call the AD routine for each point to get the whole Jacobian
 * 31.01.2020 T. Dick
 */
void GetParameterizationJacobianReverse(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian) {

 /// new version

  unsigned short nDim, nMarker, nDV, nDV_Value, nDV_Total, nPoint, nTotal_Index, nVertex;
  unsigned short iDV, iDV_Value, iDV_index, iPoint, iDim, total_index, iMarker, iVertex, subDim;
  su2double* VarCoord;
  su2double DV_Value;
  double my_Gradient, localGradient;

  // get some numbers from config
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nPoint  = geometry->GetnPoint();
  nTotal_Index = nPoint*nDim;
  nDV     = config->GetnDV();
  nDV_Total = config->GetnDV_Total();

  // structure to calculate and manage the total number of marked points
  unsigned* visitedPoints = new unsigned[nPoint];
  for (iPoint=0; iPoint<nPoint; iPoint++){
    visitedPoints[iPoint]=0;
  }
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        if (visitedPoints[iPoint]==0){
          visitedPoints[iPoint] = 1;
        } else {
          visitedPoints[iPoint] = 1;
        }
      }
    }
  }

  /*--- Discrete adjoint tape recording ---*/

  if (SU2_MPI::GetRank() == MASTER_NODE)
    cout  << endl << "Taping surface parameterization using Algorithmic Differentiation (ZONE " << config->GetiZone() << ")." << endl;

  /*--- Start recording of operations ---*/

  AD::Reset();

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  for (iDV = 0; iDV < nDV; iDV++){

    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization with su2double resets the index ---*/

      DV_Value = 0.0;

      AD::RegisterInput(DV_Value);

      config->SetDV_Value(iDV, iDV_Value, DV_Value);
    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Stop the recording --- */

  AD::StopRecording();

  /*--- Loop over all visited points and calculate the derivative with respect to that point --*/

  for (iPoint=0; iPoint<nPoint; iPoint++) {

    /*--- Go over all visited points and dimensions ---*/
    if (visitedPoints[iPoint]==1) {
    for (subDim=0; subDim<nDim; subDim++) {

      /*--- Set the seeding appropriately ---*/
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        if (config->GetMarker_All_DV(iMarker) == YES) {
          nVertex = geometry->nVertex[iMarker];
          for (iVertex = 0; iVertex <nVertex; iVertex++) {
            VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
            if (iPoint == geometry->vertex[iMarker][iVertex]->GetNode()) {
              for (iDim = 0; iDim < nDim; iDim++){
                SU2_TYPE::SetDerivative(VarCoord[iDim], 0.0);
              }
              SU2_TYPE::SetDerivative(VarCoord[subDim], 1.0);
              total_index = iPoint*nDim+subDim;
            } else {
              for (iDim = 0; iDim < nDim; iDim++){
                SU2_TYPE::SetDerivative(VarCoord[iDim], 0.0);
              }
            }
          }
        }
      }

      /*--- Compute derivatives and extract gradient ---*/

      AD::ComputeAdjoint();

      iDV_index = 0;

      for (iDV = 0; iDV  < nDV; iDV++){
        nDV_Value =  config->GetnDV_Value(iDV);
        for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

          DV_Value = config->GetDV_Value(iDV, iDV_Value);
          my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
          #ifdef HAVE_MPI
            SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          #else
            localGradient = my_Gradient;
          #endif

          Jacobian[iDV_index*nPoint*nDim+total_index] = localGradient;

          iDV_index++;
        }
      }

      AD::ClearAdjoints();

    }
    }
  }

  AD::Reset();
}


MatrixType Cast2Eigenmatrix(CGeometry *geometry, CConfig *config, su2double *Jacobian) {

  unsigned nDV_Total = config->GetnDV_Total();

  unsigned short nDim    = geometry->GetnDim();
  unsigned long nPoint  = geometry->GetnPoint();

  Eigen::Map<Eigen::Matrix<su2double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Jacobi_eigen(&Jacobian[0], nDV_Total, nPoint*nDim);

  MatrixType Jacobi_eigen_trans = Jacobi_eigen.transpose();

  return Jacobi_eigen_trans;
}
