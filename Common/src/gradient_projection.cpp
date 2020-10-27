/*!
 * \file gradient_projection.cpp
 * \brief Functions for projections between mesh coordinate sensitivities and design values
 *        similar to the functions in SU2_DOT
 * \author T.Dick
 * \version 7.0.5 "Blackbird"
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
 * 31.01.2020 T. Dick
 */
void GetParameterizationJacobianReverse(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian) {

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
        visitedPoints[iPoint] = 1;
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


/* development comments:
 * Call the AD routine with the help of tape forward evalution to get the Jacobian
 * 13.05.2020 T. Dick
 */
void GetParameterizationJacobianForward(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian) {

  unsigned short nDim, nMarker, nDV, nDV_Value, nDV_Total, nPoint, nVertex;
  unsigned short iDV, iDV_Value, iDV_index, iPoint, iMarker, iVertex, iDim;
  su2double* VarCoord;
  su2double** DV_Value;

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nPoint  = geometry->GetnPoint();
  nDV     = config->GetnDV();
  nDV_Total = config->GetnDV_Total();

  /*--- Start recording of operations ---*/

  AD::Reset();

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  DV_Value = config->GetDV_Pointer();
  for (iDV = 0; iDV < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization of su2double with 0.0 resets the index ---*/
      DV_Value[iDV][iDV_Value] = 0.0;
      AD::RegisterInput(DV_Value[iDV][iDV_Value]);

    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Register Outputs --- */
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim=0; iDim<nDim; iDim++) {
          AD::RegisterOutput(VarCoord[iDim]);
        }
      }
    }
  }


  /*--- Stop the recording --- */
  AD::StopRecording();

  /*--- Loop over all design variables and calculate the derivative with respect to that variable --*/
  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++) {

      SU2_TYPE::SetDerivative(DV_Value[iDV][iDV_Value], 1.0);

      AD::ComputeAdjointForward();

      /*--- Extract sensitivities ---*/
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
        if (config->GetMarker_All_DV(iMarker) == YES) {
          nVertex = geometry->nVertex[iMarker];
          for (iVertex = 0; iVertex <nVertex; iVertex++) {
            VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
            iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
            for (iDim=0; iDim<nDim; iDim++) {
              Jacobian[iDV_index*nPoint*nDim+iPoint*nDim+iDim] = AD::GetDerivativeValue(VarCoord[iDim]);
            }
          }
        }
      }

      iDV_index++;
      AD::ClearAdjoints();

    }
  }

  AD::Reset();
}


/* development comments:
 * Call the AD routine for each point to get the whole Jacobian
 * makes use of preaccumulation
 * 31.01.2020 T. Dick
 */
void GetParameterizationJacobianPreaccumulation(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, su2double *Jacobian) {

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
        visitedPoints[iPoint] = 1;
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

  /*--- Mark designe variables as preaccumulation input ---*/
  AD::StartPreacc();
  su2double** DV_Values = config->GetDV_Pointer();
  for (iDV = 0; iDV < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      AD::SetPreaccIn(DV_Values[iDV][iDV_Value]);
    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Mark the preaccumulation outputs ---*/

  for (iPoint=0; iPoint<nPoint; iPoint++) {

    /*--- Go over all visited points and dimensions ---*/
    if (visitedPoints[iPoint]==1) {
      for (subDim=0; subDim<nDim; subDim++) {
        for (iMarker = 0; iMarker < nMarker; iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            nVertex = geometry->nVertex[iMarker];
            for (iVertex = 0; iVertex <nVertex; iVertex++) {
              VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
              if (iPoint == geometry->vertex[iMarker][iVertex]->GetNode()) {
                AD::SetPreaccOut(VarCoord[subDim]);
              }
            }
          }
        }
      }
    }
  }

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


/* development comments:
 * Record a tape containing one evaluation of the parameterization.
 * \date 19.10.2020
 * \name T. Dick
 */
void RecordParameterizationJacobian(CGeometry *geometry, CConfig *config, CSurfaceMovement *surface_movement, CSysVector<su2double>& registeredCoord) {

  unsigned short nDim, nMarker, nDV, nDV_Value, nDV_Total, nPoint, nVertex;
  unsigned short iDV, iDV_Value, iMarker, iPoint, iVertex, iDim;
  su2double* VarCoord;
  su2double** DV_Value;

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nPoint  = geometry->GetnPoint();
  nDV     = config->GetnDV();
  nDV_Total = config->GetnDV_Total();

  /*--- Start recording of operations ---*/

  AD::Reset();

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  DV_Value = config->GetDV_Pointer();
  for (iDV = 0; iDV < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization of su2double with 0.0 resets the index ---*/
      DV_Value[iDV][iDV_Value] = 0.0;
      AD::RegisterInput(DV_Value[iDV][iDV_Value]);

    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Register Outputs --- */
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();
        for (iDim=0; iDim<nDim; iDim++) {
          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          registeredCoord(iPoint,iDim) = VarCoord[iDim];
        }
      }
    }
  }
  for (iPoint = 0; iPoint<nPoint; iPoint++) {
    for (iDim=0; iDim<nDim; iDim++) {
      AD::RegisterOutput(registeredCoord(iPoint,iDim));
    }
  }

  /*--- Stop the recording --- */
  AD::StopRecording();

}


/* development comments:
 * Evaluate the derivatives of the parameterization in tape forward mode.
 * \date 19.10.2020
 * \name T. Dick
 */
void ProjectDVtoMesh(CGeometry *geometry, CConfig *config, std::vector<su2double>& seeding, CSysVector<su2mixedfloat>& result, CSysVector<su2double>& registeredCoord) {

  unsigned short nDim, nMarker, nDV, nDV_Value, nVertex;
  unsigned short iDV, iDV_Value, iDV_index, iMarker, iVertex, iPoint, iDim;
  su2double** DV_Value = config->GetDV_Pointer();

  /*--- get information from config ---*/
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Seeding for the DV. --*/
  iDV_index = 0;
  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);
    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      SU2_TYPE::SetDerivative(DV_Value[iDV][iDV_Value], SU2_TYPE::GetValue(seeding[iDV_index]));
      iDV_index++;
    }
  }

  AD::ComputeAdjointForward();

  /*--- Extract sensitivities ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim=0; iDim<nDim; iDim++) {
          result(iPoint,iDim) = AD::GetDerivativeValue(registeredCoord(iPoint,iDim));
        }
      }
    }
  }

  AD::ClearAdjoints();

}

/* development comments:
 * Evaluate the derivatives of the parameterization in reverse mode.
 * \date 19.10.2020
 * \name T. Dick
 */
void ProjectMeshToDV(CGeometry *geometry, CConfig *config, CSysVector<su2mixedfloat>& sensitivity, std::vector<su2double>& output, CSysVector<su2double>& registeredCoord) {

  /*--- Part 1: adjoint volumetric mesh deformation, if needed ---*/
  if (true) {

  }

  /*--- Part 2: adjoint surface deformation ---*/

  unsigned short nDim, nMarker, nDV, nDV_Value, nVertex;
  unsigned short iDV, iDV_Value, iDV_index, iPoint, iDim, iMarker, iVertex;
  su2double** DV_Value = config->GetDV_Pointer();
  double my_Gradient, localGradient;

  // get some numbers from config
  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  /*--- Set the seeding appropriately ---*/
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        for (iDim = 0; iDim < nDim; iDim++){
          SU2_TYPE::SetDerivative(registeredCoord(iPoint,iDim), SU2_TYPE::GetValue(sensitivity(iPoint,iDim)));
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
      my_Gradient = AD::GetDerivativeValue(DV_Value[iDV][iDV_Value]);
      #ifdef HAVE_MPI
        SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #else
        localGradient = my_Gradient;
      #endif
      output[iDV_index] = localGradient;
      iDV_index++;
    }
  }

  AD::ClearAdjoints();

}



