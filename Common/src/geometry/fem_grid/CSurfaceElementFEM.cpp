/*!
 * \file CSurfaceElementFEM.cpp
 * \brief Implementations of the member functions of CSurfaceElementFEM.
 * \author E. van der Weide
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/geometry/fem_grid/CSurfaceElementFEM.hpp"
#include "../../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../../include/geometry/fem_grid/CVolumeElementFEM_DG.hpp"

/*---------------------------------------------------------------------*/
/*---        Public member functions of CSurfaceElementFEM.         ---*/
/*---------------------------------------------------------------------*/

void CSurfaceElementFEM::AllocateResiduals(CConfig        *config,
                                           unsigned short nVar) {

  /*--- Determine the number of padded solution DOFs of the adjacent
        element and allocate the memory for the residual. ---*/
  const unsigned short nDOFsPadElem = standardElemFlow->GetNSolDOFsPad();
  resDOFsElem.resize(nDOFsPadElem,nVar);
}

ColMajorMatrix<su2double> &CSurfaceElementFEM::ComputeSolSide0IntPoints(
                                                CVolumeElementFEM_DG *volElem) {

  /*--- Perform the gemm call to compute the solution in the
        integration points and return the matrix. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->SolIntPoints(volElem[volElemID].solDOFsWork,
                                 standardElemFlow->workSolInt[thread]);
  return standardElemFlow->workSolInt[thread];
}

void CSurfaceElementFEM::ComputeWallDistance(CADTElemClass        *WallADT,
                                             const unsigned short nDim) {

  /*--- Determine the number of integration points. ---*/
  const unsigned short nInt = standardElemGrid->GetNIntegration();

  /*--- Loop over the integration points and compute the
        distance to the wall. ---*/
  for(unsigned short i=0; i<nInt; ++i) {

    su2double coor[3] = {0.0};
    for(unsigned short k=0; k<nDim; ++k)
      coor[k] = coorIntegrationPoints(i,k);

    unsigned short markerID;
    unsigned long  elemID;
    int            rankID;
    su2double      dist;
    WallADT->DetermineNearestElement(coor, dist, markerID, elemID, rankID);

    wallDistance(i) = dist;
  }
}

void CSurfaceElementFEM::GetCornerPointsFace(unsigned short &nPointsPerFace,
                                             unsigned long  faceConn[]) {

  /*--- Get the corner connectivities of the face, local to the element. ---*/
  CPrimalGridBoundFEM::GetLocalCornerPointsFace(VTK_Type, nPolyGrid, nDOFsGrid,
                                                nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/
  for(unsigned short j=0; j<nPointsPerFace; ++j) {
    unsigned long nn = faceConn[j];
    faceConn[j] = nodeIDsGrid[nn];
  }
}

void CSurfaceElementFEM::InitGridVelocities(const unsigned short nDim) {

  /*--- Determine the padded number of integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();

  /*--- Allocate the memory for the grid velocities and initialize
        them to zero. ---*/
  gridVelocities.resize(nIntPad, nDim);
  gridVelocities.setConstant(0.0);
}

void CSurfaceElementFEM::MetricTermsIntegrationPoints(const unsigned short         nDim,
                                                      vector<CVolumeElementFEM_DG> &volElem) {

  /*--- Determine the padded number of integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();

  /*--- Allocate the memory for the coordinates of the integration
        points and determine them. The first argument in the function
        CoorIntPoints is a dummy to be consistent with the declaration. ---*/
  coorIntegrationPoints.resize(nIntPad, nDim);

  standardElemGrid->CoorIntPoints(true, volElem[volElemID].coorGridDOFs, coorIntegrationPoints);

  /*--- Set the coordinates of the padded points to avoid problems. ---*/
  for(unsigned short j=0; j<nDim; ++j) 
    for(unsigned short i=standardElemGrid->GetNIntegration(); i<nIntPad; ++i)
      coorIntegrationPoints(i,j) = coorIntegrationPoints(0,j);

  /*--- Allocate the memory for the metric terms of the boundary face. ---*/
  JacobiansFace.resize(nIntPad);
  metricNormalsFace.resize(nIntPad, nDim);

  metricCoorDerivFace.resize(nDim);
  for(unsigned short k=0; k<nDim; ++k) {
    metricCoorDerivFace[k].resize(nIntPad, nDim);
    metricCoorDerivFace[k].setConstant(0.0);      // To avoid uninitialized data
  }                                               // for the padded points.

  /*--- Compute the metric terms in the surface integration points. ---*/
  standardElemGrid->MetricTermsSurfaceIntPoints(volElem[volElemID].coorGridDOFs,
                                                JacobiansFace, metricNormalsFace,
                                                metricCoorDerivFace);
}

void CSurfaceElementFEM::ResidualBasisFunctions(ColMajorMatrix<su2double> &scalarDataInt) {
  standardElemFlow->ResidualBasisFunctions(scalarDataInt, resDOFsElem);
}

void CSurfaceElementFEM::SetWallDistance(su2double val) {

  /*--- Determine the padded number of integration points, allocate the
        memory for the wall distances and set them. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();
  wallDistance.resize(nIntPad);
  wallDistance.setConstant(val);
}
