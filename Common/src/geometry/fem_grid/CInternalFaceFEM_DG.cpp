/*!
 * \file CInternalFaceFEM_DG.cpp
 * \brief Implementations of the member functions of CInternalFaceFEM_DG.
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

#include "../../../include/geometry/fem_grid/CInternalFaceFEM_DG.hpp"
#include "../../../include/geometry/fem_grid/CVolumeElementFEM_DG.hpp"

/*---------------------------------------------------------------------*/
/*---        Public member functions of CInternalFaceFEM_DG.        ---*/
/*---------------------------------------------------------------------*/

void CInternalFaceFEM_DG::AllocateResiduals(CConfig        *config,
                                            unsigned short nVar) {

  /*--- Determine the number of padded solution DOFs of the adjacent
        elements and allocate the memory for the residual. ---*/
  const unsigned short nDOFsPadSide0 = standardElemFlow->elem0->GetNSolDOFsPad();
  const unsigned short nDOFsPadSide1 = standardElemFlow->elem1->GetNSolDOFsPad();

  resDOFsSide0.resize(nDOFsPadSide0,nVar);
  resDOFsSide1.resize(nDOFsPadSide1,nVar);
}

vector<ColMajorMatrix<su2double> > &CInternalFaceFEM_DG::ComputeGradSolSide0IntPoints(
                                                    CVolumeElementFEM_DG *volElem) {

  /*--- Perform the gemm calls to compute the gradients in the integration
        points and return the vector of matrices. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->elem0->GradSolIntPoints(volElem[elemID0].solDOFsWork, 
                                            standardElemFlow->elem0->workGradSolInt[thread]);
  return standardElemFlow->elem0->workGradSolInt[thread];
}

vector<ColMajorMatrix<su2double> > &CInternalFaceFEM_DG::ComputeGradSolSide1IntPoints(
                                                    CVolumeElementFEM_DG *volElem) {

  /*--- Perform the gemm calls to compute the gradients in the integration
        points and return the vector of matrices. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->elem1->GradSolIntPoints(volElem[elemID1].solDOFsWork,
                                            standardElemFlow->elem1->workGradSolInt[thread]);
  return standardElemFlow->elem1->workGradSolInt[thread];
}

ColMajorMatrix<su2double> &CInternalFaceFEM_DG::ComputeSolSide0IntPoints(
                                                CVolumeElementFEM_DG *volElem) {

  /*--- Perform the gemm call to compute the solution in the
        integration points and return the matrix. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->elem0->SolIntPoints(volElem[elemID0].solDOFsWork,
                                        standardElemFlow->elem0->workSolInt[thread]);
  return standardElemFlow->elem0->workSolInt[thread];
}

ColMajorMatrix<su2double> &CInternalFaceFEM_DG::ComputeSolSide1IntPoints(
                                                CVolumeElementFEM_DG *volElem) {

  /*--- Perform the gemm call to compute the solution in the
        integration points and return the matrix. ---*/
  const int thread = omp_get_thread_num();
  standardElemFlow->elem1->SolIntPoints(volElem[elemID1].solDOFsWork,
                                        standardElemFlow->elem1->workSolInt[thread]);
  return standardElemFlow->elem1->workSolInt[thread];
}

void CInternalFaceFEM_DG::ComputeWallDistance(CADTElemClass        *WallADT,
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

void CInternalFaceFEM_DG::InitGridVelocities(const unsigned short nDim) {

  /*--- Determine the padded number of integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();

  /*--- Allocate the memory for the grid velocities and initialize
        them to zero. ---*/
  gridVelocities.resize(nIntPad, nDim);
  gridVelocities.setConstant(0.0);
}

void CInternalFaceFEM_DG::MetricTermsIntegrationPoints(const unsigned short         nDim,
                                                       vector<CVolumeElementFEM_DG> &volElem) {

  /*--- Determine the padded number of integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();

  /*--- Allocate the memory for the coordinates of the integration
        points and determine them. Note that if the polynomial order
        of the grid of both elements differ, the element with the
        highest polynomial order is always on side 0. Hence, it
        suffices to determine the coordinates from the data of the
        element on side 0 of the face.  ---*/
  coorIntegrationPoints.resize(nIntPad, nDim);

  standardElemGrid->CoorIntPoints(volElem[elemID0].coorGridDOFs, coorIntegrationPoints);

#ifndef NDEBUG
  /*--- DEBUG: Compute the coordinates of the integration points from the data
               of th element on side 1 and compare the results. ---*/
  ColMajorMatrix<su2double> coorSide1(nIntPad, nDim);
  standardElemGrid->CoorIntPointsFromSide1(volElem[elemID1].coorGridDOFs, coorSide1);

  for(unsigned short i=0; i<standardElemGrid->GetNIntegration(); ++i) {
    su2double dist2 = 0.0;
    for(unsigned short j=0; j<nDim; ++j) {
      const su2double ds = coorIntegrationPoints(i,j) - coorSide1(i,j);
      dist2 += ds*ds;
    }
    const su2double dist = sqrt(dist2);
    assert(dist <= 1.e-6);
  }

  /*--- END DEBUG. ---*/
#endif

  /*--- Set the coordinates of the padded points to avoid problems. ---*/
  for(unsigned short j=0; j<nDim; ++j)
    for(unsigned short i=standardElemGrid->GetNIntegration(); i<nIntPad; ++i)
      coorIntegrationPoints(i,j) = coorIntegrationPoints(0,j);

  /*--- Allocate the memory for the metric terms of the internal face. ---*/
  JacobiansFace.resize(nIntPad);
  metricNormalsFace.resize(nIntPad, nDim);

  metricCoorDerivFace0.resize(nDim);
  metricCoorDerivFace1.resize(nDim);

  for(unsigned short k=0; k<nDim; ++k) {
    metricCoorDerivFace0[k].resize(nIntPad, nDim);
    metricCoorDerivFace1[k].resize(nIntPad, nDim);

    metricCoorDerivFace0[k].setConstant(0.0); // To avoid uninitialized data for
    metricCoorDerivFace1[k].setConstant(0.0); // the padded points.
  }

  /*--- Compute the metric terms in the surface integration points. ---*/
  standardElemGrid->MetricTermsSurfaceIntPoints(volElem[elemID0].coorGridDOFs,
                                                volElem[elemID1].coorGridDOFs,
                                                JacobiansFace, metricNormalsFace,
                                                metricCoorDerivFace0,
                                                metricCoorDerivFace1);
}

void CInternalFaceFEM_DG::SetWallDistance(su2double val) {

  /*--- Determine the padded number of integration points, allocate the
        memory for the wall distances and set them. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();
  wallDistance.resize(nIntPad);
  wallDistance.setConstant(val);
}

void CInternalFaceFEM_DG::ResidualBasisFunctionsSide0(ColMajorMatrix<su2double> &scalarDataInt) {
  standardElemFlow->elem0->ResidualBasisFunctions(scalarDataInt, resDOFsSide0);
}

void CInternalFaceFEM_DG::ResidualBasisFunctionsSide1(ColMajorMatrix<su2double> &scalarDataInt) {
  standardElemFlow->elem1->ResidualBasisFunctions(scalarDataInt, resDOFsSide1);
}