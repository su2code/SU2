/*!
 * \file CVolumeElementFEM_Base.cpp
 * \brief Implementations of the member functions of CVolumeElementFEM_Base.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#include "../../../include/geometry/fem_grid/CVolumeElementFEM_Base.hpp"
#include "../../../include/geometry/fem_grid/CPointFEM.hpp"
#include "../../../include/geometry/primal_grid/CPrimalGridFEM.hpp"

/*---------------------------------------------------------------------*/
/*---      Public member functions of CVolumeElementFEM_Base.       ---*/
/*---------------------------------------------------------------------*/

void CVolumeElementFEM_Base::DerMetricTermsIntegrationPoints(const unsigned short nDim) {

  /*--- Allocate the memory for the metric terms to compute the second derivatives
        in the integration points. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();
  const unsigned short n2ndDer = nDim + nDim*(nDim-1)/2;

  metricTerms2ndDerInt.resize(n2ndDer);
  for(unsigned short k=0; k<n2ndDer; ++k)
    metricTerms2ndDerInt[k].resize(nIntPad, nDim);

  /*--- Compute the metric terms for the 2nd derivatives in
        the volume integration points. ---*/
  standardElemGrid->MetricTerms2ndDerVolumeIntPoints(coorGridDOFs, metricTermsInt,
                                                     JacobiansInt, metricTerms2ndDerInt);
}

void CVolumeElementFEM_Base::GetCornerPointsAllFaces(unsigned short &numFaces,
                                                     unsigned short nPointsPerFace[],
                                                     unsigned long  faceConn[6][4]) const {

  /*--- Get the corner connectivities of the faces, local to the element. ---*/
  CPrimalGridFEM::GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid,
                                               numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values, local to the element, of faceConn to
        global values, i.e. the numbering used in this partition. ---*/
  for(unsigned short i=0; i<numFaces; ++i) {
    for(unsigned short j=0; j<nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = nodeIDsGrid[nn];
    }
  }
}

void CVolumeElementFEM_Base::InitGridVelocities(const unsigned short nDim) {

  /*--- Determine the padded number of integration points and solution DOFs. ---*/
  const unsigned short nIntPad     = standardElemGrid->GetNIntegrationPad();
  const unsigned short nDOFsSolPad = standardElemGrid->GetNSolDOFsPad();

  /*--- Allocate the memory for the grid velocities and initialize
        them to zero. ---*/
  gridVelocitiesInt.resize(nIntPad, nDim);
  gridVelocitiesInt.setConstant(0.0);

  gridVelocitiesSolDOFs.resize(nDOFsSolPad, nDim);
  gridVelocitiesSolDOFs.setConstant(0.0);
}

bool CVolumeElementFEM_Base::MetricTermsIntegrationPoints(const bool           LGLDistribution,
                                                          const unsigned short nDim) {

  /*--- Allocate the memory for the metric terms. ---*/
  const unsigned short nIntPad = standardElemGrid->GetNIntegrationPad();
  JacobiansInt.resize(nIntPad);

  metricTermsInt.resize(nDim);
  for(unsigned short k=0; k<nDim; ++k)
    metricTermsInt[k].resize(nIntPad, nDim);

  /*--- Compute the metric terms in the volume integration points. ---*/
  standardElemGrid->MetricTermsVolumeIntPoints(LGLDistribution, coorGridDOFs,
                                               metricTermsInt, JacobiansInt);

  /*--- Determine the volume of the element and check for negative Jacobians
        in the integration points. ---*/
  const unsigned short nInt = standardElemGrid->GetNIntegration();
  const passivedouble *wInt = standardElemGrid->GetIntegrationWeights();
  bool elemIsGood = true;
  volume = 0.0;

  for(unsigned short i=0; i<nInt; ++i) {
    volume += wInt[i]*JacobiansInt(i);
    if(JacobiansInt(i) <= 0.0) elemIsGood = false;
  }

  /*--- Make sure that the Jacobian for the padded data is
        not zero, to avoid problems later on. ---*/
  for(unsigned short i=nInt; i<nIntPad; ++i) JacobiansInt(i) = 1.0;

  /*--- Return the value of elemIsGood. ---*/
  return elemIsGood;
}

bool CVolumeElementFEM_Base::MetricTermsSolDOFs(const unsigned short nDim) {

  /*--- Allocate the memory for the metric terms of the solution DOFs. ---*/
  const unsigned short nDOFsSol    = standardElemGrid->GetNSolDOFs();
  const unsigned short nDOFsSolPad = standardElemGrid->GetNSolDOFsPad();

  JacobiansSolDOFs.resize(nDOFsSolPad);

  metricTermsSolDOFs.resize(nDim);
  for(unsigned short k=0; k<nDim; ++k)
    metricTermsSolDOFs[k].resize(nDOFsSolPad, nDim);

  /*--- Compute the metric terms in the solution DOFs. ---*/
  standardElemGrid->MetricTermsSolDOFs(coorGridDOFs, metricTermsSolDOFs,
                                       JacobiansSolDOFs);

  /*--- Check for negative Jacobians in the solution DOFs. ---*/
  bool elemIsGood = true;
  for(unsigned short i=0; i<nDOFsSol; ++i)
    if(JacobiansSolDOFs(i) <= 0.0) elemIsGood = false;

  /*--- Make sure that the Jacobian for the padded data is
        not zero, to avoid problems later on. ---*/
  for(unsigned short i=nDOFsSol; i<nDOFsSolPad; ++i) JacobiansSolDOFs(i) = 1.0;

  /*--- Return the value of elemIsGood. ---*/
  return elemIsGood;
}

void CVolumeElementFEM_Base::SetCoorGridDOFs(const unsigned short    nDim,
                                             const vector<CPointFEM> &meshPoints) {

  /*--- Allocate the memory for the coordinates of the grid DOFs. ---*/
  coorGridDOFs.resize(nodeIDsGrid.size(), nDim);

  /*--- Loop over the dimensions and number of grid DOFs. ---*/
  for(unsigned short k=0; k<nDim; ++k)
    for(unsigned long i=0; i<nodeIDsGrid.size(); ++i)
      coorGridDOFs(i,k) = meshPoints[nodeIDsGrid[i]].coor[k];
}
