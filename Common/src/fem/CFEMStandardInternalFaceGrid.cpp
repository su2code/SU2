/*!
 * \file CFEMStandardInternalFaceGrid.cpp
 * \brief Functions for the class CFEMStandardInternalFaceGrid.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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

#include "../../include/fem/CFEMStandardInternalFaceGrid.hpp"

/*----------------------------------------------------------------------------------*/
/*            Public member functions of CFEMStandardInternalFaceGrid.              */
/*----------------------------------------------------------------------------------*/

CFEMStandardInternalFaceGrid::CFEMStandardInternalFaceGrid(CFEMStandardElementBase *val_elem0,
		                                           CFEMStandardElementBase *val_elem1) {
  elem0 = val_elem0;
  elem1 = val_elem1;
}

void CFEMStandardInternalFaceGrid::CoorIntPoints(ColMajorMatrix<su2double> &coorGridDOFsVol,
                                                 ColMajorMatrix<su2double> &coorIntPointsFace) {

  /*--- Per definition the polynomial degree of the element on side 0 is higher or equal
        to the polynomial degree of the element on side 1. Hence the element on side 0
        is used to compute the coordinates of the integration points. The first argument
        is a dummy to be consistent with other classes. ---*/
  elem0->CoorIntPoints(true, coorGridDOFsVol, coorIntPointsFace);
}

void CFEMStandardInternalFaceGrid::CoorIntPointsFromSide1(ColMajorMatrix<su2double> &coorGridDOFsVol,
                                                          ColMajorMatrix<su2double> &coorIntPointsFace) {

  /*--- Use the data of the element on side 1 to compute the coordinates
        of the integration points. This data is only used for debugging.
        The first argument is a dummy to be consistent with other classes. ---*/
  elem1->CoorIntPoints(true, coorGridDOFsVol, coorIntPointsFace);
}

void CFEMStandardInternalFaceGrid::MetricTermsSurfaceIntPoints(ColMajorMatrix<su2double>          &matCoorElem0,
                                                               ColMajorMatrix<su2double>          &matCoorElem1,
                                                               su2activevector                    &JacobiansFace,
                                                               ColMajorMatrix<su2double>          &normalsFace,
                                                               vector<ColMajorMatrix<su2double> > &metricTermsSide0,
                                                               vector<ColMajorMatrix<su2double> > &metricTermsSide1) {

  /*--- In order to compute drdx, drdy, etc. on side 1, two calls must be made to the functions
        of elem1, because the normals are not computed here. Note that the first argument to
        DerivativesCoorIntPoints is a dummy and JacobiansFace is used as temporary storage. ---*/
  elem1->DerivativesCoorIntPoints(true, matCoorElem1, metricTermsSide1);
  elem1->MetricTermsVolume(metricTermsSide1, JacobiansFace);

  /*--- Determine the inverse of the Jacobians. First set the added values
        to 1 to avoid divisions by zero. ---*/
  const unsigned short nItems = JacobiansFace.rows();
  for(unsigned short i=GetNIntegration(); i<nItems; ++i) JacobiansFace[i] = 1.0;

  SU2_OMP_SIMD_IF_NOT_AD
  for(unsigned short i=0; i<nItems; ++i)
    JacobiansFace[i] = 1.0/JacobiansFace[i];
  
  /*--- Compute the final metric terms on side 1. ---*/
  for(unsigned short k=0; k<metricTermsSide1.size(); ++k) {
    for(unsigned short j=0; j<metricTermsSide1[k].cols(); ++j) {
      SU2_OMP_SIMD_IF_NOT_AD
      for(unsigned short i=0; i<nItems; ++i)
        metricTermsSide1[k](i,j) *= JacobiansFace[i];
    }
  }

  /*--- Call MetricTermsSurfaceIntPoints of elem0 to compute the full metric terms
        in the integration points on side 0. ---*/
  elem0->MetricTermsSurfaceIntPoints(matCoorElem0, JacobiansFace, normalsFace, metricTermsSide0);
}
