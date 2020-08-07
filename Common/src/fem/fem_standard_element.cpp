/*!
 * \file fem_standard_element.cpp
 * \brief Functions for the FEM standard elements.
 * \author E. van der Weide
 * \version 7.0.6 "Blackbird"
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

#include "../../include/fem/fem_standard_element.hpp"

/*----------------------------------------------------------------------------------*/
/*          Public member functions of CFEMStandardElementBase.                     */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElementBase::GetNDOFsStatic(unsigned short VTK_Type,
                                                       unsigned short nPoly) {
  unsigned short nDOFsEdge = nPoly + 1;
  unsigned short nDOFs;

  switch(VTK_Type) {

    case LINE:
      nDOFs = nDOFsEdge;
      break;

    case TRIANGLE:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case QUADRILATERAL:
      nDOFs = nDOFsEdge*nDOFsEdge;
      break;

    case TETRAHEDRON:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(nDOFsEdge+2)/6;
      break;

    case HEXAHEDRON:
      nDOFs = nDOFsEdge*nDOFsEdge*nDOFsEdge;
      break;

    case PRISM:
      nDOFs = nDOFsEdge*nDOFsEdge*(nDOFsEdge+1)/2;
      break;

    case PYRAMID:
      nDOFs = nDOFsEdge*(nDOFsEdge+1)*(2*nDOFsEdge+1)/6;
      break;

    default:
      nDOFs = 0;  // Indicates an invalid element.
      break;
  }

  return nDOFs;
}


/*----------------------------------------------------------------------------------*/
/*         Protected member functions of CFEMStandardElementBase.                   */
/*----------------------------------------------------------------------------------*/

void CFEMStandardElementBase::SubConnForPlottingLine(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}

void CFEMStandardElementBase::SubConnForPlottingQuadrilateral(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}

void CFEMStandardElementBase::SubConnForPlottingTriangle(
                                         const unsigned short   nPoly,
                                         vector<unsigned short> &subConn) {

}


/*----------------------------------------------------------------------------------*/
/*           Private member functions of CFEMStandardElement.                       */
/*----------------------------------------------------------------------------------*/

unsigned short CFEMStandardElement::GetNDOFsPerSubElem(unsigned short val_VTK_Type) const {

  return 0;
}
