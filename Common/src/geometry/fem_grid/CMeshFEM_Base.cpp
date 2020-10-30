/*!
 * \file CMeshFEM_Base.cpp
 * \brief Implementations of the member functions of CMeshFEM_Base.
 * \author E. van der Weide
 * \version 7.0.7 "Blackbird"
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

#include "../../../include/geometry/fem_grid/CMeshFEM_Base.hpp"
#include "../../../include/fem/CFEMStandardHexGrid.hpp"
#include "../../../include/fem/CFEMStandardPrismGrid.hpp"
#include "../../../include/fem/CFEMStandardPyraGrid.hpp"
#include "../../../include/fem/CFEMStandardQuadGrid.hpp"
#include "../../../include/fem/CFEMStandardTetGrid.hpp"
#include "../../../include/fem/CFEMStandardTriGrid.hpp"

/*---------------------------------------------------------------------*/
/*---          Public member functions of CMeshFEM_Base.            ---*/
/*---------------------------------------------------------------------*/

CMeshFEM_Base::CMeshFEM_Base(CGeometry *geometry, CConfig *config) {

  /*--- The new FEM mesh class has the same problem dimension/zone. ---*/
  nDim         = geometry->GetnDim();
  nZone        = geometry->GetnZone();
  Global_nElem = geometry->GetGlobal_nElem();
}

CMeshFEM_Base::~CMeshFEM_Base(void) {

  for(unsigned long i=0; i<standardVolumeElementsGrid.size(); ++i) {
    if( standardVolumeElementsGrid[i] ) delete standardVolumeElementsGrid[i];
    standardVolumeElementsGrid[i] = nullptr;
  }
}

/*---------------------------------------------------------------------*/
/*---        Protected member functions of CMeshFEM_Base.           ---*/
/*---------------------------------------------------------------------*/

void CMeshFEM_Base::CreateStandardVolumeElementsGrid(const vector<CUnsignedShort3T> &elemTypes) {

  /*--- Allocate the memory for the pointers. ---*/
  standardVolumeElementsGrid.resize(elemTypes.size(), nullptr);

  /*--- Loop over the different element types for the grid. ---*/
  for(unsigned long i=0; i<elemTypes.size(); ++i) {

    /*--- Abbreviate the element type, polynomial degree and polynomial order that
          must be integrated exactly for readability. ---*/
    const unsigned short VTK_Type   = elemTypes[i].short0;
    const unsigned short nPoly      = elemTypes[i].short1;
    const unsigned short orderExact = elemTypes[i].short2;

    /*--- Determine the element type and allocate the appropriate object. ---*/
    switch( VTK_Type ) {
      case TRIANGLE:
        standardVolumeElementsGrid[i] = new CFEMStandardTriGrid(nPoly, orderExact, false);
        break;
      case QUADRILATERAL:
        standardVolumeElementsGrid[i] = new CFEMStandardQuadGrid(nPoly, orderExact, false);
        break;
      case TETRAHEDRON:
        standardVolumeElementsGrid[i] = new CFEMStandardTetGrid(nPoly, orderExact);
        break;
      case PYRAMID:
        standardVolumeElementsGrid[i] = new CFEMStandardPyraGrid(nPoly, orderExact);
        break;
      case PRISM:
        standardVolumeElementsGrid[i] = new CFEMStandardPrismGrid(nPoly, orderExact);
        break;
      case HEXAHEDRON:
        standardVolumeElementsGrid[i] = new CFEMStandardHexGrid(nPoly, orderExact);
        break;
      default:  /*--- To avoid a compiler warning. ---*/
        SU2_MPI::Error(string("Unknown volume element. This should not happen"),
                         CURRENT_FUNCTION);
    }
  }
}
