/*!
 * \file CMeshFEM_Base.cpp
 * \brief Implementations of the member functions of CMeshFEM_Base.
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

  for(unsigned long i=0; i<standardSurfaceElementsGrid.size(); ++i) {
    if( standardSurfaceElementsGrid[i] ) delete standardSurfaceElementsGrid[i];
    standardSurfaceElementsGrid[i] = nullptr;
  }

  for(unsigned long i=0; i<gemmTypesFaces.size(); ++i) {
    if( gemmTypesFaces[i] ) delete gemmTypesFaces[i];
    gemmTypesFaces[i] = nullptr;
  }
}

/*---------------------------------------------------------------------*/
/*---        Protected member functions of CMeshFEM_Base.           ---*/
/*---------------------------------------------------------------------*/

void CMeshFEM_Base::CreateStandardVolumeElementsGrid(const vector<CUnsignedShort4T> &elemTypes,
                                                     const unsigned short           locGridDOFs) {

  /*--- Allocate the memory for the pointers. ---*/
  standardVolumeElementsGrid.resize(elemTypes.size(), nullptr);

  /*--- Loop over the different element types for the grid. ---*/
  for(unsigned long i=0; i<elemTypes.size(); ++i) {

    /*--- Abbreviate the element type, polynomial degree and polynomial order that
          must be integrated exactly for readability. ---*/
    const unsigned short VTK_Type   = elemTypes[i].short0;
    const unsigned short nPolyGrid  = elemTypes[i].short1;
    const unsigned short nPolySol   = elemTypes[i].short2;
    const unsigned short orderExact = elemTypes[i].short3;

    /*--- Determine the element type and allocate the appropriate object. ---*/
    switch( VTK_Type ) {
      case TRIANGLE:
        standardVolumeElementsGrid[i] = new CFEMStandardTriGrid(nPolyGrid,  nPolySol,
                                                                orderExact, locGridDOFs);
        break;
      case QUADRILATERAL:
        standardVolumeElementsGrid[i] = new CFEMStandardQuadGrid(nPolyGrid,  nPolySol,
                                                                 orderExact, locGridDOFs);
        break;
      case TETRAHEDRON:
        standardVolumeElementsGrid[i] = new CFEMStandardTetGrid(nPolyGrid,  nPolySol,
                                                                orderExact, locGridDOFs); 
        break;
      case PYRAMID:
        standardVolumeElementsGrid[i] = new CFEMStandardPyraGrid(nPolyGrid,  nPolySol,
                                                                 orderExact, locGridDOFs); 
        break;
      case PRISM:
        standardVolumeElementsGrid[i] = new CFEMStandardPrismGrid(nPolyGrid,  nPolySol,
                                                                  orderExact, locGridDOFs); 
        break;
      case HEXAHEDRON:
        standardVolumeElementsGrid[i] = new CFEMStandardHexGrid(nPolyGrid,  nPolySol,
                                                                orderExact, locGridDOFs); 
        break;
      default:  /*--- To avoid a compiler warning. ---*/
        SU2_MPI::Error(string("Unknown volume element. This should not happen"),
                       CURRENT_FUNCTION);
    }
  }
}

std::unique_ptr<CADTElemClass> CMeshFEM_Base::ComputeViscousWallADT(const CConfig *config) const {

  /*--- Initialize an array for the mesh points, which eventually contains the
        mapping from the local nodes to the number used in the connectivity of the
        local boundary faces. However, in a first pass it is an indicator whether
        or not a mesh point is on a local wall boundary. ---*/
  vector<unsigned long> meshToSurface(meshPoints.size(), 0);

  /*--- Define the vectors for the connectivity of the local linear subelements,
        the element ID's, the element type and marker ID's. ---*/
  vector<unsigned long> surfaceConn;
  vector<unsigned long> elemIDs;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned short> markerIDs;

  /*--- Loop over the boundary markers. ---*/
  for(unsigned short iMarker=0; iMarker<boundaries.size(); ++iMarker) {
    if( !boundaries[iMarker].periodicBoundary ) {

      /*--- Check for a viscous wall. ---*/
      if( (config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL) ) {

        /*--- Loop over the surface elements of this marker. ---*/
        const vector<CSurfaceElementFEM> &surfElem = boundaries[iMarker].surfElem;
        for(unsigned long i=0; i<surfElem.size(); ++i) {

          /* Set the flag of the mesh points on this surface to true. */
          for(unsigned short j=0; j<surfElem[i].nDOFsGrid; ++j)
            meshToSurface[surfElem[i].nodeIDsGrid[j]] = 1;

          /*--- Determine the necessary data from the corresponding standard
                element for the grid. Note that the linear sub-elements of a
                a surface element are of only one type, so no need to consider
                SubType2. ---*/
          const unsigned short VTK_Type      = surfElem[i].standardElemGrid->GetVTK_SubType1();
          const unsigned short nSubFaces     = surfElem[i].standardElemGrid->GetNSubElemsType1();
          const unsigned short nDOFsPerFace  = surfElem[i].standardElemGrid->GetNDOFsPerSubElem(VTK_Type);
          const unsigned short *connSubFaces = surfElem[i].standardElemGrid->GetSubConnType1();

          /*--- Loop over the number of subfaces and store the required data. ---*/
          unsigned short ii = 0;
          for(unsigned short j=0; j<nSubFaces; ++j) {
            markerIDs.push_back(iMarker);
            VTK_TypeElem.push_back(VTK_Type);
            elemIDs.push_back(i);

            for(unsigned short k=0; k<nDOFsPerFace; ++k, ++ii)
              surfaceConn.push_back(surfElem[i].nodeIDsGrid[connSubFaces[ii]]);
          }
        }
      }
    }
  }

  /*--- Create the coordinates of the local points on the viscous surfaces and
        create the final version of the mapping from all volume points to the
        points on the viscous surfaces. ---*/
  vector<su2double> surfaceCoor;
  unsigned long nVertex_SolidWall = 0;

  for(unsigned long i=0; i<meshPoints.size(); ++i) {
    if( meshToSurface[i] ) {
      meshToSurface[i] = nVertex_SolidWall++;

      for(unsigned short k=0; k<nDim; ++k)
        surfaceCoor.push_back(meshPoints[i].coor[k]);
    }
  }

  /*--- Change the surface connectivity, such that it corresponds to
        the entries in surfaceCoor rather than in meshPoints. ---*/
  for(unsigned long i=0; i<surfaceConn.size(); ++i)
    surfaceConn[i] = meshToSurface[surfaceConn[i]];

  /* Build the ADT and return it. */
  std::unique_ptr<CADTElemClass> WallADT(new CADTElemClass(nDim, surfaceCoor, surfaceConn, VTK_TypeElem,
                                                           markerIDs, elemIDs, true));
  return WallADT;
}
