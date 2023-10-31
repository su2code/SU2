/*!
 * \file fem_wall_distance.cpp
 * \brief Main subroutines for computing the wall distance for the FEM solver.
 * \author E. van der Weide
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/fem/fem_geometry_structure.hpp"
#include "../../include/adt/CADTElemClass.hpp"

std::unique_ptr<CADTElemClass> CMeshFEM_DG::ComputeViscousWallADT(const CConfig* config) const {
  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Create the coordinates and connectivity of the linear      ---*/
  /*---         subelements of the local boundaries that must be taken     ---*/
  /*---         into account in the wall distance computation.             ---*/
  /*--------------------------------------------------------------------------*/

  /* Initialize an array for the mesh points, which eventually contains the
     mapping from the local nodes to the number used in the connectivity of the
     local boundary faces. However, in a first pass it is an indicator whether
     or not a mesh point is on a local wall boundary. */
  vector<unsigned long> meshToSurface(meshPoints.size(), 0);

  /* Define the vectors for the connectivity of the local linear subelements,
     the element ID's, the element type and marker ID's. */
  vector<unsigned long> surfaceConn;
  vector<unsigned long> elemIDs;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned short> markerIDs;

  /* Loop over the boundary markers. */
  for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
    if (!boundaries[iMarker].periodicBoundary) {
      /* Check for a viscous wall. */
      if ((config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX) ||
          (config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL)) {
        /* Loop over the surface elements of this marker. */
        const vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;
        for (unsigned long i = 0; i < surfElem.size(); ++i) {
          /* Set the flag of the mesh points on this surface to true. */
          for (unsigned short j = 0; j < surfElem[i].nDOFsGrid; ++j) meshToSurface[surfElem[i].DOFsGridFace[j]] = 1;

          /* Determine the necessary data from the corresponding standard face,
             such as the number of linear subfaces, the number of DOFs per
             linear subface and the corresponding local connectivity. */
          const unsigned short ind = surfElem[i].indStandardElement;
          const unsigned short VTK_Type = standardBoundaryFacesGrid[ind].GetVTK_Type();
          const unsigned short nSubFaces = standardBoundaryFacesGrid[ind].GetNSubFaces();
          const unsigned short nDOFsPerFace = standardBoundaryFacesGrid[ind].GetNDOFsPerSubFace();
          const unsigned short* connSubFaces = standardBoundaryFacesGrid[ind].GetSubFaceConn();

          /* Loop over the number of subfaces and store the required data. */
          unsigned short ii = 0;
          for (unsigned short j = 0; j < nSubFaces; ++j) {
            markerIDs.push_back(iMarker);
            VTK_TypeElem.push_back(VTK_Type);
            elemIDs.push_back(i);

            for (unsigned short k = 0; k < nDOFsPerFace; ++k, ++ii)
              surfaceConn.push_back(surfElem[i].DOFsGridFace[connSubFaces[ii]]);
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

  for (unsigned long i = 0; i < meshPoints.size(); ++i) {
    if (meshToSurface[i]) {
      meshToSurface[i] = nVertex_SolidWall++;

      for (unsigned short k = 0; k < nDim; ++k) surfaceCoor.push_back(meshPoints[i].coor[k]);
    }
  }

  /*--- Change the surface connectivity, such that it corresponds to
        the entries in surfaceCoor rather than in meshPoints. ---*/
  for (unsigned long i = 0; i < surfaceConn.size(); ++i) surfaceConn[i] = meshToSurface[surfaceConn[i]];

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Build the ADT, which is an ADT of bounding boxes of the    ---*/
  /*---         surface elements. A nearest point search does not give     ---*/
  /*---         accurate results, especially not for the integration       ---*/
  /*---         points of the elements close to a wall boundary.           ---*/
  /*--------------------------------------------------------------------------*/

  /* Build the ADT. */
  std::unique_ptr<CADTElemClass> WallADT(
      new CADTElemClass(nDim, surfaceCoor, surfaceConn, VTK_TypeElem, markerIDs, elemIDs, true));

  return WallADT;
}

/*!
 * \brief Set wall distances a specific value
 */
void CMeshFEM_DG::SetWallDistance(su2double val) {
  for (unsigned long l = 0; l < nVolElemOwned; ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = volElem[l].indStandardElement;
    const unsigned short nInt = standardElementsGrid[ind].GetNIntegration();

    /* Store the number of solDOFS a bit easier. */
    const unsigned short nDOFsSol = volElem[l].nDOFsSol;

    /* Allocate the memory for the wall distance of this element. */
    volElem[l].wallDistance.resize(nInt, val);
    volElem[l].wallDistanceSolDOFs.resize(nDOFsSol, val);
  }

  for (unsigned long l = 0; l < matchingFaces.size(); ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = matchingFaces[l].indStandardElement;
    const unsigned short nInt = standardMatchingFacesGrid[ind].GetNIntegration();

    /* Allocate the memory for the wall distance for this matching face. */
    matchingFaces[l].wallDistance.resize(nInt, val);
  }

  for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
    if (!boundaries[iMarker].periodicBoundary) {
      /* Loop over the boundary faces and determine the wall distances
         in the integration points. */
      vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;
      for (unsigned long l = 0; l < surfElem.size(); ++l) {
        /* Get the required data from the corresponding standard element. */
        const unsigned short ind = surfElem[l].indStandardElement;
        const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();

        /* Allocate the memory for the wall distance for this boundary face. */
        surfElem[l].wallDistance.resize(nInt, val);
      }
    }
  }
}

void CMeshFEM_DG::SetWallDistance(CADTElemClass* WallADT, const CConfig* config, unsigned short iZone) {
  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Determine the wall distance of the integration points of   ---*/
  /*---         locally owned volume elements.                             ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the owned elements to compute the wall distance
        in the integration points. ---*/
  for (unsigned long l = 0; l < nVolElemOwned; ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = volElem[l].indStandardElement;
    const unsigned short nInt = standardElementsGrid[ind].GetNIntegration();

    /* Allocate the memory for the wall distance of this element. */
    volElem[l].wallDistance.resize(nInt);

    if (!WallADT->IsEmpty()) {
      /*--- The tree is not empty. Loop over the integration points
            and determine the wall distance. ---*/
      for (unsigned short i = 0; i < nInt; ++i) {
        const su2double* coor = volElem[l].coorIntegrationPoints.data() + i * nDim;
        unsigned short markerID;
        unsigned long elemID;
        int rankID;
        su2double dist;
        WallADT->DetermineNearestElement(coor, dist, markerID, elemID, rankID);

        volElem[l].wallDistance[i] = dist;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 4: Determine the wall distance of the solution DOFs of        ---*/
  /*---         locally owned volume elements.                             ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the owned elements to compute the wall distance
        in the integration points. ---*/
  for (unsigned long l = 0; l < nVolElemOwned; ++l) {
    /* Store the number of solDOFS a bit easier. */
    const unsigned short nDOFsSol = volElem[l].nDOFsSol;

    /* Allocate the memory for the wall distance of the solution DOFs. */
    volElem[l].wallDistanceSolDOFs.resize(nDOFsSol);

    if (!WallADT->IsEmpty()) {
      /*--- The tree is not empty. Loop over the solution DOFs
            and determine the wall distance. ---*/
      for (unsigned short i = 0; i < nDOFsSol; ++i) {
        const su2double* coor = volElem[l].coorSolDOFs.data() + i * nDim;
        unsigned short markerID;
        unsigned long elemID;
        int rankID;
        su2double dist;
        WallADT->DetermineNearestElement(coor, dist, markerID, elemID, rankID);

        volElem[l].wallDistanceSolDOFs[i] = dist;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 5: Determine the wall distance of the integration points of   ---*/
  /*---         the internal matching faces.                               ---*/
  /*--------------------------------------------------------------------------*/

  for (unsigned long l = 0; l < matchingFaces.size(); ++l) {
    /* Get the required data from the corresponding standard element. */
    const unsigned short ind = matchingFaces[l].indStandardElement;
    const unsigned short nInt = standardMatchingFacesGrid[ind].GetNIntegration();

    /* Allocate the memory for the wall distance for this matching face. */
    matchingFaces[l].wallDistance.resize(nInt);

    if (!WallADT->IsEmpty()) {
      /*--- The tree is not empty. Loop over the integration points
            and determine the wall distance. */
      for (unsigned short i = 0; i < nInt; ++i) {
        const su2double* coor = matchingFaces[l].coorIntegrationPoints.data() + i * nDim;
        unsigned short markerID;
        unsigned long elemID;
        int rankID;
        su2double dist;
        WallADT->DetermineNearestElement(coor, dist, markerID, elemID, rankID);

        matchingFaces[l].wallDistance[i] = dist;
      }
    }
  }

  /*--------------------------------------------------------------------------*/
  /*--- Step 6: Determine the wall distance of the integration points of   ---*/
  /*---         locally owned boundary faces.                              ---*/
  /*--------------------------------------------------------------------------*/

  /*--- Loop over the boundary markers. Make sure to exclude the periodic
        boundaries, because these are not physical. ---*/
  for (unsigned short iMarker = 0; iMarker < boundaries.size(); ++iMarker) {
    if (!boundaries[iMarker].periodicBoundary) {
      /* Determine whether or not this is a viscous wall boundary. */
      const bool viscousWall =
          config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL;

      /* Loop over the boundary faces and determine the wall distances
         in the integration points. */
      vector<CSurfaceElementFEM>& surfElem = boundaries[iMarker].surfElem;
      for (unsigned long l = 0; l < surfElem.size(); ++l) {
        /* Get the required data from the corresponding standard element. */
        const unsigned short ind = surfElem[l].indStandardElement;
        const unsigned short nInt = standardBoundaryFacesGrid[ind].GetNIntegration();

        /* Allocate the memory for the wall distance for this boundary face. */
        surfElem[l].wallDistance.resize(nInt);

        /* Check for an empty tree or a viscous wall.
           In those case the wall distance is set to zero. */
        if (viscousWall) {
          /* Wall distance must be set to zero. */
          for (unsigned short i = 0; i < nInt; ++i) surfElem[l].wallDistance[i] = 0.0;
        } else if (!WallADT->IsEmpty()) {
          /*--- Not a viscous wall boundary, while viscous walls are present.
                The distance must be computed. Loop over the integration points
                and do so. ---*/
          for (unsigned short i = 0; i < nInt; ++i) {
            const su2double* coor = surfElem[l].coorIntegrationPoints.data() + i * nDim;
            unsigned short markerID;
            unsigned long elemID;
            int rankID;
            su2double dist;
            WallADT->DetermineNearestElement(coor, dist, markerID, elemID, rankID);

            surfElem[l].wallDistance[i] = dist;
          }
        }
      }
    }
  }
}
