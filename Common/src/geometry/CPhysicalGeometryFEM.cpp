/*!
 * \file CPhysicalGeometry.cpp
 * \brief Implementation of the FEM part physical geometry class.
 * \author F. Palacios, T. Economon
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

#include "../../include/geometry/primal_grid/CPrimalGridFEM.hpp"
#include "../../include/geometry/primal_grid/CPrimalGridBoundFEM.hpp"
#include "../../include/geometry/CPhysicalGeometry.hpp"

void CPhysicalGeometry::LoadLinearlyPartitionedPointsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Get the partitioned coordinates and their global IDs from the mesh object. ---*/
  const auto &gridCoords     = mesh->GetLocalPointCoordinates();
  const auto &globalPointIDs = mesh->GetGlobalPointIDs();

  /*--- Initialize point counts and the grid node data structure. ---*/
  nPointNode = nPoint;
  nodes = new CPoint(nPoint, nDim);

  /*--- Loop over the points and set the coordinates and global index. ---*/
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    for (unsigned short iDim = 0; iDim < nDim; ++iDim)
      nodes->SetCoord(iPoint, iDim, gridCoords[iDim][iPoint]);
    nodes->SetGlobalIndex(iPoint, globalPointIDs[iPoint]);
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedVolumeElementsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Reset the global to local element mapping. ---*/
  Global_to_Local_Elem.clear();

  /*--- Get the volume connectivity from the mesh object. ---*/
  const auto &dataElems = mesh->GetLocalVolumeElementConnectivity();

  /*--- Allocate space for the interior elements in our SU2 data
        structure. Note that we only instantiate our rank's local set. ---*/
  elem = new CPrimalGrid*[nElem] ();

  /*--- Loop over all of the internal, local volumetric elements. ---*/
  unsigned long ind = 0;
  for (unsigned long jElem=0; jElem<nElem; ++jElem) {

    /*--- Create a FEM element from the data dataElems. ---*/
    const auto *dataElem = dataElems.data() + ind;
    elem[jElem] = new CPrimalGridFEM(dataElem);

    /*--- Store the global to local mapping in Global_to_Local_Elem. ---*/
    Global_to_Local_Elem[dataElem[4]] = jElem;

    /*--- Update ind for the next element. ---*/
    ind += dataElem[3] + 5;
  }
}

void CPhysicalGeometry::LoadLinearlyPartitionedSurfaceElementsFEM(CConfig *config, CMeshReader *mesh) {

  /*--- Store the number of markers and print to the screen. ---*/
  nMarker = mesh->GetNumberOfMarkers();
  config->SetnMarker_All(nMarker);
  if (rank == MASTER_NODE)
    cout << nMarker << " surface markers." << endl;

  /*--- Create the data structure for boundary elements. ---*/
  bound         = new CPrimalGrid**[nMarker];
  nElem_Bound   = new unsigned long [nMarker];
  Tag_to_Marker = new string [config->GetnMarker_Max()];

  /*--- Retrieve the name of the surface markers as well as
        the number of surface elements for every marker. ---*/
  const auto &sectionNames       = mesh->GetMarkerNames();
  const auto &nSurfElemPerMarker = mesh->GetNumberOfSurfaceElementsAllMarkers();

  /*--- Loop over all sections that we extracted from the grid file
        that were identified as boundary element sections so that we can
        store those elements into our SU2 data structures. ---*/
  for (int iMarker = 0; iMarker < nMarker; ++iMarker) {

    /*--- Get the string name and set the number of surface elements
          for this marker. ---*/
    string Marker_Tag    = sectionNames[iMarker];
    nElem_Bound[iMarker] = nSurfElemPerMarker[iMarker];

    /*--- Allocate the memory of the pointers for the surface
          elements for this marker. ---*/
    bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

    /*--- Retrieve the boundary element data for this marker. ---*/
    const auto &dataElems = mesh->GetSurfaceElementConnectivityForMarker(iMarker);

    /*--- Loop over the number of boundary elements for this marker. ---*/
    unsigned long ind = 0;
    for (unsigned long jElem=0; jElem<nElem_Bound[iMarker]; ++jElem) {

      /*--- Create a boundary FEM element from the data dataElems. ---*/
      const auto *dataElem = dataElems.data() + ind;
      bound[iMarker][jElem] = new CPrimalGridBoundFEM(dataElem);

      /*--- Update ind for the next element. ---*/
      ind += dataElem[2] + 5;
    }

    /*--- Update config file lists in order to store the boundary
          information for this marker in the correct place. ---*/
    Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
    config->SetMarker_All_TagBound(iMarker, Marker_Tag);
    config->SetMarker_All_KindBC(iMarker, config->GetMarker_CfgFile_KindBC(Marker_Tag));
    config->SetMarker_All_Monitoring(iMarker, config->GetMarker_CfgFile_Monitoring(Marker_Tag));
    config->SetMarker_All_GeoEval(iMarker, config->GetMarker_CfgFile_GeoEval(Marker_Tag));
    config->SetMarker_All_Designing(iMarker, config->GetMarker_CfgFile_Designing(Marker_Tag));
    config->SetMarker_All_Plotting(iMarker, config->GetMarker_CfgFile_Plotting(Marker_Tag));
    config->SetMarker_All_Analyze(iMarker, config->GetMarker_CfgFile_Analyze(Marker_Tag));
    config->SetMarker_All_ZoneInterface(iMarker, config->GetMarker_CfgFile_ZoneInterface(Marker_Tag));
    config->SetMarker_All_DV(iMarker, config->GetMarker_CfgFile_DV(Marker_Tag));
    config->SetMarker_All_Moving(iMarker, config->GetMarker_CfgFile_Moving(Marker_Tag));
    config->SetMarker_All_Deform_Mesh(iMarker, config->GetMarker_CfgFile_Deform_Mesh(Marker_Tag));
    config->SetMarker_All_Fluid_Load(iMarker, config->GetMarker_CfgFile_Fluid_Load(Marker_Tag));
    config->SetMarker_All_PyCustom(iMarker, config->GetMarker_CfgFile_PyCustom(Marker_Tag));
    config->SetMarker_All_PerBound(iMarker, config->GetMarker_CfgFile_PerBound(Marker_Tag));
    config->SetMarker_All_SendRecv(iMarker, NONE);
    config->SetMarker_All_Turbomachinery(iMarker, config->GetMarker_CfgFile_Turbomachinery(Marker_Tag));
    config->SetMarker_All_TurbomachineryFlag(iMarker, config->GetMarker_CfgFile_TurbomachineryFlag(Marker_Tag));
    config->SetMarker_All_MixingPlaneInterface(iMarker, config->GetMarker_CfgFile_MixingPlaneInterface(Marker_Tag));
  }
}

void CPhysicalGeometry::SetColorFEMGrid_Parallel(CConfig *config) {

  SU2_MPI::Error(string("Not implemented yet"), CURRENT_FUNCTION);
}

