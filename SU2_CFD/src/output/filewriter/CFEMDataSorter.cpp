/*!
 * \file CFEMDataSorter.cpp
 * \brief Datasorter class for FEM solvers.
 * \author T. Albring
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

#include "../../../include/output/filewriter/CFEMDataSorter.hpp"
#include "../../../../Common/include/fem/fem_geometry_structure.hpp"
#include <numeric>

CFEMDataSorter::CFEMDataSorter(CConfig *config, CGeometry *geometry, const vector<string> &valFieldNames) :
  CParallelDataSorter(config, valFieldNames){

  nDim = geometry->GetnDim();

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
   geometrical information for the FEM DG solver. ---*/

  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  /*--- Create the map from the global DOF ID to the local index. ---*/

  vector<unsigned long> globalID;

  /*--- Update the solution by looping over the owned volume elements. ---*/

  for(unsigned long l=0; l<nVolElemOwned; ++l) {

    /* Count up the number of local points we have for allocating storage. */

    for(unsigned short j=0; j<volElem[l].nDOFsSol; ++j) {

      const unsigned long globalIndex = volElem[l].offsetDOFsSolGlobal + j;
      globalID.push_back(globalIndex);

      nLocalPointsBeforeSort++;
    }
  }

  SU2_MPI::Allreduce(&nLocalPointsBeforeSort, &nGlobalPointBeforeSort, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, SU2_MPI::GetComm());

  /*--- Create a linear partition --- */

  linearPartitioner.Initialize(nGlobalPointBeforeSort, 0);

  /*--- Prepare the send buffers ---*/

  PrepareSendBuffers(globalID);

}

void CFEMDataSorter::SortConnectivity(CConfig *config, CGeometry *geometry, bool val_sort) {

  /*--- Sort connectivity for each type of element (excluding halos). Note
   In these routines, we sort the connectivity into a linear partitioning
   across all processors based on the global index of the grid nodes. ---*/

  /*--- Sort volumetric grid connectivity. ---*/

  nElemPerType.fill(0);

  SortVolumetricConnectivity(config, geometry, TRIANGLE     );
  SortVolumetricConnectivity(config, geometry, QUADRILATERAL);
  SortVolumetricConnectivity(config, geometry, TETRAHEDRON  );
  SortVolumetricConnectivity(config, geometry, HEXAHEDRON   );
  SortVolumetricConnectivity(config, geometry, PRISM        );
  SortVolumetricConnectivity(config, geometry, PYRAMID      );

  SetTotalElements();

  connectivitySorted = true;

}

void CFEMDataSorter::SortVolumetricConnectivity(CConfig *config, CGeometry *geometry, unsigned short Elem_Type) {

  /* Determine the number of nodes for this element type. */
  unsigned short NODES_PER_ELEMENT = 0;
  switch (Elem_Type) {
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    case TETRAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
      break;
    case HEXAHEDRON:
      NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
      break;
    case PRISM:
      NODES_PER_ELEMENT = N_POINTS_PRISM;
      break;
    case PYRAMID:
      NODES_PER_ELEMENT = N_POINTS_PYRAMID;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
  }

  /*--- Create an object of the class CMeshFEM_DG and retrieve the necessary
        geometrical information for the FEM DG solver. ---*/
  auto *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  unsigned long nVolElemOwned = DGGeometry->GetNVolElemOwned();
  CVolumeElementFEM *volElem  = DGGeometry->GetVolElem();

  const CFEMStandardElement *standardElementsSol = DGGeometry->GetStandardElementsSol();

  /*--- Determine the number of sub-elements on this rank. ---*/
  unsigned long nSubElem_Local = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

     /* Only store the linear sub elements if they are of
        the current type that we are storing. */
     if(Elem_Type == VTK_Type1) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType1();
     if(Elem_Type == VTK_Type2) nSubElem_Local += standardElementsSol[ind].GetNSubElemsType2();
  }

  /* Allocate the memory to store the connectivity if the size is
     larger than zero. */
  int *Conn_SubElem = nullptr;
  if(nSubElem_Local > 0) Conn_SubElem = new int[nSubElem_Local*NODES_PER_ELEMENT];

  /*--- Loop again over the local volume elements and store the global
        connectivities of the sub-elements. Note one is added to the
        index value, because visualization softwares typically use
        1-based indexing. ---*/
  unsigned long kNode = 0;
  for(unsigned long i=0; i<nVolElemOwned; ++i) {

    /* Determine the necessary data from the corresponding standard elem. */
    const unsigned short ind       = volElem[i].indStandardElement;
    const unsigned short VTK_Type1 = standardElementsSol[ind].GetVTK_Type1();
    const unsigned short VTK_Type2 = standardElementsSol[ind].GetVTK_Type2();

    /* Check if the first sub-element is of the required type. */
    if(Elem_Type == VTK_Type1) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType1();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType1();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }

    /* Check if the second sub-element is of the required type. */
    if(Elem_Type == VTK_Type2) {

      /* Get the number of sub-elements and the local connectivity of
         the sub-elements. */
      const unsigned short nSubElems     = standardElementsSol[ind].GetNSubElemsType2();
      const unsigned short *connSubElems = standardElementsSol[ind].GetSubConnType2();

      /* Store the global connectivities. */
      const unsigned short kk = NODES_PER_ELEMENT*nSubElems;
      for(unsigned short k=0; k<kk; ++k, ++kNode)
        Conn_SubElem[kNode] = connSubElems[k] + volElem[i].offsetDOFsSolGlobal + 1;
    }
  }

  nElemPerType[TypeMap.at(Elem_Type)] = nSubElem_Local;

  /*--- Store the particular global element count in the class data,
        and set the class data pointer to the connectivity array. ---*/
  switch (Elem_Type) {
    case TRIANGLE:
      delete [] Conn_Tria_Par;
      Conn_Tria_Par = Conn_SubElem;
      break;
    case QUADRILATERAL:
      delete [] Conn_Quad_Par;
      Conn_Quad_Par = Conn_SubElem;
      break;
    case TETRAHEDRON:
      delete [] Conn_Tetr_Par;
      Conn_Tetr_Par = Conn_SubElem;
      break;
    case HEXAHEDRON:
      delete [] Conn_Hexa_Par;
      Conn_Hexa_Par = Conn_SubElem;
      break;
    case PRISM:
      delete [] Conn_Pris_Par;
      Conn_Pris_Par = Conn_SubElem;
      break;
    case PYRAMID:
      delete [] Conn_Pyra_Par;
      Conn_Pyra_Par = Conn_SubElem;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type", CURRENT_FUNCTION);
      break;
  }
}
