/*!
 * \file CFEMDataSorter.cpp
 * \brief Datasorter class for FEM solvers.
 * \author T. Albring
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

#include "../../../include/output/filewriter/CFEMDataSorter.hpp"
#include "../../../../Common/include/geometry/fem_grid/CMeshFEM_DG.hpp"
#include <numeric>

CFEMDataSorter::CFEMDataSorter(CConfig *config, CGeometry *geometry, const vector<string> &valFieldNames) :
  CParallelDataSorter(config, valFieldNames){

  nDim = geometry->GetnDim();

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}

CFEMDataSorter::~CFEMDataSorter(){

        delete [] Index;
       delete [] idSend;
  delete linearPartitioner;

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

  SU2_MPI::Error("Not implemented yet", CURRENT_FUNCTION);
}
