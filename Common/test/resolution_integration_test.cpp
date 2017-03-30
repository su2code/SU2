/*!
 * \file resolution_integration_test.cpp
 * \brief This test checks whether the resolution tensor is correctly set for a grid
 * of quadrilateral cells.
 * \author C. Pederson
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"

void WriteQuadMeshFile () {
  /*--- Local variables ---*/
    int KindElem, KindBound, nDim;
    int iElem, iDim, jDim;
    int iNode, jNode, kNode;
    int iPoint;
    int num_Nodes;
    double xSpacing, ySpacing;

    std::ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    nDim = 2;
    KindElem  = 9; // Quadrilateral
    KindBound = 3; // Line

    /*--- Store the number of nodes in each direction ---*/
    iDim = 4;
    jDim = 4;

    /*--- The grid spacing in each direction ---*/
    xSpacing = 4.0;
    ySpacing = 2.0;

    /*--- Open .su2 grid file ---*/
    Mesh_File.open("quadtestgrid.su2", ios::out);
    Mesh_File << fixed << setprecision(15);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Problem dimension" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NDIME= 2" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Inner element connectivity" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NELEM= " <<  (iDim-1)*(jDim-1) << std::endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
      for (jNode = 0; jNode < jDim-1; jNode++) {
        for (iNode = 0; iNode < iDim-1; iNode++) {
          Mesh_File << KindElem << "\t";
          Mesh_File << iNode     + (jNode*jDim)     << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     << "\t";
          // NOTE: Reverse ordering here is essential
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) << "\t";
          Mesh_File << iElem << std::endl;
          iElem++;
        }
      }


    /*--- Compute the number of nodes and write the node coordinates ---*/
    num_Nodes = iDim*jDim;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Node coordinates" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NPOIN= " << num_Nodes << std::endl;
    iPoint = 0;
      for (jNode = 0; jNode < jDim; jNode++) {
        for (iNode = 0; iNode < iDim; iNode++) {
          Mesh_File << iNode*xSpacing << "\t";
          Mesh_File << jNode*ySpacing << "\t";
          Mesh_File << iPoint << std::endl;
          iPoint++;
        }
      }



    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Boundary elements" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NMARK= 4" << std::endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= lower" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1) << std::endl;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      Mesh_File << KindBound << "\t";
      Mesh_File << iNode       << "\t";
      Mesh_File << (iNode + 1) << std::endl;
    }
    Mesh_File << "MARKER_TAG= right" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1) << std::endl;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      Mesh_File << KindBound << "\t";
      Mesh_File << (jNode+1)*iDim - 1 << "\t";
      Mesh_File << (jNode+2)*iDim - 1 << std::endl;
    }
    Mesh_File << "MARKER_TAG= upper" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1) << std::endl;
    for (iNode = jDim*(iDim-1); iNode < iDim*jDim-1; ++iNode) {
      Mesh_File << KindBound << "\t";
      Mesh_File << iNode       << "\t";
      Mesh_File << (iNode + 1) << std::endl;
    }
    Mesh_File << "MARKER_TAG= left" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1) << std::endl;
    for (jNode = 0; jNode < jDim-1; ++jNode) {
      Mesh_File << KindBound << "\t";
      Mesh_File << (jNode  )*iDim << "\t";
      Mesh_File << (jNode+1)*iDim << std::endl;
    }

    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteCfgFile() {
  std::ofstream cfg_file;

  cfg_file.open("quadtest.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "MARKER_FAR= ( lower upper left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= quadtestgrid.su2" << std::endl;
  cfg_file << "MESH_FORMAT= SU2" << std::endl;

  cfg_file.close();

}

int main() {

  //---------------------------------------------------------------------------
  // Setup
  //---------------------------------------------------------------------------
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  // Write out the mesh and configuration files.
  WriteQuadMeshFile();
  WriteCfgFile();

  // The use of "geometry_aux" is necessary to mock a multigrid configuration
  CConfig* config = new CConfig("quadtest.cfg", SU2_CFD, 0, 1, 2, VERB_NONE);
  CGeometry *geometry_aux = NULL;
  geometry_aux = new CPhysicalGeometry(config, 0, 1);
  CGeometry* geometry = new CGeometry();
  geometry = new CPhysicalGeometry(geometry_aux, config);
  delete geometry_aux;

  // Initialize the geometry
  geometry->SetBoundaries(config);
  geometry->SetPoint_Connectivity();
  geometry->SetElement_Connectivity();
  geometry->SetBoundVolume();
  geometry->Check_IntElem_Orientation(config);
  geometry->Check_BoundElem_Orientation(config);
  geometry->SetEdges();
  geometry->SetVertex(config);
  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, ALLOCATE);
  geometry->SetBoundControlVolume(config, ALLOCATE);

  //---------------------------------------------------------------------------
  // Tests
  //---------------------------------------------------------------------------
  unsigned short iDim, nDim = 2;
  unsigned short iPoint;

  geometry->SetResolutionTensor();

  bool entries_correct = true;

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    su2double** Mij = geometry->node[iPoint]->GetResolutionTensor();

    // ---------------------------------------------------------------------------
    // Check that the values of Mij are correct
    bool entries_correct = true;
    if (Mij[0][0] != 4.0 ) entries_correct = false;
    if (Mij[0][1] != 0.0 ) entries_correct = false;
    if (Mij[1][0] != 0.0 ) entries_correct = false;
    if (Mij[1][1] != 2.0 ) entries_correct = false;

    if (not(entries_correct)) {
      std::cout << "ERROR: The resolution tensor for a quadrilateral was not correct."
          << std::endl;
      std::cout << "The elements of the array were incorrect." << std::endl;
      std::cout << "Array elements: [[" << Mij[0][0] << "," << Mij[0][1] << "],[";
      std::cout << Mij[1][0] << "," << Mij[1][1] << "]]" << std::endl;
      return_flag = 1;
    }
  }

  //---------------------------------------------------------------------------
  // Teardown
  //---------------------------------------------------------------------------
  delete config;
  delete geometry;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





