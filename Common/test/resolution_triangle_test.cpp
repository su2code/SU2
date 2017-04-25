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

void WriteTriangleMeshFile () {
  /*--- Local variables ---*/
    int KindElem, KindBound;
    int iElem, nNode, mNode;
    int iNode, jNode, nPoint, iPoint, jPoint, kPoint;
    ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    KindElem  = 5; // Triangle
    KindBound = 3; // Line

    /*--- Store the number of nodes and output mesh filename ---*/
    nNode     = 3;
    mNode     = 3;

    /*--- Open .su2 grid file ---*/
    Mesh_File.precision(15);
    Mesh_File.open("triangletestgrid.su2", ios::out);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << endl;
    Mesh_File << "% Problem dimension" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NDIME= 2" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "% Inner element connectivity" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NELEM= " <<  2*(nNode-1)*(mNode-1) << endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
    for (jNode = 0; jNode < mNode-1; jNode++) {
      for (iNode = 0; iNode < nNode-1; iNode++) {
        iPoint = jNode*nNode + iNode;
        jPoint = jNode*nNode + iNode + 1;
        kPoint = (jNode + 1)*nNode + iNode;
        Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
        iElem ++;
        iPoint = jNode*nNode + (iNode + 1);
        jPoint = (jNode + 1)*nNode + (iNode + 1);
        kPoint = (jNode + 1)*nNode + iNode;
        Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
        iElem++;
      }
    }

    /*--- Compute the number of nodes and write the node coordinates ---*/
    nPoint = (nNode)*(mNode);
    Mesh_File << "%" << endl;
    Mesh_File << "% Node coordinates" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NPOIN= " << nNode*mNode << endl;
    iPoint = 0;
    for (jNode = 0; jNode < mNode; jNode++) {
      for (iNode = 0; iNode < nNode; iNode++) {
        Mesh_File << ((double)iNode)/((double)(nNode-1)) << "\t" << ((double)jNode)/((double)(mNode-1)) << "\t" << iPoint << endl;
        iPoint++;
      }
    }

    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << endl;
    Mesh_File << "% Boundary elements" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NMARK= 4" << endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= lower" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
    for (iNode = 0; iNode < nNode-1; iNode++) {
      Mesh_File << KindBound << "\t" << iNode << "\t" << (iNode + 1) << endl;
    }
    Mesh_File << "MARKER_TAG= right" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
    for (jNode = 0; jNode < mNode-1; jNode++) {
      Mesh_File << KindBound << "\t" << jNode*nNode + (nNode - 1) << "\t" << (jNode + 1)*nNode + (nNode - 1) << endl;
    }
    Mesh_File << "MARKER_TAG= upper" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
    for (iNode = 0; iNode < nNode-1; iNode++) {
      Mesh_File << KindBound << "\t" << (nNode*mNode - 1) - iNode << "\t" << (nNode*mNode - 1) - (iNode + 1) << endl;
    }
    Mesh_File << "MARKER_TAG= left" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
    for (jNode = mNode-2; jNode > mNode-4; jNode--) {
      Mesh_File << KindBound << "\t" << (jNode + 1)*nNode << "\t" << jNode*nNode << endl;
    }

    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteCfgFile() {
  std::ofstream cfg_file;

  cfg_file.open("triangletestgrid.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "MARKER_FAR= ( lower upper left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= triangletestgrid.su2" << std::endl;
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
  WriteTriangleMeshFile();
  WriteCfgFile();

  // The use of "geometry_aux" is necessary to mock a multigrid configuration
  CConfig* config = new CConfig("triangletestgrid.cfg", SU2_CFD, 0, 1, 2, VERB_NONE);
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
  su2double machine_eps = std::numeric_limits<su2double>::epsilon();

  geometry->SetResolutionTensor();

  bool entries_correct = true;

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    su2double** Mij = geometry->node[iPoint]->GetResolutionTensor();

    // ---------------------------------------------------------------------------
    // Check that the values of Mij are correct
    if (std::abs(Mij[0][0] - 0.5) > machine_eps) entries_correct = false;
    if (std::abs(Mij[1][0] - 0.0) > machine_eps) entries_correct = false;
    if (std::abs(Mij[0][1] - 0.0) > machine_eps) entries_correct = false;
    if (std::abs(Mij[1][1] - 0.5) > machine_eps) entries_correct = false;

    if (not(entries_correct)) {
      std::cout << "ERROR: The resolution tensor for a triangle was not correct."
          << std::endl;
      std::cout << "The elements of the array were incorrect." << std::endl;
      std::cout << "Array elements: [[" << Mij[0][0] << "," << Mij[0][1] << "],[";
      std::cout << Mij[1][0] << "," << Mij[1][1] << "]]" << std::endl;
      return_flag = 1;
      break;
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





