/*!
 * \file resolution_gradient_test.cpp
 * \brief This test checks whether the gradient of the resolution tensor is
 *        calculated correctly for 3D hexahedra.
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
    int iElem, iDim, jDim, kDim;
    int iNode, jNode, kNode;
    int iPoint;
    int num_Nodes;
    double xSpacing, ySpacing, zSpacing;
    double xFactor, yFactor, zFactor;
    double xCoord[5] = {0, 1, 2, 5, 8}; // Variable gradient
    double yCoord[5] = {0, 2, 4, 6, 8}; // Gradient of 0
    double zCoord[5] = {0, 1, 2, 3, 4}; // Gradient of 0

    std::ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    nDim = 3;
    KindElem  = 12; // Hexahedra
    KindBound = 9; // Quadrilateral

    /*--- Store the number of nodes in each direction ---*/
    iDim = 5;
    jDim = 5;
    kDim = 5;

    /*--- Open .su2 grid file ---*/
    Mesh_File.open("gradtestgrid.su2", ios::out);
    Mesh_File << fixed << setprecision(15);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Problem dimension" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NDIME= 3" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Inner element connectivity" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NELEM= " <<  (iDim-1)*(jDim-1)*(kDim-1) << std::endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
    for (kNode = 0; kNode < kDim-1; kNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        for (iNode = 0; iNode < iDim-1; iNode++) {
          Mesh_File << KindElem << "\t";
          // Proper ordering here is essential.
          // See VTK documentation for hexahedral cells for ordering.
          Mesh_File << iNode     + (jNode*jDim)     + (kNode*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     + (kNode*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) + (kNode*jDim*iDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) + (kNode*jDim*iDim) << "\t";
          Mesh_File << iNode     + (jNode*jDim)     + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << iElem << std::endl;
          iElem++;
        }
      }
    }


    /*--- Compute the number of nodes and write the node coordinates ---*/
    num_Nodes = iDim*jDim*kDim;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Node coordinates" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NPOIN= " << num_Nodes << std::endl;
    iPoint = 0;
    for (kNode = 0; kNode < kDim; kNode++) {
      for (jNode = 0; jNode < jDim; jNode++) {
        for (iNode = 0; iNode < iDim; iNode++) {
          Mesh_File << xCoord[iNode] << "\t";
          Mesh_File << yCoord[jNode] << "\t";
          Mesh_File << zCoord[kNode] << "\t";
          Mesh_File << iPoint << std::endl;
          iPoint++;
        }
      }
    }

    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Boundary elements" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NMARK= 6" << std::endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= bottom" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(jDim-1) << std::endl;
    kNode = 0;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + (jNode+1)*jDim + kNode*(iDim*jDim) << "\t";
        Mesh_File << iNode     + (jNode+1)*jDim + kNode*(iDim*jDim) << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= top" << std::endl;
    Mesh_File << "MARKER_ELEMS= " << (iDim-1)*(jDim-1) << std::endl;
    kNode = kDim-1;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + (jNode+1)*jDim + kNode*(iDim*jDim) << "\t";
        Mesh_File << iNode     + (jNode+1)*jDim + kNode*(iDim*jDim) << std::endl;
      }
    }

    Mesh_File << "MARKER_TAG= left" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1)*(kDim-1) << std::endl;
    iNode = 0;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode + jNode*jDim     + kNode*(iDim*jDim)     << "\t";
        Mesh_File << iNode + jNode*jDim     + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + kNode*(iDim*jDim)     << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= right" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1)*(kDim-1) << std::endl;
    iNode = iDim-1;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode + jNode*jDim     + kNode*(iDim*jDim)     << "\t";
        Mesh_File << iNode + jNode*jDim     + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + kNode*(iDim*jDim)     << std::endl;
      }
    }

    Mesh_File << "MARKER_TAG= front" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(kDim-1) << std::endl;
    jNode = 0;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode     + jNode*jDim + (kNode+1)*(iDim*jDim) << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= back" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(kDim-1) << std::endl;
    jNode = jDim-1;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode     + jNode*jDim + (kNode+1)*(iDim*jDim) << std::endl;
      }
    }
    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteCfgFile() {
  std::ofstream cfg_file;

  cfg_file.open("gradtest.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "MARKER_FAR= ( top bottom back front left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= gradtestgrid.su2" << std::endl;
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
  CConfig* config = new CConfig("gradtest.cfg", SU2_CFD, 0, 1, 2, VERB_NONE);
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

  /**
   * Due to the way that the dual mesh is set up, the edge control volumes
   * (the ones lying on the boundary) have a different size than the interior
   * control volumes.
   *
   * The gradient correct should ignore these boundary artifacts.
   */

  unsigned short iDim, nDim = 3;
  unsigned short iPoint;

  geometry->SetResolutionTensor();

  bool entries_correct = true;

  // These are hand-calculated numerical derivatives.
  su2double correct_grads[5] = {0, 1.5, 7.0/3, 5.0/6, 0};

/*  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
     su2double** Mij = geometry->node[iPoint]->GetResolutionTensor();
     su2double* coord;
     coord = geometry->node[iPoint]->GetCoord();

     std::cout << "The resolution tensor for point: ";
     std::cout << iPoint << std::endl;
     std::cout << "    at location: [";
     std::cout << coord[0] << ", " << coord[1] << ", " << coord[2];
     std::cout << "]" << std::endl;
     std::cout << "    and direction: " << iDim << std::endl;

     std::cout << "Is:" << std::endl;
     std::cout << "    [[";
     std::cout << Mij[0][0] << "," << Mij[0][1] << "," << Mij[0][2];
     std::cout << "]" << std::endl;
     std::cout << "     [";
     std::cout << Mij[1][0] << "," << Mij[1][1] << "," << Mij[1][2];
     std::cout << "]" << std::endl;
     std::cout << "     [";
     std::cout << Mij[2][0] << "," << Mij[2][1] << "," << Mij[2][2];
     std::cout << "]]" << std::endl;
  }*/

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double** dMsqdx = geometry->node[iPoint]->GetResolutionGradient(iDim);
      su2double* coord;
      coord = geometry->node[iPoint]->GetCoord();

      // Entries of dMdx are indexed as d(M_ij)/dx = dMdx[i][j]
      for (int i = 0; i<nDim; ++i) {
        for (int j = 0; j<nDim; ++j) {
          if (iDim == 0 and i==0 and j == 0) {
            int index;
            switch (int(coord[0])) {
              case 0:
                index = 0;
                break;
              case 1:
                index = 1;
                break;
              case 2:
                index = 2;
                break;
              case 5:
                index = 3;
                break;
              case 8:
                index = 4;
                break;
              default:
                std::cout << "An error occurred in the test program at: ";
                std::cout << __FILE__ << ", line # : " << __LINE__ << std::endl;
                entries_correct = false;
            }
            if (std::abs(dMsqdx[i][j] - correct_grads[index]) > 1e-6)
              entries_correct = false;
          } else {
            if (std::abs(dMsqdx[i][j]) > 1e-6) entries_correct = false;
          }
        }
      }

      if (not(entries_correct)) {
        std::cout << "ERROR: The gradient of the resolution tensor was not correct."
            << std::endl;
        std::cout << "The elements of the array were incorrect for point: ";
        std::cout << iPoint << std::endl;
        std::cout << "    at location: [";
        std::cout << coord[0] << ", " << coord[1] << ", " << coord[2];
        std::cout << "]" << std::endl;
        std::cout << "    and direction: " << iDim << std::endl;

        std::cout << "    Found:" << std::endl;
        std::cout << "    [[";
        std::cout << dMsqdx[0][0] << "," << dMsqdx[0][1] << "," << dMsqdx[0][2];
        std::cout << "]" << std::endl;
        std::cout << "     [";
        std::cout << dMsqdx[1][0] << "," << dMsqdx[1][1] << "," << dMsqdx[1][2];
        std::cout << "]" << std::endl;
        std::cout << "     [";
        std::cout << dMsqdx[2][0] << "," << dMsqdx[2][1] << "," << dMsqdx[2][2];
        std::cout << "]]" << std::endl;

        return_flag = 1;
        break;
      }
    }
    if (not(entries_correct)) break;
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





