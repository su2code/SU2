/*!
 * \file CPhysicalGeometry.cpp
 * \brief Implementation of the physical geometry class.
 * \author F. Palacios, T. Economon
 * \version 7.0.4 "Blackbird"
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

#include "../../include/geometry/CPhysicalGeometry.hpp"
#include "../../include/adt_structure.hpp"
#include "../../include/toolboxes/printing_toolbox.hpp"
#include "../../include/toolboxes/CLinearPartitioner.hpp"
#include "../../include/geometry/meshreader/CSU2ASCIIMeshReaderFVM.hpp"
#include "../../include/geometry/meshreader/CCGNSMeshReaderFVM.hpp"
#include "../../include/geometry/meshreader/CRectangularMeshReaderFVM.hpp"
#include "../../include/geometry/meshreader/CBoxMeshReaderFVM.hpp"

#include "../../include/geometry/primal_grid/CPrimalGrid.hpp"
#include "../../include/geometry/primal_grid/CLine.hpp"
#include "../../include/geometry/primal_grid/CTriangle.hpp"
#include "../../include/geometry/primal_grid/CQuadrilateral.hpp"
#include "../../include/geometry/primal_grid/CTetrahedron.hpp"
#include "../../include/geometry/primal_grid/CHexahedron.hpp"
#include "../../include/geometry/primal_grid/CPyramid.hpp"
#include "../../include/geometry/primal_grid/CPrism.hpp"
#include "../../include/geometry/primal_grid/CVertexMPI.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <iterator>
#ifdef _MSC_VER
#include <direct.h>
#endif


CPhysicalGeometry::CPhysicalGeometry() : CGeometry() {

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  Local_to_Global_Point  = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;

  beg_node = NULL;
  end_node = NULL;

  nPointLinear      = NULL;
  nPointCumulative  = NULL;

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  adjacency = NULL;
  xadj      = NULL;
#endif
#endif

  /*--- Arrays for defining the turbomachinery structure ---*/

  nSpanWiseSections       = NULL;
  nSpanSectionsByMarker   = NULL;
  SpanWiseValue           = NULL;
  nVertexSpan             = NULL;
  nTotVertexSpan          = NULL;
  turbovertex             = NULL;
  AverageTurboNormal      = NULL;
  AverageNormal           = NULL;
  AverageGridVel          = NULL;
  AverageTangGridVel      = NULL;
  SpanArea                = NULL;
  TurboRadius             = NULL;
  MaxAngularCoord         = NULL;
  MinAngularCoord         = NULL;
  MinRelAngularCoord      = NULL;

  TangGridVelIn           = NULL;
  SpanAreaIn              = NULL;
  TurboRadiusIn           = NULL;
  TangGridVelOut          = NULL;
  SpanAreaOut             = NULL;
  TurboRadiusOut          = NULL;
  SpanWiseValue           = NULL;

}

CPhysicalGeometry::CPhysicalGeometry(CConfig *config, unsigned short val_iZone, unsigned short val_nZone) : CGeometry() {

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  Local_to_Global_Point = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;

  beg_node = NULL;
  end_node = NULL;

  nPointLinear      = NULL;
  nPointCumulative  = NULL;

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  adjacency = NULL;
  xadj      = NULL;
#endif
#endif

  edgeColorGroupSize = config->GetEdgeColoringGroupSize();

  /*--- Arrays for defining the turbomachinery structure ---*/

  nSpanWiseSections       = NULL;
  nSpanSectionsByMarker   = NULL;
  SpanWiseValue           = NULL;
  nVertexSpan             = NULL;
  nTotVertexSpan          = NULL;
  turbovertex             = NULL;
  AverageTurboNormal      = NULL;
  AverageNormal           = NULL;
  AverageGridVel          = NULL;
  AverageTangGridVel      = NULL;
  SpanArea                = NULL;
  TurboRadius             = NULL;
  MaxAngularCoord         = NULL;
  MinAngularCoord         = NULL;
  MinRelAngularCoord      = NULL;

  TangGridVelIn           = NULL;
  SpanAreaIn              = NULL;
  TurboRadiusIn           = NULL;
  TangGridVelOut          = NULL;
  SpanAreaOut             = NULL;
  TurboRadiusOut          = NULL;
  SpanWiseValue           = NULL;

  string text_line, Marker_Tag;
  ifstream mesh_file;
  unsigned short iDim, iMarker, iNodes;
  unsigned long iPoint, iElem_Bound;
  su2double *NewCoord;
  nZone = val_nZone;
  ofstream boundary_file;
  string Grid_Marker;

  string val_mesh_filename  = config->GetMesh_FileName();
  unsigned short val_format = config->GetMesh_FileFormat();

  /*--- Determine whether or not a FEM discretization is used ---*/

  const bool fem_solver = config->GetFEMSolver();

  /*--- Initialize counters for local/global points & elements ---*/

  if( fem_solver ) {
    switch (val_format) {
      case SU2:
        Read_SU2_Format_Parallel_FEM(config, val_mesh_filename, val_iZone, val_nZone);
        break;

      case CGNS_GRID:
        Read_CGNS_Format_Parallel_FEM(config, val_mesh_filename, val_iZone, val_nZone);
        break;

      default:
        SU2_MPI::Error("Unrecognized mesh format specified for the FEM solver!", CURRENT_FUNCTION);
        break;
    }
  }
  else {

    switch (val_format) {
      case SU2: case CGNS_GRID: case RECTANGLE: case BOX:
        Read_Mesh_FVM(config, val_mesh_filename, val_iZone, val_nZone);
        break;
      default:
        SU2_MPI::Error("Unrecognized mesh format specified!", CURRENT_FUNCTION);
        break;
    }
  }

  /*--- After reading the mesh, assert that the dimension is equal to 2 or 3. ---*/

  assert((nDim == 2) || (nDim == 3));

  /*--- Loop over the points element to re-scale the mesh, and plot it (only SU2_CFD) ---*/

  if (config->GetKind_SU2() == SU2_CFD) {

    NewCoord = new su2double [nDim];

    /*--- The US system uses feet, but SU2 assumes that the grid is in inches ---*/

    if (config->GetSystemMeasurements() == US) {
      for (iPoint = 0; iPoint < nPoint; iPoint++) {
        for (iDim = 0; iDim < nDim; iDim++) {
          NewCoord[iDim] = node[iPoint]->GetCoord(iDim)/12.0;
        }
        node[iPoint]->SetCoord(NewCoord);
      }
    }

    delete [] NewCoord;

  }

  /*--- If SU2_DEF then write a file with the boundary information ---*/

  if ((config->GetKind_SU2() == SU2_DEF) && (rank == MASTER_NODE)) {

    string str = "boundary.dat";

    str = config->GetMultizone_FileName(str, val_iZone, ".dat");

    /*--- Open .su2 grid file ---*/

    boundary_file.open(str.c_str(), ios::out);

    /*--- Loop through and write the boundary info ---*/

    boundary_file << "NMARK= " << nMarker << endl;

    for (iMarker = 0; iMarker < nMarker; iMarker++) {

      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      boundary_file << "MARKER_TAG= " << Grid_Marker << endl;
      boundary_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      boundary_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;

          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == VERTEX) {
            boundary_file << bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t";
          }
          boundary_file	<< iElem_Bound << endl;
        }
      }

      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          boundary_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            boundary_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;

          if (bound[iMarker][iElem_Bound]->GetVTK_Type() == VERTEX) {
            boundary_file << bound[iMarker][iElem_Bound]->GetRotation_Type() << "\t";
          }
          boundary_file	<< iElem_Bound << endl;
        }
      }

    }

    boundary_file.close();

  }

}

CPhysicalGeometry::CPhysicalGeometry(CGeometry *geometry,
                                     CConfig *config) {

  /*--- Get rank and size. ---*/

  size = SU2_MPI::GetSize();
  rank = SU2_MPI::GetRank();

  /*--- Initialize several class data members for later. ---*/

  Local_to_Global_Point  = NULL;
  Local_to_Global_Marker = NULL;
  Global_to_Local_Marker = NULL;

  /*--- Arrays for defining the linear partitioning. ---*/

  beg_node = NULL;
  end_node = NULL;

  nPointLinear      = NULL;
  nPointCumulative  = NULL;

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS
  adjacency = NULL;
  xadj      = NULL;
#endif
#endif

  edgeColorGroupSize = config->GetEdgeColoringGroupSize();

  /*--- Arrays for defining the turbomachinery structure ---*/

  nSpanWiseSections       = NULL;
  nSpanSectionsByMarker   = NULL;
  SpanWiseValue           = NULL;
  nVertexSpan             = NULL;
  nTotVertexSpan          = NULL;
  turbovertex             = NULL;
  AverageTurboNormal      = NULL;
  AverageNormal           = NULL;
  AverageGridVel          = NULL;
  AverageTangGridVel      = NULL;
  SpanArea                = NULL;
  TurboRadius             = NULL;
  MaxAngularCoord         = NULL;
  MinAngularCoord         = NULL;
  MinRelAngularCoord      = NULL;

  TangGridVelIn           = NULL;
  SpanAreaIn              = NULL;
  TurboRadiusIn           = NULL;
  TangGridVelOut          = NULL;
  SpanAreaOut             = NULL;
  TurboRadiusOut          = NULL;

  /*--- Initialize counters for the points/elements local to a rank. ---*/

  nLocal_Point         = 0;
  nLocal_PointDomain   = 0;
  nLocal_PointGhost    = 0;
  nLocal_PointPeriodic = 0;
  nLocal_Line          = 0;
  nLocal_BoundTria     = 0;
  nLocal_BoundQuad     = 0;
  nLocal_Tria          = 0;
  nLocal_Quad          = 0;
  nLocal_Tetr          = 0;
  nLocal_Hexa          = 0;
  nLocal_Pris          = 0;
  nLocal_Pyra          = 0;

  Local_Coords = NULL;
  Local_Points = NULL;
  Local_Colors = NULL;

  /*--- Arrays for holding the element connectivity. ---*/

  Conn_Line      = NULL;
  Conn_BoundTria = NULL;
  Conn_BoundQuad = NULL;

  Conn_Line_Linear      = NULL;
  Conn_BoundTria_Linear = NULL;
  Conn_BoundQuad_Linear = NULL;

  Conn_Tria = NULL;
  Conn_Quad = NULL;
  Conn_Tetr = NULL;
  Conn_Hexa = NULL;
  Conn_Pris = NULL;
  Conn_Pyra = NULL;

  /*--- Arrays for holding the element IDs. ---*/

  ID_Line             = NULL;
  ID_BoundTria        = NULL;
  ID_BoundQuad        = NULL;
  ID_Line_Linear      = NULL;
  ID_BoundTria_Linear = NULL;
  ID_BoundQuad_Linear = NULL;

  ID_Tria = NULL;
  ID_Quad = NULL;
  ID_Tetr = NULL;
  ID_Hexa = NULL;
  ID_Pris = NULL;
  ID_Pyra = NULL;

  Elem_ID_Line             = NULL;
  Elem_ID_BoundTria        = NULL;
  Elem_ID_BoundQuad        = NULL;
  Elem_ID_Line_Linear      = NULL;
  Elem_ID_BoundTria_Linear = NULL;
  Elem_ID_BoundQuad_Linear = NULL;

  /*--- The new geometry class has the same problem dimension/zone. ---*/

  nDim  = geometry->GetnDim();
  nZone = geometry->GetnZone();

  /*--- Recompute the linear partitioning offsets. ---*/

  PrepareOffsets(geometry->GetGlobal_nPoint());

  /*--- Communicate the coloring data so that each rank has a complete set
   of colors for all points that reside on it, including repeats. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Distributing ParMETIS coloring." << endl;

  DistributeColoring(config, geometry);

  /*--- Redistribute the points to all ranks based on the coloring. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Rebalancing vertices." << endl;

  DistributePoints(config, geometry);

  /*--- Distribute the element information to all ranks based on coloring. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Rebalancing volume element connectivity." << endl;

  DistributeVolumeConnectivity(config, geometry, TRIANGLE     );
  DistributeVolumeConnectivity(config, geometry, QUADRILATERAL);
  DistributeVolumeConnectivity(config, geometry, TETRAHEDRON  );
  DistributeVolumeConnectivity(config, geometry, HEXAHEDRON   );
  DistributeVolumeConnectivity(config, geometry, PRISM        );
  DistributeVolumeConnectivity(config, geometry, PYRAMID      );

  /*--- Distribute the marker information to all ranks based on coloring. ---*/

  if ((rank == MASTER_NODE) && (size != SINGLE_NODE))
    cout <<"Rebalancing markers and surface elements." << endl;

  /*--- First, perform a linear partitioning of the marker information, as
   the grid readers currently store all boundary information on the master
   rank. In the future, this process can be moved directly into the grid
   reader to avoid reading the markers to the master rank alone at first. ---*/

  DistributeMarkerTags(config, geometry);
  PartitionSurfaceConnectivity(config, geometry, LINE         );
  PartitionSurfaceConnectivity(config, geometry, TRIANGLE     );
  PartitionSurfaceConnectivity(config, geometry, QUADRILATERAL);

  /*--- Once the markers are distributed according to the linear partitioning
   of the grid points, we can use similar techniques as above for distributing
   the surface element connectivity. ---*/

  DistributeSurfaceConnectivity(config, geometry, LINE         );
  DistributeSurfaceConnectivity(config, geometry, TRIANGLE     );
  DistributeSurfaceConnectivity(config, geometry, QUADRILATERAL);

  /*--- Reduce the total number of elements that we have on each rank. ---*/

  nLocal_Elem = (nLocal_Tria +
                 nLocal_Quad +
                 nLocal_Tetr +
                 nLocal_Hexa +
                 nLocal_Pris +
                 nLocal_Pyra);
  nLocal_Bound_Elem = nLocal_Line + nLocal_BoundTria + nLocal_BoundQuad;

  SU2_MPI::Allreduce(&nLocal_Elem, &nGlobal_Elem, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&nLocal_Bound_Elem, &nGlobal_Bound_Elem, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  /*--- With the distribution of all points, elements, and markers based
   on the ParMETIS coloring complete, as a final step, load this data into
   our geometry class data structures. ---*/

  LoadPoints(config, geometry);
  LoadVolumeElements(config, geometry);
  LoadSurfaceElements(config, geometry);

  /*--- Free memory associated with the partitioning of points and elems. ---*/

  Neighbors.clear();
  Color_List.clear();

  if (Local_Points != NULL) delete [] Local_Points;
  if (Local_Colors != NULL) delete [] Local_Colors;
  if (Local_Coords != NULL) delete [] Local_Coords;

  if (nLinear_Line > 0      && Conn_Line_Linear      != NULL)
    delete [] Conn_Line_Linear;
  if (nLinear_BoundTria > 0 && Conn_BoundTria_Linear != NULL)
    delete [] Conn_BoundTria_Linear;
  if (nLinear_BoundQuad > 0 && Conn_BoundQuad_Linear != NULL)
    delete [] Conn_BoundQuad_Linear;

  if (nLocal_Line > 0      && Conn_Line      != NULL) delete [] Conn_Line;
  if (nLocal_BoundTria > 0 && Conn_BoundTria != NULL) delete [] Conn_BoundTria;
  if (nLocal_BoundQuad > 0 && Conn_BoundQuad != NULL) delete [] Conn_BoundQuad;
  if (nLocal_Tria > 0      && Conn_Tria      != NULL) delete [] Conn_Tria;
  if (nLocal_Quad > 0      && Conn_Quad      != NULL) delete [] Conn_Quad;
  if (nLocal_Tetr > 0      && Conn_Tetr      != NULL) delete [] Conn_Tetr;
  if (nLocal_Hexa > 0      && Conn_Hexa      != NULL) delete [] Conn_Hexa;
  if (nLocal_Pris > 0      && Conn_Pris      != NULL) delete [] Conn_Pris;
  if (nLocal_Pyra > 0      && Conn_Pyra      != NULL) delete [] Conn_Pyra;

  if (ID_Line             != NULL) delete [] ID_Line;
  if (ID_BoundTria        != NULL) delete [] ID_BoundTria;
  if (ID_BoundQuad        != NULL) delete [] ID_BoundQuad;
  if (ID_Line_Linear      != NULL) delete [] ID_Line_Linear;
  if (ID_BoundTria_Linear != NULL) delete [] ID_BoundTria_Linear;
  if (ID_BoundQuad_Linear != NULL) delete [] ID_BoundQuad_Linear;

  if (ID_Tria != NULL) delete [] ID_Tria;
  if (ID_Quad != NULL) delete [] ID_Quad;
  if (ID_Tetr != NULL) delete [] ID_Tetr;
  if (ID_Hexa != NULL) delete [] ID_Hexa;
  if (ID_Pris != NULL) delete [] ID_Pris;
  if (ID_Pyra != NULL) delete [] ID_Pyra;

  if (Elem_ID_Line             != NULL) delete [] Elem_ID_Line;
  if (Elem_ID_BoundTria        != NULL) delete [] Elem_ID_BoundTria;
  if (Elem_ID_BoundQuad        != NULL) delete [] Elem_ID_BoundQuad;
  if (Elem_ID_Line_Linear      != NULL) delete [] Elem_ID_Line_Linear;
  if (Elem_ID_BoundTria_Linear != NULL) delete [] Elem_ID_BoundTria_Linear;
  if (Elem_ID_BoundQuad_Linear != NULL) delete [] Elem_ID_BoundQuad_Linear;

}

CPhysicalGeometry::~CPhysicalGeometry(void) {

  if (Local_to_Global_Point  != NULL) delete [] Local_to_Global_Point;
  if (Global_to_Local_Marker != NULL) delete [] Global_to_Local_Marker;
  if (Local_to_Global_Marker != NULL) delete [] Local_to_Global_Marker;

  /*--- Free up memory from turbomachinery performance computation  ---*/

  unsigned short iMarker;
  if (TangGridVelIn != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (TangGridVelIn[iMarker] != NULL) delete [] TangGridVelIn[iMarker];
    delete [] TangGridVelIn;
  }
  if (SpanAreaIn != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (SpanAreaIn[iMarker] != NULL) delete [] SpanAreaIn[iMarker];
    delete [] SpanAreaIn;
  }
  if (TurboRadiusIn != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (TurboRadiusIn[iMarker] != NULL) delete [] TurboRadiusIn[iMarker];
    delete [] TurboRadiusIn;
  }
  if (TangGridVelOut != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (TangGridVelOut[iMarker] != NULL) delete [] TangGridVelOut[iMarker];
    delete [] TangGridVelOut;
  }
  if (SpanAreaOut != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (SpanAreaOut[iMarker] != NULL) delete [] SpanAreaOut[iMarker];
    delete [] SpanAreaOut;
  }
  if (TurboRadiusOut != NULL) {
    for (iMarker = 0; iMarker < nTurboPerf; iMarker++)
      if (TurboRadiusOut[iMarker] != NULL) delete [] TurboRadiusOut[iMarker];
    delete [] TurboRadiusOut;
  }

  /*--- Free up memory from turbomachinery computations
   * If there are send/receive boundaries, nMarker isn't the same number
   * as in the constructor. There must be an explicit check to ensure
   * that iMarker doesn't point us to memory that was never allocated. ---*/

  unsigned short iSpan, iVertex;
  if (turbovertex != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Marker_All_SendRecv[iMarker] == 0 && turbovertex[iMarker] != NULL) {
        for (iSpan= 0; iSpan < nSpanSectionsByMarker[iMarker]; iSpan++) {
          if (turbovertex[iMarker][iSpan] != NULL) {
            for (iVertex = 0; iVertex < nVertexSpan[iMarker][iSpan]; iVertex++)
              if (turbovertex[iMarker][iSpan][iVertex] != NULL)
                delete turbovertex[iMarker][iSpan][iVertex];
            delete [] turbovertex[iMarker][iSpan];
          }
        }
        delete [] turbovertex[iMarker];
      }
    }
    delete [] turbovertex;
  }
  if (AverageTurboNormal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Marker_All_SendRecv[iMarker] == 0 && AverageTurboNormal[iMarker] != NULL) {
        for (iSpan= 0; iSpan < nSpanSectionsByMarker[iMarker]+1; iSpan++)
          delete [] AverageTurboNormal[iMarker][iSpan];
        delete [] AverageTurboNormal[iMarker];
      }
    }
    delete [] AverageTurboNormal;
  }
  if (AverageNormal != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Marker_All_SendRecv[iMarker] == 0 && AverageNormal[iMarker] != NULL) {
        for (iSpan= 0; iSpan < nSpanSectionsByMarker[iMarker]+1; iSpan++)
          delete [] AverageNormal[iMarker][iSpan];
        delete [] AverageNormal[iMarker];
      }
    }
    delete [] AverageNormal;
  }
  if (AverageGridVel != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (Marker_All_SendRecv[iMarker] == 0 && AverageGridVel[iMarker] != NULL) {
        for (iSpan= 0; iSpan < nSpanSectionsByMarker[iMarker]+1; iSpan++)
          delete [] AverageGridVel[iMarker][iSpan];
        delete [] AverageGridVel[iMarker];
      }
    }
    delete [] AverageGridVel;
  }

  if (AverageTangGridVel != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && AverageTangGridVel[iMarker] != NULL)
        delete [] AverageTangGridVel[iMarker];
    delete [] AverageTangGridVel;
  }
  if (SpanArea != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && SpanArea[iMarker] != NULL)
        delete [] SpanArea[iMarker];
    delete [] SpanArea;
  }
  if (TurboRadius != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && TurboRadius[iMarker] != NULL)
        delete [] TurboRadius[iMarker];
    delete [] TurboRadius;
  }
  if (MaxAngularCoord != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && MaxAngularCoord[iMarker] != NULL)
        delete [] MaxAngularCoord[iMarker];
    delete [] MaxAngularCoord;
  }
  if (MinAngularCoord != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && MinAngularCoord[iMarker] != NULL)
        delete [] MinAngularCoord[iMarker];
    delete [] MinAngularCoord;
  }
  if (MinRelAngularCoord != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && MinRelAngularCoord[iMarker] != NULL)
        delete [] MinRelAngularCoord[iMarker];
    delete [] MinRelAngularCoord;
  }

  if (nSpanWiseSections != NULL) delete [] nSpanWiseSections;
  if (nSpanSectionsByMarker != NULL) delete [] nSpanSectionsByMarker;
  if (SpanWiseValue != NULL) {
    for (iMarker = 0; iMarker < 2; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && SpanWiseValue[iMarker] != NULL)
        delete [] SpanWiseValue[iMarker];
    delete [] SpanWiseValue;
  }
  if (nVertexSpan != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && nVertexSpan[iMarker] != NULL)
        delete [] nVertexSpan[iMarker];
    delete [] nVertexSpan;
  }
  if (nTotVertexSpan != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      if (Marker_All_SendRecv[iMarker] == 0 && nTotVertexSpan[iMarker] != NULL)
        delete [] nTotVertexSpan[iMarker];
    delete [] nTotVertexSpan;
  }

}

void CPhysicalGeometry::SetGlobal_to_Local_Point(void) {
  Global_to_Local_Point.clear();
  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {
    Global_to_Local_Point[node[iPoint]->GetGlobalIndex()] = iPoint;
  }
}

void CPhysicalGeometry::DistributeColoring(CConfig *config,
                                           CGeometry *geometry) {

  /*--- To start, each linear partition carries the color only for the
   owned nodes (nPoint), but we have repeated elems on each linear partition.
   We need to complete the coloring information such that the repeated
   points on each rank also have their color values. ---*/

  unsigned short iNode, jNode;
  unsigned long iPoint, iNeighbor, jPoint, iElem, iProcessor;

  map<unsigned long, unsigned long> Point_Map;
  map<unsigned long, unsigned long>::iterator MI;

  vector<unsigned long>::iterator it;

  SU2_MPI::Request *colorSendReq = NULL, *idSendReq = NULL;
  SU2_MPI::Request *colorRecvReq = NULL, *idRecvReq = NULL;
  int iProc, iSend, iRecv, myStart, myFinal;

  /*--- Get a linear partitioner to track the partition counts. ---*/

  CLinearPartitioner pointPartitioner(geometry->GetGlobal_nPoint(),0);

  /*--- First, create a complete map of the points on this rank (excluding
   repeats) and their neighbors so that we can efficiently loop through the
   points and decide how to distribute the colors. ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
      iPoint = geometry->elem[iElem]->GetNode(iNode);
      Point_Map[iPoint] = iPoint;
    }
  }

  /*--- Error check to ensure that the number of points found for this
   rank matches the number in the mesh file (in serial). ---*/

  if ((size == SINGLE_NODE) && (Point_Map.size() < geometry->GetnPoint())) {
    SU2_MPI::Error( string("Mismatch between NPOIN and number of points")
                   +string(" listed in mesh file.\n")
                   +string("Please check the mesh file for correctness.\n"),
                   CURRENT_FUNCTION);
  }

  /*--- Create a global to local mapping that includes the unowned points. ---*/

  map<unsigned long, unsigned long> Global2Local;
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
  }

  /*--- Find extra points that carry an index higher than nPoint. ---*/

  jPoint = geometry->GetnPoint();
  for (MI = Point_Map.begin(); MI != Point_Map.end(); MI++) {
    iPoint = MI->first;
    if ((Point_Map[iPoint] <  pointPartitioner.GetFirstIndexOnRank(rank)) ||
        (Point_Map[iPoint] >= pointPartitioner.GetLastIndexOnRank(rank))){
      Global2Local[Point_Map[iPoint]] = jPoint;
      jPoint++;
    }
  }

  /*--- Now create the neighbor list for each owned node (self-inclusive). ---*/

  Neighbors.clear();
  Neighbors.resize(Point_Map.size());
  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    for (iNode = 0; iNode < geometry->elem[iElem]->GetnNodes(); iNode++) {
      iPoint = Global2Local[geometry->elem[iElem]->GetNode(iNode)];
      for (jNode = 0; jNode < geometry->elem[iElem]->GetnNodes(); jNode++) {
        jPoint = geometry->elem[iElem]->GetNode(jNode);
        Neighbors[iPoint].push_back(jPoint);
      }
    }
  }

  /*--- Post-process the neighbor lists. ---*/

  for (iPoint = 0; iPoint < Point_Map.size(); iPoint++) {
    sort(Neighbors[iPoint].begin(), Neighbors[iPoint].end());
    it = unique(Neighbors[iPoint].begin(), Neighbors[iPoint].end());
    Neighbors[iPoint].resize(it - Neighbors[iPoint].begin());
  }

  /*--- Prepare structures for communication. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (iProc = 0; iProc < size; iProc++) {
    nPoint_Send[iProc] = 0; nPoint_Recv[iProc] = 0; nPoint_Flag[iProc]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  /*--- Loop over the owned points and check all the neighbors for unowned
   points. The colors of all owned points will be communicated to any ranks
   that will require them, which is due to the repeated points/elements
   that were needed to perform the coloring. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iNeighbor = 0; iNeighbor < Neighbors[iPoint].size(); iNeighbor++) {

      /*--- Global ID of the neighbor ---*/

      jPoint = Neighbors[iPoint][iNeighbor];

      /*--- Search for the processor that owns this neighbor. ---*/

      iProcessor = pointPartitioner.GetRankContainingIndex(jPoint);

      /*--- If we have not visited this node yet, increment our
       number of points that must be sent to a particular proc. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many points it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Prepare to send colors. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (iProc = 0; iProc < size; iProc++) nPoint_Flag[iProc] = -1;

  for (iProc = 0; iProc < size; iProc++) {
    if ((iProc != rank) && (nPoint_Send[iProc+1] > 0)) nSends++;
    if ((iProc != rank) && (nPoint_Recv[iProc+1] > 0)) nRecvs++;

    nPoint_Send[iProc+1] += nPoint_Send[iProc];
    nPoint_Recv[iProc+1] += nPoint_Recv[iProc];
  }

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (iSend = 0; iSend < nPoint_Send[size]; iSend++) idSend[iSend] = 0;

  /*--- Allocate memory to hold the colors that we are sending. ---*/

  unsigned long *colorSend = new unsigned long[nPoint_Send[size]];
  for (iSend = 0; iSend < nPoint_Send[size]; iSend++) colorSend[iSend] = 0;

  /*--- Create an index variable to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++) index[iProc] = nPoint_Send[iProc];

  /*--- Now load up our buffers with the Global IDs and colors. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iNeighbor = 0; iNeighbor < Neighbors[iPoint].size(); iNeighbor++) {

      /*--- Global ID of the neighbor ---*/

      jPoint = Neighbors[iPoint][iNeighbor];

      /*--- Search for the processor that owns this neighbor ---*/

      iProcessor = pointPartitioner.GetRankContainingIndex(jPoint);

      /*--- If we have not visited this node yet, increment our
       counters and load up the global ID and color. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {

        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];

        /*--- Load the data values. ---*/

        idSend[nn]    = geometry->node[iPoint]->GetGlobalIndex();
        colorSend[nn] = geometry->node[iPoint]->GetColor();

        /*--- Increment the index by the message length ---*/

        index[iProcessor]++;

      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *colorRecv = new unsigned long[nPoint_Recv[size]];
  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++)
    colorRecv[iRecv] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++)
    idRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nSends > 0) {
    colorSendReq = new SU2_MPI::Request[nSends];
    idSendReq    = new SU2_MPI::Request[nSends];
  }
  if (nRecvs > 0) {
    colorRecvReq = new SU2_MPI::Request[nRecvs];
    idRecvReq    = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(colorSend, nPoint_Send, colorSendReq,
                   colorRecv, nPoint_Recv, colorRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(idSend, nPoint_Send, idSendReq,
                   idRecv, nPoint_Recv, idRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  iRecv   = nPoint_Recv[rank];
  myStart = nPoint_Send[rank];
  myFinal = nPoint_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    colorRecv[iRecv] = colorSend[iSend];
    idRecv[iRecv]    = idSend[iSend];
    iRecv++;
  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends, colorSendReq, nRecvs, colorRecvReq);
  CompleteCommsAll(nSends,    idSendReq, nRecvs,    idRecvReq);

  /*--- Store the complete color map for this rank in class data. Now,
   each rank has a color value for all owned nodes as well as any repeated
   grid points on the rank. Note that there may be repeats that are
   communicated in the routine above, but since we are storing in a map,
   it will simply overwrite the same entries. ---*/

  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++) {
    Color_List[idRecv[iRecv]] = colorRecv[iRecv];
  }

  /*--- Free temporary memory from communications ---*/

  if (colorSendReq != NULL) delete [] colorSendReq;
  if (idSendReq    != NULL) delete [] idSendReq;

  if (colorRecvReq != NULL) delete [] colorRecvReq;
  if (idRecvReq    != NULL) delete [] idRecvReq;

  delete [] colorSend;
  delete [] colorRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

}

void CPhysicalGeometry::DistributeVolumeConnectivity(CConfig *config,
                                                     CGeometry *geometry,
                                                     unsigned short Elem_Type) {

  unsigned short NODES_PER_ELEMENT = 0;

  unsigned long iProcessor;
  unsigned long iElem, iNode, jNode, nElem_Total = 0, Global_Index;
  unsigned long *Conn_Elem  = NULL;
  unsigned long *ID_Elems   = NULL;

  SU2_MPI::Request *connSendReq = NULL, *idSendReq = NULL;
  SU2_MPI::Request *connRecvReq = NULL, *idRecvReq = NULL;
  int iProc, iSend, iRecv, myStart, myFinal;

  /*--- Store the number of nodes per this element type. ---*/

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
      NODES_PER_ELEMENT = 0;
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  /*--- Prepare a mapping for local to global element index. ---*/

  map<unsigned long, unsigned long> Local2GlobalElem;
  map<unsigned long, unsigned long>::iterator MI;

  for (MI = geometry->Global_to_Local_Elem.begin();
       MI != geometry->Global_to_Local_Elem.end(); MI++) {
    Local2GlobalElem[MI->second] = MI->first;
  }

  /*--- We start with the connectivity distributed across all procs in a
   linear partitioning. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements distributed according to the ParMETIS coloring. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (iProc = 0; iProc < size; iProc++) {
    nElem_Send[iProc] = 0; nElem_Recv[iProc] = 0; nElem_Flag[iProc]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++ ) {
    if (geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++ ) {

        /*--- Get the index of the current point. ---*/

        Global_Index = geometry->elem[iElem]->GetNode(iNode);

        /*--- We have the color stored in a map for all local points. ---*/

        iProcessor = Color_List[Global_Index];

        /*--- If we have not visited this element yet, increment our
         number of elements that must be sent to a particular proc. ---*/

        if ((nElem_Flag[iProcessor] != (int)iElem)) {
          nElem_Flag[iProcessor] = (int)iElem;
          nElem_Send[iProcessor+1]++;
        }

      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (iProc = 0; iProc < size; iProc++) nElem_Flag[iProc] = -1;

  for (iProc = 0; iProc < size; iProc++) {
    if ((iProc != rank) && (nElem_Send[iProc+1] > 0)) nSends++;
    if ((iProc != rank) && (nElem_Recv[iProc+1] > 0)) nRecvs++;

    nElem_Send[iProc+1] += nElem_Send[iProc];
    nElem_Recv[iProc+1] += nElem_Recv[iProc];
  }

  /*--- Allocate memory to hold the connectivity and element IDs
   that we are sending. ---*/

  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (iSend = 0; iSend < NODES_PER_ELEMENT*nElem_Send[size]; iSend++)
    connSend[iSend] = 0;

  /*--- Allocate arrays for storing element global index. ---*/

  unsigned long *idSend = new unsigned long[nElem_Send[size]];
  for (iSend = 0; iSend < nElem_Send[size]; iSend++) idSend[iSend] = 0;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    index[iProc] = NODES_PER_ELEMENT*nElem_Send[iProc];

  unsigned long *idIndex = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    idIndex[iProc] = nElem_Send[iProc];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
    if (geometry->elem[iElem]->GetVTK_Type() == Elem_Type) {
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++ ) {

        /*--- Get the index of the current point. ---*/

        Global_Index = geometry->elem[iElem]->GetNode(iNode);

        /*--- We have the color stored in a map for all local points. ---*/

        iProcessor = Color_List[Global_Index];

        /*--- Load connectivity and IDs into the buffer for sending ---*/

        if (nElem_Flag[iProcessor] != (int)iElem) {

          nElem_Flag[iProcessor] = (int)iElem;
          unsigned long nn = index[iProcessor];
          unsigned long mm = idIndex[iProcessor];

          /*--- Load the connectivity values. Note that elements are already
          stored directly based on their global index for the nodes.---*/

          for (jNode = 0; jNode < NODES_PER_ELEMENT; jNode++) {
            connSend[nn] = geometry->elem[iElem]->GetNode(jNode); nn++;
          }

          /*--- Global ID for this element. ---*/

          idSend[mm] = Local2GlobalElem[iElem];

          /*--- Increment the index by the message length ---*/

          index[iProcessor] += NODES_PER_ELEMENT;
          idIndex[iProcessor]++;

        }
      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] idIndex;

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (iRecv = 0; iRecv < NODES_PER_ELEMENT*nElem_Recv[size]; iRecv++)
    connRecv[iRecv] = 0;

  unsigned long *idRecv = new unsigned long[nElem_Recv[size]];
  for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) idRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nSends > 0) {
    connSendReq = new SU2_MPI::Request[nSends];
    idSendReq   = new SU2_MPI::Request[nSends];
  }
  if (nRecvs > 0) {
    connRecvReq = new SU2_MPI::Request[nRecvs];
    idRecvReq   = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(connSend, nElem_Send, connSendReq,
                   connRecv, nElem_Recv, connRecvReq,
                   NODES_PER_ELEMENT, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(idSend, nElem_Send, idSendReq,
                   idRecv, nElem_Recv, idRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  iRecv   = NODES_PER_ELEMENT*nElem_Recv[rank];
  myStart = NODES_PER_ELEMENT*nElem_Send[rank];
  myFinal = NODES_PER_ELEMENT*nElem_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    connRecv[iRecv] = connSend[iSend];
    iRecv++;
  }

  iRecv   = nElem_Recv[rank];
  myStart = nElem_Send[rank];
  myFinal = nElem_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    idRecv[iRecv] = idSend[iSend];
    iRecv++;
  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends, connSendReq, nRecvs, connRecvReq);
  CompleteCommsAll(nSends,   idSendReq, nRecvs,   idRecvReq);

  /*--- Store the connectivity for this rank in the proper structure
   It will be loaded into the geometry objects in a later step. ---*/

  if (nElem_Recv[size] > 0) {
    Conn_Elem = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
    int count = 0; nElem_Total = 0;
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      nElem_Total++;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        Conn_Elem[count] = connRecv[iRecv*NODES_PER_ELEMENT+iNode];
        count++;
      }
    }
  }

  /*--- Store the global element IDs too. ---*/

  if (nElem_Recv[size] > 0) {
    ID_Elems = new unsigned long[nElem_Recv[size]];
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      ID_Elems[iRecv] = idRecv[iRecv];
    }
  }

  /*--- Store the particular element count, IDs, & conn. in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case TRIANGLE:
      nLocal_Tria = nElem_Total;
      if (nLocal_Tria > 0) {
        Conn_Tria = Conn_Elem;
        ID_Tria   = ID_Elems;
      }
      break;
    case QUADRILATERAL:
      nLocal_Quad = nElem_Total;
      if (nLocal_Quad > 0) {
        Conn_Quad = Conn_Elem;
        ID_Quad   = ID_Elems;
      }
      break;
    case TETRAHEDRON:
      nLocal_Tetr = nElem_Total;
      if (nLocal_Tetr > 0) {
        Conn_Tetr = Conn_Elem;
        ID_Tetr   = ID_Elems;
      }
      break;
    case HEXAHEDRON:
      nLocal_Hexa = nElem_Total;
      if (nLocal_Hexa > 0) {
        Conn_Hexa = Conn_Elem;
        ID_Hexa   = ID_Elems;
      }
      break;
    case PRISM:
      nLocal_Pris = nElem_Total;
      if (nLocal_Pris > 0) {
        Conn_Pris = Conn_Elem;
        ID_Pris   = ID_Elems;
      }
      break;
    case PYRAMID:
      nLocal_Pyra = nElem_Total;
      if (nLocal_Pyra > 0) {
        Conn_Pyra = Conn_Elem;
        ID_Pyra   = ID_Elems;
      }
      break;
    default:
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  Local2GlobalElem.clear();

  if (connSendReq != NULL) delete [] connSendReq;
  if (idSendReq   != NULL) delete [] idSendReq;

  if (connRecvReq != NULL) delete [] connRecvReq;
  if (idRecvReq   != NULL) delete [] idRecvReq;

  delete [] connSend;
  delete [] connRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;

}

void CPhysicalGeometry::DistributePoints(CConfig *config, CGeometry *geometry) {

  /*--- We now know all of the coloring for our local points and neighbors.
   From this, we can communicate the owned nodes in our linear partitioning
   to all other ranks, including coordinates and coloring info, so that the
   receivers will be able to sort the data. ---*/

  unsigned short iDim;
  unsigned long iPoint, iNeighbor, jPoint, iProcessor;
  vector<unsigned long>::iterator it;

  SU2_MPI::Request *colorSendReq = NULL, *idSendReq = NULL, *coordSendReq = NULL;
  SU2_MPI::Request *colorRecvReq = NULL, *idRecvReq = NULL, *coordRecvReq = NULL;
  int iProc, iSend, iRecv, myStart, myFinal;

  /*--- Prepare structures for communication. ---*/

  int *nPoint_Send = new int[size+1]; nPoint_Send[0] = 0;
  int *nPoint_Recv = new int[size+1]; nPoint_Recv[0] = 0;
  int *nPoint_Flag = new int[size];

  for (iProc = 0; iProc < size; iProc++) {
    nPoint_Send[iProc] = 0; nPoint_Recv[iProc] = 0; nPoint_Flag[iProc]= -1;
  }
  nPoint_Send[size] = 0; nPoint_Recv[size] = 0;

  /*--- Loop over the owned points and check all the neighbors for unowned
   points. The colors of all owned points will be communicated to any ranks
   that will require them, which is due to the repeated points/elements
   that were needed to perform the coloring. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iNeighbor = 0; iNeighbor < Neighbors[iPoint].size(); iNeighbor++) {

      /*--- Global ID of the neighbor ---*/

      jPoint = Neighbors[iPoint][iNeighbor];

      /*--- We have the color stored in a map for all local points. ---*/

      iProcessor = Color_List[jPoint];

      /*--- If we have not visited this node yet, increment our
       number of points that must be sent to a particular proc. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {
        nPoint_Flag[iProcessor] = (int)iPoint;
        nPoint_Send[iProcessor+1]++;
      }
    }
  }

  /*--- Communicate the number of nodes to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many points it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nPoint_Send[1]), 1, MPI_INT,
                    &(nPoint_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Prepare to send colors, ids, and coords. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (iProc = 0; iProc < size; iProc++) nPoint_Flag[iProc] = -1;

  for (iProc = 0; iProc < size; iProc++) {
    if ((iProc != rank) && (nPoint_Send[iProc+1] > 0)) nSends++;
    if ((iProc != rank) && (nPoint_Recv[iProc+1] > 0)) nRecvs++;

    nPoint_Send[iProc+1] += nPoint_Send[iProc];
    nPoint_Recv[iProc+1] += nPoint_Recv[iProc];
  }

  /*--- Allocate arrays for sending the global ID. ---*/

  unsigned long *idSend = new unsigned long[nPoint_Send[size]];
  for (iSend = 0; iSend < nPoint_Send[size]; iSend++) idSend[iSend] = 0;

  /*--- Allocate memory to hold the colors that we are sending. ---*/

  unsigned long *colorSend = new unsigned long[nPoint_Send[size]];
  for (iSend = 0; iSend < nPoint_Send[size]; iSend++) colorSend[iSend] = 0;

  /*--- Allocate memory to hold the coordinates that we are sending. ---*/

  su2double *coordSend = NULL;
  coordSend = new su2double[nDim*nPoint_Send[size]];
  for (iSend = 0; iSend < nDim*nPoint_Send[size]; iSend++)
    coordSend[iSend] = 0;

  /*--- Create index variables to keep track of our index
   positions as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    index[iProc] = nPoint_Send[iProc];

  unsigned long *coordIndex = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    coordIndex[iProc] = nDim*nPoint_Send[iProc];

  /*--- Now load up our buffers with the colors, ids, and coords. ---*/

  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++) {
    for (iNeighbor = 0; iNeighbor < Neighbors[iPoint].size(); iNeighbor++) {

      /*--- Global ID of the neighbor ---*/

      jPoint = Neighbors[iPoint][iNeighbor];

      /*--- We have the color stored in a map for all local points. ---*/

      iProcessor = Color_List[jPoint];

      /*--- If we have not visited this node yet, increment our
       counters and load up the colors, ids, and coords. ---*/

      if (nPoint_Flag[iProcessor] != (int)iPoint) {

        nPoint_Flag[iProcessor] = (int)iPoint;
        unsigned long nn = index[iProcessor];

        /*--- Load the global ID, color, and coordinate values. ---*/

        idSend[nn]    = geometry->node[iPoint]->GetGlobalIndex();
        colorSend[nn] = geometry->node[iPoint]->GetColor();

        nn = coordIndex[iProcessor];
        for (iDim  = 0; iDim < nDim; iDim++) {
          coordSend[nn] = geometry->node[iPoint]->GetCoord(iDim); nn++;
        }

        /*--- Increment the index by the message length ---*/

        coordIndex[iProcessor] += nDim;
        index[iProcessor]++;

      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] coordIndex;

  /*--- Allocate the memory that we need for receiving the
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *colorRecv = new unsigned long[nPoint_Recv[size]];
  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++)
    colorRecv[iRecv] = 0;

  unsigned long *idRecv = new unsigned long[nPoint_Recv[size]];
  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++)
    idRecv[iRecv] = 0;

  su2double *coordRecv = NULL;
  coordRecv = new su2double[nDim*nPoint_Recv[size]];
  for (iRecv = 0; iRecv < nDim*nPoint_Recv[size]; iRecv++)
    coordRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nSends > 0) {
    colorSendReq = new SU2_MPI::Request[nSends];
    idSendReq    = new SU2_MPI::Request[nSends];
    coordSendReq = new SU2_MPI::Request[nSends];

  }
  if (nRecvs > 0) {
    colorRecvReq = new SU2_MPI::Request[nRecvs];
    idRecvReq    = new SU2_MPI::Request[nRecvs];
    coordRecvReq = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(colorSend, nPoint_Send, colorSendReq,
                   colorRecv, nPoint_Recv, colorRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(idSend, nPoint_Send, idSendReq,
                   idRecv, nPoint_Recv, idRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(coordSend, nPoint_Send, coordSendReq,
                   coordRecv, nPoint_Recv, coordRecvReq,
                   nDim, COMM_TYPE_DOUBLE);

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  iRecv   = nPoint_Recv[rank];
  myStart = nPoint_Send[rank];
  myFinal = nPoint_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    colorRecv[iRecv] = colorSend[iSend];
    idRecv[iRecv]    = idSend[iSend];
    iRecv++;
  }

  iRecv   = nDim*nPoint_Recv[rank];
  myStart = nDim*nPoint_Send[rank];
  myFinal = nDim*nPoint_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    coordRecv[iRecv] = coordSend[iSend];
    iRecv++;
  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends, colorSendReq, nRecvs, colorRecvReq);
  CompleteCommsAll(nSends,    idSendReq, nRecvs,    idRecvReq);
  CompleteCommsAll(nSends, coordSendReq, nRecvs, coordRecvReq);

  /*--- Store the total number of local points my rank has for
   the current section after completing the communications. ---*/

  nLocal_Point = nPoint_Recv[size];

  /*--- Store the proper local IDs, colors, and coordinates. We will load
   all of this information into our geometry classes in a later step. ---*/

  Local_Points = new unsigned long[nPoint_Recv[size]];
  Local_Colors = new unsigned long[nPoint_Recv[size]];
  Local_Coords = new su2double[nDim*nPoint_Recv[size]];

  nLocal_PointDomain = 0; nLocal_PointGhost = 0;
  for (iRecv = 0; iRecv < nPoint_Recv[size]; iRecv++) {
    Local_Points[iRecv] = idRecv[iRecv];
    Local_Colors[iRecv] = colorRecv[iRecv];
    for (iDim = 0; iDim < nDim; iDim++)
      Local_Coords[iRecv*nDim+iDim] = coordRecv[iRecv*nDim+iDim];
    if (Local_Colors[iRecv] == (unsigned long)rank) nLocal_PointDomain++;
    else nLocal_PointGhost++;
  }

  /*--- Free temporary memory from communications ---*/

  if (colorSendReq != NULL) delete [] colorSendReq;
  if (idSendReq    != NULL) delete [] idSendReq;
  if (coordSendReq != NULL) delete [] coordSendReq;

  if (colorRecvReq != NULL) delete [] colorRecvReq;
  if (idRecvReq    != NULL) delete [] idRecvReq;
  if (coordRecvReq != NULL) delete [] coordRecvReq;

  delete [] colorSend;
  delete [] colorRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] coordSend;
  delete [] coordRecv;
  delete [] nPoint_Recv;
  delete [] nPoint_Send;
  delete [] nPoint_Flag;

}

void CPhysicalGeometry::PartitionSurfaceConnectivity(CConfig *config,
                                                     CGeometry *geometry,
                                                     unsigned short Elem_Type) {

  /*--- We begin with all marker information residing on the master rank,
   as the master currently stores all marker info when reading the grid.
   We first check and communicate basic information that each rank will
   need to hold its portion of the linearly partitioned markers. In a
   later step, we will distribute the markers according to the ParMETIS
   coloring. This intermediate step is necessary since we already have the
   correct coloring distributed by the linear partitions, which we would
   like to reuse when partitioning the markers. Plus, the markers should
   truly be linearly partitioned upon reading the mesh, which we will
   change eventually. ---*/

  unsigned short NODES_PER_ELEMENT = 0;

  unsigned long iMarker, iProcessor, iElem, iNode, jNode;
  unsigned long nElem_Total = 0, Global_Index, Global_Elem_Index;

  unsigned long *Conn_Elem      = NULL;
  unsigned long *Linear_Markers = NULL;
  unsigned long *ID_SurfElem    = NULL;

  SU2_MPI::Request *connSendReq = NULL, *markerSendReq = NULL, *idSendReq = NULL;
  SU2_MPI::Request *connRecvReq = NULL, *markerRecvReq = NULL, *idRecvReq = NULL;
  int iProc, iSend, iRecv, myStart, myFinal;

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case LINE:
      NODES_PER_ELEMENT = N_POINTS_LINE;
      break;
    case TRIANGLE:
      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      break;
    case QUADRILATERAL:
      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (iProc = 0; iProc < size; iProc++) {
    nElem_Send[iProc] = 0; nElem_Recv[iProc] = 0; nElem_Flag[iProc]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  /*--- We know that the master owns all of the info and will be the only
   rank sending anything, although all ranks might receive something. ---*/

  if (rank == MASTER_NODE) {
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {

      /*--- Reset the flag in between markers, just to ensure that we
       don't miss some elements on different markers with the same local
       index. ---*/

      for (iProc = 0; iProc < size; iProc++) nElem_Flag[iProc]= -1;

      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {

        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {

          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++ ) {

            /*--- Get the index of the current point (stored as global). ---*/

            Global_Index = geometry->bound[iMarker][iElem]->GetNode(iNode);

            /*--- Search for the processor that owns this point ---*/

            iProcessor = GetLinearPartition(Global_Index);

            /*--- If we have not visited this element yet, increment our
             number of elements that must be sent to a particular proc. ---*/

            if ((nElem_Flag[iProcessor] != (int)iElem)) {
              nElem_Flag[iProcessor] = (int)iElem;
              nElem_Send[iProcessor+1]++;
            }
          }
        }
      }
    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Scatter(&(nElem_Send[1]), 1, MPI_INT,
                   &(nElem_Recv[1]), 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (iProc = 0; iProc < size; iProc++) nElem_Flag[iProc] = -1;

  for (iProc = 0; iProc < size; iProc++) {
    if ((iProc != rank) && (nElem_Send[iProc+1] > 0)) nSends++;
    if ((iProc != rank) && (nElem_Recv[iProc+1] > 0)) nRecvs++;

    nElem_Send[iProc+1] += nElem_Send[iProc];
    nElem_Recv[iProc+1] += nElem_Recv[iProc];
  }

  /*--- Allocate memory to hold the connectivity that we are sending. ---*/

  unsigned long *connSend   = NULL;
  unsigned long *markerSend = NULL;
  unsigned long *idSend     = NULL;

  if (rank == MASTER_NODE) {

    connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
    for (iSend = 0; iSend < NODES_PER_ELEMENT*nElem_Send[size]; iSend++)
      connSend[iSend] = 0;

    markerSend = new unsigned long[nElem_Send[size]];
    for (iSend = 0; iSend < nElem_Send[size]; iSend++)
      markerSend[iSend] = 0;

    idSend = new unsigned long[nElem_Send[size]];
    for (iSend = 0; iSend < nElem_Send[size]; iSend++)
      idSend[iSend] = 0;

    /*--- Create an index variable to keep track of our index
     position as we load up the send buffer. ---*/

    unsigned long *index = new unsigned long[size];
    for (iProc = 0; iProc < size; iProc++)
      index[iProc] = NODES_PER_ELEMENT*nElem_Send[iProc];

    unsigned long *markerIndex = new unsigned long[size];
    for (iProc = 0; iProc < size; iProc++)
      markerIndex[iProc] = nElem_Send[iProc];

    /*--- Loop through our elements and load the elems and their
     additional data that we will send to the other procs. ---*/

    Global_Elem_Index = 0;
    for (iMarker = 0; iMarker < geometry->GetnMarker(); iMarker++) {

      /*--- Reset the flag in between markers, just to ensure that we
       don't miss some elements on different markers with the same local
       index. ---*/

      for (iProc = 0; iProc < size; iProc++) nElem_Flag[iProc]= -1;

      for (iElem = 0; iElem < geometry->GetnElem_Bound(iMarker); iElem++) {
        if (geometry->bound[iMarker][iElem]->GetVTK_Type() == Elem_Type) {
          for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++ ) {

            /*--- Get the index of the current point. ---*/

            Global_Index = geometry->bound[iMarker][iElem]->GetNode(iNode);

            /*--- Search for the processor that owns this point ---*/

            iProcessor = GetLinearPartition(Global_Index);

            /*--- Load connectivity into the buffer for sending ---*/

            if ((nElem_Flag[iProcessor] != (int)iElem)) {

              nElem_Flag[iProcessor] = (int)iElem;
              unsigned long nn = index[iProcessor];
              unsigned long mm = markerIndex[iProcessor];

              /*--- Load the connectivity values. ---*/

              for (jNode = 0; jNode < NODES_PER_ELEMENT; jNode++) {
                connSend[nn] = geometry->bound[iMarker][iElem]->GetNode(jNode);
                nn++;
              }

              /*--- Store the marker index and surface elem global ID ---*/

              markerSend[mm] = iMarker;
              idSend[mm]     = Global_Elem_Index;

              /*--- Increment the index by the message length ---*/

              index[iProcessor] += NODES_PER_ELEMENT;
              markerIndex[iProcessor]++;
            }

          }
        }

        Global_Elem_Index++;

      }
    }

    /*--- Free memory after loading up the send buffer. ---*/

    delete [] index;
    delete [] markerIndex;

  }

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (iRecv = 0; iRecv < NODES_PER_ELEMENT*nElem_Recv[size]; iRecv++)
    connRecv[iRecv] = 0;

  unsigned long *markerRecv = new unsigned long[nElem_Recv[size]];
  for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++)
    markerRecv[iRecv] = 0;

  unsigned long *idRecv = new unsigned long[nElem_Recv[size]];
  for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++)
    idRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nSends > 0) {
    connSendReq   = new SU2_MPI::Request[nSends];
    markerSendReq = new SU2_MPI::Request[nSends];
    idSendReq     = new SU2_MPI::Request[nSends];
  }
  if (nRecvs > 0) {
    connRecvReq   = new SU2_MPI::Request[nRecvs];
    markerRecvReq = new SU2_MPI::Request[nRecvs];
    idRecvReq     = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(connSend, nElem_Send, connSendReq,
                   connRecv, nElem_Recv, connRecvReq,
                   NODES_PER_ELEMENT, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(markerSend, nElem_Send, markerSendReq,
                   markerRecv, nElem_Recv, markerRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(idSend, nElem_Send, idSendReq,
                   idRecv, nElem_Recv, idRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  if (rank == MASTER_NODE) {

    iRecv   = NODES_PER_ELEMENT*nElem_Recv[rank];
    myStart = NODES_PER_ELEMENT*nElem_Send[rank];
    myFinal = NODES_PER_ELEMENT*nElem_Send[rank+1];
    for (iSend = myStart; iSend < myFinal; iSend++) {
      connRecv[iRecv] = connSend[iSend];
      iRecv++;
    }

    iRecv   = nElem_Recv[rank];
    myStart = nElem_Send[rank];
    myFinal = nElem_Send[rank+1];
    for (iSend = myStart; iSend < myFinal; iSend++) {
      markerRecv[iRecv] = markerSend[iSend];
      idRecv[iRecv]     = idSend[iSend];
      iRecv++;
    }

  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends,   connSendReq, nRecvs,   connRecvReq);
  CompleteCommsAll(nSends, markerSendReq, nRecvs, markerRecvReq);
  CompleteCommsAll(nSends,     idSendReq, nRecvs,     idRecvReq);

  /*--- Store the connectivity for this rank in the proper data
   structure before post-processing below. First, allocate
   appropriate amount of memory for this section. ---*/

  if (nElem_Recv[size] > 0) {
    Conn_Elem = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
    int count = 0; nElem_Total = 0;
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      nElem_Total++;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        Conn_Elem[count] = connRecv[iRecv*NODES_PER_ELEMENT+iNode];
        count++;
      }
    }
  }

  /*--- Store the global marker ID for each element. ---*/

  if (nElem_Recv[size] > 0) {
    Linear_Markers = new unsigned long[nElem_Recv[size]];
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      Linear_Markers[iRecv] = markerRecv[iRecv];
    }
  }

  /*--- Store the global surface elem ID for each element. ---*/

  if (nElem_Recv[size] > 0) {
    ID_SurfElem = new unsigned long[nElem_Recv[size]];
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      ID_SurfElem[iRecv] = idRecv[iRecv];
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case LINE:
      nLinear_Line = nElem_Total;
      if (nLinear_Line > 0) {
        Conn_Line_Linear    = Conn_Elem;
        ID_Line_Linear      = Linear_Markers;
        Elem_ID_Line_Linear = ID_SurfElem;
      }
      break;
    case TRIANGLE:
      nLinear_BoundTria = nElem_Total;
      if (nLinear_BoundTria > 0) {
        Conn_BoundTria_Linear    = Conn_Elem;
        ID_BoundTria_Linear      = Linear_Markers;
        Elem_ID_BoundTria_Linear = ID_SurfElem;
      }
      break;
    case QUADRILATERAL:
      nLinear_BoundQuad = nElem_Total;
      if (nLinear_BoundQuad > 0) {
        Conn_BoundQuad_Linear    = Conn_Elem;
        ID_BoundQuad_Linear      = Linear_Markers;
        Elem_ID_BoundQuad_Linear = ID_SurfElem;
      }
      break;
    default:
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  if (connSendReq   != NULL) delete [] connSendReq;
  if (markerSendReq != NULL) delete [] markerSendReq;
  if (idSendReq     != NULL) delete [] idSendReq;

  if (connRecvReq   != NULL) delete [] connRecvReq;
  if (markerRecvReq != NULL) delete [] markerRecvReq;
  if (idRecvReq     != NULL) delete [] idRecvReq;

  if (connSend   != NULL) delete [] connSend;
  if (markerSend != NULL) delete [] markerSend;
  if (idSend     != NULL) delete [] idSend;

  delete [] connRecv;
  delete [] markerRecv;
  delete [] idRecv;

  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;

}

void CPhysicalGeometry::DistributeSurfaceConnectivity(CConfig *config,
                                                      CGeometry *geometry,
                                                      unsigned short Elem_Type) {

  unsigned short NODES_PER_ELEMENT = 0;

  unsigned long iProcessor, NELEM = 0;
  unsigned long iElem, iNode, jNode, nElem_Total = 0, Global_Index;

  unsigned long *Conn_Linear        = NULL;
  unsigned long *Conn_Elem          = NULL;
  unsigned long *Linear_Markers     = NULL;
  unsigned long *ID_SurfElem_Linear = NULL;
  unsigned long *Local_Markers      = NULL;
  unsigned long *ID_SurfElem        = NULL;

  SU2_MPI::Request *connSendReq = NULL,*markerSendReq = NULL,*idSendReq = NULL;
  SU2_MPI::Request *connRecvReq = NULL,*markerRecvReq = NULL,*idRecvReq = NULL;
  int iProc, iSend, iRecv, myStart, myFinal;

  /*--- Store the local number of this element type and the number of nodes
   per this element type. In serial, this will be the total number of this
   element type in the entire mesh. In parallel, it is the number on only
   the current partition. ---*/

  switch (Elem_Type) {
    case LINE:
      NELEM              = nLinear_Line;
      NODES_PER_ELEMENT  = N_POINTS_LINE;
      Conn_Linear        = Conn_Line_Linear;
      Linear_Markers     = ID_Line_Linear;
      ID_SurfElem_Linear = Elem_ID_Line_Linear;
      break;
    case TRIANGLE:
      NELEM              = nLinear_BoundTria;
      NODES_PER_ELEMENT  = N_POINTS_TRIANGLE;
      Conn_Linear        = Conn_BoundTria_Linear;
      Linear_Markers     = ID_BoundTria_Linear;
      ID_SurfElem_Linear = Elem_ID_BoundTria_Linear;
      break;
    case QUADRILATERAL:
      NELEM              = nLinear_BoundQuad;
      NODES_PER_ELEMENT  = N_POINTS_QUADRILATERAL;
      Conn_Linear        = Conn_BoundQuad_Linear;
      Linear_Markers     = ID_BoundQuad_Linear;
      ID_SurfElem_Linear = Elem_ID_BoundQuad_Linear;
      break;
    default:
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  /*--- We start with the connectivity distributed across all procs in a
   linear partitioning. We need to loop through our local partition
   and decide how many elements we must send to each other rank in order to
   have all elements distributed according to the ParMETIS coloring. ---*/

  int *nElem_Send = new int[size+1]; nElem_Send[0] = 0;
  int *nElem_Recv = new int[size+1]; nElem_Recv[0] = 0;
  int *nElem_Flag = new int[size];

  for (iProc = 0; iProc < size; iProc++) {
    nElem_Send[iProc] = 0; nElem_Recv[iProc] = 0; nElem_Flag[iProc]= -1;
  }
  nElem_Send[size] = 0; nElem_Recv[size] = 0;

  for (iElem = 0; iElem < NELEM; iElem++) {
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Get the index of the current point. ---*/

      Global_Index = Conn_Linear[iElem*NODES_PER_ELEMENT+iNode];

      /*--- We have the color stored in a map for all local points. ---*/

      iProcessor = Color_List[Global_Index];

      /*--- If we have not visited this element yet, increment our
       number of elements that must be sent to a particular proc. ---*/

      if ((nElem_Flag[iProcessor] != (int)iElem)) {
        nElem_Flag[iProcessor] = (int)iElem;
        nElem_Send[iProcessor+1]++;
      }

    }
  }

  /*--- Communicate the number of cells to be sent/recv'd amongst
   all processors. After this communication, each proc knows how
   many cells it will receive from each other processor. ---*/

  SU2_MPI::Alltoall(&(nElem_Send[1]), 1, MPI_INT,
                    &(nElem_Recv[1]), 1, MPI_INT, MPI_COMM_WORLD);

  /*--- Prepare to send connectivities. First check how many
   messages we will be sending and receiving. Here we also put
   the counters into cumulative storage format to make the
   communications simpler. ---*/

  int nSends = 0, nRecvs = 0;
  for (iProc = 0; iProc < size; iProc++) nElem_Flag[iProc] = -1;

  for (iProc = 0; iProc < size; iProc++) {
    if ((iProc != rank) && (nElem_Send[iProc+1] > 0)) nSends++;
    if ((iProc != rank) && (nElem_Recv[iProc+1] > 0)) nRecvs++;

    nElem_Send[iProc+1] += nElem_Send[iProc];
    nElem_Recv[iProc+1] += nElem_Recv[iProc];
  }

  /*--- Allocate memory to hold the connectivity that we are
   sending. ---*/

  unsigned long *connSend = NULL;
  connSend = new unsigned long[NODES_PER_ELEMENT*nElem_Send[size]];
  for (iSend = 0; iSend < NODES_PER_ELEMENT*nElem_Send[size]; iSend++)
    connSend[iSend] = 0;

  /*--- Allocate arrays for storing the marker global index. ---*/

  unsigned long *markerSend = new unsigned long[nElem_Send[size]];
  for (iSend = 0; iSend < nElem_Send[size]; iSend++) markerSend[iSend] = 0;

  unsigned long *idSend = new unsigned long[nElem_Send[size]];
  for (iSend = 0; iSend < nElem_Send[size]; iSend++) idSend[iSend] = 0;

  /*--- Create an index variable to keep track of our index
   position as we load up the send buffer. ---*/

  unsigned long *index = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    index[iProc] = NODES_PER_ELEMENT*nElem_Send[iProc];

  unsigned long *markerIndex = new unsigned long[size];
  for (iProc = 0; iProc < size; iProc++)
    markerIndex[iProc] = nElem_Send[iProc];

  /*--- Loop through our elements and load the elems and their
   additional data that we will send to the other procs. ---*/

  for (iElem = 0; iElem < NELEM; iElem++) {
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {

      /*--- Get the index of the current point. ---*/

      Global_Index = Conn_Linear[iElem*NODES_PER_ELEMENT+iNode];

      /*--- We have the color stored in a map for all local points. ---*/

      iProcessor = Color_List[Global_Index];

      /*--- If we have not visited this element yet, load up the data
       for sending. ---*/

      if (nElem_Flag[iProcessor] != (int)iElem) {

        nElem_Flag[iProcessor] = (int)iElem;
        unsigned long nn = index[iProcessor];
        unsigned long mm = markerIndex[iProcessor];

        /*--- Load the connectivity values. ---*/

        for (jNode = 0; jNode < NODES_PER_ELEMENT; jNode++) {

          /*--- Note that elements are already stored directly based on
           their global index for the nodes. ---*/

          connSend[nn] = Conn_Linear[iElem*NODES_PER_ELEMENT+jNode]; nn++;

        }

        /*--- Global marker ID for this element. ---*/

        markerSend[mm] = Linear_Markers[iElem];
        idSend[mm]     = ID_SurfElem_Linear[iElem];

        /*--- Increment the index by the message length ---*/

        index[iProcessor] += NODES_PER_ELEMENT;
        markerIndex[iProcessor]++;

      }
    }
  }

  /*--- Free memory after loading up the send buffer. ---*/

  delete [] index;
  delete [] markerIndex;

  /*--- Allocate the memory that we need for receiving the conn
   values and then cue up the non-blocking receives. Note that
   we do not include our own rank in the communications. We will
   directly copy our own data later. ---*/

  unsigned long *connRecv = NULL;
  connRecv = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
  for (iRecv = 0; iRecv < NODES_PER_ELEMENT*nElem_Recv[size]; iRecv++)
    connRecv[iRecv] = 0;

  unsigned long *markerRecv = new unsigned long[nElem_Recv[size]];
  for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) markerRecv[iRecv] = 0;

  unsigned long *idRecv = new unsigned long[nElem_Recv[size]];
  for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) idRecv[iRecv] = 0;

  /*--- Allocate memory for the MPI requests if we need to communicate. ---*/

  if (nSends > 0) {
    connSendReq   = new SU2_MPI::Request[nSends];
    markerSendReq = new SU2_MPI::Request[nSends];
    idSendReq     = new SU2_MPI::Request[nSends];
  }
  if (nRecvs > 0) {
    connRecvReq   = new SU2_MPI::Request[nRecvs];
    markerRecvReq = new SU2_MPI::Request[nRecvs];
    idRecvReq     = new SU2_MPI::Request[nRecvs];
  }

  /*--- Launch the non-blocking sends and receives. ---*/

  InitiateCommsAll(connSend, nElem_Send, connSendReq,
                   connRecv, nElem_Recv, connRecvReq,
                   NODES_PER_ELEMENT, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(markerSend, nElem_Send, markerSendReq,
                   markerRecv, nElem_Recv, markerRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  InitiateCommsAll(idSend, nElem_Send, idSendReq,
                   idRecv, nElem_Recv, idRecvReq,
                   1, COMM_TYPE_UNSIGNED_LONG);

  /*--- Copy my own rank's data into the recv buffer directly. ---*/

  iRecv   = NODES_PER_ELEMENT*nElem_Recv[rank];
  myStart = NODES_PER_ELEMENT*nElem_Send[rank];
  myFinal = NODES_PER_ELEMENT*nElem_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    connRecv[iRecv] = connSend[iSend];
    iRecv++;
  }

  iRecv   = nElem_Recv[rank];
  myStart = nElem_Send[rank];
  myFinal = nElem_Send[rank+1];
  for (iSend = myStart; iSend < myFinal; iSend++) {
    markerRecv[iRecv] = markerSend[iSend];
    idRecv[iRecv]     = idSend[iSend];
    iRecv++;
  }

  /*--- Complete the non-blocking communications. ---*/

  CompleteCommsAll(nSends,   connSendReq, nRecvs,   connRecvReq);
  CompleteCommsAll(nSends, markerSendReq, nRecvs, markerRecvReq);
  CompleteCommsAll(nSends,     idSendReq, nRecvs,     idRecvReq);

  /*--- Store the connectivity for this rank in the proper data
   structure. It will be loaded into the geometry objects in a later step. ---*/

  if (nElem_Recv[size] > 0) {
    Conn_Elem = new unsigned long[NODES_PER_ELEMENT*nElem_Recv[size]];
    int count = 0; nElem_Total = 0;
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      nElem_Total++;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        Conn_Elem[count] = connRecv[iRecv*NODES_PER_ELEMENT+iNode];
        count++;
      }
    }
  }

  /*--- Store the global marker IDs too. ---*/

  if (nElem_Recv[size] > 0) {
    Local_Markers = new unsigned long[nElem_Recv[size]];
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      Local_Markers[iRecv] = markerRecv[iRecv];
    }
  }

  /*--- Store the global surface elem IDs too. ---*/

  if (nElem_Recv[size] > 0) {
    ID_SurfElem = new unsigned long[nElem_Recv[size]];
    for (iRecv = 0; iRecv < nElem_Recv[size]; iRecv++) {
      ID_SurfElem[iRecv] = idRecv[iRecv];
    }
  }

  /*--- Store the particular global element count in the class data,
   and set the class data pointer to the connectivity array. ---*/

  switch (Elem_Type) {
    case LINE:
      nLocal_Line = nElem_Total;
      if (nLocal_Line > 0) {
        Conn_Line    = Conn_Elem;
        ID_Line      = Local_Markers;
        Elem_ID_Line = ID_SurfElem;
      }
      break;
    case TRIANGLE:
      nLocal_BoundTria = nElem_Total;
      if (nLocal_BoundTria > 0) {
        Conn_BoundTria    = Conn_Elem;
        ID_BoundTria      = Local_Markers;
        Elem_ID_BoundTria = ID_SurfElem;
      }
      break;
    case QUADRILATERAL:
      nLocal_BoundQuad = nElem_Total;
      if (nLocal_BoundQuad > 0) {
        Conn_BoundQuad    = Conn_Elem;
        ID_BoundQuad      = Local_Markers;
        Elem_ID_BoundQuad = ID_SurfElem;
      }
      break;
    default:
      SU2_MPI::Error("Unrecognized element type.", CURRENT_FUNCTION);
      break;
  }

  /*--- Free temporary memory from communications ---*/

  if (connSendReq   != NULL) delete [] connSendReq;
  if (markerSendReq != NULL) delete [] markerSendReq;
  if (idSendReq     != NULL) delete [] idSendReq;

  if (connRecvReq   != NULL) delete [] connRecvReq;
  if (markerRecvReq != NULL) delete [] markerRecvReq;
  if (idRecvReq     != NULL) delete [] idRecvReq;

  delete [] connSend;
  delete [] connRecv;
  delete [] markerSend;
  delete [] markerRecv;
  delete [] idSend;
  delete [] idRecv;
  delete [] nElem_Recv;
  delete [] nElem_Send;
  delete [] nElem_Flag;

}

void CPhysicalGeometry::DistributeMarkerTags(CConfig *config, CGeometry *geometry) {

  unsigned long iMarker, index, iChar;

  char str_buf[MAX_STRING_SIZE];

  /*--- The master node will communicate the entire list of marker tags
   (in global ordering) so that it will be simple for each rank to grab
   the string name for each marker. ---*/

  nMarker_Global = 0;
  if (rank == MASTER_NODE) nMarker_Global = config->GetnMarker_All();

  /*--- Broadcast the global number of markers in the mesh. ---*/

  SU2_MPI::Bcast(&nMarker_Global, 1, MPI_UNSIGNED_LONG,
                 MASTER_NODE, MPI_COMM_WORLD);

  char *mpi_str_buf = new char[nMarker_Global*MAX_STRING_SIZE]();
  if (rank == MASTER_NODE) {
    for (iMarker = 0; iMarker < nMarker_Global; iMarker++) {
      SPRINTF(&mpi_str_buf[iMarker*MAX_STRING_SIZE], "%s",
              config->GetMarker_All_TagBound(iMarker).c_str());
    }
  }

  /*--- Broadcast the string names of the variables. ---*/

  SU2_MPI::Bcast(mpi_str_buf, (int)nMarker_Global*MAX_STRING_SIZE, MPI_CHAR,
                 MASTER_NODE, MPI_COMM_WORLD);

  /*--- Now parse the string names and load into our marker tag vector.
   We also need to set the values of all markers into the config. ---*/

  for (iMarker = 0; iMarker < nMarker_Global; iMarker++) {
    index = iMarker*MAX_STRING_SIZE;
    for (iChar = 0; iChar < MAX_STRING_SIZE; iChar++) {
      str_buf[iChar] = mpi_str_buf[index + iChar];
    }
    Marker_Tags.push_back(str_buf);
    config->SetMarker_All_TagBound(iMarker,str_buf);
    config->SetMarker_All_SendRecv(iMarker,NO);
  }

  /*--- Free string buffer memory. ---*/

  delete [] mpi_str_buf;

}

void CPhysicalGeometry::LoadPoints(CConfig *config, CGeometry *geometry) {

  unsigned long iPoint, jPoint, iOwned, iPeriodic, iGhost;

  /*--- Create the basic point structures before storing the points. ---*/

  nPoint       = nLocal_Point;
  nPointDomain = nLocal_PointDomain;
  nPointNode   = nPoint;

  node = new CPoint*[nPoint];

  Local_to_Global_Point = new long[nPoint];

  /*--- Array initialization ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    Local_to_Global_Point[iPoint] = -1;
  }

  /*--- Set our counters correctly based on the number of owned and ghost
   nodes that we counted during the partitioning. ---*/

  jPoint    = 0;
  iOwned    = 0;
  iPeriodic = nLocal_PointDomain;
  iGhost    = nLocal_PointDomain + nLocal_PointPeriodic;

  /*--- Loop over all of the points that we have recv'd and store the
   coordinates, global index, and colors ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Set the starting point to the correct counter for this point. ---*/

    if (Local_Colors[iPoint] == (unsigned long)rank) {
      if (Local_Points[iPoint] < geometry->GetGlobal_nPointDomain())
        jPoint = iOwned;
      else jPoint = iPeriodic;
    } else {
      jPoint = iGhost;
    }

    /*--- Get the global index ---*/

    Local_to_Global_Point[jPoint] = Local_Points[iPoint];

    /*--- Allocating the Point object ---*/

    if ( nDim == 2 )
      node[jPoint] = new CPoint(Local_Coords[iPoint*nDim+0],
                                Local_Coords[iPoint*nDim+1],
                                Local_to_Global_Point[jPoint], config);
    if ( nDim == 3 )
      node[jPoint] = new CPoint(Local_Coords[iPoint*nDim+0],
                                Local_Coords[iPoint*nDim+1],
                                Local_Coords[iPoint*nDim+2],
                                Local_to_Global_Point[jPoint], config);

    /*--- Set the color ---*/

    node[jPoint]->SetColor(Local_Colors[iPoint]);

    /*--- Increment the correct counter before moving to the next point. ---*/

    if (Local_Colors[iPoint] == (unsigned long)rank) {
      if (Local_Points[iPoint] < geometry->GetGlobal_nPointDomain())
        iOwned++;
      else iPeriodic++;
    } else {
      iGhost++;
    }
  }

  /*--- Create the global to local mapping, which will be useful for loading
   the elements and boundaries in subsequent steps. ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    Global_to_Local_Point[Local_to_Global_Point[iPoint]] = iPoint;

  /*--- Set the value of Global_nPoint and Global_nPointDomain ---*/

  unsigned long Local_nPoint = nPoint;
  unsigned long Local_nPointDomain = nPointDomain;

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Local_nPoint, &Global_nPoint, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nPointDomain, &Global_nPointDomain, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nPoint = Local_nPoint;
  Global_nPointDomain = Local_nPointDomain;
#endif

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << Global_nPoint << " vertices including ghost points. " << endl;

}

void CPhysicalGeometry::LoadVolumeElements(CConfig *config, CGeometry *geometry) {

  unsigned short NODES_PER_ELEMENT;

  unsigned long iElem, jElem, kElem, iNode, Local_Elem, iGlobal_Index;
  unsigned long Local_Nodes[N_POINTS_HEXAHEDRON];

  unsigned long iElemTria = 0;
  unsigned long iElemQuad = 0;
  unsigned long iElemTetr = 0;
  unsigned long iElemHexa = 0;
  unsigned long iElemPris = 0;
  unsigned long iElemPyra = 0;

  unsigned long nTria, nQuad, nTetr, nHexa, nPris, nPyra;

  map<unsigned long, unsigned long> Tria_List;
  map<unsigned long, unsigned long> Quad_List;
  map<unsigned long, unsigned long> Tetr_List;
  map<unsigned long, unsigned long> Hexa_List;
  map<unsigned long, unsigned long> Pris_List;
  map<unsigned long, unsigned long> Pyra_List;
  map<unsigned long, unsigned long>::iterator it;

  /*--- It is possible that we have repeated elements during the previous
   communications, as we mostly focus on the grid points and their colors.
   First, loop through our local elements and build a mapping by simply
   overwriting the duplicate entries. ---*/

  jElem = 0;
  for (iElem=0; iElem < nLocal_Tria; iElem++) {
    Tria_List[ID_Tria[iElem]] = iElem;
  }
  nTria = Tria_List.size();

  jElem = 0;
  for (iElem=0; iElem < nLocal_Quad; iElem++) {
    Quad_List[ID_Quad[iElem]] = iElem;
  }
  nQuad = Quad_List.size();

  jElem = 0;
  for (iElem=0; iElem < nLocal_Tetr; iElem++) {
    Tetr_List[ID_Tetr[iElem]] = iElem;
  }
  nTetr = Tetr_List.size();

  jElem = 0;
  for (iElem=0; iElem < nLocal_Hexa; iElem++) {
    Hexa_List[ID_Hexa[iElem]] = iElem;
  }
  nHexa = Hexa_List.size();

  jElem = 0;
  for (iElem=0; iElem < nLocal_Pris; iElem++) {
    Pris_List[ID_Pris[iElem]] = iElem;
  }
  nPris = Pris_List.size();

  jElem = 0;
  for (iElem=0; iElem < nLocal_Pyra; iElem++) {
    Pyra_List[ID_Pyra[iElem]] = iElem;
  }
  nPyra = Pyra_List.size();

  /*--- Reduce the final count of non-repeated elements on this rank. ---*/

  Local_Elem = nTria + nQuad + nTetr + nHexa + nPris + nPyra;

  /*--- Create the basic structures for holding the grid elements. ---*/

  jElem = 0;
  nElem = Local_Elem;
  elem  = new CPrimalGrid*[nElem];

  /*--- Store the elements of each type in the proper containers. ---*/

  for (it = Tria_List.begin(); it != Tria_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Tria[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CTriangle(Local_Nodes[0],
                                Local_Nodes[1],
                                Local_Nodes[2], 2);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemTria++;

  }

  /*--- Free memory as we go. ---*/

  Tria_List.clear();

  for (it = Quad_List.begin(); it != Quad_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Quad[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CQuadrilateral(Local_Nodes[0],
                                     Local_Nodes[1],
                                     Local_Nodes[2],
                                     Local_Nodes[3], 2);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemQuad++;

  }

  /*--- Free memory as we go. ---*/

  Quad_List.clear();

  for (it = Tetr_List.begin(); it != Tetr_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_TETRAHEDRON;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Tetr[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CTetrahedron(Local_Nodes[0],
                                   Local_Nodes[1],
                                   Local_Nodes[2],
                                   Local_Nodes[3]);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemTetr++;

  }

  /*--- Free memory as we go. ---*/

  Tetr_List.clear();

  for (it = Hexa_List.begin(); it != Hexa_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_HEXAHEDRON;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Hexa[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CHexahedron(Local_Nodes[0],
                                  Local_Nodes[1],
                                  Local_Nodes[2],
                                  Local_Nodes[3],
                                  Local_Nodes[4],
                                  Local_Nodes[5],
                                  Local_Nodes[6],
                                  Local_Nodes[7]);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemHexa++;

  }

  /*--- Free memory as we go. ---*/

  Hexa_List.clear();

  for (it = Pris_List.begin(); it != Pris_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_PRISM;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Pris[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CPrism(Local_Nodes[0],
                             Local_Nodes[1],
                             Local_Nodes[2],
                             Local_Nodes[3],
                             Local_Nodes[4],
                             Local_Nodes[5]);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemPris++;

  }

  /*--- Free memory as we go. ---*/

  Pris_List.clear();

  for (it = Pyra_List.begin(); it != Pyra_List.end(); it++) {

    kElem = it->first;
    iElem = it->second;

    /*--- Transform the stored connectivity for this element from global
     to local values on this rank. ---*/

    NODES_PER_ELEMENT = N_POINTS_PYRAMID;
    for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
      iGlobal_Index      = Conn_Pyra[iElem*NODES_PER_ELEMENT+iNode];
      Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
    }

    /*--- Create the element object. ---*/

    elem[jElem] = new CPyramid(Local_Nodes[0],
                               Local_Nodes[1],
                               Local_Nodes[2],
                               Local_Nodes[3],
                               Local_Nodes[4]);

    elem[jElem]->SetGlobalIndex(kElem);

    /*--- Increment our local counters. ---*/

    jElem++; iElemPyra++;

  }

  /*--- Free memory as we go. ---*/

  Pyra_List.clear();

  /*--- Communicate the number of each element type to all processors. These
   values are important for merging and writing output later. ---*/

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Local_Elem, &Global_nElem, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nElem = nElem;
#endif

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << Global_nElem << " interior elements including halo cells. " << endl;

  /*--- Set the value of Global_nElemDomain (stored in the geometry
   container that is passed in). ---*/

  Global_nElemDomain = geometry->GetGlobal_nElemDomain();

  /*--- Store total number of each element type after incrementing the
   counters in the recv loop above (to make sure there aren't repeats). ---*/

  nelem_triangle = iElemTria;
  nelem_quad     = iElemQuad;
  nelem_tetra    = iElemTetr;
  nelem_hexa     = iElemHexa;
  nelem_prism    = iElemPris;
  nelem_pyramid  = iElemPyra;

#ifdef HAVE_MPI
  unsigned long Local_nElemTri     = nelem_triangle;
  unsigned long Local_nElemQuad    = nelem_quad;
  unsigned long Local_nElemTet     = nelem_tetra;
  unsigned long Local_nElemHex     = nelem_hexa;
  unsigned long Local_nElemPrism   = nelem_prism;
  unsigned long Local_nElemPyramid = nelem_pyramid;

  SU2_MPI::Allreduce(&Local_nElemTri, &Global_nelem_triangle, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemQuad, &Global_nelem_quad, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemTet, &Global_nelem_tetra, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemHex, &Global_nelem_hexa, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPrism, &Global_nelem_prism, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid, 1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nelem_triangle = nelem_triangle;
  Global_nelem_quad     = nelem_quad;
  Global_nelem_tetra    = nelem_tetra;
  Global_nelem_hexa     = nelem_hexa;
  Global_nelem_prism    = nelem_prism;
  Global_nelem_pyramid  = nelem_pyramid;
#endif

  /*--- Print information about the elements to the console ---*/

  if (rank == MASTER_NODE) {
    if (Global_nelem_triangle > 0)
      cout << Global_nelem_triangle << " triangles."      << endl;
    if (Global_nelem_quad     > 0)
      cout << Global_nelem_quad     << " quadrilaterals." << endl;
    if (Global_nelem_tetra    > 0)
      cout << Global_nelem_tetra    << " tetrahedra."     << endl;
    if (Global_nelem_hexa     > 0)
      cout << Global_nelem_hexa     << " hexahedra."      << endl;
    if (Global_nelem_prism    > 0)
      cout << Global_nelem_prism    << " prisms."         << endl;
    if (Global_nelem_pyramid  > 0)
      cout << Global_nelem_pyramid  << " pyramids."       << endl;
  }

}

void CPhysicalGeometry::LoadSurfaceElements(CConfig *config, CGeometry *geometry) {

  unsigned short NODES_PER_ELEMENT;
  unsigned short iNode, nMarker_Max = config->GetnMarker_Max();

  unsigned long iElem, iMarker, Global_Marker, iGlobal_Index;

  unsigned long iElem_Line = 0;
  unsigned long iElem_Tria = 0;
  unsigned long iElem_Quad = 0;

  unsigned long Local_Nodes[N_POINTS_HEXAHEDRON];

  vector<vector<unsigned long> > Line_List;
  vector<vector<unsigned long> > BoundTria_List;
  vector<vector<unsigned long> > BoundQuad_List;

  vector<unsigned long> Marker_Local;
  vector<unsigned long>::iterator it;

  /*--- Compute how many markers we have local to this rank by looping
   through the global marker numbers of each local surface element and
   counting the unique set. ---*/

  for (iElem = 0; iElem < nLocal_Line; iElem++) {
    if (find(Marker_Local.begin(), Marker_Local.end(),
             ID_Line[iElem]) == Marker_Local.end()) {
      Marker_Local.push_back(ID_Line[iElem]);
    }
  }

  for (iElem = 0; iElem < nLocal_BoundTria; iElem++) {
    if (find(Marker_Local.begin(), Marker_Local.end(),
             ID_BoundTria[iElem]) == Marker_Local.end()) {
      Marker_Local.push_back(ID_BoundTria[iElem]);
    }
  }

  for (iElem = 0; iElem < nLocal_BoundQuad; iElem++) {
    if (find(Marker_Local.begin(), Marker_Local.end(),
             ID_BoundQuad[iElem]) == Marker_Local.end()) {
      Marker_Local.push_back(ID_BoundQuad[iElem]);
    }
  }

  /*--- Create a mapping from global to local marker ID (and vice-versa). ---*/

  map<unsigned long, unsigned long> Marker_Global_to_Local;
  map<unsigned long, unsigned long> Marker_Local_to_Global;

  for (iMarker = 0; iMarker < Marker_Local.size(); iMarker++) {
    Marker_Global_to_Local[Marker_Local[iMarker]] = iMarker;
    Marker_Local_to_Global[iMarker] = Marker_Local[iMarker];
  }

  /*--- Set up our element counters on each marker so that we can avoid
   duplicating any elements from the previous communications. ---*/

  Line_List.resize(Marker_Local.size());
  BoundTria_List.resize(Marker_Local.size());
  BoundQuad_List.resize(Marker_Local.size());

  /*--- Count the number of elements on each marker and store in a
   vector by marker. ---*/

  vector<unsigned long> nElemBound_Local;
  nElemBound_Local.resize(Marker_Local.size());
  for (iMarker = 0; iMarker < Marker_Local.size(); iMarker++)
    nElemBound_Local[iMarker] = 0;

  for (iElem = 0; iElem < nLocal_Line; iElem++) {
    iMarker = Marker_Global_to_Local[ID_Line[iElem]];
    if (find(Line_List[iMarker].begin(), Line_List[iMarker].end(),
             Elem_ID_Line[iElem]) == Line_List[iMarker].end()) {
      nElemBound_Local[iMarker]++;
      Line_List[iMarker].push_back(Elem_ID_Line[iElem]);
    }
  }

  for (iElem = 0; iElem < nLocal_BoundTria; iElem++) {
    iMarker = Marker_Global_to_Local[ID_BoundTria[iElem]];
    if (find(BoundTria_List[iMarker].begin(), BoundTria_List[iMarker].end(),
             Elem_ID_BoundTria[iElem]) == BoundTria_List[iMarker].end()) {
      nElemBound_Local[iMarker]++;
      BoundTria_List[iMarker].push_back(Elem_ID_BoundTria[iElem]);
    }
  }

  for (iElem = 0; iElem < nLocal_BoundQuad; iElem++) {
    iMarker = Marker_Global_to_Local[ID_BoundQuad[iElem]];
    if (find(BoundQuad_List[iMarker].begin(), BoundQuad_List[iMarker].end(),
             Elem_ID_BoundQuad[iElem]) == BoundQuad_List[iMarker].end()) {
      nElemBound_Local[iMarker]++;
      BoundQuad_List[iMarker].push_back(Elem_ID_BoundQuad[iElem]);
    }
  }

  /*--- Create the domain structures for the boundaries. Initially, stick
   with nMarkerMax here, but come back and compute size we need. Same for
   OVERHEAD - this can precomputed. ---*/

  nMarker             = Marker_Local.size();
  nElem_Bound         = new unsigned long[nMarker_Max];
  Tag_to_Marker       = new string[nMarker_Max];
  Marker_All_SendRecv = new short[nMarker_Max];

  /*--- Allocate space for the elements on each marker ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nElem_Bound[iMarker] = nElemBound_Local[iMarker];

  bound = new CPrimalGrid**[nMarker+(OVERHEAD*size)];
  for (iMarker = 0; iMarker < nMarker+(OVERHEAD*size); iMarker++)
    bound[iMarker] = NULL;

  for (iMarker = 0; iMarker < nMarker; iMarker++)
    bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

  /*--- Initialize boundary element counters ---*/

  iElem_Line = 0;
  iElem_Tria = 0;
  iElem_Quad = 0;

  Line_List.clear();      Line_List.resize(Marker_Local.size());
  BoundTria_List.clear(); BoundTria_List.resize(Marker_Local.size());
  BoundQuad_List.clear(); BoundQuad_List.resize(Marker_Local.size());

  /*--- Reset our element counter on a marker-basis. ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++)
    nElemBound_Local[iMarker] = 0;

  /*--- Store the boundary element connectivity. Note here that we have
   communicated the global index values for the elements, so we need to
   convert this to the local index when instantiating the element. ---*/

  for (iElem = 0; iElem < nLocal_Line; iElem++) {

    iMarker = Marker_Global_to_Local[ID_Line[iElem]];

    /*--- Avoid duplicates on this marker. ---*/

    if (find(Line_List[iMarker].begin(), Line_List[iMarker].end(),
             Elem_ID_Line[iElem]) == Line_List[iMarker].end()) {

      /*--- Transform the stored connectivity for this element from global
       to local values on this rank. ---*/

      NODES_PER_ELEMENT = N_POINTS_LINE;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iGlobal_Index      = Conn_Line[iElem*NODES_PER_ELEMENT+iNode];
        Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
      }

      /*--- Create the geometry object for this element. ---*/

      bound[iMarker][nElemBound_Local[iMarker]] = new CLine(Local_Nodes[0],
                                                            Local_Nodes[1], 2);

      /*--- Increment our counters for this marker and element type. ---*/

      nElemBound_Local[iMarker]++; iElem_Line++;

      Line_List[iMarker].push_back(Elem_ID_Line[iElem]);

    }
  }

  for (iElem = 0; iElem < nLocal_BoundTria; iElem++) {

    iMarker = Marker_Global_to_Local[ID_BoundTria[iElem]];

    /*--- Avoid duplicates on this marker. ---*/

    if (find(BoundTria_List[iMarker].begin(), BoundTria_List[iMarker].end(),
             Elem_ID_BoundTria[iElem]) == BoundTria_List[iMarker].end()) {

      /*--- Transform the stored connectivity for this element from global
       to local values on this rank. ---*/

      NODES_PER_ELEMENT = N_POINTS_TRIANGLE;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iGlobal_Index      = Conn_BoundTria[iElem*NODES_PER_ELEMENT+iNode];
        Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
      }

      /*--- Create the geometry object for this element. ---*/

      bound[iMarker][nElemBound_Local[iMarker]] = new CTriangle(Local_Nodes[0],
                                                                Local_Nodes[1],
                                                                Local_Nodes[2], 3);

      /*--- Increment our counters for this marker and element type. ---*/

      nElemBound_Local[iMarker]++; iElem_Tria++;

      BoundTria_List[iMarker].push_back(Elem_ID_BoundTria[iElem]);

    }
  }

  for (iElem = 0; iElem < nLocal_BoundQuad; iElem++) {

    iMarker = Marker_Global_to_Local[ID_BoundQuad[iElem]];

    /*--- Avoid duplicates on this marker. ---*/

    if (find(BoundQuad_List[iMarker].begin(), BoundQuad_List[iMarker].end(),
             Elem_ID_BoundQuad[iElem]) == BoundQuad_List[iMarker].end()) {

      /*--- Transform the stored connectivity for this element from global
       to local values on this rank. ---*/

      NODES_PER_ELEMENT = N_POINTS_QUADRILATERAL;
      for (iNode = 0; iNode < NODES_PER_ELEMENT; iNode++) {
        iGlobal_Index      = Conn_BoundQuad[iElem*NODES_PER_ELEMENT+iNode];
        Local_Nodes[iNode] = Global_to_Local_Point[iGlobal_Index];
      }

      /*--- Create the geometry object for this element. ---*/

      bound[iMarker][nElemBound_Local[iMarker]] = new CQuadrilateral(Local_Nodes[0],
                                                                     Local_Nodes[1],
                                                                     Local_Nodes[2],
                                                                     Local_Nodes[3], 3);

      /*--- Increment our counters for this marker and element type. ---*/

      nElemBound_Local[iMarker]++; iElem_Quad++;

      BoundQuad_List[iMarker].push_back(Elem_ID_BoundQuad[iElem]);

    }
  }

  /*--- Store total number of each boundary element type ---*/

  nelem_edge_bound     = iElem_Line;
  nelem_triangle_bound = iElem_Tria;
  nelem_quad_bound     = iElem_Quad;

  /*--- Set some auxiliary information on a per-marker basis. ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    Global_Marker = Marker_Local_to_Global[iMarker];

    /*--- Now each domain has the right information ---*/

    string Grid_Marker = config->GetMarker_All_TagBound(Marker_Local_to_Global[iMarker]);
    short SendRecv     = config->GetMarker_All_SendRecv(Marker_Local_to_Global[iMarker]);

    Tag_to_Marker[iMarker] = Marker_Tags[Global_Marker];
    Marker_All_SendRecv[iMarker] = SendRecv;

    /*--- Set the marker tags correctly to match the values in config. ---*/

    config->SetMarker_All_TagBound(iMarker, Tag_to_Marker[iMarker]);
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);

  }

  /*--- Initialize pointers for turbomachinery computations  ---*/

  nSpanWiseSections       = new unsigned short[2] ();
  SpanWiseValue           = new su2double*[2] ();

  nSpanSectionsByMarker             = new unsigned short[nMarker] ();
  nVertexSpan                       = new long* [nMarker] ();
  nTotVertexSpan                    = new unsigned long* [nMarker] ();
  turbovertex                       = new CTurboVertex***[nMarker] ();
  AverageTurboNormal                = new su2double**[nMarker] ();
  AverageNormal                     = new su2double**[nMarker] ();
  AverageGridVel                    = new su2double**[nMarker] ();
  AverageTangGridVel                = new su2double*[nMarker] ();
  SpanArea                          = new su2double*[nMarker] ();
  TurboRadius                       = new su2double*[nMarker] ();
  MaxAngularCoord                   = new su2double*[nMarker] ();
  MinAngularCoord                   = new su2double*[nMarker] ();
  MinRelAngularCoord                = new su2double*[nMarker] ();

  /*--- Initialize pointers for turbomachinery performance computation  ---*/

  nTurboPerf     = config->GetnMarker_TurboPerformance();
  TangGridVelIn  = new su2double*[nTurboPerf] ();
  SpanAreaIn     = new su2double*[nTurboPerf] ();
  TurboRadiusIn  = new su2double*[nTurboPerf] ();
  TangGridVelOut = new su2double*[nTurboPerf] ();
  SpanAreaOut    = new su2double*[nTurboPerf] ();
  TurboRadiusOut = new su2double*[nTurboPerf] ();

}

void CPhysicalGeometry::InitiateCommsAll(void *bufSend,
                                         int *nElemSend,
                                         SU2_MPI::Request *sendReq,
                                         void *bufRecv,
                                         int *nElemRecv,
                                         SU2_MPI::Request *recvReq,
                                         unsigned short countPerElem,
                                         unsigned short commType) {

  /*--- Local variables ---*/

  int iMessage, iProc, offset, nElem, count, source, dest, tag;

  /*--- Launch the non-blocking recv's first. ---*/

  iMessage = 0;
  for (iProc = 0; iProc < size; iProc++) {

    /*--- Post recv's only if another proc is sending us data. We do
     not communicate with ourselves or post recv's for zero length
     messages to keep overhead down. ---*/

    if ((nElemRecv[iProc+1] > nElemRecv[iProc]) && (iProc != rank)) {

      /*--- Compute our location in the recv buffer. ---*/

      offset = countPerElem*nElemRecv[iProc];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to recv. ---*/

      nElem = nElemRecv[iProc+1] - nElemRecv[iProc];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerElem*nElem;

      /*--- Post non-blocking recv for this proc. ---*/

      source = iProc; tag = iProc + 1;

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Irecv(&(static_cast<su2double*>(bufRecv)[offset]),
                         count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_LONG:
          SU2_MPI::Irecv(&(static_cast<unsigned long*>(bufRecv)[offset]),
                         count, MPI_UNSIGNED_LONG, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_LONG:
          SU2_MPI::Irecv(&(static_cast<long*>(bufRecv)[offset]),
                         count, MPI_LONG, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Irecv(&(static_cast<unsigned short*>(bufRecv)[offset]),
                         count, MPI_UNSIGNED_SHORT, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_CHAR:
          SU2_MPI::Irecv(&(static_cast<char*>(bufRecv)[offset]),
                         count, MPI_CHAR, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_SHORT:
          SU2_MPI::Irecv(&(static_cast<short*>(bufRecv)[offset]),
                         count, MPI_SHORT, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        case COMM_TYPE_INT:
          SU2_MPI::Irecv(&(static_cast<int*>(bufRecv)[offset]),
                         count, MPI_INT, source, tag, MPI_COMM_WORLD,
                         &(recvReq[iMessage]));
          break;
        default:
          break;
      }

      /*--- Increment message counter. ---*/

      iMessage++;

    }
  }

  /*--- Launch the non-blocking sends next. ---*/

  iMessage = 0;
  for (iProc = 0; iProc < size; iProc++) {

    /*--- Post sends only if we are sending another proc data. We do
     not communicate with ourselves or post sends for zero length
     messages to keep overhead down. ---*/

    if ((nElemSend[iProc+1] > nElemSend[iProc]) && (iProc != rank)) {

      /*--- Compute our location in the send buffer. ---*/

      offset = countPerElem*nElemSend[iProc];

      /*--- Take advantage of cumulative storage format to get the number
       of elems that we need to send. ---*/

      nElem = nElemSend[iProc+1] - nElemSend[iProc];

      /*--- Total count can include multiple pieces of data per element. ---*/

      count = countPerElem*nElem;

      /*--- Post non-blocking send for this proc. ---*/

      dest = iProc; tag = rank + 1;

      switch (commType) {
        case COMM_TYPE_DOUBLE:
          SU2_MPI::Isend(&(static_cast<su2double*>(bufSend)[offset]),
                         count, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_LONG:
          SU2_MPI::Isend(&(static_cast<unsigned long*>(bufSend)[offset]),
                         count, MPI_UNSIGNED_LONG, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_LONG:
          SU2_MPI::Isend(&(static_cast<long*>(bufSend)[offset]),
                         count, MPI_LONG, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_UNSIGNED_SHORT:
          SU2_MPI::Isend(&(static_cast<unsigned short*>(bufSend)[offset]),
                         count, MPI_UNSIGNED_SHORT, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_CHAR:
          SU2_MPI::Isend(&(static_cast<char*>(bufSend)[offset]),
                         count, MPI_CHAR, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_SHORT:
          SU2_MPI::Isend(&(static_cast<short*>(bufSend)[offset]),
                         count, MPI_SHORT, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        case COMM_TYPE_INT:
          SU2_MPI::Isend(&(static_cast<int*>(bufSend)[offset]),
                         count, MPI_INT, dest, tag, MPI_COMM_WORLD,
                         &(sendReq[iMessage]));
          break;
        default:
          break;
      }

      /*--- Increment message counter. ---*/

      iMessage++;

    }
  }

}

void CPhysicalGeometry::CompleteCommsAll(int nSends,
                                         SU2_MPI::Request *sendReq,
                                         int nRecvs,
                                         SU2_MPI::Request *recvReq) {

  /*--- Local variables ---*/

  int ind, iSend, iRecv;
  SU2_MPI::Status status;

  /*--- Wait for the non-blocking sends to complete. ---*/

  for (iSend = 0; iSend < nSends; iSend++)
    SU2_MPI::Waitany(nSends, sendReq, &ind, &status);

  /*--- Wait for the non-blocking recvs to complete. ---*/

  for (iRecv = 0; iRecv < nRecvs; iRecv++)
    SU2_MPI::Waitany(nRecvs, recvReq, &ind, &status);

}

void CPhysicalGeometry::PrepareOffsets(unsigned long val_npoint_global) {

  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/

  if (beg_node == NULL) beg_node = new unsigned long[size];
  if (end_node == NULL) end_node = new unsigned long[size];

  if (nPointLinear     == NULL) nPointLinear     = new unsigned long[size];
  if (nPointCumulative == NULL) nPointCumulative = new unsigned long[size+1];

  unsigned long quotient = val_npoint_global/size;
  int remainder = int(val_npoint_global%size);
  for (int ii = 0; ii < size; ii++) {
    nPointLinear[ii] = quotient + int(ii < remainder);
  }

  /*--- Store the local number of nodes on each proc in the linear
   partitioning, the beginning/end index, and the linear partitioning
   within an array in cumulative storage format. ---*/

  beg_node[0] = 0;
  end_node[0] = beg_node[0] + nPointLinear[0];
  nPointCumulative[0] = 0;
  for (int iProc = 1; iProc < size; iProc++) {
    beg_node[iProc] = end_node[iProc-1];
    end_node[iProc] = beg_node[iProc] + nPointLinear[iProc];
    nPointCumulative[iProc] = nPointCumulative[iProc-1] + nPointLinear[iProc-1];
  }
  nPointCumulative[size] = val_npoint_global;

}

unsigned long CPhysicalGeometry::GetLinearPartition(unsigned long val_global_index) {

  unsigned long iProcessor = 0;

  /*--- Initial guess ---*/

  iProcessor = val_global_index/nPointLinear[0];

  /*--- Guard against going over size. ---*/

  if (iProcessor >= (unsigned long)size)
    iProcessor = (unsigned long)size-1;

  /*--- Move up or down until we find the processor. ---*/

  if (val_global_index >= nPointCumulative[iProcessor])
    while(val_global_index >= nPointCumulative[iProcessor+1])
      iProcessor++;
  else
    while(val_global_index < nPointCumulative[iProcessor])
      iProcessor--;

  return iProcessor;

}

void CPhysicalGeometry::SortAdjacency(CConfig *config) {

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Executing the partitioning functions." << endl;

  /*--- Post process the adjacency information in order to get it into the
   CSR format before sending the data to ParMETIS. We need to remove
   repeats and adjust the size of the array for each local node. ---*/

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Building the graph adjacency structure." << endl;

  /*--- Create a partitioner object so we can transform the global
   index values stored in the elements to a local index. ---*/

  CLinearPartitioner pointPartitioner(Global_nPointDomain,0);

  /*--- We can already create the array that indexes the adjacency. ---*/

  if (xadj == NULL) xadj = new idx_t[pointPartitioner.GetSizeOnRank(rank)+1];
  xadj[0] = 0;

  /*--- Here, we transfer the adjacency information from a multi-dim vector
   on a node-by-node basis into a single vector container. First, we sort
   the entries and remove the duplicates we find for each node, then we
   copy it into the single vect and clear memory from the multi-dim vec. ---*/

  unsigned long total_adj_size = 0;
  vector<unsigned long>::iterator it;
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- For each point, sort the adjacency in ascending order
     so that we can remove duplicates and complete the size for
     unique set of adjacent nodes for that point. ---*/

    sort(adj_nodes[iPoint].begin(), adj_nodes[iPoint].end());
    it = unique(adj_nodes[iPoint].begin(), adj_nodes[iPoint].end());
    const unsigned long local_size = it - adj_nodes[iPoint].begin();
    adj_nodes[iPoint].resize(local_size);
    total_adj_size += local_size;

  }

  /*--- Now that we know the size, create the final adjacency array. This
   is the array that we will feed to ParMETIS for partitioning. ---*/

  if (adjacency == NULL) adjacency = new idx_t[total_adj_size];

  unsigned long kPoint = 0;
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Local size of the adjacency for the current point. ---*/

    const unsigned long local_size = adj_nodes[iPoint].size();

    /*--- Move the sorted adjacency into a 1-D vector for all
     points for loading into ParMETIS for partitioning next. ---*/

    for (unsigned long jPoint = 0; jPoint < local_size; jPoint++) {
      adjacency[kPoint] = (idx_t)adj_nodes[iPoint][jPoint]; kPoint++;
    }

    /*--- Increment the starting index for the next point (CSR). ---*/

    xadj[iPoint+1] = xadj[iPoint] + local_size;

    /*--- Free vector memory as we go. ---*/

    vector<unsigned long>().swap(adj_nodes[iPoint]);

  }

  /*--- Force free the entire old multi-dim. adjacency vector. ---*/

  vector< vector<unsigned long> >().swap(adj_nodes);

#endif
#endif

}

void CPhysicalGeometry::SetSendReceive(CConfig *config) {

  unsigned short Counter_Send, Counter_Receive, iMarkerSend, iMarkerReceive;
  unsigned long iVertex, LocalNode;
  unsigned short nMarker_Max = config->GetnMarker_Max();
  unsigned long  iPoint, jPoint, iElem, nDomain, iDomain, jDomain;
  unsigned long *nVertexDomain = new unsigned long[nMarker_Max];
  unsigned short iNode, jNode;
  vector<unsigned long>::iterator it;

  vector<vector<unsigned long> > SendTransfLocal;       /*!< \brief Vector to store the type of transformation for this send point. */
  vector<vector<unsigned long> > ReceivedTransfLocal;   /*!< \brief Vector to store the type of transformation for this received point. */
  vector<vector<unsigned long> > SendDomainLocal;       /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */
  vector<vector<unsigned long> > ReceivedDomainLocal;   /*!< \brief SendDomain[from domain][to domain] and return the point index of the node that must me sended. */

  map<unsigned long, unsigned long>::const_iterator MI;

  if (rank == MASTER_NODE && size > SINGLE_NODE)
    cout << "Establishing MPI communication patterns." << endl;

  nDomain = size;

  SendTransfLocal.resize(nDomain);
  ReceivedTransfLocal.resize(nDomain);
  SendDomainLocal.resize(nDomain);
  ReceivedDomainLocal.resize(nDomain);

  /*--- Loop over the all the points of the elements on this rank in order
   to find the points with different colors. Create the send/received lists
   from this information. ---*/

  for (iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {

      iPoint  = elem[iElem]->GetNode(iNode);
      iDomain = node[iPoint]->GetColor();

      if (iDomain == (unsigned long) rank) {
        for (jNode = 0; jNode < elem[iElem]->GetnNodes(); jNode++) {

          jPoint  = elem[iElem]->GetNode(jNode);
          jDomain = node[jPoint]->GetColor();

          /*--- If one of the neighbors is a different color and connected
           by an edge, then we add them to the list. ---*/

          if (iDomain != jDomain) {

            /*--- We send from iDomain to jDomain the value of iPoint,
             we save the global value becuase we need to sort the lists. ---*/

            SendDomainLocal[jDomain].push_back(Local_to_Global_Point[iPoint]);

            /*--- We send from jDomain to iDomain the value of jPoint,
             we save the global value becuase we need to sort the lists. ---*/

            ReceivedDomainLocal[jDomain].push_back(Local_to_Global_Point[jPoint]);

          }
        }
      }
    }
  }

  /*--- Sort the points that must be sent and delete repeated points, note
   that the sorting should be done with the global index (not the local). ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort(SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    it = unique(SendDomainLocal[iDomain].begin(), SendDomainLocal[iDomain].end());
    SendDomainLocal[iDomain].resize(it - SendDomainLocal[iDomain].begin());
  }

  /*--- Sort the points that must be received and delete repeated points, note
   that the sorting should be done with the global point (not the local). ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    sort(ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    it = unique( ReceivedDomainLocal[iDomain].begin(), ReceivedDomainLocal[iDomain].end());
    ReceivedDomainLocal[iDomain].resize(it - ReceivedDomainLocal[iDomain].begin());
  }

  /*--- Create Global to Local Point array, note that the array is smaller (Max_GlobalPoint) than the total
   number of points in the simulation  ---*/
  Max_GlobalPoint = 0;
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    if (Local_to_Global_Point[iPoint] > (long)Max_GlobalPoint)
      Max_GlobalPoint = Local_to_Global_Point[iPoint];
  }

  /*--- Set the value of some of the points ---*/
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    Global_to_Local_Point[Local_to_Global_Point[iPoint]] = iPoint;

  /*--- Add the new MPI send boundaries, reset the transformation,
   and save the local value. ---*/

  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = SendDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {

        MI = Global_to_Local_Point.find(SendDomainLocal[iDomain][iVertex]);
        if (MI != Global_to_Local_Point.end())
          iPoint = Global_to_Local_Point[SendDomainLocal[iDomain][iVertex]];
        else iPoint = -1;

        SendDomainLocal[iDomain][iVertex] = iPoint;
        SendTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }

  /*--- Add the new MPI receive boundaries, reset the transformation, and save the local value ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      nVertexDomain[nMarker] = ReceivedDomainLocal[iDomain].size();
      for (iVertex = 0; iVertex < nVertexDomain[nMarker]; iVertex++) {

        MI = Global_to_Local_Point.find(ReceivedDomainLocal[iDomain][iVertex]);
        if (MI != Global_to_Local_Point.end())
          iPoint = Global_to_Local_Point[ReceivedDomainLocal[iDomain][iVertex]];
        else iPoint = -1;

        ReceivedDomainLocal[iDomain][iVertex] = iPoint;
        ReceivedTransfLocal[iDomain].push_back(0);
      }
      nElem_Bound[nMarker] = nVertexDomain[nMarker];
      bound[nMarker] = new CPrimalGrid*[nElem_Bound[nMarker]];
      nMarker++;
    }
  }

  /*--- First compute the Send/Receive boundaries ---*/
  Counter_Send = 0; 	Counter_Receive = 0;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (SendDomainLocal[iDomain].size() != 0) Counter_Send++;

  for (iDomain = 0; iDomain < nDomain; iDomain++)
    if (ReceivedDomainLocal[iDomain].size() != 0) Counter_Receive++;

  iMarkerSend    = nMarker - Counter_Send - Counter_Receive;
  iMarkerReceive = nMarker - Counter_Receive;

  /*--- First we do the send ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (SendDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerSend); iVertex++) {
        LocalNode = SendDomainLocal[iDomain][iVertex];
        bound[iMarkerSend][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerSend][iVertex]->SetRotation_Type(SendTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerSend] = iDomain+1;
      iMarkerSend++;
    }
  }

  /*--- Second we do the receive ---*/
  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    if (ReceivedDomainLocal[iDomain].size() != 0) {
      for (iVertex = 0; iVertex < GetnElem_Bound(iMarkerReceive); iVertex++) {
        LocalNode = ReceivedDomainLocal[iDomain][iVertex];
        bound[iMarkerReceive][iVertex] = new CVertexMPI(LocalNode, nDim);
        bound[iMarkerReceive][iVertex]->SetRotation_Type(ReceivedTransfLocal[iDomain][iVertex]);
      }
      Marker_All_SendRecv[iMarkerReceive] = -(iDomain+1);
      iMarkerReceive++;
    }
  }

  /*--- Free memory ---*/

  delete [] nVertexDomain;

}

void CPhysicalGeometry::SetBoundaries(CConfig *config) {

  unsigned long iElem_Bound, TotalElem, *nElem_Bound_Copy, iVertex_;
  string Grid_Marker;
  unsigned short iDomain, nDomain, iMarkersDomain, iLoop, *DomainCount, nMarker_Physical, Duplicate_SendReceive, *DomainSendCount,
                 **DomainSendMarkers, *DomainReceiveCount, **DomainReceiveMarkers, nMarker_SendRecv, iMarker, iMarker_;
  CPrimalGrid*** bound_Copy;
  short *Marker_All_SendRecv_Copy;
  bool CheckStart;

  nDomain = size+1;

  /*--- Count the number of physical markers
   in the boundaries ---*/

  nMarker_Physical = 0;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {
      nMarker_Physical++;
    }
  }

  /*--- Identify if there are markers that send/received with the same domain,
   they should be together---*/

  Duplicate_SendReceive = 0;
  for (iLoop = 0; iLoop < 2; iLoop++) {

    DomainCount = new unsigned short [nDomain];

    for (iDomain = 0; iDomain < nDomain; iDomain++)
      DomainCount[iDomain] = 0;

    if (iLoop == 0) {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == iDomain) DomainCount[iDomain]++;
    }
    else {
      for (iDomain = 0; iDomain < nDomain; iDomain++)
        for (iMarker = 0; iMarker < nMarker; iMarker++)
          if (bound[iMarker][0]->GetVTK_Type() == VERTEX)
            if (Marker_All_SendRecv[iMarker] == -iDomain) DomainCount[iDomain]++;
    }

    for (iDomain = 0; iDomain < nDomain; iDomain++)
      if (DomainCount[iDomain] > 1) Duplicate_SendReceive++;

    delete [] DomainCount;

  }

  DomainSendCount = new unsigned short [nDomain];
  DomainSendMarkers = new unsigned short *[nDomain];
  DomainReceiveCount = new unsigned short [nDomain];
  DomainReceiveMarkers = new unsigned short *[nDomain];

  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    DomainSendCount[iDomain] = 0;
    DomainSendMarkers[iDomain] = new unsigned short [nMarker];

    DomainReceiveCount[iDomain] = 0;
    DomainReceiveMarkers[iDomain] = new unsigned short [nMarker];
  }

  for (iDomain = 0; iDomain < nDomain; iDomain++) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
        if (Marker_All_SendRecv[iMarker] == iDomain) {
          DomainSendMarkers[iDomain][DomainSendCount[iDomain]] = iMarker;
          DomainSendCount[iDomain]++;
        }
        if (Marker_All_SendRecv[iMarker] == -iDomain) {
          DomainReceiveMarkers[iDomain][DomainReceiveCount[iDomain]] = iMarker;
          DomainReceiveCount[iDomain]++;
        }
      }
    }
  }

  /*--- Create an structure to store the Send/Receive
   boundaries, because they require some reorganization ---*/

  nMarker_SendRecv = nMarker - nMarker_Physical - Duplicate_SendReceive;
  bound_Copy = new CPrimalGrid**[nMarker_Physical + nMarker_SendRecv];
  nElem_Bound_Copy = new unsigned long [nMarker_Physical + nMarker_SendRecv];
  Marker_All_SendRecv_Copy = new short [nMarker_Physical + nMarker_SendRecv];
  iMarker_ = nMarker_Physical;
  iVertex_ = 0;
  CheckStart = false;

  /*--- Copy and allocate the physical markers in the data structure ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {

      nElem_Bound_Copy[iMarker] = nElem_Bound[iMarker];
      bound_Copy[iMarker] = new CPrimalGrid* [nElem_Bound[iMarker]];

      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == LINE)
          bound_Copy[iMarker][iElem_Bound] = new CLine(bound[iMarker][iElem_Bound]->GetNode(0),
                                                       bound[iMarker][iElem_Bound]->GetNode(1), 2);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == TRIANGLE)

          bound_Copy[iMarker][iElem_Bound] = new CTriangle(bound[iMarker][iElem_Bound]->GetNode(0),
                                                           bound[iMarker][iElem_Bound]->GetNode(1),
                                                           bound[iMarker][iElem_Bound]->GetNode(2), 3);
        if (bound[iMarker][iElem_Bound]->GetVTK_Type() == QUADRILATERAL)
          bound_Copy[iMarker][iElem_Bound] = new CQuadrilateral(bound[iMarker][iElem_Bound]->GetNode(0),
                                                            bound[iMarker][iElem_Bound]->GetNode(1),
                                                            bound[iMarker][iElem_Bound]->GetNode(2),
                                                            bound[iMarker][iElem_Bound]->GetNode(3), 3);
      }
    }
  }


  for (iDomain = 0; iDomain < nDomain; iDomain++) {

    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/

    if (DomainSendCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];
    }

    for (iMarkersDomain = 0; iMarkersDomain < DomainSendCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainSendMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];

      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }

    }

    /*--- Compute the total number of elements (adding all the
     boundaries with the same Send/Receive ---*/

    if (DomainReceiveCount[iDomain] != 0) {
      TotalElem = 0;
      for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
        iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
        TotalElem += nElem_Bound[iMarker];
      }
      if (CheckStart) iMarker_++;
      CheckStart = true;
      iVertex_ = 0;
      nElem_Bound_Copy[iMarker_] = TotalElem;
      bound_Copy[iMarker_] = new CPrimalGrid*[TotalElem];

    }

    for (iMarkersDomain = 0; iMarkersDomain < DomainReceiveCount[iDomain]; iMarkersDomain++) {
      iMarker = DomainReceiveMarkers[iDomain][iMarkersDomain];
      Marker_All_SendRecv_Copy[iMarker_] = Marker_All_SendRecv[iMarker];

      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        bound_Copy[iMarker_][iVertex_] = new CVertexMPI(bound[iMarker][iElem_Bound]->GetNode(0), nDim);
        bound_Copy[iMarker_][iVertex_]->SetRotation_Type(bound[iMarker][iElem_Bound]->GetRotation_Type());
        iVertex_++;
      }

    }

  }

  delete [] DomainSendCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete [] DomainSendMarkers[iDomain];
  delete[] DomainSendMarkers;

  delete [] DomainReceiveCount;
  for (iDomain = 0; iDomain < nDomain; iDomain++)
    delete [] DomainReceiveMarkers[iDomain];
  delete[] DomainReceiveMarkers;

   /*--- Deallocate the bound variables ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
   for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++)
     if (bound[iMarker][iElem_Bound] != NULL) delete bound[iMarker][iElem_Bound];
    if (bound[iMarker] != NULL) delete [] bound[iMarker];
  }
  if (bound != NULL) delete [] bound;

  /*--- Allocate the new bound variables, and set the number of markers ---*/

  bound = bound_Copy;
  nMarker = nMarker_Physical + nMarker_SendRecv;

  config->SetnMarker_All(nMarker);

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    nElem_Bound[iMarker] = nElem_Bound_Copy[iMarker];
  }

  for (iMarker = nMarker_Physical; iMarker < nMarker; iMarker++) {
    Marker_All_SendRecv[iMarker] = Marker_All_SendRecv_Copy[iMarker];
    config->SetMarker_All_SendRecv(iMarker, Marker_All_SendRecv[iMarker]);
    config->SetMarker_All_TagBound(iMarker, "SEND_RECEIVE");
  }

  /*--- Update config information storing the boundary information in the right place ---*/

  for (iMarker = 0 ; iMarker < nMarker; iMarker++) {

    string Marker_Tag = config->GetMarker_All_TagBound(iMarker);

    if (Marker_Tag != "SEND_RECEIVE") {

      /*--- Update config information storing the boundary information in the right place ---*/

      Tag_to_Marker[config->GetMarker_CfgFile_TagBound(Marker_Tag)] = Marker_Tag;
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
      config->SetMarker_All_Turbomachinery(iMarker, config->GetMarker_CfgFile_Turbomachinery(Marker_Tag));
      config->SetMarker_All_TurbomachineryFlag(iMarker, config->GetMarker_CfgFile_TurbomachineryFlag(Marker_Tag));
      config->SetMarker_All_MixingPlaneInterface(iMarker, config->GetMarker_CfgFile_MixingPlaneInterface(Marker_Tag));
    }

    /*--- Send-Receive boundaries definition ---*/

    else {

      config->SetMarker_All_KindBC(iMarker, SEND_RECEIVE);
      config->SetMarker_All_Monitoring(iMarker, NO);
      config->SetMarker_All_GeoEval(iMarker, NO);
      config->SetMarker_All_Designing(iMarker, NO);
      config->SetMarker_All_Plotting(iMarker, NO);
      config->SetMarker_All_Analyze(iMarker, NO);
      config->SetMarker_All_ZoneInterface(iMarker, NO);
      config->SetMarker_All_DV(iMarker, NO);
      config->SetMarker_All_Moving(iMarker, NO);
      config->SetMarker_All_Deform_Mesh(iMarker, NO);
      config->SetMarker_All_Fluid_Load(iMarker, NO);
      config->SetMarker_All_PyCustom(iMarker, NO);
      config->SetMarker_All_PerBound(iMarker, NO);
      config->SetMarker_All_Turbomachinery(iMarker, NO);
      config->SetMarker_All_TurbomachineryFlag(iMarker, NO);
      config->SetMarker_All_MixingPlaneInterface(iMarker, NO);

      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        if (config->GetMarker_All_SendRecv(iMarker) < 0)
          node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
      }

    }

    /*--- Loop over the surface element to set the boundaries ---*/

    unsigned long Point_Surface, iElem_Surface;
    unsigned short iNode_Surface;

    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {
      for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
        Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
        node[Point_Surface]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPhysicalBoundary(true);

        if (config->GetSolid_Wall(iMarker))
          node[Point_Surface]->SetSolidBoundary(true);

        if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
          node[Point_Surface]->SetPeriodicBoundary(true);
      }
    }

  }

  delete [] Marker_All_SendRecv_Copy;
  delete [] nElem_Bound_Copy;
}

void CPhysicalGeometry::Read_Mesh_FVM(CConfig        *config,
                                      string         val_mesh_filename,
                                      unsigned short val_iZone,
                                      unsigned short val_nZone) {

  /*--- Initialize counters for local/global points & elements ---*/

  Global_nPoint  = 0; Global_nPointDomain   = 0;
  Global_nElem   = 0; Global_nElemDomain    = 0;
  nelem_edge     = 0; Global_nelem_edge     = 0;
  nelem_triangle = 0; Global_nelem_triangle = 0;
  nelem_quad     = 0; Global_nelem_quad     = 0;
  nelem_tetra    = 0; Global_nelem_tetra    = 0;
  nelem_hexa     = 0; Global_nelem_hexa     = 0;
  nelem_prism    = 0; Global_nelem_prism    = 0;
  nelem_pyramid  = 0; Global_nelem_pyramid  = 0;

  /*--- Set the zone number from the input value. ---*/

  nZone = val_nZone;

  /*--- Create a mesh reader to read a CGNS grid into linear partitions. ---*/

  unsigned short val_format = config->GetMesh_FileFormat();

  CMeshReaderFVM *MeshFVM = NULL;
  switch (val_format) {
    case SU2:
      MeshFVM = new CSU2ASCIIMeshReaderFVM(config, val_iZone, val_nZone);
      break;
    case CGNS_GRID:
      MeshFVM = new CCGNSMeshReaderFVM(config, val_iZone, val_nZone);
      break;
    case RECTANGLE:
      MeshFVM = new CRectangularMeshReaderFVM(config, val_iZone, val_nZone);
      break;
    case BOX:
      MeshFVM = new CBoxMeshReaderFVM(config, val_iZone, val_nZone);
      break;
    default:
      SU2_MPI::Error("Unrecognized mesh format specified!", CURRENT_FUNCTION);
      break;
  }

  /*--- Store the dimension of the problem ---*/

  nDim = MeshFVM->GetDimension();
  if (rank == MASTER_NODE) {
    if (nDim == 2) cout << "Two dimensional problem." << endl;
    if (nDim == 3) cout << "Three dimensional problem." << endl;
  }

  /*--- Store the local and global number of nodes for this rank. ---*/

  nPoint              = MeshFVM->GetNumberOfLocalPoints();
  nPointDomain        = MeshFVM->GetNumberOfLocalPoints();
  Global_nPoint       = MeshFVM->GetNumberOfGlobalPoints();
  Global_nPointDomain = MeshFVM->GetNumberOfGlobalPoints();

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
    cout << Global_nPoint << " grid points before partitioning." << endl;
  } else if (rank == MASTER_NODE) {
    cout << Global_nPoint << " grid points." << endl;
  }

  /*--- Store the local and global number of interior elements. ---*/

  nElem              = MeshFVM->GetNumberOfLocalElements();
  Global_nElem       = MeshFVM->GetNumberOfGlobalElements();
  Global_nElemDomain = MeshFVM->GetNumberOfGlobalElements();

  if ((rank == MASTER_NODE) && (size > SINGLE_NODE)) {
    cout << Global_nElem << " volume elements before partitioning." << endl;
  } else if (rank == MASTER_NODE) {
    cout << Global_nElem << " volume elements." << endl;
  }

  /*--- Load the grid points, volume elements, and surface elements
   from the mesh object into the proper SU2 data structures. ---*/

  LoadLinearlyPartitionedPoints(config,         MeshFVM);
  LoadLinearlyPartitionedVolumeElements(config, MeshFVM);
  LoadUnpartitionedSurfaceElements(config,      MeshFVM);

  /*--- Prepare the nodal adjacency structures for ParMETIS. ---*/

  PrepareAdjacency(config);

  /*--- Now that we have loaded all information from the mesh,
   delete the mesh reader object. ---*/

  if (MeshFVM != NULL) delete MeshFVM;

}

void CPhysicalGeometry::LoadLinearlyPartitionedPoints(CConfig        *config,
                                                      CMeshReaderFVM *mesh) {

  /*--- Get the linearly partitioned coordinates from the mesh object. ---*/

  const vector<vector<passivedouble> > &gridCoords =
  mesh->GetLocalPointCoordinates();

  /*--- Initialize point counts and the grid node data structure. ---*/

  nPointNode = nPoint;
  node       = new CPoint*[nPoint];

  /*--- Loop over the CGNS grid nodes and load into the SU2 data
   structure. Note that since we have performed a linear partitioning
   of the grid nodes, we can simply initialize the global index to
   the first node that lies on our rank and increment. ---*/

  CLinearPartitioner pointPartitioner(Global_nPointDomain,0);
  unsigned long GlobalIndex = pointPartitioner.GetFirstIndexOnRank(rank);
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    switch(nDim) {
      case 2:
        node[iPoint] = new CPoint(gridCoords[0][iPoint],
                                  gridCoords[1][iPoint],
                                  GlobalIndex,
                                  config);
        GlobalIndex++;
        break;
      case 3:
        node[iPoint] = new CPoint(gridCoords[0][iPoint],
                                  gridCoords[1][iPoint],
                                  gridCoords[2][iPoint],
                                  GlobalIndex,
                                  config);
        GlobalIndex++;
        break;
    }
  }

}

void CPhysicalGeometry::LoadLinearlyPartitionedVolumeElements(CConfig        *config,
                                                              CMeshReaderFVM *mesh) {

  /*--- Reset the global to local element mapping. ---*/

  Global_to_Local_Elem.clear();

  /*--- Get the volume connectivity from the mesh object. ---*/

  const vector<unsigned long> &connElems =
  mesh->GetLocalVolumeElementConnectivity();

  /*--- Allocate space for the CGNS interior elements in our SU2 data
   structure. Note that we only instantiate our rank's local set. ---*/

  elem = new CPrimalGrid*[nElem];

  /*--- Some temporaries for the loop below. ---*/

  unsigned long Global_Index_Elem = 0, iElem = 0;
  vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON);

  /*--- Loop over all of the internal, local volumetric elements. ---*/

  for (unsigned long jElem = 0; jElem < nElem; jElem++) {

    /*--- Get the global ID for this element. This is stored in
     the first entry of our connectivity stucture. ---*/

    Global_Index_Elem = connElems[jElem*SU2_CONN_SIZE + 0];

    /*--- Get the VTK type for this element. This is stored in the
     second entry of the connectivity structure. ---*/

    int vtk_type = (int)connElems[jElem*SU2_CONN_SIZE + 1];

    /*--- Instantiate this element in the proper SU2 data structure.
     During this loop, we also set the global to local element map
     for later use and increment the element counts for all types. ---*/

    switch(vtk_type) {

      case TRIANGLE:

        for (unsigned long j = 0; j < N_POINTS_TRIANGLE; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CTriangle(connectivity[0],
                                    connectivity[1],
                                    connectivity[2], nDim);
        iElem++; nelem_triangle++;
        break;

      case QUADRILATERAL:

        for (unsigned long j = 0; j < N_POINTS_QUADRILATERAL; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CQuadrilateral(connectivity[0],
                                         connectivity[1],
                                         connectivity[2],
                                         connectivity[3], nDim);
        iElem++; nelem_quad++;
        break;

      case TETRAHEDRON:

        for (unsigned long j = 0; j < N_POINTS_TETRAHEDRON; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CTetrahedron(connectivity[0],
                                       connectivity[1],
                                       connectivity[2],
                                       connectivity[3]);
        iElem++; nelem_tetra++;
        break;

      case HEXAHEDRON:

        for (unsigned long j = 0; j < N_POINTS_HEXAHEDRON; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CHexahedron(connectivity[0],
                                      connectivity[1],
                                      connectivity[2],
                                      connectivity[3],
                                      connectivity[4],
                                      connectivity[5],
                                      connectivity[6],
                                      connectivity[7]);
        iElem++; nelem_hexa++;
        break;

      case PRISM:

        for (unsigned long j = 0; j < N_POINTS_PRISM; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CPrism(connectivity[0],
                                 connectivity[1],
                                 connectivity[2],
                                 connectivity[3],
                                 connectivity[4],
                                 connectivity[5]);
        iElem++; nelem_prism++;
        break;

      case PYRAMID:

        for (unsigned long j = 0; j < N_POINTS_PYRAMID; j++) {
          connectivity[j] = connElems[jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j];
        }
        Global_to_Local_Elem[Global_Index_Elem] = iElem;
        elem[iElem] = new CPyramid(connectivity[0],
                                   connectivity[1],
                                   connectivity[2],
                                   connectivity[3],
                                   connectivity[4]);
        iElem++; nelem_pyramid++;
        break;

      default:
        SU2_MPI::Error("Element type not supported!", CURRENT_FUNCTION);
        break;
    }
  }


  /*--- Reduce the global counts of all element types found in
   the CGNS grid with all ranks. ---*/

#ifdef HAVE_MPI
  unsigned long Local_nElemTri     = nelem_triangle;
  unsigned long Local_nElemQuad    = nelem_quad;
  unsigned long Local_nElemTet     = nelem_tetra;
  unsigned long Local_nElemHex     = nelem_hexa;
  unsigned long Local_nElemPrism   = nelem_prism;
  unsigned long Local_nElemPyramid = nelem_pyramid;
  SU2_MPI::Allreduce(&Local_nElemTri,     &Global_nelem_triangle,  1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemQuad,    &Global_nelem_quad,      1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemTet,     &Global_nelem_tetra,     1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemHex,     &Global_nelem_hexa,      1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPrism,   &Global_nelem_prism,     1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Local_nElemPyramid, &Global_nelem_pyramid,   1,
                     MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  Global_nelem_triangle = nelem_triangle;
  Global_nelem_quad     = nelem_quad;
  Global_nelem_tetra    = nelem_tetra;
  Global_nelem_hexa     = nelem_hexa;
  Global_nelem_prism    = nelem_prism;
  Global_nelem_pyramid  = nelem_pyramid;
#endif

}

void CPhysicalGeometry::LoadUnpartitionedSurfaceElements(CConfig        *config,
                                                         CMeshReaderFVM *mesh) {

  /*--- The master node takes care of loading all markers and
   surface elements from the file. This information is later
   put into linear partitions to make its redistribution easier
   after we call ParMETIS. ---*/

  if (rank == MASTER_NODE) {

    const vector<string> &sectionNames = mesh->GetMarkerNames();

    /*--- Store the number of markers and print to the screen. ---*/

    nMarker = mesh->GetNumberOfMarkers();
    config->SetnMarker_All(nMarker);
    cout << nMarker << " surface markers." << endl;

    /*--- Create the data structure for boundary elements. ---*/

    bound         = new CPrimalGrid**[nMarker];
    nElem_Bound   = new unsigned long [nMarker];
    Tag_to_Marker = new string [config->GetnMarker_Max()];

    /*--- Set some temporaries for the loop below. ---*/

    int npe, vtk_type;
    unsigned long iElem = 0;
    vector<unsigned long> connectivity(N_POINTS_HEXAHEDRON);

    /*--- Loop over all sections that we extracted from the CGNS file
     that were identified as boundary element sections so that we can
     store those elements into our SU2 data structures. ---*/

    for (int iMarker = 0; iMarker < nMarker; iMarker++) {

      /*--- Initialize some counter variables ---*/

      nelem_edge_bound = 0; nelem_triangle_bound = 0;
      nelem_quad_bound = 0; iElem = 0;

      /*--- Get the string name for this marker. ---*/

      string Marker_Tag = sectionNames[iMarker];

      /* Get the marker info and surface connectivity from the mesh object. */

      const unsigned long surfElems =
      mesh->GetNumberOfSurfaceElementsForMarker(iMarker);

      const vector<unsigned long> &connElems =
      mesh->GetSurfaceElementConnectivityForMarker(iMarker);

      /*--- Set the number of boundary elements in this marker. ---*/

      nElem_Bound[iMarker] = surfElems;

      /*--- Report the number and name of the marker to the console. ---*/

      cout << nElem_Bound[iMarker]  << " boundary elements in index ";
      cout << iMarker <<" (Marker = " <<Marker_Tag<< ")." << endl;

      /*--- Instantiate the list of elements in the data structure. ---*/

      bound[iMarker] = new CPrimalGrid*[nElem_Bound[iMarker]];

      for (unsigned long jElem = 0; jElem < nElem_Bound[iMarker]; jElem++ ) {

        /*--- Not a mixed section. We already know the element type,
         which is stored ---*/

        vtk_type = (int)connElems[jElem*SU2_CONN_SIZE + 1];

        /*--- Store the loop size more easily. ---*/

        npe = (int)(SU2_CONN_SIZE-SU2_CONN_SKIP);

        /*--- Store the nodes for this element more clearly. ---*/

        for (int j = 0; j < npe; j++) {
          unsigned long nn = jElem*SU2_CONN_SIZE + SU2_CONN_SKIP + j;
          connectivity[j] = connElems[nn];
        }

        /*--- Instantiate the boundary element object. ---*/

        switch(vtk_type) {
          case LINE:
            bound[iMarker][iElem] = new CLine(connectivity[0],
                                              connectivity[1],2);
            iElem++; nelem_edge_bound++; break;
          case TRIANGLE:
            bound[iMarker][iElem] = new CTriangle(connectivity[0],
                                                  connectivity[1],
                                                  connectivity[2],3);
            iElem++; nelem_triangle_bound++; break;
          case QUADRILATERAL:
            bound[iMarker][iElem] = new CQuadrilateral(connectivity[0],
                                                       connectivity[1],
                                                       connectivity[2],
                                                       connectivity[3],3);
            iElem++; nelem_quad_bound++; break;
        }

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

}

void CPhysicalGeometry::PrepareAdjacency(CConfig *config) {

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS

  /*--- Resize the vector for the adjacency information (ParMETIS). ---*/

  adj_nodes.clear();
  adj_nodes.resize(nPoint);
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
  adj_nodes[iPoint].resize(0);

  /*--- Create a partitioner object so we can transform the global
   index values stored in the elements to a local index. ---*/

  CLinearPartitioner pointPartitioner(Global_nPointDomain,0);
  const unsigned long firstIndex = pointPartitioner.GetFirstIndexOnRank(rank);

  /*--- Loop over all elements that are now loaded and store adjacency. ---*/

  unsigned long connectivity[8] = {0,0,0,0,0,0,0,0};
  for (unsigned long iElem = 0; iElem < nElem; iElem++) {

    /*--- Get the VTK type for this element. This is stored in the
     first entry of the connectivity structure. ---*/

    unsigned short VTK_Type = elem[iElem]->GetVTK_Type();

    /*--- Instantiate this element and build adjacency structure. ---*/

    switch(VTK_Type) {

      case TRIANGLE:

        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_TRIANGLE; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_TRIANGLE; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            for (unsigned long jNode = 0; jNode < N_POINTS_TRIANGLE; jNode++) {

              /*--- Build adjacency assuming the VTK connectivity ---*/

              if (iNode != jNode)
              adj_nodes[local_index].push_back(connectivity[jNode]);

            }

          }
        }

        break;

      case QUADRILATERAL:

        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_QUADRILATERAL; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_QUADRILATERAL; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            /*--- Build adjacency assuming the VTK connectivity ---*/

            adj_nodes[local_index].push_back(connectivity[(iNode+1)%4]);
            adj_nodes[local_index].push_back(connectivity[(iNode+3)%4]);

          }
        }

        break;

      case TETRAHEDRON:

        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_TETRAHEDRON; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_TETRAHEDRON; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            for (unsigned long jNode = 0; jNode < N_POINTS_TETRAHEDRON; jNode++) {

              /*--- Build adjacency assuming the VTK connectivity ---*/

              if (iNode != jNode)
              adj_nodes[local_index].push_back(connectivity[jNode]);

            }

          }
        }

        break;

      case HEXAHEDRON:

        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_HEXAHEDRON; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_HEXAHEDRON; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            /*--- Build adjacency assuming the VTK connectivity ---*/

            if (iNode < 4) {
              adj_nodes[local_index].push_back(connectivity[(iNode+1)%4]);
              adj_nodes[local_index].push_back(connectivity[(iNode+3)%4]);
            } else {
              adj_nodes[local_index].push_back(connectivity[(iNode-3)%4+4]);
              adj_nodes[local_index].push_back(connectivity[(iNode-1)%4+4]);
            }
            adj_nodes[local_index].push_back(connectivity[(iNode+4)%8]);

          }
        }

        break;

      case PRISM:

        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_PRISM; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_PRISM; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            /*--- Build adjacency assuming the VTK connectivity ---*/

            if (iNode < 3) {
              adj_nodes[local_index].push_back(connectivity[(iNode+1)%3]);
              adj_nodes[local_index].push_back(connectivity[(iNode+2)%3]);
            } else {
              adj_nodes[local_index].push_back(connectivity[(iNode-2)%3+3]);
              adj_nodes[local_index].push_back(connectivity[(iNode-1)%3+3]);
            }
            adj_nodes[local_index].push_back(connectivity[(iNode+3)%6]);

          }
        }

        break;

      case PYRAMID:


        /*--- Store the connectivity for this element more easily. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_PYRAMID; iNode++) {
          connectivity[iNode] = elem[iElem]->GetNode(iNode);
        }

        /*--- Decide whether we need to store the adjacency for any nodes
         in the current element, i.e., check if any of the nodes have a
         global index value within the range of our linear partitioning. ---*/

        for (unsigned long iNode = 0; iNode < N_POINTS_PYRAMID; iNode++) {

          const long local_index = connectivity[iNode]-firstIndex;

          if ((local_index >= 0) && (local_index < (long)nPoint)) {

            /*--- This node is within our linear partition.
             Add the neighboring nodes to this nodes' adjacency list. ---*/

            /*--- Build adjacency assuming the VTK connectivity ---*/

            if (iNode < 4) {
              adj_nodes[local_index].push_back(connectivity[(iNode+1)%4]);
              adj_nodes[local_index].push_back(connectivity[(iNode+3)%4]);
              adj_nodes[local_index].push_back(connectivity[4]);
            } else {
              adj_nodes[local_index].push_back(connectivity[0]);
              adj_nodes[local_index].push_back(connectivity[1]);
              adj_nodes[local_index].push_back(connectivity[2]);
              adj_nodes[local_index].push_back(connectivity[3]);
            }

          }
        }

        break;

      default:
        SU2_MPI::Error("Element type not supported!", CURRENT_FUNCTION);
        break;
    }
  }

  /*--- Prepare the adjacency information that ParMETIS will need for
   completing the graph partitioning in parallel. ---*/

  SortAdjacency(config);

#endif
#endif

}

void CPhysicalGeometry::Check_IntElem_Orientation(CConfig *config) {

  unsigned long Point_1, Point_2, Point_3, Point_4, Point_5, Point_6,
  iElem, triangle_flip = 0, quad_flip = 0, tet_flip = 0, prism_flip = 0,
  hexa_flip = 0, pyram_flip = 0;
  su2double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, *Coord_6, a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, test;
  unsigned short iDim;

  /*--- Loop over all the elements ---*/

  for (iElem = 0; iElem < nElem; iElem++) {

    /*--- 2D grid, triangle case ---*/

    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test = a[0]*b[1]-b[0]*a[1];

      if (test < 0.0) {
    	  elem[iElem]->Change_Orientation();
    	  triangle_flip++;
      }
    }

    /*--- 2D grid, quadrilateral case ---*/

    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]); }
      test_1 = a[0]*b[1]-b[0]*a[1];

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
        b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]); }
      test_2 = a[0]*b[1]-b[0]*a[1];

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
        b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]); }
      test_3 = a[0]*b[1]-b[0]*a[1];

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]); }
      test_4 = a[0]*b[1]-b[0]*a[1];

      if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0)) {
        elem[iElem]->Change_Orientation();
        quad_flip++;
      }
    }

    /*--- 3D grid, tetrahedron case ---*/

    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
      if (test < 0.0) {
    	  elem[iElem]->Change_Orientation();
    	  tet_flip++;
      }

    }

    /*--- 3D grid, prism case ---*/

    if (elem[iElem]->GetVTK_Type() == PRISM) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(3); Coord_4 = node[Point_4]->GetCoord();
      Point_5 = elem[iElem]->GetNode(4); Coord_5 = node[Point_5]->GetCoord();
      Point_6 = elem[iElem]->GetNode(5); Coord_6 = node[Point_6]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        c[iDim] = (Coord_4[iDim]-Coord_1[iDim])+
        (Coord_5[iDim]-Coord_2[iDim])+
        (Coord_6[iDim]-Coord_3[iDim]); }

      /*--- The normal vector should point to the interior of the element ---*/

      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_5[iDim]-Coord_4[iDim]);
        b[iDim] = 0.5*(Coord_6[iDim]-Coord_4[iDim]);
        c[iDim] = (Coord_1[iDim]-Coord_4[iDim])+
        (Coord_2[iDim]-Coord_5[iDim])+
        (Coord_3[iDim]-Coord_6[iDim]); }

      /*--- The normal vector should point to the interior of the element ---*/

      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      if ((test_1 < 0.0) || (test_2 < 0.0)) {
          elem[iElem]->Change_Orientation();
          prism_flip++;
      }

    }

    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(5); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(7); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      Point_1 = elem[iElem]->GetNode(1); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(2); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(3); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(6); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      Point_1 = elem[iElem]->GetNode(3); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(0); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(1); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      if ((test_1 < 0.0) || (test_2 < 0.0) || (test_3 < 0.0)
          || (test_4 < 0.0)) {
    	  elem[iElem]->Change_Orientation();
      	  hexa_flip++;
      }

    }

    if (elem[iElem]->GetVTK_Type() == PYRAMID) {

      Point_1 = elem[iElem]->GetNode(0); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(1); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(2); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      Point_1 = elem[iElem]->GetNode(2); Coord_1 = node[Point_1]->GetCoord();
      Point_2 = elem[iElem]->GetNode(3); Coord_2 = node[Point_2]->GetCoord();
      Point_3 = elem[iElem]->GetNode(0); Coord_3 = node[Point_3]->GetCoord();
      Point_4 = elem[iElem]->GetNode(4); Coord_4 = node[Point_4]->GetCoord();

      for (iDim = 0; iDim < nDim; iDim++) {
        a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
        b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
        c[iDim] = Coord_4[iDim]-Coord_1[iDim]; }
      n[0] = a[1]*b[2]-b[1]*a[2];
      n[1] = -(a[0]*b[2]-b[0]*a[2]);
      n[2] = a[0]*b[1]-b[0]*a[1];

      test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

      if ((test_1 < 0.0) || (test_2 < 0.0)) {
          elem[iElem]->Change_Orientation();
      	  pyram_flip++;
      }

    }

  }

#ifdef HAVE_MPI
  unsigned long Mytriangle_flip  = triangle_flip;
  unsigned long Myquad_flip      = quad_flip;
  unsigned long Mytet_flip       = tet_flip;
  unsigned long Myprism_flip     = prism_flip;
  unsigned long Myhexa_flip      = hexa_flip;
  unsigned long Mypyram_flip     = pyram_flip;

  SU2_MPI::Allreduce(&Mytriangle_flip, &triangle_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myquad_flip, &quad_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mytet_flip, &tet_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myprism_flip, &prism_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myhexa_flip, &hexa_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mypyram_flip, &pyram_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE) {
    if (triangle_flip > 0) cout << "There has been a re-orientation of the TRIANGLE volume elements." << endl;
    if (quad_flip > 0) cout << "There has been a re-orientation of the QUADRILATERAL volume elements." << endl;
    if (tet_flip > 0) cout << "There has been a re-orientation of the TETRAHEDRON volume elements." << endl;
    if (prism_flip > 0) cout << "There has been a re-orientation of the PRISM volume elements." << endl;
    if (hexa_flip > 0) cout << "There has been a re-orientation of the HEXAHEDRON volume elements." << endl;
    if (pyram_flip > 0) cout << "There has been a re-orientation of the PYRAMID volume elements." << endl;
  }

}

void CPhysicalGeometry::Check_BoundElem_Orientation(CConfig *config) {

  unsigned long Point_1_Surface, Point_2_Surface, Point_3_Surface, Point_4_Surface,
  iElem_Domain, Point_Domain = 0, Point_Surface, iElem_Surface,
  line_flip = 0, triangle_flip = 0, quad_flip = 0;
  su2double test_1, test_2, test_3, test_4, *Coord_1, *Coord_2, *Coord_3, *Coord_4,
  *Coord_5, a[3] = {0.0,0.0,0.0}, b[3] = {0.0,0.0,0.0}, c[3] = {0.0,0.0,0.0}, n[3] = {0.0,0.0,0.0}, test;
  unsigned short iDim, iMarker, iNode_Domain, iNode_Surface;
  bool find;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) {

      for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {

        iElem_Domain = bound[iMarker][iElem_Surface]->GetDomainElement();
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          find = false;
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) {find = true; break;}
          }
          if (!find) break;
        }

        /*--- 2D grid, line case ---*/

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == LINE) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Coord_3 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
          }
          test = a[0]*b[1]-b[0]*a[1];

          if (test < 0.0) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            line_flip++;
          }

        }

        /*--- 3D grid, triangle case ---*/

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == TRIANGLE) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
          Coord_4 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
            c[iDim] = Coord_4[iDim]-Coord_1[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];

          test = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];
          if (test < 0.0) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            node[Point_3_Surface]->SetFlip_Orientation();
            triangle_flip++;
          }

        }

        /*--- 3D grid, rectangle case ---*/

        if (bound[iMarker][iElem_Surface]->GetVTK_Type() == QUADRILATERAL) {

          Point_1_Surface = bound[iMarker][iElem_Surface]->GetNode(0); Coord_1 = node[Point_1_Surface]->GetCoord();
          Point_2_Surface = bound[iMarker][iElem_Surface]->GetNode(1); Coord_2 = node[Point_2_Surface]->GetCoord();
          Point_3_Surface = bound[iMarker][iElem_Surface]->GetNode(2); Coord_3 = node[Point_3_Surface]->GetCoord();
          Point_4_Surface = bound[iMarker][iElem_Surface]->GetNode(3); Coord_4 = node[Point_4_Surface]->GetCoord();
          Coord_5 = node[Point_Domain]->GetCoord();

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_2[iDim]-Coord_1[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_1[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_1[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_1 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_3[iDim]-Coord_2[iDim]);
            b[iDim] = 0.5*(Coord_4[iDim]-Coord_2[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_2[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_2 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_4[iDim]-Coord_3[iDim]);
            b[iDim] = 0.5*(Coord_1[iDim]-Coord_3[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_3[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_3 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          for (iDim = 0; iDim < nDim; iDim++) {
            a[iDim] = 0.5*(Coord_1[iDim]-Coord_4[iDim]);
            b[iDim] = 0.5*(Coord_3[iDim]-Coord_4[iDim]);
            c[iDim] = Coord_5[iDim]-Coord_4[iDim];
          }
          n[0] = a[1]*b[2]-b[1]*a[2];
          n[1] = -(a[0]*b[2]-b[0]*a[2]);
          n[2] = a[0]*b[1]-b[0]*a[1];
          test_4 = n[0]*c[0]+n[1]*c[1]+n[2]*c[2];

          if ((test_1 < 0.0) && (test_2 < 0.0) && (test_3 < 0.0) && (test_4 < 0.0)) {
            bound[iMarker][iElem_Surface]->Change_Orientation();
            node[Point_1_Surface]->SetFlip_Orientation();
            node[Point_2_Surface]->SetFlip_Orientation();
            node[Point_3_Surface]->SetFlip_Orientation();
            node[Point_4_Surface]->SetFlip_Orientation();
            quad_flip++;
          }

        }
      }
    }
  }

#ifdef HAVE_MPI
  unsigned long Myline_flip   = line_flip;
  unsigned long Mytriangle_flip  = triangle_flip;
  unsigned long Myquad_flip   = quad_flip;
  SU2_MPI::Allreduce(&Myline_flip, &line_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Mytriangle_flip, &triangle_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&Myquad_flip, &quad_flip, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (rank == MASTER_NODE) {
    if (line_flip > 0) cout << "There has been a re-orientation of the LINE surface elements." << endl;
    if (triangle_flip > 0) cout << "There has been a re-orientation of the TRIANGLE surface elements." << endl;
    if (quad_flip > 0) cout << "There has been a re-orientation of the QUADRILATERAL surface elements." << endl;
  }

}

void CPhysicalGeometry::SetPositive_ZArea(CConfig *config) {
  unsigned short iMarker, Boundary, Monitoring;
  unsigned long iVertex, iPoint;
  su2double *Normal, PositiveXArea, PositiveYArea, PositiveZArea, WettedArea, CoordX = 0.0, CoordY = 0.0, CoordZ = 0.0, MinCoordX = 1E10, MinCoordY = 1E10,
  MinCoordZ = 1E10, MaxCoordX = -1E10, MaxCoordY = -1E10, MaxCoordZ = -1E10, TotalMinCoordX = 1E10, TotalMinCoordY = 1E10,
  TotalMinCoordZ = 1E10, TotalMaxCoordX = -1E10, TotalMaxCoordY = -1E10, TotalMaxCoordZ = -1E10;
  su2double TotalPositiveXArea = 0.0, TotalPositiveYArea = 0.0, TotalPositiveZArea = 0.0, TotalWettedArea = 0.0, AxiFactor;

  bool axisymmetric = config->GetAxisymmetric();
  bool fea = ((config->GetKind_Solver() == FEM_ELASTICITY) || (config->GetKind_Solver() == DISC_ADJ_FEM));

  PositiveXArea = 0.0;
  PositiveYArea = 0.0;
  PositiveZArea = 0.0;
  WettedArea = 0.0;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    Boundary = config->GetMarker_All_KindBC(iMarker);
    Monitoring = config->GetMarker_All_Monitoring(iMarker);

    if ((((Boundary == EULER_WALL)              ||
          (Boundary == HEAT_FLUX)               ||
          (Boundary == ISOTHERMAL)              ||
          (Boundary == LOAD_BOUNDARY)           ||
          (Boundary == DISPLACEMENT_BOUNDARY)) && (Monitoring == YES))
        || (fea))

      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();

        if (node[iPoint]->GetDomain()) {
          Normal = vertex[iMarker][iVertex]->GetNormal();
          CoordX = node[iPoint]->GetCoord(0);
          CoordY = node[iPoint]->GetCoord(1);
          if (nDim == 3) CoordZ = node[iPoint]->GetCoord(2);

          if (axisymmetric) AxiFactor = 2.0*PI_NUMBER*node[iPoint]->GetCoord(1);
          else AxiFactor = 1.0;

          if (nDim == 2) WettedArea += AxiFactor * sqrt (Normal[0]*Normal[0] + Normal[1]*Normal[1]);
          if (nDim == 3) WettedArea += sqrt (Normal[0]*Normal[0] + Normal[1]*Normal[1] + Normal[2]*Normal[2]);

          if (Normal[0] < 0) PositiveXArea -= Normal[0];
          if (Normal[1] < 0) PositiveYArea -= Normal[1];
          if ((nDim == 3) && (Normal[2] < 0)) PositiveZArea -= Normal[2];

          if (CoordX < MinCoordX) MinCoordX = CoordX;
          if (CoordX > MaxCoordX) MaxCoordX = CoordX;

          if (CoordY < MinCoordY) MinCoordY = CoordY;
          if (CoordY > MaxCoordY) MaxCoordY = CoordY;

          if (nDim == 3) {
            if (CoordZ < MinCoordZ) MinCoordZ = CoordZ;
            if (CoordZ > MaxCoordZ) MaxCoordZ = CoordZ;
          }

        }
      }

  }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&PositiveXArea, &TotalPositiveXArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&PositiveYArea, &TotalPositiveYArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&PositiveZArea, &TotalPositiveZArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  SU2_MPI::Allreduce(&MinCoordX, &TotalMinCoordX, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MinCoordY, &TotalMinCoordY, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MinCoordZ, &TotalMinCoordZ, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  SU2_MPI::Allreduce(&MaxCoordX, &TotalMaxCoordX, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MaxCoordY, &TotalMaxCoordY, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&MaxCoordZ, &TotalMaxCoordZ, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  SU2_MPI::Allreduce(&WettedArea, &TotalWettedArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  TotalPositiveXArea = PositiveXArea;
  TotalPositiveYArea = PositiveYArea;
  TotalPositiveZArea = PositiveZArea;

  TotalMinCoordX = MinCoordX;
  TotalMinCoordY = MinCoordY;
  TotalMinCoordZ = MinCoordZ;

  TotalMaxCoordX = MaxCoordX;
  TotalMaxCoordY = MaxCoordY;
  TotalMaxCoordZ = MaxCoordZ;

  TotalWettedArea    = WettedArea;
#endif

  /*--- Set a reference area if no value is provided ---*/

  if (config->GetRefArea() == 0.0) {

    if (nDim == 3) config->SetRefArea(TotalPositiveZArea);
    else config->SetRefArea(TotalPositiveYArea);

    if (rank == MASTER_NODE) {
      if (nDim == 3) {
        cout << "Reference area = "<< TotalPositiveZArea;
        if (config->GetSystemMeasurements() == SI) cout <<" m^2." << endl; else cout <<" ft^2." << endl;
      }
      else {
        cout << "Reference length = "<< TotalPositiveYArea;
        if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl;
      }
    }

  }

  /*--- Set a semi-span value if no value is provided ---*/

  if (config->GetSemiSpan() == 0.0) {

    if (nDim == 3) config->SetSemiSpan(fabs(TotalMaxCoordY));
    else config->SetSemiSpan(1.0);

    if ((nDim == 3) && (rank == MASTER_NODE)) {
      cout << "Semi-span length = "<< TotalMaxCoordY;
      if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl;
    }

  }

  if (rank == MASTER_NODE) {

    if (fea) cout << "Surface area = "<< TotalWettedArea;
    else cout << "Wetted area = "<< TotalWettedArea;

    if ((nDim == 3) || (axisymmetric)) { if (config->GetSystemMeasurements() == SI) cout <<" m^2." << endl; else cout <<" ft^2." << endl; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl; }

    cout << "Area projection in the x-plane = "<< TotalPositiveXArea;
    if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout <<" m^2,"; else cout <<" ft^2,"; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m,"; else cout <<" ft,"; }

    cout << " y-plane = "<< TotalPositiveYArea;
    if (nDim == 3) { if (config->GetSystemMeasurements() == SI) cout <<" m^2,"; else cout <<" ft^2,"; }
    else { if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft." << endl; }

    if (nDim == 3) { cout << " z-plane = "<< TotalPositiveZArea;
      if (config->GetSystemMeasurements() == SI) cout <<" m^2." << endl; else cout <<" ft^2."<< endl; }

    cout << "Max. coordinate in the x-direction = "<< TotalMaxCoordX;
    if (config->GetSystemMeasurements() == SI) cout <<" m,"; else cout <<" ft,";

    cout << " y-direction = "<< TotalMaxCoordY;
    if (config->GetSystemMeasurements() == SI) cout <<" m"; else cout <<" ft";

    if (nDim == 3) {
    	cout << ", z-direction = "<< TotalMaxCoordZ;
      if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft."<< endl;
    }
    else cout << "." << endl;

    cout << "Min. coordinate in the x-direction = "<< TotalMinCoordX;
    if (config->GetSystemMeasurements() == SI) cout <<" m,"; else cout <<" ft";

    cout << " y-direction = "<< TotalMinCoordY;
    if (config->GetSystemMeasurements() == SI) cout <<" m"; else cout <<" ft";

    if (nDim == 3) {
    	cout << ", z-direction = "<< TotalMinCoordZ;
      if (config->GetSystemMeasurements() == SI) cout <<" m." << endl; else cout <<" ft."<< endl;
    }
    else cout << "." << endl;

  }

}

void CPhysicalGeometry::SetPoint_Connectivity(void) {

  unsigned short Node_Neighbor, iNode, iNeighbor;
  unsigned long jElem, Point_Neighbor, iPoint, iElem;

  /*--- Loop over all the elements ---*/

  for (iElem = 0; iElem < nElem; iElem++)

  /*--- Loop over all the nodes of an element ---*/

    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);

      /*--- Store the element into the point ---*/

      node[iPoint]->SetElem(iElem);
    }

  /*--- Loop over all the points ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)

  /*--- Loop over all elements shared by the point ---*/

    for (iElem = 0; iElem < node[iPoint]->GetnElem(); iElem++) {

      jElem = node[iPoint]->GetElem(iElem);

      /*--- If we find the point iPoint in the surronding element ---*/

      for (iNode = 0; iNode < elem[jElem]->GetnNodes(); iNode++)

        if (elem[jElem]->GetNode(iNode) == iPoint)

        /*--- Localize the local index of the neighbor of iPoint in the element ---*/

          for (iNeighbor = 0; iNeighbor < elem[jElem]->GetnNeighbor_Nodes(iNode); iNeighbor++) {
            Node_Neighbor = elem[jElem]->GetNeighbor_Nodes(iNode, iNeighbor);
            Point_Neighbor = elem[jElem]->GetNode(Node_Neighbor);

            /*--- Store the point into the point ---*/

            node[iPoint]->SetPoint(Point_Neighbor);
          }
    }

  /*--- Set the number of neighbors variable, this is
   important for JST and multigrid in parallel ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetnNeighbor(node[iPoint]->GetnPoint());

}

void CPhysicalGeometry::SetRCM_Ordering(CConfig *config) {
  unsigned long iPoint, AdjPoint, AuxPoint, AddPoint, iElem, iNode, jNode;
  vector<unsigned long> Queue, AuxQueue, Result;
  unsigned short Degree, MinDegree, iDim, iMarker;
  bool *inQueue;

  inQueue = new bool [nPoint];

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    inQueue[iPoint] = false;

  /*--- Select the node with the lowest degree in the grid. ---*/

  MinDegree = node[0]->GetnNeighbor(); AddPoint = 0;
  for (iPoint = 1; iPoint < nPointDomain; iPoint++) {
    Degree = node[iPoint]->GetnPoint();
    if (Degree < MinDegree) { MinDegree = Degree; AddPoint = iPoint; }
  }

  /*--- Add the node in the first free position. ---*/

  Result.push_back(AddPoint); inQueue[AddPoint] = true;

  /*--- Loop until reorganize all the nodes ---*/

  do {

    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/

    AuxQueue.clear();
    for (iNode = 0; iNode < node[AddPoint]->GetnPoint(); iNode++) {
      AdjPoint = node[AddPoint]->GetPoint(iNode);
      if ((!inQueue[AdjPoint]) && (AdjPoint < nPointDomain)) {
        AuxQueue.push_back(AdjPoint);
      }
    }

    if (AuxQueue.size() != 0) {

      /*--- Sort the auxiliar queue based on the number of neighbors ---*/

      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        for (jNode = 0; jNode < AuxQueue.size() - 1 - iNode; jNode++) {
          if (node[AuxQueue[jNode]]->GetnPoint() > node[AuxQueue[jNode+1]]->GetnPoint()) {
            AuxPoint = AuxQueue[jNode];
            AuxQueue[jNode] = AuxQueue[jNode+1];
            AuxQueue[jNode+1] = AuxPoint;
          }
        }
      }

      Queue.insert(Queue.end(), AuxQueue.begin(), AuxQueue.end());
      for (iNode = 0; iNode < AuxQueue.size(); iNode++) {
        inQueue[AuxQueue[iNode]] = true;
      }

    }

    /*--- Extract the first node from the queue and add it in the first free
     position. ---*/

    if (Queue.size() != 0) {
      AddPoint = Queue[0];
      Result.push_back(Queue[0]);
      Queue.erase (Queue.begin(), Queue.begin()+1);
    }

    /*--- Add to the queue all the nodes adjacent in the increasing
     order of their degree, checking if the element is already
     in the Queue. ---*/

  } while (Queue.size() != 0);

  /*--- Check that all the points have been added ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    if (inQueue[iPoint] == false) Result.push_back(iPoint);
  }

  delete[] inQueue;

  reverse(Result.begin(), Result.end());

  /*--- Add the MPI points ---*/

  for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
    Result.push_back(iPoint);
  }

  /*--- Reset old data structures ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->ResetElem();
    node[iPoint]->ResetPoint();
    node[iPoint]->ResetBoundary();
    node[iPoint]->SetPhysicalBoundary(false);
    node[iPoint]->SetSolidBoundary(false);
    node[iPoint]->SetPeriodicBoundary(false);
    node[iPoint]->SetDomain(true);
  }

  /*--- Set the new coordinates ---*/

  su2double **AuxCoord;
  unsigned long *AuxGlobalIndex;

  AuxGlobalIndex = new unsigned long [nPoint];
  AuxCoord = new su2double* [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    AuxCoord[iPoint] = new su2double [nDim];

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    AuxGlobalIndex[iPoint] = node[iPoint]->GetGlobalIndex();
    for (iDim = 0; iDim < nDim; iDim++) {
      AuxCoord[iPoint][iDim] = node[iPoint]->GetCoord(iDim);
    }
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint]->SetGlobalIndex(AuxGlobalIndex[Result[iPoint]]);
    for (iDim = 0; iDim < nDim; iDim++)
      node[iPoint]->SetCoord(iDim, AuxCoord[Result[iPoint]][iDim]);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    delete[] AuxCoord[iPoint];
  delete[] AuxCoord;
  delete[] AuxGlobalIndex;

  /*--- Set the new conectivities ---*/

  unsigned long *InvResult;
  InvResult = new unsigned long [nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    InvResult[Result[iPoint]] = iPoint;

  for (iElem = 0; iElem < nElem; iElem++) {
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      iPoint = elem[iElem]->GetNode(iNode);
      elem[iElem]->SetNode(iNode, InvResult[iPoint]);
    }
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {

      string Marker_Tag = config->GetMarker_All_TagBound(iMarker);
      if (Marker_Tag == "SEND_RECEIVE") {
        for (unsigned long iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          if (config->GetMarker_All_SendRecv(iMarker) < 0)
            node[bound[iMarker][iElem_Bound]->GetNode(0)]->SetDomain(false);
        }
      }

      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        bound[iMarker][iElem]->SetNode(iNode, InvResult[iPoint]);
        node[InvResult[iPoint]]->SetBoundary(nMarker);
        if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
            config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY &&
            config->GetMarker_All_KindBC(iMarker) != PERIODIC_BOUNDARY)
          node[InvResult[iPoint]]->SetPhysicalBoundary(true);

        if (config->GetSolid_Wall(iMarker))
          node[InvResult[iPoint]]->SetSolidBoundary(true);

        if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY)
          node[InvResult[iPoint]]->SetPeriodicBoundary(true);
      }
    }
  }


  delete[] InvResult;

}

void CPhysicalGeometry::SetElement_Connectivity(void) {
  unsigned short first_elem_face, second_elem_face, iFace, iNode, jElem;
  unsigned long face_point, Test_Elem, iElem;

  /*--- Loop over all the elements, faces and nodes ---*/

  for (iElem = 0; iElem < nElem; iElem++)
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++)
      for (iNode = 0; iNode < elem[iElem]->GetnNodesFace(iFace); iNode++) {
        face_point = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iNode));

        /*--- Loop over all elements sharing the face point ---*/

        for (jElem = 0; jElem < node[face_point]->GetnElem(); jElem++) {
          Test_Elem = node[face_point]->GetElem(jElem);

          /*--- If it is a new element in this face ---*/

          if ((elem[iElem]->GetNeighbor_Elements(iFace) == -1) && (iElem < Test_Elem) &&
              (FindFace(iElem, Test_Elem, first_elem_face, second_elem_face))) {

            /*--- Localice which faces are sharing both elements ---*/

            elem[iElem]->SetNeighbor_Elements(Test_Elem, first_elem_face);

            /*--- Store the element for both elements ---*/

            elem[Test_Elem]->SetNeighbor_Elements(iElem, second_elem_face);

          }
        }
      }
}

void CPhysicalGeometry::SetBoundVolume(void) {
  unsigned short cont, iMarker, iElem, iNode_Domain, iNode_Surface;
  unsigned long Point_Domain, Point_Surface, Point, iElem_Surface, iElem_Domain;
  bool CheckVol;

  for (iMarker = 0; iMarker < nMarker; iMarker++)
    for (iElem_Surface = 0; iElem_Surface < nElem_Bound[iMarker]; iElem_Surface++) {

      /*--- Choose and arbitrary point from the surface --*/
      Point = bound[iMarker][iElem_Surface]->GetNode(0);
      CheckVol = false;

      for (iElem = 0; iElem < node[Point]->GetnElem(); iElem++) {
        /*--- Look for elements surronding that point --*/
        cont = 0; iElem_Domain = node[Point]->GetElem(iElem);
        for (iNode_Domain = 0; iNode_Domain < elem[iElem_Domain]->GetnNodes(); iNode_Domain++) {
          Point_Domain = elem[iElem_Domain]->GetNode(iNode_Domain);
          for (iNode_Surface = 0; iNode_Surface < bound[iMarker][iElem_Surface]->GetnNodes(); iNode_Surface++) {
            Point_Surface = bound[iMarker][iElem_Surface]->GetNode(iNode_Surface);
            if (Point_Surface == Point_Domain) cont++;
            if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
          }
          if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) break;
        }

        if (cont == bound[iMarker][iElem_Surface]->GetnNodes()) {
          bound[iMarker][iElem_Surface]->SetDomainElement(iElem_Domain);
          CheckVol = true;
          break;
        }
      }
      if (!CheckVol) {
        char buf[100];
        SPRINTF(buf,"The surface element (%u, %lu) doesn't have an associated volume element", iMarker, iElem_Surface );
        SU2_MPI::Error(buf, CURRENT_FUNCTION);
      }
    }
}

void CPhysicalGeometry::SetVertex(CConfig *config) {
  unsigned long  iPoint, iVertex, iElem;
  unsigned short iMarker, iNode;

  /*--- Initialize the Vertex vector for each node of the grid ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1, iMarker);

  /*--- Create and compute the vector with the number of vertex per marker ---*/

  nVertex = new unsigned long [nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Initialize the number of Bound Vertex for each Marker ---*/

    nVertex[iMarker] = 0;
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);

        /*--- Set the vertex in the node information ---*/

        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)) {
          node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
          nVertex[iMarker]++;
        }
      }
  }

  /*--- Initialize the Vertex vector for each node, the previous result is deleted ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      node[iPoint]->SetVertex(-1, iMarker);

  /*--- Create the bound vertex structure, note that the order
   is the same as in the input file, this is important for Send/Receive part ---*/

  vertex = new CVertex**[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    vertex[iMarker] = new CVertex* [nVertex[iMarker]];
    nVertex[iMarker] = 0;

    /*--- Initialize the number of Bound Vertex for each Marker ---*/

    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);

        /*--- Set the vertex in the node information ---*/

        if ((node[iPoint]->GetVertex(iMarker) == -1) || (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE)) {
          iVertex = nVertex[iMarker];
          vertex[iMarker][iVertex] = new CVertex(iPoint, nDim);

          if (config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) {
            vertex[iMarker][iVertex]->SetRotation_Type(bound[iMarker][iElem]->GetRotation_Type());
          }
          node[iPoint]->SetVertex(nVertex[iMarker], iMarker);
          nVertex[iMarker]++;
        }
      }
  }
}

void CPhysicalGeometry::ComputeNSpan(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {
  unsigned short iMarker, jMarker, iMarkerTP, iSpan, jSpan, kSpan = 0;
  unsigned long iPoint, iVertex;
  long jVertex;
  int nSpan, nSpan_loc;
  su2double *coord, *valueSpan, min, max, radius, delta;
  short PeriodicBoundary;
  unsigned short SpanWise_Kind = config->GetKind_SpanWise();

#ifdef HAVE_MPI
  unsigned short iSize;
  int nSpan_max;
  int My_nSpan, My_MaxnSpan, *My_nSpan_loc = NULL;
  su2double MyMin, MyMax, *MyTotValueSpan =NULL,*MyValueSpan =NULL;
#endif

  nSpan = 0;
  nSpan_loc = 0;
  if (nDim == 2){
    nSpanWiseSections[marker_flag-1] = 1;
    //TODO (turbo) make it more genral
    if(marker_flag == OUTFLOW)	config->SetnSpanWiseSections(1);

    /*---Initilize the vector of span-wise values that will be ordered ---*/
    SpanWiseValue[marker_flag -1] = new su2double[1];
    for (iSpan = 0; iSpan < 1; iSpan++){
      SpanWiseValue[marker_flag -1][iSpan] = 0;
    }
  }
  else{
    if(SpanWise_Kind == AUTOMATIC){
      /*--- loop to find inflow of outflow marker---*/
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
              for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint = vertex[iMarker][iVertex]->GetNode();

                /*--- loop to find the vertex that ar both of inflow or outflow marker and on the periodic
                 * in order to caount the number of Span ---*/
                for (jMarker = 0; jMarker < nMarker; jMarker++){
                  if (config->GetMarker_All_KindBC(jMarker) == PERIODIC_BOUNDARY) {
                    PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
                    jVertex = node[iPoint]->GetVertex(jMarker);
                    if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){
                        nSpan++;
                    }
                  }
                }
              }
            }
          }
        }
      }

      /*--- storing the local number of span---*/
      nSpan_loc = nSpan;
      /*--- if parallel computing the global number of span---*/
#ifdef HAVE_MPI
      nSpan_max = nSpan;
      My_nSpan  = nSpan; nSpan = 0;
      My_MaxnSpan = nSpan_max; nSpan_max = 0;
      SU2_MPI::Allreduce(&My_nSpan, &nSpan, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&My_MaxnSpan, &nSpan_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif

      /*--- initialize the vector that will contain the disordered values span-wise ---*/
      nSpanWiseSections[marker_flag -1] = nSpan;
      valueSpan = new su2double[nSpan];

      for (iSpan = 0; iSpan < nSpan; iSpan ++ ){
        valueSpan[iSpan] = -1001.0;
      }


      /*--- store the local span-wise value for each processor ---*/
      nSpan_loc = 0;
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
              for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint = vertex[iMarker][iVertex]->GetNode();
                for (jMarker = 0; jMarker < nMarker; jMarker++){
                  if (config->GetMarker_All_KindBC(jMarker) == PERIODIC_BOUNDARY) {
                    PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
                  	jVertex = node[iPoint]->GetVertex(jMarker);
                    if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){
                      coord = node[iPoint]->GetCoord();
                      switch (config->GetKind_TurboMachinery(val_iZone)){
                      case CENTRIFUGAL:
                        valueSpan[nSpan_loc] = coord[2];
                        break;
                      case CENTRIPETAL:
                        valueSpan[nSpan_loc] = coord[2];
                        break;
                      case AXIAL:
                        valueSpan[nSpan_loc] = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                        break;
                      case CENTRIPETAL_AXIAL:
                        if (marker_flag == OUTFLOW){
                          valueSpan[nSpan_loc] = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                        }
                        else{
                          valueSpan[nSpan_loc] = coord[2];
                        }
                        break;
                      case AXIAL_CENTRIFUGAL:
                        if (marker_flag == INFLOW){
                          valueSpan[nSpan_loc] = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                        }
                        else{
                          valueSpan[nSpan_loc] = coord[2];
                        }
                        break;

                      }
                      nSpan_loc++;
                    }
                  }
                }
              }
            }
          }
        }
      }

      /*--- Gather the span-wise values on all the processor ---*/

#ifdef HAVE_MPI
      MyTotValueSpan  = new su2double[nSpan_max*size];
      MyValueSpan     = new su2double[nSpan_max];
      My_nSpan_loc    = new int[size];
      for(iSpan = 0; iSpan < nSpan_max; iSpan++){
        MyValueSpan[iSpan] = -1001.0;
        for (iSize = 0; iSize< size; iSize++){
          MyTotValueSpan[iSize*nSpan_max + iSpan] = -1001.0;
        }
      }

      for(iSpan = 0; iSpan <nSpan_loc; iSpan++){
        MyValueSpan[iSpan] = valueSpan[iSpan];
      }

      for(iSpan = 0; iSpan <nSpan; iSpan++){
        valueSpan[iSpan] = -1001.0;
      }

      SU2_MPI::Allgather(MyValueSpan, nSpan_max , MPI_DOUBLE, MyTotValueSpan, nSpan_max, MPI_DOUBLE, MPI_COMM_WORLD);
      SU2_MPI::Allgather(&nSpan_loc, 1 , MPI_INT, My_nSpan_loc, 1, MPI_INT, MPI_COMM_WORLD);

      jSpan = 0;
      for (iSize = 0; iSize< size; iSize++){
        for(iSpan = 0; iSpan < My_nSpan_loc[iSize]; iSpan++){
          valueSpan[jSpan] = MyTotValueSpan[iSize*nSpan_max + iSpan];
          jSpan++;
        }
      }

      delete [] MyTotValueSpan; delete [] MyValueSpan; delete [] My_nSpan_loc;

#endif

      // check if the value are gathered correctly
      //
      //  for (iSpan = 0; iSpan < nSpan; iSpan++){
      //  	if(rank == MASTER_NODE){
      //  		cout << setprecision(16)<<  iSpan +1 << " with a value of " <<valueSpan[iSpan]<< " at flag " << marker_flag <<endl;
      //  	}
      //  }


      /*--- Find the minimum value among the span-wise values  ---*/
      min = 10.0E+06;
      for (iSpan = 0; iSpan < nSpan; iSpan++){
        if(valueSpan[iSpan]< min) min = valueSpan[iSpan];
      }

      /*---Initilize the vector of span-wise values that will be ordered ---*/
      SpanWiseValue[marker_flag -1] = new su2double[nSpan];
      for (iSpan = 0; iSpan < nSpan; iSpan++){
        SpanWiseValue[marker_flag -1][iSpan] = 0;
      }

      /*---Ordering the vector of span-wise values---*/
      SpanWiseValue[marker_flag -1][0] = min;
      for (iSpan = 1; iSpan < nSpan; iSpan++){
        min = 10.0E+06;
        for (jSpan = 0; jSpan < nSpan; jSpan++){
          if((valueSpan[jSpan] - SpanWiseValue[marker_flag -1][iSpan-1]) < min && (valueSpan[jSpan] - SpanWiseValue[marker_flag -1][iSpan-1]) > EPS){
            min    = valueSpan[jSpan] - SpanWiseValue[marker_flag -1][iSpan-1];
            kSpan = jSpan;
          }
        }
        SpanWiseValue[marker_flag -1][iSpan] = valueSpan[kSpan];
      }

      delete [] valueSpan;
    }
    /*--- Compute equispaced Span-wise sections using number of section specified by the User---*/
    else{
      /*--- Initialize number of span---*/
      nSpanWiseSections[marker_flag-1] = config->Get_nSpanWiseSections_User();
      SpanWiseValue[marker_flag -1] = new su2double[config->Get_nSpanWiseSections_User()];
      for (iSpan = 0; iSpan < config->Get_nSpanWiseSections_User(); iSpan++){
        SpanWiseValue[marker_flag -1][iSpan] = 0;
      }
      /*--- Compute maximum and minimum value span-wise---*/
      min = 10.0E+06;
      max = -10.0E+06;
      for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
          if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
            if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
              for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
                iPoint = vertex[iMarker][iVertex]->GetNode();
                for (jMarker = 0; jMarker < nMarker; jMarker++){
                  if (config->GetMarker_All_KindBC(jMarker) == PERIODIC_BOUNDARY) {
                    PeriodicBoundary = config->GetMarker_All_PerBound(jMarker);
                    jVertex = node[iPoint]->GetVertex(jMarker);
                    if ((jVertex != -1) && (PeriodicBoundary == (val_iZone + 1))){
                      coord = node[iPoint]->GetCoord();
                      switch (config->GetKind_TurboMachinery(val_iZone)){
                      case CENTRIFUGAL: case CENTRIPETAL:
                        if (coord[2] < min) min = coord[2];
                        if (coord[2] > max) max = coord[2];
                        break;
                      case AXIAL:
                        radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                        if (radius < min) min = radius;
                        if (radius > max) max = radius;
                        break;
                      case CENTRIPETAL_AXIAL:
                        if (marker_flag == OUTFLOW){
                          radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                          if (radius < min) min = radius;
                          if (radius > max) max = radius;
                        }
                        else{
                          if (coord[2] < min) min = coord[2];
                          if (coord[2] > max) max = coord[2];
                        }
                        break;

                      case AXIAL_CENTRIFUGAL:
                        if (marker_flag == INFLOW){
                          radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                          if (radius < min) min = radius;
                          if (radius > max) max = radius;
                        }
                        else{
                          if (coord[2] < min) min = coord[2];
                          if (coord[2] > max) max = coord[2];
                        }
                        break;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      /*--- compute global minimum and maximum value on span-wise ---*/
#ifdef HAVE_MPI
      MyMin= min;  min = 0;
      MyMax= max;  max = 0;
      SU2_MPI::Allreduce(&MyMin, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      SU2_MPI::Allreduce(&MyMax, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

      //  cout <<"min  " <<  min << endl;
      //  cout <<"max  " << max << endl;
      /*--- compute height value for each spanwise section---*/
      delta = (max - min)/(nSpanWiseSections[marker_flag-1] -1);
      for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
        SpanWiseValue[marker_flag - 1][iSpan]= min + delta*iSpan;
      }
    }


    if(marker_flag == OUTFLOW){
      if(nSpanWiseSections[INFLOW -1] != nSpanWiseSections[OUTFLOW - 1]){
        char buf[100];
        SPRINTF(buf, "nSpan inflow %u, nSpan outflow %u", nSpanWiseSections[INFLOW], nSpanWiseSections[OUTFLOW]);
        SU2_MPI::Error(string(" At the moment only turbomachinery with the same amount of span-wise section can be simulated\n") + buf, CURRENT_FUNCTION);
      }
      else{
        config->SetnSpanWiseSections(nSpanWiseSections[OUTFLOW -1]);
      }
    }

  }

}
void CPhysicalGeometry::SetTurboVertex(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {
  unsigned long  iPoint, **ordered, **disordered, **oldVertex3D, iInternalVertex;
  unsigned long nVert, nVertMax;
  unsigned short iMarker, iMarkerTP, iSpan, jSpan, iDim;
  su2double min, minInt, max, *coord, dist, Normal2, *TurboNormal, *NormalArea, target = 0.0, **area, ***unitnormal, Area = 0.0;
  bool **checkAssign;
  min    =  10.0E+06;
  minInt =  10.0E+06;
  max    = -10.0E+06;

  su2double radius;
  long iVertex, iSpanVertex, jSpanVertex, kSpanVertex = 0;
  int *nTotVertex_gb, *nVertexSpanHalo;
  su2double **x_loc, **y_loc, **z_loc, **angCoord_loc, **deltaAngCoord_loc, **angPitch, **deltaAngPitch, *minIntAngPitch,
  *minAngPitch, *maxAngPitch;
  int       **rank_loc;
#ifdef HAVE_MPI
  unsigned short iSize, kSize = 0, jSize;
  su2double MyMin,MyIntMin, MyMax;
  su2double *x_gb = NULL, *y_gb = NULL, *z_gb = NULL, *angCoord_gb = NULL, *deltaAngCoord_gb = NULL;
  bool *checkAssign_gb =NULL;
  unsigned long My_nVert;

#endif
  string multizone_filename;

  x_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  y_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  z_loc              = new su2double*[nSpanWiseSections[marker_flag-1]];
  angCoord_loc       = new su2double*[nSpanWiseSections[marker_flag-1]];
  deltaAngCoord_loc  = new su2double*[nSpanWiseSections[marker_flag-1]];
  angPitch           = new su2double*[nSpanWiseSections[marker_flag-1]];
  deltaAngPitch      = new su2double*[nSpanWiseSections[marker_flag-1]];
  rank_loc           = new int*[nSpanWiseSections[marker_flag-1]];
  minAngPitch        = new su2double[nSpanWiseSections[marker_flag-1]];
  minIntAngPitch     = new su2double[nSpanWiseSections[marker_flag-1]];
  maxAngPitch        = new su2double[nSpanWiseSections[marker_flag-1]];

  nTotVertex_gb      = new int[nSpanWiseSections[marker_flag-1]];
  nVertexSpanHalo    = new int[nSpanWiseSections[marker_flag-1]];
  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    nTotVertex_gb[iSpan]   = -1;
    nVertexSpanHalo[iSpan] = 0;
    minAngPitch[iSpan]     = 10.0E+06;
    minIntAngPitch[iSpan]  = 10.0E+06;
    maxAngPitch[iSpan]     = -10.0E+06;
  }

  /*--- Initialize auxiliary pointers ---*/
  TurboNormal        = new su2double[3];
  NormalArea         = new su2double[3];
  ordered            = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  disordered         = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  oldVertex3D        = new unsigned long* [nSpanWiseSections[marker_flag-1]];
  area               = new su2double* [nSpanWiseSections[marker_flag-1]];
  unitnormal         = new su2double** [nSpanWiseSections[marker_flag-1]];
  checkAssign        = new bool* [nSpanWiseSections[marker_flag-1]];

  /*--- Initialize the new Vertex structure. The if statement ensures that these vectors are initialized
   * only once even if the routine is called more than once.---*/

  if (allocate){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nSpanSectionsByMarker[iMarker]       = nSpanWiseSections[marker_flag-1];
            nVertexSpan[iMarker]                 = new long[nSpanWiseSections[marker_flag-1]];
            turbovertex[iMarker]                 = new CTurboVertex** [nSpanWiseSections[marker_flag-1]];
            nTotVertexSpan[iMarker]              = new unsigned long [nSpanWiseSections[marker_flag-1] +1];
            MaxAngularCoord[iMarker]             = new su2double [nSpanWiseSections[marker_flag-1]];
            MinAngularCoord[iMarker]             = new su2double [nSpanWiseSections[marker_flag-1]];
            MinRelAngularCoord[iMarker]          = new su2double [nSpanWiseSections[marker_flag-1]];
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
              nVertexSpan[iMarker][iSpan]        = 0;
              turbovertex[iMarker][iSpan]        = NULL;
              MinAngularCoord[iMarker][iSpan]    = 10.0E+06;
              MaxAngularCoord[iMarker][iSpan]    = -10.0E+06;
              MinRelAngularCoord[iMarker][iSpan] = 10.0E+06;
            }
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1] +1; iSpan++){
              nTotVertexSpan[iMarker][iSpan]     = 0;
            }
          }
        }
      }
    }
  }

  //this works only for turbomachinery rotating around the Z-Axes.
  // the reordering algorithm pitch-wise assumes that X-coordinate of each boundary vertex is positive so that reordering can be based on the Y-coordinate.
    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){

            /*--- compute the amount of vertexes for each span-wise section to initialize the CTurboVertex pointers and auxiliary pointers  ---*/
            for (iVertex = 0; (unsigned long)iVertex  < nVertex[iMarker]; iVertex++) {
              iPoint = vertex[iMarker][iVertex]->GetNode();
              if (nDim == 3){
                dist = 10E+06;
                jSpan = -1;
                coord = node[iPoint]->GetCoord();

                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL: case CENTRIPETAL:
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case AXIAL:
                  radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case CENTRIPETAL_AXIAL:
                  if (marker_flag == OUTFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;

                case AXIAL_CENTRIFUGAL:
                  if (marker_flag == INFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;
                }
              }

              /*--- 2D problem do not need span-wise separation---*/
              else{
                jSpan = 0;
              }

              if(node[iPoint]->GetDomain()){
                nVertexSpan[iMarker][jSpan]++;
              }
              nVertexSpanHalo[jSpan]++;
            }

            /*--- initialize the CTurboVertex pointers and auxiliary pointers  ---*/
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
              if (allocate){
                turbovertex[iMarker][iSpan] = new CTurboVertex* [nVertexSpan[iMarker][iSpan]];
                for (iVertex = 0; iVertex < nVertexSpan[iMarker][iSpan]; iVertex++){
                  turbovertex[iMarker][iSpan][iVertex] = NULL;
                }
              }
              ordered[iSpan]                           = new unsigned long [nVertexSpanHalo[iSpan]];
              disordered[iSpan]                        = new unsigned long [nVertexSpanHalo[iSpan]];
              oldVertex3D[iSpan]                       = new unsigned long [nVertexSpanHalo[iSpan]];
              checkAssign[iSpan]                       = new bool [nVertexSpanHalo[iSpan]];
              area[iSpan]                              = new su2double [nVertexSpanHalo[iSpan]];
              unitnormal[iSpan]                        = new su2double* [nVertexSpanHalo[iSpan]];
              for (iVertex = 0; iVertex < nVertexSpanHalo[iSpan]; iVertex++){
                unitnormal[iSpan][iVertex]             = new su2double [nDim];
              }
              angPitch[iSpan]                          = new su2double [nVertexSpanHalo[iSpan]];
              deltaAngPitch[iSpan]                     = new su2double [nVertexSpanHalo[iSpan]];
              nVertexSpanHalo[iSpan]                   = 0;
            }

            /*--- store the vertexes in a ordered manner in span-wise directions but not yet ordered pitch-wise ---*/
            for (iVertex = 0; (unsigned long)iVertex < nVertex[iMarker]; iVertex++) {
              iPoint = vertex[iMarker][iVertex]->GetNode();
              if(nDim == 3){
                dist  = 10E+06;
                jSpan = -1;

                coord = node[iPoint]->GetCoord();
                switch (config->GetKind_TurboMachinery(val_iZone)){
                case CENTRIFUGAL: case CENTRIPETAL:
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case AXIAL:
                  radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                    if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                      dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                      jSpan=iSpan;
                    }
                  }
                  break;
                case CENTRIPETAL_AXIAL:
                  if(marker_flag == OUTFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;

                case AXIAL_CENTRIFUGAL:
                  if(marker_flag == INFLOW){
                    radius = sqrt(coord[0]*coord[0]+coord[1]*coord[1]);
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(radius - SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(radius-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }else{
                    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
                      if (dist > (abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]))){
                        dist= abs(coord[2]-SpanWiseValue[marker_flag-1][iSpan]);
                        jSpan=iSpan;
                      }
                    }
                  }
                  break;
                }
              }
              /*--- 2D problem do not need span-wise separation---*/
              else{
                jSpan = 0;
              }
              /*--- compute the face area associated with the vertex ---*/
              vertex[iMarker][iVertex]->GetNormal(NormalArea);
              for (iDim = 0; iDim < nDim; iDim++) NormalArea[iDim] = -NormalArea[iDim];
              Area = 0.0;
              for (iDim = 0; iDim < nDim; iDim++)
                Area += NormalArea[iDim]*NormalArea[iDim];
              Area = sqrt(Area);
              for (iDim = 0; iDim < nDim; iDim++) NormalArea[iDim] /= Area;
              /*--- store all the all the info into the auxiliary containers ---*/
              disordered[jSpan][nVertexSpanHalo[jSpan]]  = iPoint;
              oldVertex3D[jSpan][nVertexSpanHalo[jSpan]] = iVertex;
              area[jSpan][nVertexSpanHalo[jSpan]]        = Area;
              for (iDim = 0; iDim < nDim; iDim++){
                unitnormal[jSpan][nVertexSpanHalo[jSpan]][iDim] = NormalArea[iDim];
              }
              checkAssign[jSpan][nVertexSpanHalo[jSpan]] = false;
              nVertexSpanHalo[jSpan]++;
            }

            /*--- using the auxiliary container reordered the vertexes pitch-wise direction at each span ---*/
            // the reordering algorithm can be based on the Y-coordinate.
            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

              /*--- find the local minimum and maximum pitch-wise for each processor---*/
              min    = 10E+06;
              minInt = 10E+06;
              max    = -10E+06;
              for(iSpanVertex = 0; iSpanVertex < nVertexSpanHalo[iSpan]; iSpanVertex++){
                iPoint = disordered[iSpan][iSpanVertex];
                coord = node[iPoint]->GetCoord();
                /*--- find nodes at minimum pitch among all nodes---*/
                if (coord[1]<min){
                  min = coord[1];
                  if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                    MinAngularCoord[iMarker][iSpan] = coord[1];
                  }
                  else{
                  MinAngularCoord[iMarker][iSpan] = atan(coord[1]/coord[0]);
                  }
                  minAngPitch[iSpan]= MinAngularCoord[iMarker][iSpan];
                  kSpanVertex =iSpanVertex;
                }

                /*--- find nodes at minimum pitch among the internal nodes---*/
                if (coord[1]<minInt){
                  if(node[iPoint]->GetDomain()){
                    minInt = coord[1];
                    if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                      minIntAngPitch[iSpan] = coord[1];
                    }
                    else{
                      minIntAngPitch[iSpan] = atan(coord[1]/coord[0]);
                    }
                  }
                }

                /*--- find nodes at maximum pitch among the internal nodes---*/
                if (coord[1]>max){
                  if(node[iPoint]->GetDomain()){
                    max =coord[1];
                    if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                      MaxAngularCoord[iMarker][iSpan] = coord[1];
                    }
                    else{
                      MaxAngularCoord[iMarker][iSpan] = atan(coord[1]/coord[0]);
                    }
                    maxAngPitch[iSpan]= MaxAngularCoord[iMarker][iSpan];
                  }
                }
              }

              iInternalVertex = 0;

              /*--- reordering the vertex pitch-wise, store the ordered vertexes span-wise and pitch-wise---*/
              for(iSpanVertex = 0; iSpanVertex<nVertexSpanHalo[iSpan]; iSpanVertex++){
                dist = 10E+06;
                ordered[iSpan][iSpanVertex] = disordered[iSpan][kSpanVertex];
                checkAssign[iSpan][kSpanVertex] = true;
                coord = node[ordered[iSpan][iSpanVertex]]->GetCoord();
                target = coord[1];
                if (nDim == 2 && config->GetKind_TurboMachinery(val_iZone) == AXIAL){
                   angPitch[iSpan][iSpanVertex]=coord[1];
                }
                else{
                  angPitch[iSpan][iSpanVertex]=atan(coord[1]/coord[0]);
                }
                if(iSpanVertex == 0){
                  deltaAngPitch[iSpan][iSpanVertex]=0.0;
                }
                else{
                  deltaAngPitch[iSpan][iSpanVertex]= angPitch[iSpan][iSpanVertex] - angPitch[iSpan][iSpanVertex - 1];
                }
                /*---create turbovertex structure only for the internal nodes---*/
                if(node[ordered[iSpan][iSpanVertex]]->GetDomain()){
                  if (allocate){
                    turbovertex[iMarker][iSpan][iInternalVertex] = new CTurboVertex(ordered[iSpan][iSpanVertex], nDim);
                  }
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetArea(area[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetNormal(unitnormal[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetOldVertex(oldVertex3D[iSpan][kSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetAngularCoord(angPitch[iSpan][iSpanVertex]);
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetDeltaAngularCoord(deltaAngPitch[iSpan][iSpanVertex]);
                  switch (config->GetKind_TurboMachinery(val_iZone)){
                  case CENTRIFUGAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = -coord[0]/sqrt(Normal2);
                      TurboNormal[1] = -coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;
                  case CENTRIPETAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == OUTFLOW){
                      TurboNormal[0] = -coord[0]/sqrt(Normal2);
                      TurboNormal[1] = -coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;
                  case AXIAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if(nDim == 3){
                      if (marker_flag == INFLOW){
                        TurboNormal[0] = coord[0]/sqrt(Normal2);
                        TurboNormal[1] = coord[1]/sqrt(Normal2);
                        TurboNormal[2] = 0.0;
                      }else{
                        TurboNormal[0] = coord[0]/sqrt(Normal2);
                        TurboNormal[1] = coord[1]/sqrt(Normal2);
                        TurboNormal[2] = 0.0;
                      }
                    }
                    else{
                      if (marker_flag == INFLOW){
                        TurboNormal[0] = -1.0;
                        TurboNormal[1] = 0.0;
                        TurboNormal[2] = 0.0;
                      }else{
                        TurboNormal[0] = 1.0;
                        TurboNormal[1] = 0.0;
                        TurboNormal[2] = 0.0;
                      }
                    }

                    break;
                  case CENTRIPETAL_AXIAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;

                  case AXIAL_CENTRIFUGAL:
                    Normal2 = 0.0;
                    for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                    if (marker_flag == INFLOW){
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }else{
                      TurboNormal[0] = coord[0]/sqrt(Normal2);
                      TurboNormal[1] = coord[1]/sqrt(Normal2);
                      TurboNormal[2] = 0.0;
                    }
                    break;

                  }
                  turbovertex[iMarker][iSpan][iInternalVertex]->SetTurboNormal(TurboNormal);
                  iInternalVertex++;
                }


                for(jSpanVertex = 0; jSpanVertex<nVertexSpanHalo[iSpan]; jSpanVertex++){
                  coord = node[disordered[iSpan][jSpanVertex]]->GetCoord();
                  if(dist >= (coord[1] - target) && !checkAssign[iSpan][jSpanVertex] && (coord[1] - target) >= 0.0){
                    dist= coord[1] - target;
                    kSpanVertex =jSpanVertex;
                  }
                }
              }
            }

            for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

              delete [] ordered[iSpan];
              delete [] disordered[iSpan];
              delete [] oldVertex3D[iSpan];
              delete [] checkAssign[iSpan];
              delete [] area[iSpan];
              delete [] angPitch[iSpan];
              delete [] deltaAngPitch[iSpan];

              for(iVertex=0; iVertex < nVertexSpanHalo[iSpan]; iVertex++){
                delete [] unitnormal[iSpan][iVertex];
              }
              delete [] unitnormal[iSpan];
            }
          }
        }
      }
    }

  /*--- to be set for all the processor to initialize an appropriate number of frequency for the NR BC ---*/
  nVertMax = 0;

  /*--- compute global max and min pitch per span ---*/
  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    nVert    = 0;

#ifdef HAVE_MPI
    MyMin     = minAngPitch[iSpan];      minAngPitch[iSpan]    = 10.0E+6;
    MyIntMin  = minIntAngPitch[iSpan];   minIntAngPitch[iSpan] = 10.0E+6;
    MyMax     = maxAngPitch[iSpan];      maxAngPitch[iSpan]    = -10.0E+6;

    SU2_MPI::Allreduce(&MyMin, &minAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyIntMin, &minIntAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyMax, &maxAngPitch[iSpan], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif


    /*--- compute the relative angular pitch with respect to the minimum value ---*/

    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nVert = nVertexSpan[iMarker][iSpan];
            MinAngularCoord[iMarker][iSpan]    = minAngPitch[iSpan];
            MaxAngularCoord[iMarker][iSpan]    = maxAngPitch[iSpan];
            MinRelAngularCoord[iMarker][iSpan] = minIntAngPitch[iSpan] - minAngPitch[iSpan];
            for(iSpanVertex = 0; iSpanVertex< nVertexSpan[iMarker][iSpan]; iSpanVertex++){
             turbovertex[iMarker][iSpan][iSpanVertex]->SetRelAngularCoord(MinAngularCoord[iMarker][iSpan]);
            }
          }
        }
      }
    }


#ifdef HAVE_MPI
    My_nVert = nVert;nVert = 0;
    SU2_MPI::Allreduce(&My_nVert, &nVert, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif

    /*--- to be set for all the processor to initialize an appropriate number of frequency for the NR BC ---*/
    if(nVert > nVertMax){
      SetnVertexSpanMax(marker_flag,nVert);
    }
    /*--- for all the processor should be known the amount of total turbovertex per span  ---*/
    nTotVertex_gb[iSpan]= (int)nVert;

    for (iMarker = 0; iMarker < nMarker; iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            nTotVertexSpan[iMarker][iSpan]= nVert;
            nTotVertexSpan[iMarker][nSpanWiseSections[marker_flag-1]]+= nVert;
          }
        }
      }
    }
  }


  /*--- Printing Tec file to check the global ordering of the turbovertex pitch-wise ---*/
  /*--- Send all the info to the MASTERNODE ---*/

  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    x_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    y_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    z_loc[iSpan]             = new su2double[nTotVertex_gb[iSpan]];
    angCoord_loc[iSpan]      = new su2double[nTotVertex_gb[iSpan]];
    deltaAngCoord_loc[iSpan] = new su2double[nTotVertex_gb[iSpan]];
    rank_loc[iSpan]          = new int[nTotVertex_gb[iSpan]];
    for(iSpanVertex = 0; iSpanVertex<nTotVertex_gb[iSpan]; iSpanVertex++){
      x_loc[iSpan][iSpanVertex]             = -1.0;
      y_loc[iSpan][iSpanVertex]             = -1.0;
      z_loc[iSpan][iSpanVertex]             = -1.0;
      angCoord_loc[iSpan][iSpanVertex]      = -1.0;
      deltaAngCoord_loc[iSpan][iSpanVertex] = -1.0;
      rank_loc[iSpan][iSpanVertex]          = -1;
    }
  }

  for (iMarker = 0; iMarker < nMarker; iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
            for(iSpanVertex = 0; iSpanVertex<nVertexSpan[iMarker][iSpan]; iSpanVertex++){
              iPoint = turbovertex[iMarker][iSpan][iSpanVertex]->GetNode();
              coord  = node[iPoint]->GetCoord();
              x_loc[iSpan][iSpanVertex]   = coord[0];
              y_loc[iSpan][iSpanVertex]   = coord[1];
              if (nDim == 3){
                z_loc[iSpan][iSpanVertex] = coord[2];
              }
              else{
                z_loc[iSpan][iSpanVertex] = 0.0;
              }
              angCoord_loc[iSpan][iSpanVertex]      = turbovertex[iMarker][iSpan][iSpanVertex]->GetRelAngularCoord();
              deltaAngCoord_loc[iSpan][iSpanVertex] = turbovertex[iMarker][iSpan][iSpanVertex]->GetDeltaAngularCoord();
            }
          }
        }
      }
    }
  }

#ifdef HAVE_MPI

  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    if (rank == MASTER_NODE){
      x_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      y_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      z_gb                = new su2double[nTotVertex_gb[iSpan]*size];
      angCoord_gb         = new su2double[nTotVertex_gb[iSpan]*size];
      deltaAngCoord_gb    = new su2double[nTotVertex_gb[iSpan]*size];
      checkAssign_gb      = new bool[nTotVertex_gb[iSpan]*size];

     for(iSize= 0; iSize < size; iSize++){
       for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
         checkAssign_gb[iSize*nTotVertex_gb[iSpan] + iSpanVertex] = false;
       }
     }
    }
    SU2_MPI::Gather(y_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, y_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(x_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, x_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(z_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, z_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(angCoord_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, angCoord_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Gather(deltaAngCoord_loc[iSpan], nTotVertex_gb[iSpan] , MPI_DOUBLE, deltaAngCoord_gb, nTotVertex_gb[iSpan], MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE){
      for(iSpanVertex = 0; iSpanVertex<nTotVertex_gb[iSpan]; iSpanVertex++){
        x_loc[iSpan][iSpanVertex]             = -1.0;
        y_loc[iSpan][iSpanVertex]             = -1.0;
        z_loc[iSpan][iSpanVertex]             = -1.0;
        angCoord_loc[iSpan][iSpanVertex]      = -1.0;
        deltaAngCoord_loc[iSpan][iSpanVertex] = -1.0;
      }



      min = 10.0E+06;
      for(iSize= 0; iSize < size; iSize++){
        if (angCoord_gb[iSize*nTotVertex_gb[iSpan]] < min && angCoord_gb[iSize*nTotVertex_gb[iSpan]] >= 0.0){
          kSize = iSize;
          min = angCoord_gb[iSize*nTotVertex_gb[iSpan]];
        }
      }

      kSpanVertex = 0;
      for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
        x_loc[iSpan][iSpanVertex]              = x_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        y_loc[iSpan][iSpanVertex]              = y_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        z_loc[iSpan][iSpanVertex]              = z_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        angCoord_loc[iSpan][iSpanVertex]       = angCoord_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        deltaAngCoord_loc[iSpan][iSpanVertex]  = deltaAngCoord_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex];
        rank_loc[iSpan][iSpanVertex]           = kSize;
        target = angCoord_loc[iSpan][iSpanVertex];
        checkAssign_gb[kSize*nTotVertex_gb[iSpan] + kSpanVertex] = true;
        min = 10.0E+06;
        for(jSize= 0; jSize < size; jSize++){
          for(jSpanVertex = 0; jSpanVertex < nTotVertex_gb[iSpan]; jSpanVertex++){
            if ((angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex] < min) &&
                (angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex] >= target) &&
                !checkAssign_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex]) {
              kSize = jSize;
              kSpanVertex = jSpanVertex;
              min = angCoord_gb[jSize*nTotVertex_gb[iSpan] + jSpanVertex];
            }
          }
        }
      }

      delete [] x_gb;	delete [] y_gb; delete [] z_gb;	 delete [] angCoord_gb; delete [] deltaAngCoord_gb; delete[] checkAssign_gb;

    }
  }

#endif

  if (rank == MASTER_NODE){
    if (marker_flag == INFLOW && val_iZone ==0){
      std::string sPath = "TURBOMACHINERY";
      int nError = 0;
#if defined(_WIN32)
#ifdef __MINGW32__
      nError = mkdir(sPath.c_str());  // MINGW on Windows
#else
      nError = _mkdir(sPath.c_str()); // can be used on Windows
#endif
#else
      mode_t nMode = 0733; // UNIX style permissions
      nError = mkdir(sPath.c_str(),nMode); // can be used on non-Windows
#endif
      if (nError != 0) {
        cout << "TURBOMACHINERY folder creation failed." <<endl;
      }
    }
    if (marker_flag == INFLOW){
      multizone_filename = "TURBOMACHINERY/spanwise_division_inflow.dat";
    }
    else{
      multizone_filename = "TURBOMACHINERY/spanwise_division_outflow.dat";
    }
    char buffer[50];

    if (GetnZone() > 1){
      unsigned short lastindex = multizone_filename.find_last_of(".");
      multizone_filename = multizone_filename.substr(0, lastindex);
      SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(val_iZone));
      multizone_filename.append(string(buffer));
    }

    // File to print the vector x_loc, y_loc, z_loc, globIdx_loc to check vertex ordering
    ofstream myfile;
    myfile.open (multizone_filename.data(), ios::out | ios::trunc);
    myfile.setf(ios::uppercase | ios::scientific);
    myfile.precision(8);

    myfile << "TITLE = \"Global index visualization file\"" << endl;
    myfile << "VARIABLES =" << endl;
    myfile.width(10); myfile << "\"iSpan\"";
    myfile.width(20); myfile << "\"x_coord\"" ;
    myfile.width(20); myfile << "\"y_coord\"" ;
    myfile.width(20); myfile << "\"z_coord\"" ;
    myfile.width(20); myfile << "\"radius\"" ;
    myfile.width(20); myfile << "\"Relative Angular Coord \"" ;
    myfile.width(20); myfile << "\"Delta Angular Coord \"" ;
    myfile.width(20); myfile << "\"processor\"" <<endl;
    for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
      for(iSpanVertex = 0; iSpanVertex < nTotVertex_gb[iSpan]; iSpanVertex++){
        radius = sqrt(x_loc[iSpan][iSpanVertex]*x_loc[iSpan][iSpanVertex] + y_loc[iSpan][iSpanVertex]*y_loc[iSpan][iSpanVertex]);
        myfile.width(10); myfile << iSpan;
        myfile.width(20); myfile << x_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << y_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << z_loc[iSpan][iSpanVertex];
        myfile.width(20); myfile << radius;
        if (nDim ==2 && config->GetKind_TurboMachinery(val_iZone)){
          myfile.width(20); myfile << angCoord_loc[iSpan][iSpanVertex];
          myfile.width(20); myfile << deltaAngCoord_loc[iSpan][iSpanVertex];
        }
        else{
          myfile.width(20); myfile << angCoord_loc[iSpan][iSpanVertex]*180.0/PI_NUMBER;
          myfile.width(20); myfile << deltaAngCoord_loc[iSpan][iSpanVertex]*180.0/PI_NUMBER;
        }
        myfile.width(20); myfile << rank_loc[iSpan][iSpanVertex]<<endl;
      }
      myfile << endl;
    }
  }


  for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
    delete [] x_loc[iSpan];
    delete [] y_loc[iSpan];
    delete [] z_loc[iSpan];
    delete [] angCoord_loc[iSpan];
    delete [] deltaAngCoord_loc[iSpan];
    delete [] rank_loc[iSpan];

  }


  delete [] area;
  delete [] ordered;
  delete [] disordered;
  delete [] oldVertex3D;
  delete [] checkAssign;
  delete [] TurboNormal;
  delete [] unitnormal;
  delete [] NormalArea;
  delete [] x_loc;
  delete [] y_loc;
  delete [] z_loc;
  delete [] angCoord_loc;
  delete [] nTotVertex_gb;
  delete [] nVertexSpanHalo;
  delete [] angPitch;
  delete [] deltaAngPitch;
  delete [] deltaAngCoord_loc;
  delete [] rank_loc;
  delete [] minAngPitch;
  delete [] maxAngPitch;
  delete [] minIntAngPitch;

}


void CPhysicalGeometry::UpdateTurboVertex(CConfig *config, unsigned short val_iZone, unsigned short marker_flag) {
  unsigned short iMarker, iMarkerTP, iSpan, iDim;
  long iSpanVertex, iPoint;
  su2double *coord, *TurboNormal, Normal2;

  /*--- Initialize auxiliary pointers ---*/
  TurboNormal      = new su2double[3];

  for (iMarker = 0; iMarker < nMarker; iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          for(iSpan = 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){
            for(iSpanVertex = 0; iSpanVertex<nVertexSpan[iMarker][iSpan]; iSpanVertex++){
              iPoint = turbovertex[iMarker][iSpan][iSpanVertex]->GetNode();
              coord  = node[iPoint]->GetCoord();
              /*--- compute appropriate turbo normal ---*/
              switch (config->GetKind_TurboMachinery(val_iZone)){
              case CENTRIFUGAL:
                Normal2 = 0.0;
                for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                if (marker_flag == INFLOW){
                  TurboNormal[0] = -coord[0]/sqrt(Normal2);
                  TurboNormal[1] = -coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }else{
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }
                break;
              case CENTRIPETAL:
                Normal2 = 0.0;
                for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                if (marker_flag == OUTFLOW){
                  TurboNormal[0] = -coord[0]/sqrt(Normal2);
                  TurboNormal[1] = -coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }else{
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }
                break;
              case AXIAL:
                Normal2 = 0.0;
                for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                if(nDim == 3){
                  if (marker_flag == INFLOW){
                    TurboNormal[0] = coord[0]/sqrt(Normal2);
                    TurboNormal[1] = coord[1]/sqrt(Normal2);
                    TurboNormal[2] = 0.0;
                  }else{
                    TurboNormal[0] = coord[0]/sqrt(Normal2);
                    TurboNormal[1] = coord[1]/sqrt(Normal2);
                    TurboNormal[2] = 0.0;
                  }
                }
                else{
                  if (marker_flag == INFLOW){
                    TurboNormal[0] = -1.0;
                    TurboNormal[1] = 0.0;
                    TurboNormal[2] = 0.0;
                  }else{
                    TurboNormal[0] = 1.0;
                    TurboNormal[1] = 0.0;
                    TurboNormal[2] = 0.0;
                  }
                }

                break;
              case CENTRIPETAL_AXIAL:
                Normal2 = 0.0;
                for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                if (marker_flag == INFLOW){
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }else{
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }
                break;

              case AXIAL_CENTRIFUGAL:
                Normal2 = 0.0;
                for(iDim = 0; iDim < 2; iDim++) Normal2 +=coord[iDim]*coord[iDim];
                if (marker_flag == INFLOW){
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }else{
                  TurboNormal[0] = coord[0]/sqrt(Normal2);
                  TurboNormal[1] = coord[1]/sqrt(Normal2);
                  TurboNormal[2] = 0.0;
                }
                break;
              }

              /*--- store the new turbo normal ---*/
              turbovertex[iMarker][iSpan][iSpanVertex]->SetTurboNormal(TurboNormal);
            }
          }
        }
      }
    }
  }

  delete [] TurboNormal;
}

void CPhysicalGeometry::SetAvgTurboValue(CConfig *config, unsigned short val_iZone, unsigned short marker_flag, bool allocate) {

  unsigned short iMarker, iMarkerTP, iSpan, iDim;
  unsigned long iPoint;
  su2double *TurboNormal,*coord, *Normal, turboNormal2, Normal2, *gridVel, TotalArea, TotalRadius, radius;
  su2double *TotalTurboNormal,*TotalNormal, *TotalGridVel, Area;
  long iVertex;
  /*-- Variables declaration and allocation ---*/
  TotalTurboNormal = new su2double[nDim];
  TotalNormal      = new su2double[nDim];
  TurboNormal      = new su2double[nDim];
  TotalGridVel     = new su2double[nDim];
  Normal           = new su2double[nDim];

  bool grid_movement        = config->GetGrid_Movement();
#ifdef HAVE_MPI
  su2double MyTotalArea, MyTotalRadius, *MyTotalTurboNormal= NULL, *MyTotalNormal= NULL, *MyTotalGridVel= NULL;
#endif

  /*--- Intialization of the vector for the interested boundary ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          if(allocate){
            AverageTurboNormal[iMarker]               = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageNormal[iMarker]                    = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageGridVel[iMarker]                   = new su2double *[nSpanWiseSections[marker_flag-1] + 1];
            AverageTangGridVel[iMarker]               = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            SpanArea[iMarker]                         = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            TurboRadius[iMarker]                      = new su2double [nSpanWiseSections[marker_flag-1] + 1];
            for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1] + 1; iSpan++){
              AverageTurboNormal[iMarker][iSpan]      = new su2double [nDim];
              AverageNormal[iMarker][iSpan]           = new su2double [nDim];
              AverageGridVel[iMarker][iSpan]          = new su2double [nDim];
            }
          }
          for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1] + 1; iSpan++){
            AverageTangGridVel[iMarker][iSpan]          = 0.0;
            SpanArea[iMarker][iSpan]                    = 0.0;
            TurboRadius[iMarker][iSpan]                 = 0.0;
            for(iDim=0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][iSpan][iDim]  = 0.0;
              AverageNormal[iMarker][iSpan][iDim]       = 0.0;
              AverageGridVel[iMarker][iSpan][iDim]      = 0.0;
            }
          }
        }
      }
    }
  }



  /*--- start computing the average quantities span wise --- */
  for (iSpan= 0; iSpan < nSpanWiseSections[marker_flag-1]; iSpan++){

    /*--- Forces initialization for contenitors to zero ---*/
    for (iDim=0; iDim<nDim; iDim++) {
      TotalTurboNormal[iDim]	=0.0;
      TotalNormal[iDim]         =0.0;
      TotalGridVel[iDim]        =0.0;
    }
    TotalArea = 0.0;
    TotalRadius = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
            for(iVertex = 0; iVertex < nVertexSpan[iMarker][iSpan]; iVertex++){
              iPoint = turbovertex[iMarker][iSpan][iVertex]->GetNode();
              turbovertex[iMarker][iSpan][iVertex]->GetTurboNormal(TurboNormal);
              turbovertex[iMarker][iSpan][iVertex]->GetNormal(Normal);
              coord  = node[iPoint]->GetCoord();

              if (nDim == 3){
                radius = sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
              }
              else{
                radius = 0.0;
              }
              Area = turbovertex[iMarker][iSpan][iVertex]->GetArea();
              TotalArea   += Area;
              TotalRadius += radius;
              for (iDim = 0; iDim < nDim; iDim++) {
                TotalTurboNormal[iDim]  +=TurboNormal[iDim];
                TotalNormal[iDim]       +=Normal[iDim];
              }
              if (grid_movement){
                gridVel = node[iPoint]->GetGridVel();
                for (iDim = 0; iDim < nDim; iDim++) TotalGridVel[iDim] +=gridVel[iDim];
              }
            }
          }
        }
      }
    }

#ifdef HAVE_MPI

    MyTotalArea            = TotalArea;                 TotalArea            = 0;
    MyTotalRadius          = TotalRadius;               TotalRadius          = 0;
    SU2_MPI::Allreduce(&MyTotalArea, &TotalArea, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&MyTotalRadius, &TotalRadius, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    MyTotalTurboNormal     = new su2double[nDim];
    MyTotalNormal          = new su2double[nDim];
    MyTotalGridVel         = new su2double[nDim];

    for (iDim = 0; iDim < nDim; iDim++) {
      MyTotalTurboNormal[iDim]	= TotalTurboNormal[iDim];
      TotalTurboNormal[iDim]    = 0.0;
      MyTotalNormal[iDim]       = TotalNormal[iDim];
      TotalNormal[iDim]         = 0.0;
      MyTotalGridVel[iDim]      = TotalGridVel[iDim];
      TotalGridVel[iDim]        = 0.0;
    }

    SU2_MPI::Allreduce(MyTotalTurboNormal, TotalTurboNormal, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MyTotalNormal, TotalNormal, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(MyTotalGridVel, TotalGridVel, nDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    delete [] MyTotalTurboNormal;delete [] MyTotalNormal; delete [] MyTotalGridVel;

#endif

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){


            SpanArea[iMarker][iSpan]           = TotalArea;
            TurboRadius[iMarker][iSpan]        = TotalRadius/nTotVertexSpan[iMarker][iSpan];

            turboNormal2    = 0.0;
            Normal2         = 0.0;
            for (iDim = 0; iDim < nDim; iDim++){
              turboNormal2 += TotalTurboNormal[iDim]*TotalTurboNormal[iDim];
              Normal2      += TotalNormal[iDim]*TotalNormal[iDim];
            }
            for (iDim = 0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][iSpan][iDim] = TotalTurboNormal[iDim]/sqrt(turboNormal2);
              AverageNormal[iMarker][iSpan][iDim]      = TotalNormal[iDim]/sqrt(Normal2);
            }
            if (grid_movement){
              for (iDim = 0; iDim < nDim; iDim++){
                AverageGridVel[iMarker][iSpan][iDim]   =TotalGridVel[iDim]/nTotVertexSpan[iMarker][iSpan];
              }
              switch (config->GetKind_TurboMachinery(val_iZone)){
              case CENTRIFUGAL:case CENTRIPETAL:
                if (marker_flag == INFLOW ){
                  AverageTangGridVel[iMarker][iSpan]= -(AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                break;
              case AXIAL:
                if (marker_flag == INFLOW && nDim == 2){
                  AverageTangGridVel[iMarker][iSpan]= -AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1] + AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                  break;
              case CENTRIPETAL_AXIAL:
                if (marker_flag == OUTFLOW){
                  AverageTangGridVel[iMarker][iSpan]= (AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                else{
                  AverageTangGridVel[iMarker][iSpan]= -(AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0]);
                }
                break;
              case AXIAL_CENTRIFUGAL:
                if (marker_flag == INFLOW)
                {
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }else
                {
                  AverageTangGridVel[iMarker][iSpan]= AverageTurboNormal[iMarker][iSpan][0]*AverageGridVel[iMarker][iSpan][1]-AverageTurboNormal[iMarker][iSpan][1]*AverageGridVel[iMarker][iSpan][0];
                }
                break;

              default:
                  SU2_MPI::Error("Tang grid velocity NOT IMPLEMENTED YET for this configuration", CURRENT_FUNCTION);
                break;
              }
            }

            /*--- Compute the 1D average values ---*/
            AverageTangGridVel[iMarker][nSpanWiseSections[marker_flag-1]]             += AverageTangGridVel[iMarker][iSpan]/nSpanWiseSections[marker_flag-1];
            SpanArea[iMarker][nSpanWiseSections[marker_flag-1]]                       += SpanArea[iMarker][iSpan];
            for(iDim=0; iDim < nDim; iDim++){
              AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]     += AverageTurboNormal[iMarker][iSpan][iDim];
              AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]          += AverageNormal[iMarker][iSpan][iDim];
              AverageGridVel[iMarker][nSpanWiseSections[marker_flag-1]][iDim]         += AverageGridVel[iMarker][iSpan][iDim]/nSpanWiseSections[marker_flag-1];

            }
          }
        }
      }
    }
  }

  /*--- Normalize 1D normals---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    for (iMarkerTP=1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
      if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
        if (config->GetMarker_All_TurbomachineryFlag(iMarker) == marker_flag){
          turboNormal2 = 0.0;
          Normal2 		= 0.0;

          for (iDim = 0; iDim < nDim; iDim++){
            turboNormal2 += AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]*AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim];
            Normal2      += AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim]*AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim];
          }
          for (iDim = 0; iDim < nDim; iDim++){
            AverageTurboNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim] /=sqrt(turboNormal2);
            AverageNormal[iMarker][nSpanWiseSections[marker_flag-1]][iDim] /=sqrt(Normal2);
          }
        }
      }
    }
  }


  delete [] TotalTurboNormal;
  delete [] TotalNormal;
  delete [] TotalGridVel;
  delete [] TurboNormal;
  delete [] Normal;

}


void CPhysicalGeometry::GatherInOutAverageValues(CConfig *config, bool allocate){

  unsigned short iMarker, iMarkerTP;
  unsigned short iSpan, iDim;
  int markerTP;
  su2double nBlades;
  unsigned short nSpanWiseSections = config->GetnSpanWiseSections();

  su2double tangGridVelIn, tangGridVelOut;
  su2double areaIn, areaOut, pitchIn, Pitch;
  su2double radiusIn, radiusOut, *turboNormal;

  turboNormal = new su2double[nDim];
  Pitch = 0.0;

  if(allocate){
    for (iMarkerTP=0; iMarkerTP < config->GetnMarker_TurboPerformance(); iMarkerTP++){
      SpanAreaIn[iMarkerTP]       = new su2double[config->GetnSpanMaxAllZones() +1];
      TangGridVelIn[iMarkerTP]    = new su2double[config->GetnSpanMaxAllZones() +1];
      TurboRadiusIn[iMarkerTP]    = new su2double[config->GetnSpanMaxAllZones() +1];
      SpanAreaOut[iMarkerTP]      = new su2double[config->GetnSpanMaxAllZones() +1];
      TangGridVelOut[iMarkerTP]   = new su2double[config->GetnSpanMaxAllZones() +1];
      TurboRadiusOut[iMarkerTP]   = new su2double[config->GetnSpanMaxAllZones() +1];

      for (iSpan= 0; iSpan < config->GetnSpanMaxAllZones() + 1 ; iSpan++){
        SpanAreaIn[iMarkerTP][iSpan]       = 0.0;
        TangGridVelIn[iMarkerTP][iSpan]    = 0.0;
        TurboRadiusIn[iMarkerTP][iSpan]    = 0.0;
        SpanAreaOut[iMarkerTP][iSpan]      = 0.0;
        TangGridVelOut[iMarkerTP][iSpan]   = 0.0;
        TurboRadiusOut[iMarkerTP][iSpan]   = 0.0;
      }
    }
  }



  for (iSpan= 0; iSpan < nSpanWiseSections + 1 ; iSpan++){
#ifdef HAVE_MPI
    unsigned short i, n1, n2, n1t,n2t;
    su2double *TurbGeoIn= NULL,*TurbGeoOut= NULL;
    su2double *TotTurbGeoIn = NULL,*TotTurbGeoOut = NULL;
    int *TotMarkerTP;

    n1          = 6;
    n2          = 3;
    n1t         = n1*size;
    n2t         = n2*size;
    TurbGeoIn  = new su2double[n1];
    TurbGeoOut = new su2double[n2];

    for (i=0;i<n1;i++)
      TurbGeoIn[i]    = -1.0;
    for (i=0;i<n2;i++)
      TurbGeoOut[i]   = -1.0;
#endif
    pitchIn           =  0.0;
    areaIn            = -1.0;
    tangGridVelIn     = -1.0;
    radiusIn          = -1.0;
    for(iDim = 0; iDim < nDim; iDim++){
      turboNormal[iDim] = -1.0;
    }

    areaOut           = -1.0;
    tangGridVelOut    = -1.0;
    radiusOut         = -1.0;

    markerTP          = -1;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
      for (iMarkerTP = 1; iMarkerTP < config->GetnMarker_Turbomachinery()+1; iMarkerTP++){
        if (config->GetMarker_All_Turbomachinery(iMarker) == iMarkerTP){
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == INFLOW){
            markerTP          = iMarkerTP;
            if (iSpan < nSpanWiseSections){
              pitchIn         = MaxAngularCoord[iMarker][iSpan] - MinAngularCoord[iMarker][iSpan];
            }
            areaIn            = SpanArea[iMarker][iSpan];
            tangGridVelIn     = AverageTangGridVel[iMarker][iSpan];
            radiusIn          = TurboRadius[iMarker][iSpan];
            for(iDim = 0; iDim < nDim; iDim++) turboNormal[iDim] = AverageTurboNormal[iMarker][iSpan][iDim];


#ifdef HAVE_MPI
            TurbGeoIn[0]  = areaIn;
            TurbGeoIn[1]  = tangGridVelIn;
            TurbGeoIn[2]  = radiusIn;
            TurbGeoIn[3]  = turboNormal[0];
            TurbGeoIn[4]  = turboNormal[1];
            TurbGeoIn[5]  = pitchIn;
#endif
          }

          /*--- retrieve outlet information ---*/
          if (config->GetMarker_All_TurbomachineryFlag(iMarker) == OUTFLOW){
            if (iSpan < nSpanWiseSections){
              pitchIn       = MaxAngularCoord[iMarker][iSpan] - MinAngularCoord[iMarker][iSpan];
            }
            areaOut         = SpanArea[iMarker][iSpan];
            tangGridVelOut  = AverageTangGridVel[iMarker][iSpan];
            radiusOut       = TurboRadius[iMarker][iSpan];
#ifdef HAVE_MPI
            TurbGeoOut[0]  = areaOut;
            TurbGeoOut[1]  = tangGridVelOut;
            TurbGeoOut[2]  = radiusOut;
#endif
          }
        }
      }
    }

#ifdef HAVE_MPI
    TotTurbGeoIn       = new su2double[n1t];
    TotTurbGeoOut      = new su2double[n2t];
    for (i=0;i<n1t;i++)
      TotTurbGeoIn[i]  = -1.0;
    for (i=0;i<n2t;i++)
      TotTurbGeoOut[i] = -1.0;
    TotMarkerTP = new int[size];
    for(i=0; i<size; i++){
      TotMarkerTP[i]    = -1;
    }

    SU2_MPI::Allgather(TurbGeoIn, n1, MPI_DOUBLE, TotTurbGeoIn, n1, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(TurbGeoOut, n2, MPI_DOUBLE,TotTurbGeoOut, n2, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(&markerTP, 1, MPI_INT,TotMarkerTP, 1, MPI_INT, MPI_COMM_WORLD);

    delete [] TurbGeoIn, delete [] TurbGeoOut;


    for (i=0;i<size;i++){
      if(TotTurbGeoIn[n1*i] > 0.0){
        areaIn              = 0.0;
        areaIn              = TotTurbGeoIn[n1*i];
        tangGridVelIn       = 0.0;
        tangGridVelIn       = TotTurbGeoIn[n1*i+1];
        radiusIn            = 0.0;
        radiusIn            = TotTurbGeoIn[n1*i+2];
        turboNormal[0]      = 0.0;
        turboNormal[0]      = TotTurbGeoIn[n1*i+3];
        turboNormal[1]      = 0.0;
        turboNormal[1]      = TotTurbGeoIn[n1*i+4];
        pitchIn             = 0.0;
        pitchIn             = TotTurbGeoIn[n1*i+5];

        markerTP            = -1;
        markerTP            = TotMarkerTP[i];
      }

      if(TotTurbGeoOut[n2*i] > 0.0){
        areaOut             = 0.0;
        areaOut             = TotTurbGeoOut[n2*i];
        tangGridVelOut      = 0.0;
        tangGridVelOut      = TotTurbGeoOut[n2*i+1];
        radiusOut           = 0.0;
        radiusOut           = TotTurbGeoOut[n2*i+2];
      }
    }

    delete [] TotTurbGeoIn, delete [] TotTurbGeoOut; delete [] TotMarkerTP;


#endif

    Pitch +=pitchIn/nSpanWiseSections;

    if (iSpan == nSpanWiseSections) {
      config->SetFreeStreamTurboNormal(turboNormal);
      if (config->GetKind_TurboMachinery(config->GetiZone()) == AXIAL && nDim == 2){
        nBlades = 1/Pitch;
      }
      else{
        nBlades = 2*PI_NUMBER/Pitch;
      }
      config->SetnBlades(config->GetiZone(), nBlades);
    }

    if (rank == MASTER_NODE){
      /*----Quantities needed for computing the turbomachinery performance -----*/
      SpanAreaIn[markerTP -1][iSpan]       = areaIn;
      TangGridVelIn[markerTP -1][iSpan]    = tangGridVelIn;
      TurboRadiusIn[markerTP -1][iSpan]    = radiusIn;

      SpanAreaOut[markerTP -1][iSpan]      = areaOut;
      TangGridVelOut[markerTP -1][iSpan]   = tangGridVelOut;
      TurboRadiusOut[markerTP -1][iSpan]   = radiusOut;
    }
  }
  delete [] turboNormal;

}


void CPhysicalGeometry::SetCoord_CG(void) {
  unsigned short iMarker, iNode;
  unsigned long elem_poin, edge_poin, iElem, iEdge;

  /*--- Buffer of pointers to node coordinates ---*/
  array<const su2double*, N_POINTS_MAXIMUM> Coord;

  /*--- Compute the center of gravity for elements ---*/

  for (iElem = 0; iElem<nElem; iElem++) {
    assert(elem[iElem]->GetnNodes() <= N_POINTS_MAXIMUM && "Insufficient N_POINTS_MAXIMUM");

    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < elem[iElem]->GetnNodes(); iNode++) {
      elem_poin = elem[iElem]->GetNode(iNode);
      Coord[iNode] = node[elem_poin]->GetCoord();
    }

    /*--- Compute the element CG coordinates ---*/
    elem[iElem]->SetCoord_CG(Coord.data());
  }

  /*--- Center of gravity for face elements ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {

      /*--- Store the coordinates for all the element nodes ---*/
      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        elem_poin = bound[iMarker][iElem]->GetNode(iNode);
        Coord[iNode] = node[elem_poin]->GetCoord();
      }

      /*--- Compute the element CG coordinates ---*/
      bound[iMarker][iElem]->SetCoord_CG(Coord.data());
    }
  }

  /*--- Center of gravity for edges ---*/

  for (iEdge = 0; iEdge < nEdge; iEdge++) {

    /*--- Store the coordinates for all the element nodes ---*/
    for (iNode = 0; iNode < edges->GetnNodes(); iNode++) {
      edge_poin=edges->GetNode(iEdge,iNode);
      Coord[iNode] = node[edge_poin]->GetCoord();
    }

    /*--- Compute the edge CG coordinates ---*/
    edges->SetCoord_CG(iEdge, Coord.data());
  }

}

void CPhysicalGeometry::SetBoundControlVolume(CConfig *config, unsigned short action) {
  unsigned short Neighbor_Node, iMarker, iNode, iNeighbor_Nodes, iDim;
  unsigned long Neighbor_Point, iVertex, iPoint, iElem;
  long iEdge;
  su2double Area, *NormalFace = NULL;

  /*--- Update values of faces of the edge ---*/

  if (action != ALLOCATE)
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
        vertex[iMarker][iVertex]->SetZeroValues();

  su2double *Coord_Edge_CG = new su2double [nDim];
  su2double *Coord_Elem_CG = new su2double [nDim];
  su2double *Coord_Vertex = new su2double [nDim];

  /*--- Loop over all the markers ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++)

  /*--- Loop over all the boundary elements ---*/

    for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++)

    /*--- Loop over all the nodes of the boundary ---*/

      for (iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
        iPoint = bound[iMarker][iElem]->GetNode(iNode);
        iVertex = node[iPoint]->GetVertex(iMarker);

        /*--- Loop over the neighbor nodes, there is a face for each one ---*/

        for (iNeighbor_Nodes = 0; iNeighbor_Nodes < bound[iMarker][iElem]->GetnNeighbor_Nodes(iNode); iNeighbor_Nodes++) {
          Neighbor_Node = bound[iMarker][iElem]->GetNeighbor_Nodes(iNode, iNeighbor_Nodes);
          Neighbor_Point = bound[iMarker][iElem]->GetNode(Neighbor_Node);

          /*--- Shared edge by the Neighbor Point and the point ---*/

          iEdge = FindEdge(iPoint, Neighbor_Point);
          for (iDim = 0; iDim < nDim; iDim++) {
            Coord_Edge_CG[iDim] = edges->GetCG(iEdge,iDim);
            Coord_Elem_CG[iDim] = bound[iMarker][iElem]->GetCG(iDim);
            Coord_Vertex[iDim] = node[iPoint]->GetCoord(iDim);
          }
          switch (nDim) {
            case 2:

              /*--- Store the 2D face ---*/

              if (iNode == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Vertex);
              if (iNode == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Vertex, Coord_Elem_CG);
              break;
            case 3:

              /*--- Store the 3D face ---*/

              if (iNeighbor_Nodes == 0) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Elem_CG, Coord_Edge_CG, Coord_Vertex);
              if (iNeighbor_Nodes == 1) vertex[iMarker][iVertex]->SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG, Coord_Vertex);
              break;
          }
        }
      }

  delete[] Coord_Edge_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_Vertex;

  /*--- Check if there is a normal with null area ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker ++)
    for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
      NormalFace = vertex[iMarker][iVertex]->GetNormal();
      Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
      Area = sqrt(Area);
      if (Area == 0.0) for (iDim = 0; iDim < nDim; iDim++) NormalFace[iDim] = EPS*EPS;
    }

}

void CPhysicalGeometry::SetMaxLength(CConfig* config) {

  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++){
    const unsigned short nNeigh = node[iPoint]->GetnPoint();
    const su2double* Coord_i = node[iPoint]->GetCoord();

    /*--- If using AD, computing the maximum grid length can generate
     * a lot of unnecessary overhead since we would store all computations
     * of each grid length, even though we only need the maximum value.
     * We solve that by finding the neighbor that's furthest away
     * (corresponding to the maximum distance) using passive calculations,
     * then set the max using the AD datatype. ---*/

    passivedouble passive_max_delta=0;
    unsigned short max_neighbor = 0;
    for (unsigned short iNeigh = 0; iNeigh < nNeigh; iNeigh++) {

      /*-- Calculate the cell-center to cell-center length ---*/

      const unsigned long jPoint  = node[iPoint]->GetPoint(iNeigh);
      const su2double* Coord_j = node[jPoint]->GetCoord();

      passivedouble delta_aux = 0;
      for (unsigned short iDim = 0;iDim < nDim; iDim++){
        delta_aux += pow(SU2_TYPE::GetValue(Coord_j[iDim])-SU2_TYPE::GetValue(Coord_i[iDim]), 2.);
      }

      /*--- Only keep the maximum length ---*/

      if (delta_aux > passive_max_delta) {
        passive_max_delta = delta_aux;
        max_neighbor = iNeigh;
      }
    }

    /*--- Now that we know where the maximum distance is, repeat
     * calculation with the AD-friendly su2double datatype ---*/

    const unsigned long jPoint  = node[iPoint]->GetPoint(max_neighbor);
    const su2double* Coord_j = node[jPoint]->GetCoord();

    su2double max_delta = 0;
    for (unsigned short iDim = 0;iDim < nDim; iDim++) {
      max_delta += pow((Coord_j[iDim]-Coord_i[iDim]), 2.);
    }
    max_delta = sqrt(max_delta);

    node[iPoint]->SetMaxLength(max_delta);
  }

  InitiateComms(this, config, MAX_LENGTH);
  CompleteComms(this, config, MAX_LENGTH);

}

void CPhysicalGeometry::MatchNearField(CConfig *config) {

  su2double epsilon = 1e-1;

  unsigned short nMarker_NearFieldBound = config->GetnMarker_NearFieldBound();

  if (nMarker_NearFieldBound != 0) {

    unsigned short iMarker, iDim, jMarker, pMarker = 0;
    unsigned long iVertex, iPoint, pVertex = 0, pPoint = 0, jVertex, jPoint, iPointGlobal, jPointGlobal, jVertex_, pPointGlobal = 0;
    su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local, maxdist_global;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_NearField = 0, MaxLocalVertex_NearField = 0;
    int nProcessor = size;

    unsigned long *Buffer_Send_nVertex = new unsigned long [1];
    unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];

    if (rank == MASTER_NODE) cout << "Set NearField boundary conditions (if any)." << endl;

    /*--- Compute the number of vertex that have interfase boundary condition
     without including the ghost nodes ---*/

    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_NearField ++;
        }

    Buffer_Send_nVertex[0] = nLocalVertex_NearField;

    /*--- Send NearField vertex information --*/

    SU2_MPI::Allreduce(&nLocalVertex_NearField, &MaxLocalVertex_NearField, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

    su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Send_Point = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_GlobalIndex  = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_Vertex  = new unsigned long [MaxLocalVertex_NearField];
    unsigned long *Buffer_Send_Marker  = new unsigned long [MaxLocalVertex_NearField];

    su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_NearField*nDim];
    unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_Vertex = new unsigned long [nProcessor*MaxLocalVertex_NearField];
    unsigned long *Buffer_Receive_Marker = new unsigned long [nProcessor*MaxLocalVertex_NearField];

    unsigned long nBuffer_Coord = MaxLocalVertex_NearField*nDim;
    unsigned long nBuffer_Point = MaxLocalVertex_NearField;
    unsigned long nBuffer_GlobalIndex = MaxLocalVertex_NearField;
    unsigned long nBuffer_Vertex = MaxLocalVertex_NearField;
    unsigned long nBuffer_Marker = MaxLocalVertex_NearField;

    for (iVertex = 0; iVertex < MaxLocalVertex_NearField; iVertex++) {
      Buffer_Send_Point[iVertex] = 0;
      Buffer_Send_GlobalIndex[iVertex] = 0;
      Buffer_Send_Vertex[iVertex] = 0;
      Buffer_Send_Marker[iVertex] = 0;
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
    }

    /*--- Copy coordinates and point to the auxiliar vector --*/

    nLocalVertex_NearField = 0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY)
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_NearField] = iPoint;
            Buffer_Send_GlobalIndex[nLocalVertex_NearField] = iPointGlobal;
            Buffer_Send_Vertex[nLocalVertex_NearField] = iVertex;
            Buffer_Send_Marker[nLocalVertex_NearField] = iMarker;
            for (iDim = 0; iDim < nDim; iDim++)
              Buffer_Send_Coord[nLocalVertex_NearField*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_NearField++;
          }
        }

    SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, Buffer_Receive_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, Buffer_Receive_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);


    /*--- Compute the closest point to a Near-Field boundary point ---*/

    maxdist_local = 0.0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_KindBC(iMarker) == NEARFIELD_BOUNDARY) {

        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();

          if (node[iPoint]->GetDomain()) {

            /*--- Coordinates of the boundary point ---*/

            Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

            /*--- Loop over all the boundaries to find the pair ---*/
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_NearField+jVertex];
                jPointGlobal = Buffer_Receive_GlobalIndex[iProcessor*MaxLocalVertex_NearField+jVertex];
                jVertex_ = Buffer_Receive_Vertex[iProcessor*MaxLocalVertex_NearField+jVertex];
                jMarker = Buffer_Receive_Marker[iProcessor*MaxLocalVertex_NearField+jVertex];

                if (jPointGlobal != iPointGlobal) {

                  /*--- Compute the distance ---*/

                  dist = 0.0; for (iDim = 0; iDim < nDim; iDim++) {
                    Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_NearField+jVertex)*nDim+iDim];
                    dist += pow(Coord_j[iDim]-Coord_i[iDim],2.0);
                  } dist = sqrt(dist);

                  if (((dist < mindist) && (iProcessor != rank)) ||
                      ((dist < mindist) && (iProcessor == rank) && (jPoint != iPoint))) {
                    mindist = dist; pProcessor = iProcessor; pPoint = jPoint; pPointGlobal = jPointGlobal;
                    pVertex = jVertex_; pMarker = jMarker;
                    if (dist == 0.0) break;
                  }
                }
              }

            /*--- Store the value of the pair ---*/

            maxdist_local = max(maxdist_local, mindist);
            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);

            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPoint << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iPointGlobal, pVertex, pMarker, pProcessor);
              maxdist_local = min(maxdist_local, 0.0);
            }

          }
        }
      }
    }

    SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Point;

    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Point;

    delete[] Buffer_Send_nVertex;
    delete[] Buffer_Receive_nVertex;

    delete [] Buffer_Send_GlobalIndex;
    delete [] Buffer_Send_Vertex;
    delete [] Buffer_Send_Marker;

    delete [] Buffer_Receive_GlobalIndex;
    delete [] Buffer_Receive_Vertex;
    delete [] Buffer_Receive_Marker;

  }

}

void CPhysicalGeometry::MatchActuator_Disk(CConfig *config) {

  su2double epsilon = 1e-1;

  unsigned short nMarker_ActDiskInlet = config->GetnMarker_ActDiskInlet();

  if (nMarker_ActDiskInlet != 0) {

    unsigned short iMarker, iDim;
    unsigned long iVertex, iPoint, iPointGlobal, pPoint = 0, pPointGlobal = 0, pVertex = 0, pMarker = 0, jVertex, jVertex_, jPoint, jPointGlobal, jMarker;
    su2double *Coord_i, Coord_j[3], dist = 0.0, mindist, maxdist_local = 0.0, maxdist_global = 0.0;
    int iProcessor, pProcessor = 0;
    unsigned long nLocalVertex_ActDisk = 0, MaxLocalVertex_ActDisk = 0;
    int nProcessor = size;
    unsigned short Beneficiary = 0, Donor = 0, iBC;
    bool Perimeter;

    for (iBC = 0; iBC < 2; iBC++) {

      if (iBC == 0) { Beneficiary = ACTDISK_INLET; Donor = ACTDISK_OUTLET; }
      if (iBC == 1) { Beneficiary = ACTDISK_OUTLET; Donor = ACTDISK_INLET; }

      unsigned long *Buffer_Send_nVertex = new unsigned long [1];
      unsigned long *Buffer_Receive_nVertex = new unsigned long [nProcessor];

      if ((iBC == 0) && (rank == MASTER_NODE)) cout << "Set Actuator Disk inlet boundary conditions." << endl;
      if ((iBC == 1) && (rank == MASTER_NODE)) cout << "Set Actuator Disk outlet boundary conditions." << endl;

      /*--- Compute the number of vertex that have an actuator disk outlet boundary condition
       without including the ghost nodes ---*/

      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            if (node[iPoint]->GetDomain()) nLocalVertex_ActDisk ++;
          }
        }
      }

      Buffer_Send_nVertex[0] = nLocalVertex_ActDisk;

      /*--- Send actuator disk vertex information --*/

      SU2_MPI::Allreduce(&nLocalVertex_ActDisk, &MaxLocalVertex_ActDisk, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

      /*--- Array dimensionalization --*/

      su2double *Buffer_Send_Coord = new su2double [MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Send_Point  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_GlobalIndex  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_Vertex  = new unsigned long [MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Send_Marker  = new unsigned long [MaxLocalVertex_ActDisk];

      su2double *Buffer_Receive_Coord = new su2double [nProcessor*MaxLocalVertex_ActDisk*nDim];
      unsigned long *Buffer_Receive_Point = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_Vertex = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];
      unsigned long *Buffer_Receive_Marker = new unsigned long [nProcessor*MaxLocalVertex_ActDisk];

      unsigned long nBuffer_Coord = MaxLocalVertex_ActDisk*nDim;
      unsigned long nBuffer_Point = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_GlobalIndex = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_Vertex = MaxLocalVertex_ActDisk;
      unsigned long nBuffer_Marker = MaxLocalVertex_ActDisk;

      for (iVertex = 0; iVertex < MaxLocalVertex_ActDisk; iVertex++) {
        Buffer_Send_Point[iVertex] = 0;
        Buffer_Send_GlobalIndex[iVertex] = 0;
        Buffer_Send_Vertex[iVertex] = 0;
        Buffer_Send_Marker[iVertex] = 0;
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      }

      /*--- Copy coordinates and point to the auxiliar vector --*/

      nLocalVertex_ActDisk = 0;
      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Donor) {
          for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            iPointGlobal = node[iPoint]->GetGlobalIndex();
            if (node[iPoint]->GetDomain()) {
              Buffer_Send_Point[nLocalVertex_ActDisk] = iPoint;
              Buffer_Send_GlobalIndex[nLocalVertex_ActDisk] = iPointGlobal;
              Buffer_Send_Vertex[nLocalVertex_ActDisk] = iVertex;
              Buffer_Send_Marker[nLocalVertex_ActDisk] = iMarker;
              for (iDim = 0; iDim < nDim; iDim++)
                Buffer_Send_Coord[nLocalVertex_ActDisk*nDim+iDim] = node[iPoint]->GetCoord(iDim);
              nLocalVertex_ActDisk++;
            }
          }
        }
      }

      SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, Buffer_Receive_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
      SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, Buffer_Receive_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

      /*--- Compute the closest point to an actuator disk inlet point ---*/

      maxdist_local = 0.0;

      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        if (config->GetMarker_All_KindBC(iMarker) == Beneficiary) {

          for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
            iPoint = vertex[iMarker][iVertex]->GetNode();
            iPointGlobal = node[iPoint]->GetGlobalIndex();


            if (node[iPoint]->GetDomain()) {

              /*--- Coordinates of the boundary point ---*/

              Coord_i = node[iPoint]->GetCoord(); mindist = 1E6; pProcessor = 0; pPoint = 0;

              /*--- Loop over all the boundaries to find the pair ---*/

              Perimeter = false;

              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex[iProcessor]; jVertex++) {
                  jPoint = Buffer_Receive_Point[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jPointGlobal = Buffer_Receive_GlobalIndex[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jVertex_ = Buffer_Receive_Vertex[iProcessor*MaxLocalVertex_ActDisk+jVertex];
                  jMarker = Buffer_Receive_Marker[iProcessor*MaxLocalVertex_ActDisk+jVertex];

 //                 if (jPointGlobal != iPointGlobal) {
 //                 ActDisk_Perimeter

                    /*--- Compute the distance ---*/

                    dist = 0.0;
                    for (iDim = 0; iDim < nDim; iDim++) {
                      Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_ActDisk+jVertex)*nDim+iDim];
                      dist += pow(Coord_j[iDim]-Coord_i[iDim], 2.0);
                    }
                    dist = sqrt(dist);

                    if (dist < mindist) {
                      mindist = dist; pProcessor = iProcessor; pPoint = jPoint; pPointGlobal = jPointGlobal;
                      pVertex = jVertex_; pMarker = jMarker;
                      if (dist == 0.0) break;
                    }

//                  }
//                  else { Perimeter = true; mindist = 0.0; dist = 0.0; break; }
               }
              }

              /*--- Store the value of the pair ---*/

              maxdist_local = max(maxdist_local, mindist);
              vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);
              vertex[iMarker][iVertex]->SetActDisk_Perimeter(Perimeter);

              if (mindist > epsilon) {
                cout.precision(10);
                cout << endl;
                cout << "   Bad match for point " << iPoint << ".\tNearest";
                cout << " donor distance: " << scientific << mindist << ".";
                vertex[iMarker][iVertex]->SetDonorPoint(iPoint, iPointGlobal, pVertex, pMarker, pProcessor);
                maxdist_local = min(maxdist_local, 0.0);
              }

            }
          }

        }
      }

      SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE, MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);

      if (rank == MASTER_NODE) cout <<"The max distance between points is: " << maxdist_global <<"."<< endl;

      delete[] Buffer_Send_Coord;
      delete[] Buffer_Send_Point;

      delete[] Buffer_Receive_Coord;
      delete[] Buffer_Receive_Point;

      delete[] Buffer_Send_nVertex;
      delete[] Buffer_Receive_nVertex;

      delete [] Buffer_Send_GlobalIndex;
      delete [] Buffer_Send_Vertex;
      delete [] Buffer_Send_Marker;

      delete [] Buffer_Receive_GlobalIndex;
      delete [] Buffer_Receive_Vertex;
      delete [] Buffer_Receive_Marker;

    }
  }

}

void CPhysicalGeometry::MatchPeriodic(CConfig        *config,
                                      unsigned short val_periodic) {

  unsigned short iMarker, iDim, jMarker, pMarker = 0;
  unsigned short iPeriodic, nPeriodic;

  unsigned long iVertex, iPoint, iPointGlobal, index;
  unsigned long jVertex, jVertex_, jPoint, jPointGlobal;
  unsigned long pVertex = 0, pPoint = 0, pPointGlobal = 0;
  unsigned long nLocalVertex_Periodic = 0, MaxLocalVertex_Periodic = 0;
  unsigned long nPointMatch = 0;

  int iProcessor, pProcessor = 0, nProcessor = size;

  bool isBadMatch = false;

  string Marker_Tag;

  su2double *Coord_i, Coord_j[3], dist, mindist, maxdist_local, maxdist_global;
  su2double *center, *angles, translation[3]={0.0,0.0,0.0}, *trans, dx, dy, dz;
  su2double rotMatrix[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double Theta, Phi, Psi, cosTheta, sinTheta, cosPhi, sinPhi, cosPsi, sinPsi;
  su2double rotCoord[3] = {0.0, 0.0, 0.0};

  bool pointOnAxis = false;

  bool chkSamePoint = false;

  su2double distToAxis = 0.0;

  /*--- Tolerance for distance-based match to report warning. ---*/

  su2double epsilon = 1e-6;

  /*--- Evaluate the number of periodic boundary conditions ---*/

  nPeriodic = config->GetnMarker_Periodic();

  /*--- Send an initial message to the console. ---*/

  if (rank == MASTER_NODE) {
    cout << "Matching the periodic boundary points for marker pair ";
    cout << val_periodic << "." << endl;
  }

  /*--- Compute the total number of vertices that sit on a periodic
   boundary on our local rank. We only include our "owned" nodes. ---*/

  nLocalVertex_Periodic = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          if (node[iPoint]->GetDomain()) nLocalVertex_Periodic++;
        }
      }
    }
  }

  /*--- Communicate our local periodic point count globally
   and receive the counts of periodic points from all other ranks.---*/

  unsigned long *Buffer_Send_nVertex = new unsigned long [1];
  unsigned long *Buffer_Recv_nVertex = new unsigned long [nProcessor];

  Buffer_Send_nVertex[0] = nLocalVertex_Periodic;

  /*--- Copy our own count in serial or use collective comms with MPI. ---*/

  SU2_MPI::Allreduce(&nLocalVertex_Periodic, &MaxLocalVertex_Periodic, 1,
                     MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex, 1, MPI_UNSIGNED_LONG,
                     Buffer_Recv_nVertex, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Prepare buffers to send the information for each
   periodic point to all ranks so that we can match pairs. ---*/

  su2double *Buffer_Send_Coord           = new su2double [MaxLocalVertex_Periodic*nDim];
  unsigned long *Buffer_Send_Point       = new unsigned long [MaxLocalVertex_Periodic];
  unsigned long *Buffer_Send_GlobalIndex = new unsigned long [MaxLocalVertex_Periodic];
  unsigned long *Buffer_Send_Vertex      = new unsigned long [MaxLocalVertex_Periodic];
  unsigned long *Buffer_Send_Marker      = new unsigned long [MaxLocalVertex_Periodic];

  su2double *Buffer_Recv_Coord           = new su2double [nProcessor*MaxLocalVertex_Periodic*nDim];
  unsigned long *Buffer_Recv_Point       = new unsigned long [nProcessor*MaxLocalVertex_Periodic];
  unsigned long *Buffer_Recv_GlobalIndex = new unsigned long [nProcessor*MaxLocalVertex_Periodic];
  unsigned long *Buffer_Recv_Vertex      = new unsigned long [nProcessor*MaxLocalVertex_Periodic];
  unsigned long *Buffer_Recv_Marker      = new unsigned long [nProcessor*MaxLocalVertex_Periodic];

  unsigned long nBuffer_Coord       = MaxLocalVertex_Periodic*nDim;
  unsigned long nBuffer_Point       = MaxLocalVertex_Periodic;
  unsigned long nBuffer_GlobalIndex = MaxLocalVertex_Periodic;
  unsigned long nBuffer_Vertex      = MaxLocalVertex_Periodic;
  unsigned long nBuffer_Marker      = MaxLocalVertex_Periodic;

  for (iVertex = 0; iVertex < MaxLocalVertex_Periodic; iVertex++) {
    Buffer_Send_Point[iVertex]       = 0;
    Buffer_Send_GlobalIndex[iVertex] = 0;
    Buffer_Send_Vertex[iVertex]      = 0;
    Buffer_Send_Marker[iVertex]      = 0;
    for (iDim = 0; iDim < nDim; iDim++)
    Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
  }

  /*--- Store the local index, global index, local boundary index,
   marker index, and point coordinates in the buffers for sending.
   Note again that this is only for the current pair of periodic
   markers and for only the "owned" points on each rank. ---*/

  nLocalVertex_Periodic = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {
      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {
        for (iVertex = 0; iVertex < GetnVertex(iMarker); iVertex++) {
          iPoint = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();
          if (node[iPoint]->GetDomain()) {
            Buffer_Send_Point[nLocalVertex_Periodic] = iPoint;
            Buffer_Send_GlobalIndex[nLocalVertex_Periodic] = iPointGlobal;
            Buffer_Send_Vertex[nLocalVertex_Periodic] = iVertex;
            Buffer_Send_Marker[nLocalVertex_Periodic] = iMarker;
            for (iDim = 0; iDim < nDim; iDim++)
            Buffer_Send_Coord[nLocalVertex_Periodic*nDim+iDim] = node[iPoint]->GetCoord(iDim);
            nLocalVertex_Periodic++;
          }
        }
      }
    }
  }

  /*--- Copy our own data in serial or use collective comms to gather
   the data for all points on each rank with MPI. Note that, since the
   periodic point count should be small relative to the volume grid
   and we are only storing one periodic marker pair at a time,
   repeating the data for each pair on all ranks should be manageable. ---*/

  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE,
                     Buffer_Recv_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_Point, nBuffer_Point, MPI_UNSIGNED_LONG,
                     Buffer_Recv_Point, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG,
                     Buffer_Recv_GlobalIndex, nBuffer_GlobalIndex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG,
                     Buffer_Recv_Vertex, nBuffer_Vertex, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG,
                     Buffer_Recv_Marker, nBuffer_Marker, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  /*--- Now that all ranks have the data for all periodic points for
   this pair of periodic markers, we match the individual points
   based on the translation / rotation specified for the marker pair. ---*/

  maxdist_local = 0.0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == PERIODIC_BOUNDARY) {

      iPeriodic = config->GetMarker_All_PerBound(iMarker);
      if ((iPeriodic == val_periodic) ||
          (iPeriodic == val_periodic + nPeriodic/2)) {

        /*--- Retrieve the supplied periodic information. ---*/

        Marker_Tag = config->GetMarker_All_TagBound(iMarker);
        center     = config->GetPeriodicRotCenter(Marker_Tag);
        angles     = config->GetPeriodicRotAngles(Marker_Tag);
        trans      = config->GetPeriodicTranslation(Marker_Tag);

        /*--- Store (center+trans) as it is constant and will be added. ---*/

        translation[0] = center[0] + trans[0];
        translation[1] = center[1] + trans[1];
        translation[2] = center[2] + trans[2];

        /*--- Store angles separately for clarity. Compute sines/cosines. ---*/

        Theta    = angles[0];      Phi = angles[1];     Psi = angles[2];
        cosTheta = cos(Theta);  cosPhi = cos(Phi);   cosPsi = cos(Psi);
        sinTheta = sin(Theta);  sinPhi = sin(Phi);   sinPsi = sin(Psi);

        /*--- Compute the rotation matrix. Note that the implicit
         ordering is rotation about the x-axis, y-axis, then z-axis. ---*/

        rotMatrix[0][0] = cosPhi*cosPsi;
        rotMatrix[1][0] = cosPhi*sinPsi;
        rotMatrix[2][0] = -sinPhi;

        rotMatrix[0][1] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;
        rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;
        rotMatrix[2][1] = sinTheta*cosPhi;

        rotMatrix[0][2] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
        rotMatrix[1][2] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
        rotMatrix[2][2] = cosTheta*cosPhi;

        /*--- Loop over each point on the periodic marker that this rank
         holds locally and find the matching point from the donor marker. ---*/

        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

          /*--- Local and global index for the owned periodic point. ---*/

          iPoint       = vertex[iMarker][iVertex]->GetNode();
          iPointGlobal = node[iPoint]->GetGlobalIndex();

          /*--- If this is not a ghost, find the periodic match. ---*/

          if (node[iPoint]->GetDomain()) {

            /*--- Coordinates of the current boundary point ---*/

            Coord_i = node[iPoint]->GetCoord();

            /*--- Get the position vector from rotation center to point. ---*/

            dx = Coord_i[0] - center[0];
            dy = Coord_i[1] - center[1];
            if (nDim == 3) dz = Coord_i[2] - center[2];
            else           dz = 0.0;

            /*--- Compute transformed point coordinates. ---*/

            rotCoord[0] = (rotMatrix[0][0]*dx +
                           rotMatrix[0][1]*dy +
                           rotMatrix[0][2]*dz + translation[0]);

            rotCoord[1] = (rotMatrix[1][0]*dx +
                           rotMatrix[1][1]*dy +
                           rotMatrix[1][2]*dz + translation[1]);

            rotCoord[2] = (rotMatrix[2][0]*dx +
                           rotMatrix[2][1]*dy +
                           rotMatrix[2][2]*dz + translation[2]);

            /*--- Check if the point lies on the axis of rotation. If it does,
             the rotated coordinate and the original coordinate are the same. ---*/

            pointOnAxis = false;
            distToAxis = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
               distToAxis = (rotCoord[iDim] - Coord_i[iDim])*(rotCoord[iDim] - Coord_i[iDim]);
            distToAxis = sqrt(distToAxis);

            if (distToAxis < epsilon) pointOnAxis = true;

            /*--- Our search is based on the minimum distance, so we
             initialize the distance to a large value. ---*/

            mindist = 1E6; pProcessor = 0; pPoint = 0;

            /*--- Loop over all of the periodic data that was gathered from
             all ranks in order to find the matching periodic point. ---*/

            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++)
            for (jVertex = 0; jVertex < Buffer_Recv_nVertex[iProcessor]; jVertex++) {

              /*--- Store the loop index more easily. ---*/

              index = iProcessor*MaxLocalVertex_Periodic + jVertex;

              /*--- For each candidate, we have the local and global index,
               along with the boundary vertex and marker index. ---*/

              jPoint       = Buffer_Recv_Point[index];
              jPointGlobal = Buffer_Recv_GlobalIndex[index];
              jVertex_     = Buffer_Recv_Vertex[index];
              jMarker      = Buffer_Recv_Marker[index];

              /*--- The gathered data will also include the current
               "owned" periodic point that we are matching, so first make
               sure that we avoid the original point by checking that the
               global index values are not the same. ---*/

              if ((jPointGlobal != iPointGlobal) || (pointOnAxis)) {

                /*--- Compute the distance between the candidate periodic
                 point and the transformed coordinates of the owned point. ---*/

                dist = 0.0;
                for (iDim = 0; iDim < nDim; iDim++) {
                  Coord_j[iDim] = Buffer_Recv_Coord[index*nDim + iDim];
                  dist         += pow(Coord_j[iDim]-rotCoord[iDim],2.0);
                }
                dist = sqrt(dist);

                /*--- Compare the distance against the existing minimum
                 and also perform checks just to be sure that this is an
                 independent periodic point (even if on the same rank),
                  unless it lies on the axis of rotation. ---*/

                chkSamePoint = false;
                chkSamePoint = (((dist < mindist) && (iProcessor != rank)) ||
                                ((dist < mindist) && (iProcessor == rank) &&
                                (jPoint != iPoint)));

                if (chkSamePoint || ((dist < mindist) && (pointOnAxis))) {

                  /*--- We have found an intermediate match. Store the
                   data for this point before continuing the search. ---*/

                  mindist      = dist;
                  pProcessor   = iProcessor;
                  pPoint       = jPoint;
                  pPointGlobal = jPointGlobal;
                  pVertex      = jVertex_;
                  pMarker      = jMarker;

                }
              }

            }

            /*--- Store the data for the best match found for the
             owned periodic point. ---*/

            vertex[iMarker][iVertex]->SetDonorPoint(pPoint, pPointGlobal, pVertex, pMarker, pProcessor);
            maxdist_local = max(maxdist_local, mindist);
            nPointMatch++;

            /*--- If the distance to the closest point is larger than our
             tolerance, then throw a warning for this point. ---*/

            if (mindist > epsilon) {
              cout.precision(10);
              cout << endl;
              cout << "   Bad match for point " << iPointGlobal << ".\tNearest";
              cout << " donor distance: " << scientific << mindist << ".";
              maxdist_local = min(maxdist_local, 0.0);
              isBadMatch = true;
            }

          }
        }
      }
    }
  }

  /*--- Communicate the final count of number of matched points
   for the periodic boundary pair and the max distance for all
   pairs of points. ---*/

  unsigned long nPointMatch_Local = nPointMatch;
  SU2_MPI::Reduce(&nPointMatch_Local, &nPointMatch, 1, MPI_UNSIGNED_LONG,
                  MPI_SUM, MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Reduce(&maxdist_local, &maxdist_global, 1, MPI_DOUBLE,
                  MPI_MAX, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Output some information about the matching process. ---*/

  if (rank == MASTER_NODE) {
    if (nPointMatch > 0) {
      cout <<" Matched " << nPointMatch << " points with a max distance of: ";
      cout << maxdist_global <<"."<< endl;
    } else {
      cout <<" No matching points for periodic marker pair ";
      cout << val_periodic << " in current zone." << endl;
    }

    /*--- Print final warning when finding bad matches. ---*/

    if (isBadMatch) {
      cout << endl;
      cout << "\n !!! Warning !!!" << endl;
      cout << "Bad matches found. Computation will continue, but be cautious.\n";
    }
  }

  /*--- Free local memory for communications. ---*/

  delete[] Buffer_Send_Coord;
  delete[] Buffer_Send_Point;

  delete[] Buffer_Recv_Coord;
  delete[] Buffer_Recv_Point;

  delete[] Buffer_Send_nVertex;
  delete[] Buffer_Recv_nVertex;

  delete [] Buffer_Send_GlobalIndex;
  delete [] Buffer_Send_Vertex;
  delete [] Buffer_Send_Marker;

  delete [] Buffer_Recv_GlobalIndex;
  delete [] Buffer_Recv_Vertex;
  delete [] Buffer_Recv_Marker;

}

void CPhysicalGeometry::SetControlVolume(CConfig *config, unsigned short action) {
  unsigned long face_iPoint = 0, face_jPoint = 0, iPoint, iElem;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  su2double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint, *Coord_FacejPoint, Area,
  Volume, DomainVolume, my_DomainVolume;
  bool change_face_orientation;

  /*--- Update values of faces of the edge ---*/
  if (action != ALLOCATE) {
    edges->SetZeroValues();
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetVolume (0.0);
  }

  Coord_Edge_CG = new su2double [nDim];
  Coord_FaceElem_CG = new su2double [nDim];
  Coord_Elem_CG = new su2double [nDim];
  Coord_FaceiPoint = new su2double [nDim];
  Coord_FacejPoint = new su2double [nDim];

  my_DomainVolume = 0.0;
  for (iElem = 0; iElem < nElem; iElem++) {
    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {

      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {

        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }

        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }

        /*--- We define a direction (from the smalest index to the greatest) --*/
        change_face_orientation = false;
        if (face_iPoint > face_jPoint) change_face_orientation = true;
        iEdge = FindEdge(face_iPoint, face_jPoint);

        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edges->GetCG(iEdge,iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }

        switch (nDim) {
          case 2:
            /*--- Two dimensional problem ---*/
            if (change_face_orientation) edges->SetNodes_Coord(iEdge, Coord_Elem_CG, Coord_Edge_CG);
            else edges->SetNodes_Coord(iEdge, Coord_Edge_CG, Coord_Elem_CG);
            Area = CEdge::GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Area); my_DomainVolume +=Area;
            Area = CEdge::GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Area); my_DomainVolume +=Area;
            break;
          case 3:
            /*--- Three dimensional problem ---*/
            if (change_face_orientation) edges->SetNodes_Coord(iEdge, Coord_FaceElem_CG, Coord_Edge_CG, Coord_Elem_CG);
            else edges->SetNodes_Coord(iEdge, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            Volume = CEdge::GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_iPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            Volume = CEdge::GetVolume(Coord_FacejPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
            node[face_jPoint]->AddVolume(Volume); my_DomainVolume +=Volume;
            break;
        }
      }
    }
  }

  /*--- Check if there is a normal with null area ---*/
  for (iEdge = 0; iEdge < (long)nEdge; iEdge++) {
    const auto NormalFace = edges->GetNormal(iEdge);
    Area = 0.0; for (iDim = 0; iDim < nDim; iDim++) Area += NormalFace[iDim]*NormalFace[iDim];
    Area = sqrt(Area);
    if (Area == 0.0) {
      su2double DefaultArea[3] = {EPS*EPS};
      edges->SetNormal(iEdge, DefaultArea);
    }
  }

  SU2_MPI::Allreduce(&my_DomainVolume, &DomainVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if ((rank == MASTER_NODE) && (action == ALLOCATE)) {
    if (nDim == 2) cout <<"Area of the computational grid: "<< DomainVolume <<"."<< endl;
    if (nDim == 3) cout <<"Volume of the computational grid: "<< DomainVolume <<"."<< endl;
  }

  config->SetDomainVolume(DomainVolume);

  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;
}

void CPhysicalGeometry::VisualizeControlVolume(CConfig *config, unsigned short action) {

  /*--- This routine is only meant for visualization in serial currently ---*/
#ifndef HAVE_MPI

  unsigned long face_iPoint = 0, face_jPoint = 0, iElem, iPoint_Viz;
  long iEdge;
  unsigned short nEdgesFace = 1, iFace, iEdgesFace, iDim;
  su2double *Coord_Edge_CG, *Coord_FaceElem_CG, *Coord_Elem_CG, *Coord_FaceiPoint,
  *Coord_FacejPoint;
  int counter = 0;
  char cstr[MAX_STRING_SIZE], buffer[50];
  ofstream Tecplot_File;
  string mesh_filename;
  vector<su2double> X, Y, Z, X_n, Y_n, Z_n;
  su2double r1[3], r2[3], CrossProduct[3];

  /*--- Access the point number for control volume we want to vizualize ---*/

  iPoint_Viz = config->GetVisualize_CV();

  /*--- Allocate some structures for building the dual CVs ---*/

  Coord_Edge_CG     = new su2double [nDim];
  Coord_FaceElem_CG = new su2double [nDim];
  Coord_Elem_CG     = new su2double [nDim];
  Coord_FaceiPoint  = new su2double [nDim];
  Coord_FacejPoint  = new su2double [nDim];

  /*--- Loop over each face of each element ---*/

  CrossProduct[0] = 0.0; CrossProduct[1] = 0.0; CrossProduct[2] = 0.0;

  for (iElem = 0; iElem < nElem; iElem++) {

    for (iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {

      /*--- In 2D all the faces have only one edge ---*/
      if (nDim == 2) nEdgesFace = 1;
      /*--- In 3D the number of edges per face is the same as the number of point per face ---*/
      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

      /*-- Loop over the edges of a face ---*/
      for (iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {

        /*--- In 2D only one edge (two points) per edge ---*/
        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }

        /*--- In 3D there are several edges in each face ---*/
        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }

        /*--- We define a direction (from the smallest index to the greatest) --*/
        iEdge = FindEdge(face_iPoint, face_jPoint);

        for (iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim] = edges->GetCG(iEdge,iDim);
          Coord_Elem_CG[iDim] = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
          Coord_FaceiPoint[iDim] = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim] = node[face_jPoint]->GetCoord(iDim);
        }

        /*--- Print out the coordinates for a set of triangles making
         up a single dual control volume for visualization. ---*/

        if (face_iPoint == iPoint_Viz || face_jPoint == iPoint_Viz) {

          if (nDim == 2) {
            X.push_back(Coord_Elem_CG[0]); X.push_back(Coord_Edge_CG[0]);
            Y.push_back(Coord_Elem_CG[1]); Y.push_back(Coord_Edge_CG[1]);
          } else if (nDim == 3) {
            X.push_back(Coord_FaceElem_CG[0]); X.push_back(Coord_Edge_CG[0]); X.push_back(Coord_Elem_CG[0]);
            Y.push_back(Coord_FaceElem_CG[1]); Y.push_back(Coord_Edge_CG[1]); Y.push_back(Coord_Elem_CG[1]);
            Z.push_back(Coord_FaceElem_CG[2]); Z.push_back(Coord_Edge_CG[2]); Z.push_back(Coord_Elem_CG[2]);

            for (iDim = 0; iDim < nDim; iDim++) {
              r1[iDim] = Coord_FaceElem_CG[iDim]-Coord_Elem_CG[iDim];
              r2[iDim] = Coord_Edge_CG[iDim]-Coord_Elem_CG[iDim];
            }
            CrossProduct[0] += 0.5*(r1[1]*r2[2] - r1[2]*r2[1]);
            CrossProduct[1] += 0.5*(r1[2]*r2[0] - r1[0]*r2[2]);
            CrossProduct[2] += 0.5*(r1[0]*r2[1] - r1[1]*r2[0]);
          }
          counter++;
        }
      }
    }
  }

  /*--- Write a Tecplot file to visualize the CV ---*/

  strcpy(cstr,"dual_cv");
  SPRINTF (buffer, "_%d.dat", SU2_TYPE::Int(iPoint_Viz));
  strcat(cstr, buffer);

  Tecplot_File.open(cstr, ios::out);
  Tecplot_File << "TITLE= \"Visualization of the control volume\"" << endl;

  if (nDim == 2) {
    Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*2 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
  } if (nDim == 3) {
    Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
    Tecplot_File << "ZONE NODES= "<< counter*3 <<", ELEMENTS= ";
    Tecplot_File << counter <<", DATAPACKING=POINT, ZONETYPE=FEBRICK"<< endl;
  }

  /*--- Write coordinates for the nodes in the order that they were found
   for each of the edges/triangles making up a dual control volume. ---*/

  for (vector<su2double>::size_type i = 0; i != X.size(); i++) {
    Tecplot_File << X[i] << "\t" << Y[i];
    if (nDim == 3) Tecplot_File << "\t" << Z[i];
    Tecplot_File << "\n";
  }

  /*--- Create a new connectivity table in the order the faces were found ---*/

  int j;
  for (int i= 0; i < counter; i++) {
    if (nDim == 2) {
      j = i*2;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+2 <<"\t"<<j+2 << endl;
    } if (nDim == 3) {
      j = i*3;
      Tecplot_File << j+1 <<"\t"<<j+2 <<"\t"<<j+3 <<"\t"<<j+3 <<"\t";
      Tecplot_File << j+3<<"\t" <<j+3 <<"\t"<<j+3 <<"\t"<<j+3 << endl;
    }
  }

  Tecplot_File.close();
  X.clear();
  Y.clear();
  Z.clear();

  delete[] Coord_Edge_CG;
  delete[] Coord_FaceElem_CG;
  delete[] Coord_Elem_CG;
  delete[] Coord_FaceiPoint;
  delete[] Coord_FacejPoint;

#endif

}

void CPhysicalGeometry::SetMeshFile (CConfig *config, string val_mesh_out_filename) {
  unsigned long iElem, iPoint, iElem_Bound;
  unsigned short iMarker, iNodes, iDim;
  ofstream output_file;
  string Grid_Marker;
  char *cstr;

  cstr = new char [val_mesh_out_filename.size()+1];
  strcpy (cstr, val_mesh_out_filename.c_str());

  /*--- Open .su2 grid file ---*/

  output_file.precision(15);
  output_file.open(cstr, ios::out);

  /*--- Write dimension, number of elements and number of points ---*/

  output_file << "NDIME= " << nDim << endl;
  output_file << "NELEM= " << nElem << endl;
  for (iElem = 0; iElem < nElem; iElem++) {
    output_file << elem[iElem]->GetVTK_Type();
    for (iNodes = 0; iNodes < elem[iElem]->GetnNodes(); iNodes++)
      output_file << "\t" << elem[iElem]->GetNode(iNodes);
    output_file << "\t"<<iElem<< endl;
  }

  /*--- Write the node coordinates ---*/

  output_file << "NPOIN= " << nPoint << "\t" << nPointDomain << endl;
  output_file.precision(15);
  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      output_file << scientific << "\t" << node[iPoint]->GetCoord(iDim) ;
#ifndef HAVE_MPI
    output_file << "\t" << iPoint << endl;
#else
    output_file << "\t" << iPoint << "\t" << node[iPoint]->GetGlobalIndex() << endl;
#endif

  }

  /*--- Loop through and write the boundary info ---*/

  output_file << "NMARK= " << nMarker << endl;
  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- Ignore SEND_RECEIVE for the moment ---*/
    if (bound[iMarker][0]->GetVTK_Type() != VERTEX) {

      Grid_Marker = config->GetMarker_All_TagBound(iMarker);
      output_file << "MARKER_TAG= " << Grid_Marker << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;

      if (nDim == 2) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file << iElem_Bound << endl;
        }
      }

      if (nDim == 3) {
        for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
          output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" ;
          for (iNodes = 0; iNodes < bound[iMarker][iElem_Bound]->GetnNodes(); iNodes++)
            output_file << bound[iMarker][iElem_Bound]->GetNode(iNodes) << "\t" ;
          output_file << iElem_Bound << endl;
        }
      }

    } else if (bound[iMarker][0]->GetVTK_Type() == VERTEX) {
      output_file << "MARKER_TAG= SEND_RECEIVE" << endl;
      output_file << "MARKER_ELEMS= " << nElem_Bound[iMarker]<< endl;
      if (config->GetMarker_All_SendRecv(iMarker) > 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;
      if (config->GetMarker_All_SendRecv(iMarker) < 0) output_file << "SEND_TO= " << config->GetMarker_All_SendRecv(iMarker) << endl;

      for (iElem_Bound = 0; iElem_Bound < nElem_Bound[iMarker]; iElem_Bound++) {
        output_file << bound[iMarker][iElem_Bound]->GetVTK_Type() << "\t" <<
        bound[iMarker][iElem_Bound]->GetNode(0) << "\t" <<
        bound[iMarker][iElem_Bound]->GetRotation_Type() << endl;
      }

    }
  }

  output_file.close();
}

void CPhysicalGeometry::SetCoord_Smoothing (unsigned short val_nSmooth, su2double val_smooth_coeff, CConfig *config) {
  unsigned short iSmooth, nneigh, iMarker;
  su2double *Coord_Old, *Coord_Sum, *Coord, *Coord_i, *Coord_j, Position_Plane = 0.0;
  unsigned long iEdge, iPoint, jPoint, iVertex;
  su2double eps = 1E-6;
  bool NearField = false;

  Coord = new su2double [nDim];

  for (iPoint = 0; iPoint < GetnPoint(); iPoint++) {
    su2double *Coord = node[iPoint]->GetCoord();
    node[iPoint]->SetCoord_Old(Coord);
  }

  /*--- Jacobi iterations ---*/
  for (iSmooth = 0; iSmooth < val_nSmooth; iSmooth++) {

    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint]->SetCoord_SumZero();


    /*--- Loop over Interior edges ---*/
    for (iEdge = 0; iEdge < nEdge; iEdge++) {
      iPoint = edges->GetNode(iEdge,0);
      Coord_i = node[iPoint]->GetCoord();

      jPoint = edges->GetNode(iEdge,1);
      Coord_j = node[jPoint]->GetCoord();

      /*--- Accumulate nearest neighbor Coord to Res_sum for each variable ---*/
      node[iPoint]->AddCoord_Sum(Coord_j);
      node[jPoint]->AddCoord_Sum(Coord_i);

    }

    /*--- Loop over all mesh points (Update Coords with averaged sum) ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      nneigh = node[iPoint]->GetnPoint();
      Coord_Sum = node[iPoint]->GetCoord_Sum();
      Coord_Old = node[iPoint]->GetCoord_Old();

      if (nDim == 2) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        if ((NearField) && ((Coord_Old[1] > Position_Plane-eps) && (Coord_Old[1] < Position_Plane+eps)))
          Coord[1] = Coord_Old[1];
      }

      if (nDim == 3) {
        Coord[0] =(Coord_Old[0] + val_smooth_coeff*Coord_Sum[0]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[1] =(Coord_Old[1] + val_smooth_coeff*Coord_Sum[1]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        Coord[2] =(Coord_Old[2] + val_smooth_coeff*Coord_Sum[2]) /(1.0 + val_smooth_coeff*su2double(nneigh));
        if ((NearField) && ((Coord_Old[2] > Position_Plane-eps) && (Coord_Old[2] < Position_Plane+eps)))
          Coord[2] = Coord_Old[2];
      }

      node[iPoint]->SetCoord(Coord);
    }

    /*--- Copy boundary values ---*/
    for (iMarker = 0; iMarker < nMarker; iMarker++)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        iPoint = vertex[iMarker][iVertex]->GetNode();
        Coord_Old = node[iPoint]->GetCoord_Old();
        node[iPoint]->SetCoord(Coord_Old);
      }
  }

  delete[] Coord;
}

bool CPhysicalGeometry::FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                                 unsigned short &face_second_elem) {

  /*--- Find repeated nodes between two elements to identify the common face ---*/
  unsigned long iPoint = 0, jPoint = 0;
  unsigned short face_node, iFace, iNode, jNode, nNodesFace;
  vector<unsigned long> CommonPoints, PointFaceFirst, PointFaceSecond;
  vector<unsigned long>::iterator IterPoint;
  pair<vector <unsigned long>::iterator, vector <unsigned long>::iterator> mypair;
  bool face_first_found = false, face_second_found =false;

  if (first_elem == second_elem) return false;

  for (iNode = 0; iNode < elem[first_elem]->GetnNodes(); iNode++) {
    iPoint = elem[first_elem]->GetNode(iNode);
    for (jNode = 0; jNode < elem[second_elem]->GetnNodes(); jNode++) {
      jPoint = elem[second_elem]->GetNode(jNode);
      if (iPoint == jPoint) {
        CommonPoints.push_back(iPoint);
        break;
      }
    }
  }

  /*--- Sort point in face and check that the list is unique ---*/
  sort( CommonPoints.begin(), CommonPoints.end());
  IterPoint = unique( CommonPoints.begin(), CommonPoints.end());
  CommonPoints.resize( distance(CommonPoints.begin(), IterPoint) );

  /*--- In 2D, the two elements must share two points that make up
   an edge, as all "faces" are edges in 2D. In 3D, we need to find
   exactly 3 (tri) or 4 (quad) common points. Return immediately to
   avoid a memory issue due to vectors of different lengths below. ---*/

  if ((nDim == 2) && (CommonPoints.size() != 2)) return false;
  if ((nDim == 3) && ((CommonPoints.size() != 3) &&
                      (CommonPoints.size() != 4))) return false;

  /*--- Search the sequence in the first element ---*/
  for (iFace = 0; iFace < elem[first_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[first_elem]->GetnNodesFace(iFace);

    if (nNodesFace == CommonPoints.size()) {
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[first_elem]->GetFaces(iFace, iNode);
      PointFaceFirst.push_back(elem[first_elem]->GetNode(face_node));
    }

    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceFirst.begin(), PointFaceFirst.end());

    /*--- List comparison ---*/
    mypair = mismatch (PointFaceFirst.begin(), PointFaceFirst.end(), CommonPoints.begin());
    if (mypair.first == PointFaceFirst.end()) {
      face_first_elem = iFace;
      face_first_found = true;
      break;
    }

    PointFaceFirst.erase (PointFaceFirst.begin(), PointFaceFirst.end());
  }
  }

  /*--- Search the secuence in the second element ---*/
  for (iFace = 0; iFace < elem[second_elem]->GetnFaces(); iFace++) {
    nNodesFace = elem[second_elem]->GetnNodesFace(iFace);

    if (nNodesFace == CommonPoints.size()) {
    for (iNode = 0; iNode < nNodesFace; iNode++) {
      face_node = elem[second_elem]->GetFaces(iFace, iNode);
      PointFaceSecond.push_back(elem[second_elem]->GetNode(face_node));
    }

    /*--- Sort face_poin to perform comparison ---*/
    sort( PointFaceSecond.begin(), PointFaceSecond.end());

    /*--- List comparison ---*/
    mypair = mismatch (PointFaceSecond.begin(), PointFaceSecond.end(), CommonPoints.begin());
    if (mypair.first == PointFaceSecond.end()) {
      face_second_elem = iFace;
      face_second_found = true;
      break;
    }

    PointFaceSecond.erase (PointFaceSecond.begin(), PointFaceSecond.end());
  }
  }

  if (face_first_found && face_second_found) return true;
  else return false;

}

void CPhysicalGeometry::SetTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file) {

  unsigned long iElem, iPoint;
  unsigned short iDim;
  ofstream Tecplot_File;

  /*--- Open the tecplot file and write the header ---*/

  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the volumetric grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);

  Tecplot_File << "ZONE T= ";
  if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
  else Tecplot_File << "\"Deformed grid\", C=RED, ";
  Tecplot_File << "NODES= "<< nPoint <<", ELEMENTS= "<< nElem <<", DATAPACKING= POINT";
  if (nDim == 2) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;
  if (nDim == 3) Tecplot_File << ", ZONETYPE= FEBRICK"<< endl;

  /*--- Adding coordinates ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    for (iDim = 0; iDim < nDim; iDim++)
      Tecplot_File << scientific << node[iPoint]->GetCoord(iDim) << "\t";
    Tecplot_File << "\n";
  }

  /*--- Adding conectivity ---*/

  for (iElem = 0; iElem < nElem; iElem++) {
    if (elem[iElem]->GetVTK_Type() == TRIANGLE) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == QUADRILATERAL) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == TETRAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(3)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == HEXAHEDRON) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 <<" "<<
      elem[iElem]->GetNode(6)+1 <<" "<< elem[iElem]->GetNode(7)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PYRAMID) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(2)+1 <<" "<< elem[iElem]->GetNode(3)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(4)+1 << endl;
    }
    if (elem[iElem]->GetVTK_Type() == PRISM) {
      Tecplot_File <<
      elem[iElem]->GetNode(0)+1 <<" "<< elem[iElem]->GetNode(1)+1 <<" "<<
      elem[iElem]->GetNode(1)+1 <<" "<< elem[iElem]->GetNode(2)+1 <<" "<<
      elem[iElem]->GetNode(3)+1 <<" "<< elem[iElem]->GetNode(4)+1 <<" "<<
      elem[iElem]->GetNode(4)+1 <<" "<< elem[iElem]->GetNode(5)+1 << endl;
    }
  }

  Tecplot_File.close();
}

void CPhysicalGeometry::SetBoundTecPlot(char mesh_filename[MAX_STRING_SIZE], bool new_file, CConfig *config) {

  ofstream Tecplot_File;
  unsigned long iPoint, Total_nElem_Bound, iElem, *PointSurface = NULL, nPointSurface = 0;
  unsigned short Coord_i, iMarker;

  /*--- It is important to do a renumbering to don't add points
   that do not belong to the surfaces ---*/

  PointSurface = new unsigned long[nPoint];
  for (iPoint = 0; iPoint < nPoint; iPoint++)
    if (node[iPoint]->GetBoundary()) {
      PointSurface[iPoint] = nPointSurface;
      nPointSurface++;
    }

  /*--- Compute the total number of elements ---*/

  Total_nElem_Bound = 0;
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
    if (config->GetMarker_All_Plotting(iMarker) == YES) {
      Total_nElem_Bound += nElem_Bound[iMarker];
    }
  }

  /*--- Open the tecplot file and write the header ---*/

  if (new_file) {
    Tecplot_File.open(mesh_filename, ios::out);
    Tecplot_File << "TITLE= \"Visualization of the surface grid\"" << endl;
    if (nDim == 2) Tecplot_File << "VARIABLES = \"x\",\"y\" " << endl;
    if (nDim == 3) Tecplot_File << "VARIABLES = \"x\",\"y\",\"z\" " << endl;
  }
  else Tecplot_File.open(mesh_filename, ios::out | ios::app);

  if (Total_nElem_Bound != 0) {

    /*--- Write the header of the file ---*/

    Tecplot_File << "ZONE T= ";
    if (new_file) Tecplot_File << "\"Original grid\", C=BLACK, ";
    else Tecplot_File << "\"Deformed grid\", C=RED, ";
    Tecplot_File << "NODES= "<< nPointSurface <<", ELEMENTS= "<< Total_nElem_Bound <<", DATAPACKING= POINT";
    if (nDim == 2) Tecplot_File << ", ZONETYPE= FELINESEG"<< endl;
    if (nDim == 3) Tecplot_File << ", ZONETYPE= FEQUADRILATERAL"<< endl;

    /*--- Only write the coordiantes of the points that are on the surfaces ---*/

    if (nDim == 3) {
      for (iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for (Coord_i = 0; Coord_i < nDim-1; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << node[iPoint]->GetCoord(nDim-1) << "\n";
        }
    }
    else {
      for (iPoint = 0; iPoint < nPoint; iPoint++)
        if (node[iPoint]->GetBoundary()) {
          for (Coord_i = 0; Coord_i < nDim; Coord_i++)
            Tecplot_File << node[iPoint]->GetCoord(Coord_i) << " ";
          Tecplot_File << "\n";
        }
    }

    /*--- Write the cells using the new numbering ---*/

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++)
      if (config->GetMarker_All_Plotting(iMarker) == YES)
        for (iElem = 0; iElem < nElem_Bound[iMarker]; iElem++) {
          if (nDim == 2) {
            Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
            << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << endl;
          }
          if (nDim == 3) {
            if (bound[iMarker][iElem]->GetnNodes() == 3) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << endl;
            }
            if (bound[iMarker][iElem]->GetnNodes() == 4) {
              Tecplot_File << PointSurface[bound[iMarker][iElem]->GetNode(0)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(1)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(2)]+1 << " "
              << PointSurface[bound[iMarker][iElem]->GetNode(3)]+1 << endl;
            }
          }
        }
  }
  else {

    /*--- No elements in the surface ---*/

    if (nDim == 2) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FELINESEG"<< endl;
      Tecplot_File << "0.0 0.0"<< endl;
      Tecplot_File << "1 1"<< endl;
    }
    if (nDim == 3) {
      Tecplot_File << "ZONE NODES= 1, ELEMENTS= 1, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL"<< endl;
      Tecplot_File << "0.0 0.0 0.0"<< endl;
      Tecplot_File << "1 1 1 1"<< endl;
    }
  }

  /*--- Dealocate memory and close the file ---*/

  delete[] PointSurface;
  Tecplot_File.close();

}

void CPhysicalGeometry::SetColorGrid_Parallel(CConfig *config) {

  /*--- Initialize the color vector ---*/

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint]->SetColor(0);

  /*--- We need to have parallel support with MPI and have the ParMETIS
   library compiled and linked for parallel graph partitioning. ---*/

#ifdef HAVE_MPI
#ifdef HAVE_PARMETIS

  MPI_Comm comm = MPI_COMM_WORLD;

  /*--- Linear partitioner object to help prepare parmetis data. ---*/

  CLinearPartitioner pointPartitioner(Global_nPointDomain,0);

  /*--- Only call ParMETIS if we have more than one rank to avoid errors ---*/

  if (size > SINGLE_NODE) {

    /*--- Create some structures that ParMETIS needs for partitioning. ---*/

    idx_t numflag, nparts, edgecut, wgtflag, ncon;

    idx_t *vtxdist = new idx_t[size+1];
    idx_t *part    = new idx_t[nPoint];

    real_t ubvec;
    real_t *tpwgts = new real_t[size];

    /*--- Some recommended defaults for the various ParMETIS options. ---*/

    wgtflag = 0;
    numflag = 0;
    ncon    = 1;
    ubvec   = 1.05;
    nparts  = (idx_t)size;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[1] = 0;

    /*--- Fill the necessary ParMETIS data arrays. We do not apply
     any weighting during the partitioning process. ---*/

    for (int i = 0; i < size; i++) {
      tpwgts[i] = 1.0/((real_t)size);
    }

    vtxdist[0] = 0;
    for (int i = 0; i < size; i++) {
      vtxdist[i+1] = (idx_t)pointPartitioner.GetLastIndexOnRank(i);
    }

    /*--- Calling ParMETIS ---*/

    if (rank == MASTER_NODE) cout << "Calling ParMETIS...";
    ParMETIS_V3_PartKway(vtxdist, xadj, adjacency, NULL, NULL, &wgtflag,
                         &numflag, &ncon, &nparts, tpwgts, &ubvec, options,
                         &edgecut, part, &comm);
    if (rank == MASTER_NODE) {
      cout << " graph partitioning complete (";
      cout << edgecut << " edge cuts)." << endl;
    }

    /*--- Store the results of the partitioning (note that this is local
     since each processor is calling ParMETIS in parallel and storing the
     results for its initial piece of the grid. ---*/

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint]->SetColor(part[iPoint]);
    }

    /*--- Free all memory needed for the ParMETIS structures ---*/

    if (vtxdist != NULL) delete [] vtxdist;
    if (part    != NULL) delete [] part;
    if (tpwgts  != NULL) delete [] tpwgts;

  }

  /*--- Delete the memory from the geometry class that carried the
   adjacency structure. ---*/

  if (xadj      != NULL) delete [] xadj;
  if (adjacency != NULL) delete [] adjacency;

#endif
#endif

}

void CPhysicalGeometry::ComputeMeshQualityStatistics(CConfig *config) {

  /*--- Resize our vectors for the 3 metrics: orthogonality, aspect
   ratio, and volume ratio. All are vertex-based for the dual CV. ---*/

  Orthogonality.resize(nPoint,0.0);
  Aspect_Ratio.resize(nPoint,0.0);
  Volume_Ratio.resize(nPoint,0.0);

  /*--- Helper vectors for holding intermediate values. ---*/

  vector<su2double> SurfaceArea(nPoint,0.0);
  vector<su2double> Area_Max(nPoint,0.0);
  vector<su2double> Area_Min(nPoint,1.e6);
  vector<su2double> SubVolume_Max(nPoint,0.0);
  vector<su2double> SubVolume_Min(nPoint,1.e6);

  /*--- Orthogonality and aspect ratio (areas) are computed by
   looping over all edges to check the angles and the face areas. ---*/

  for (unsigned long iEdge = 0; iEdge < nEdge; iEdge++) {

    /*--- Point identification, edge normal vector and area ---*/

    const unsigned long iPoint = edges->GetNode(iEdge,0);
    const unsigned long jPoint = edges->GetNode(iEdge,1);

    const unsigned long GlobalIndex_i = node[iPoint]->GetGlobalIndex();
    const unsigned long GlobalIndex_j = node[iPoint]->GetGlobalIndex();

    /*-- Area normal for the current edge. Recall that this normal
     is computed by summing the normals of adjacent faces along
     the edge between iPoint & jPoint. ---*/

    const su2double *Normal = edges->GetNormal(iEdge);

    /*--- Get the coordinates for point i & j. ---*/

    const su2double *Coord_i = node[iPoint]->GetCoord();
    const su2double *Coord_j = node[jPoint]->GetCoord();

    /*--- Compute the vector pointing from iPoint to jPoint and
     its distance. We also compute face area (norm of the normal vector). ---*/

    su2double distance = 0.0;
    su2double area     = 0.0;
    vector<su2double> edgeVector(nDim);
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      edgeVector[iDim] = Coord_j[iDim]-Coord_i[iDim];
      distance        += edgeVector[iDim]*edgeVector[iDim];
      area            += Normal[iDim]*Normal[iDim];
    }
    distance = sqrt(distance);
    area     = sqrt(area);

    /*--- Aspect ratio is the ratio between the largest and smallest
     faces making up the boundary of the dual CV and is a measure
     of the aspect ratio of the dual control volume. Smaller
     is better (closer to isotropic). ----*/

    if (node[iPoint]->GetDomain()) {
      Area_Min[iPoint] = min(Area_Min[iPoint], area);
      Area_Max[iPoint] = max(Area_Max[iPoint], area);
    }

    if (node[jPoint]->GetDomain()) {
      Area_Min[jPoint] = min(Area_Min[jPoint], area);
      Area_Max[jPoint] = max(Area_Max[jPoint], area);
    }

    if (area <= 0.0) {
      char buf[200];
      SPRINTF(buf, "Zero-area CV face found for edge (%lu,%lu).",
              GlobalIndex_i, GlobalIndex_j);
      SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
    }

    /*--- Compute the angle between the unit normal associated
     with the edge and the unit vector pointing from iPoint to jPoint. ---*/

    su2double dotProduct = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      dotProduct += (Normal[iDim]/area)*(edgeVector[iDim]/distance);
    }

    /*--- The definition of orthogonality is an area-weighted average of
     90 degrees minus the angle between the face area unit normal and
     the vector between i & j. If the two are perfectly aligned, then
     the orthogonality is the desired max of 90 degrees. If they are
     not aligned, the orthogonality will reduce from there. Good values
     are close to 90 degress, poor values are typically below 20 degress. ---*/

    if (node[iPoint]->GetDomain()) {
      Orthogonality[iPoint] += area*(90.0 - acos(dotProduct)*180.0/PI_NUMBER);
      SurfaceArea[iPoint]   += area;
    }
    if (node[jPoint]->GetDomain()) {
      Orthogonality[jPoint] += area*(90.0 - acos(dotProduct)*180.0/PI_NUMBER);
      SurfaceArea[jPoint]   += area;
    }

    /*--- Error check for zero volume of the dual CVs. ---*/

    if (node[iPoint]->GetVolume() <= 0.0) {
      char buf[200];
      SPRINTF(buf, "Zero-area CV face found for point %lu.", GlobalIndex_i);
      SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
    }

    if (node[jPoint]->GetVolume() <= 0.0) {
      char buf[200];
      SPRINTF(buf, "Zero-area CV face found for point %lu.", GlobalIndex_j);
      SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
    }

  }

  /*--- Loop boundary edges to include the area of the boundary elements.  ---*/

  for (unsigned short iMarker = 0; iMarker < nMarker; iMarker++) {
    if ((config->GetMarker_All_KindBC(iMarker) != INTERNAL_BOUNDARY) &&
        (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE)){

      for (unsigned long iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {
        const unsigned long iPoint = vertex[iMarker][iVertex]->GetNode();
        const su2double *Normal    = vertex[iMarker][iVertex]->GetNormal();

        if (node[iPoint]->GetDomain()) {

          /*--- Face area (norm of the normal vector) ---*/

          su2double area = 0.0;
          for (unsigned short iDim = 0; iDim < nDim; iDim++)
            area += Normal[iDim]*Normal[iDim];
          area = sqrt(area);

          /*--- Check to store the area as the min or max for i or j. ---*/

          Area_Min[iPoint] = min(Area_Min[iPoint], area);
          Area_Max[iPoint] = max(Area_Max[iPoint], area);

        }
      }
    }
  }

  /*--- Volume ratio is computed by looping over all volume elements and
   computing the sub-element volume contributions. The ratio between the
   largest and smallest sub-elements making up the dual CV is a
   measure of the volume stretching ratio for the cell. Smaller
   is better (closer to isotropic). ----*/

  unsigned long face_iPoint = 0, face_jPoint = 0;
  unsigned short nEdgesFace = 1;

  vector<su2double> Coord_Edge_CG(nDim);
  vector<su2double> Coord_FaceElem_CG(nDim);
  vector<su2double> Coord_Elem_CG(nDim);
  vector<su2double> Coord_FaceiPoint(nDim);
  vector<su2double> Coord_FacejPoint(nDim);

  for (unsigned long iElem = 0; iElem < nElem; iElem++)
    for (unsigned short iFace = 0; iFace < elem[iElem]->GetnFaces(); iFace++) {

      /*--- In 2D all the faces have only one edge ---*/

      if (nDim == 2) nEdgesFace = 1;

      /*--- In 3D the number of edges per face is the same
       as the number of points per face. ---*/

      if (nDim == 3) nEdgesFace = elem[iElem]->GetnNodesFace(iFace);

      /*-- Loop over the edges of a face ---*/

      for (unsigned short iEdgesFace = 0; iEdgesFace < nEdgesFace; iEdgesFace++) {

        /*--- In 2D only one edge (two points) per edge ---*/

        if (nDim == 2) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
          face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,1));
        }

        /*--- In 3D there are several edges in each face ---*/

        if (nDim == 3) {
          face_iPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace));
          if (iEdgesFace != nEdgesFace-1)
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace, iEdgesFace+1));
          else
            face_jPoint = elem[iElem]->GetNode(elem[iElem]->GetFaces(iFace,0));
        }

        /*--- Locate the edge for the two selected points. ---*/

        const unsigned long iEdge = FindEdge(face_iPoint, face_jPoint);

        /*--- Collect the CG and coordinates for this sub-element face. ---*/

        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          Coord_Edge_CG[iDim]     = edges->GetCG(iEdge,iDim);
          Coord_Elem_CG[iDim]     = elem[iElem]->GetCG(iDim);
          Coord_FaceElem_CG[iDim] = elem[iElem]->GetFaceCG(iFace, iDim);
          Coord_FaceiPoint[iDim]  = node[face_iPoint]->GetCoord(iDim);
          Coord_FacejPoint[iDim]  = node[face_jPoint]->GetCoord(iDim);
        }

        /*--- Access the sub-volume of the element separately in 2D or 3D. ---*/

        su2double Volume_i, Volume_j;
        switch (nDim) {
          case 2:
            Volume_i = CEdge::GetVolume(Coord_FaceiPoint.data(),
                                        Coord_Edge_CG.data(),
                                        Coord_Elem_CG.data());

            Volume_j = CEdge::GetVolume(Coord_FacejPoint.data(),
                                        Coord_Edge_CG.data(),
                                        Coord_Elem_CG.data());
            break;
          case 3:
            Volume_i = CEdge::GetVolume(Coord_FaceiPoint.data(),
                                        Coord_Edge_CG.data(),
                                        Coord_FaceElem_CG.data(),
                                        Coord_Elem_CG.data());

            Volume_j = CEdge::GetVolume(Coord_FacejPoint.data(),
                                        Coord_Edge_CG.data(),
                                        Coord_FaceElem_CG.data(),
                                        Coord_Elem_CG.data());
            break;
        }

        /*--- Check if sub-elem volume is the min or max for iPoint. ---*/

        if (node[face_iPoint]->GetDomain()) {
          SubVolume_Min[face_iPoint] = min(SubVolume_Min[face_iPoint], Volume_i);
          SubVolume_Max[face_iPoint] = max(SubVolume_Max[face_iPoint], Volume_i);
        }

        /*--- Check if sub-elem volume is the min or max for jPoint. ---*/

        if (node[face_jPoint]->GetDomain()) {
          SubVolume_Min[face_jPoint] = min(SubVolume_Min[face_jPoint], Volume_j);
          SubVolume_Max[face_jPoint] = max(SubVolume_Max[face_jPoint], Volume_j);
        }

      }
    }

  /*--- Compute the metrics with a final loop over the vertices. Also
   compute the local min and max values here for reporting. ---*/

  su2double orthoMin = 1.e6, arMin = 1.e6, vrMin = 1.e6;
  su2double orthoMax = 0.0,  arMax = 0.0,  vrMax = 0.0;
  for (unsigned long iPoint= 0; iPoint < nPointDomain; iPoint++) {
    Orthogonality[iPoint] = Orthogonality[iPoint]/SurfaceArea[iPoint];
    orthoMin = min(Orthogonality[iPoint], orthoMin);
    orthoMax = max(Orthogonality[iPoint], orthoMax);

    Aspect_Ratio[iPoint] = Area_Max[iPoint]/Area_Min[iPoint];
    arMin = min(Aspect_Ratio[iPoint], arMin);
    arMax = max(Aspect_Ratio[iPoint], arMax);

    Volume_Ratio[iPoint] = SubVolume_Max[iPoint]/SubVolume_Min[iPoint];
    vrMin = min(Volume_Ratio[iPoint], vrMin);
    vrMax = max(Volume_Ratio[iPoint], vrMax);
  }

  /*--- Reduction to find the min and max values globally. ---*/

  su2double Global_Ortho_Min, Global_Ortho_Max;
  SU2_MPI::Allreduce(&orthoMin, &Global_Ortho_Min, 1,
                     MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&orthoMax, &Global_Ortho_Max, 1,
                     MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  su2double Global_AR_Min, Global_AR_Max;
  SU2_MPI::Allreduce(&arMin, &Global_AR_Min, 1,
                     MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&arMax, &Global_AR_Max, 1,
                     MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  su2double Global_VR_Min, Global_VR_Max;
  SU2_MPI::Allreduce(&vrMin, &Global_VR_Min, 1,
                     MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  SU2_MPI::Allreduce(&vrMax, &Global_VR_Max, 1,
                     MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  /*--- Print the summary to the console for the user. ---*/

  PrintingToolbox::CTablePrinter MetricsTable(&std::cout);
  MetricsTable.AddColumn("Mesh Quality Metric", 30);
  MetricsTable.AddColumn("Minimum", 15);
  MetricsTable.AddColumn("Maximum", 15);
  if (rank == MASTER_NODE){
    MetricsTable.PrintHeader();
    MetricsTable << "Orthogonality Angle (deg.)" << Global_Ortho_Min << Global_Ortho_Max;
    MetricsTable << "CV Face Area Aspect Ratio" << Global_AR_Min << Global_AR_Max;
    MetricsTable << "CV Sub-Volume Ratio" << Global_VR_Min << Global_VR_Max;
    MetricsTable.PrintFooter();
  }

  /*--- If we will not be writing the stats to the visualization files,
   force clear the memory with the swap() function. ---*/

  if (!config->GetWrt_MeshQuality()) {
    vector<su2double>().swap(Orthogonality);
    vector<su2double>().swap(Aspect_Ratio);
    vector<su2double>().swap(Volume_Ratio);
  }

}

void CPhysicalGeometry::FindNormal_Neighbor(CConfig *config) {
  su2double cos_max, scalar_prod, norm_vect, norm_Normal, cos_alpha, diff_coord, *Normal;
  unsigned long Point_Normal, jPoint;
  unsigned short iNeigh, iMarker, iDim;
  unsigned long iPoint, iVertex;

  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {

    if (config->GetMarker_All_KindBC(iMarker) != SEND_RECEIVE &&
        config->GetMarker_All_KindBC(iMarker) != INTERFACE_BOUNDARY &&
        config->GetMarker_All_KindBC(iMarker) != NEARFIELD_BOUNDARY ) {

      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

        iPoint = vertex[iMarker][iVertex]->GetNode();
        Normal = vertex[iMarker][iVertex]->GetNormal();

        /*--- Compute closest normal neighbor, note that the normal are oriented inwards ---*/
        Point_Normal = 0; cos_max = -1.0;
        for (iNeigh = 0; iNeigh < node[iPoint]->GetnPoint(); iNeigh++) {
          jPoint = node[iPoint]->GetPoint(iNeigh);
          scalar_prod = 0.0; norm_vect = 0.0; norm_Normal = 0.0;
          for (iDim = 0; iDim < nDim; iDim++) {
            diff_coord = node[jPoint]->GetCoord(iDim)-node[iPoint]->GetCoord(iDim);
            scalar_prod += diff_coord*Normal[iDim];
            norm_vect += diff_coord*diff_coord;
            norm_Normal += Normal[iDim]*Normal[iDim];
          }
          norm_vect = sqrt(norm_vect);
          norm_Normal = sqrt(norm_Normal);
          cos_alpha = scalar_prod/(norm_vect*norm_Normal);

          /*--- Get maximum cosine ---*/
          if (cos_alpha >= cos_max) {
            Point_Normal = jPoint;
            cos_max = cos_alpha;
          }
        }
        vertex[iMarker][iVertex]->SetNormal_Neighbor(Point_Normal);
      }
    }
  }
}

void CPhysicalGeometry::SetBoundSensitivity(CConfig *config) {
  unsigned short iMarker, icommas;
  unsigned long iVertex, iPoint, (*Point2Vertex)[2], nPointLocal = 0, nPointGlobal = 0;
  su2double Sensitivity;
  bool *PointInDomain;

  nPointLocal = nPoint;
  SU2_MPI::Allreduce(&nPointLocal, &nPointGlobal, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

  Point2Vertex = new unsigned long[nPointGlobal][2];
  PointInDomain = new bool[nPointGlobal];

  for (iPoint = 0; iPoint < nPointGlobal; iPoint ++)
    PointInDomain[iPoint] = false;

  for (iMarker = 0; iMarker < nMarker; iMarker++)
    if (config->GetMarker_All_DV(iMarker) == YES)
      for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++) {

        /*--- The sensitivity file uses the global numbering ---*/
        iPoint = node[vertex[iMarker][iVertex]->GetNode()]->GetGlobalIndex();

        if (vertex[iMarker][iVertex]->GetNode() < GetnPointDomain()) {
          Point2Vertex[iPoint][0] = iMarker;
          Point2Vertex[iPoint][1] = iVertex;
          PointInDomain[iPoint] = true;
          vertex[iMarker][iVertex]->SetAuxVar(0.0);
        }
      }

  /*--- Time-average any unsteady surface sensitivities ---*/

  unsigned long iTimeIter, nTimeIter;
  su2double delta_T, total_T;
  if (config->GetTime_Marching() && config->GetTime_Domain()) {
    nTimeIter = config->GetUnst_AdjointIter();
    delta_T  = config->GetTime_Step();
    total_T  = (su2double)nTimeIter*delta_T;
  } else if (config->GetTime_Marching() == HARMONIC_BALANCE) {

    /*--- Compute period of oscillation & compute time interval using nTimeInstances ---*/

    su2double period = config->GetHarmonicBalance_Period();
    nTimeIter  = config->GetnTimeInstances();
    delta_T   = period/(su2double)nTimeIter;
    total_T   = period;

  } else {
    nTimeIter = 1;
    delta_T  = 1.0;
    total_T  = 1.0;
  }

  for (iTimeIter = 0; iTimeIter < nTimeIter; iTimeIter++) {

    /*--- Prepare to read surface sensitivity files (CSV) ---*/

    string text_line;
    ifstream Surface_file;
    char buffer[50];
    char cstr[MAX_STRING_SIZE];
    string surfadj_filename = config->GetSurfAdjCoeff_FileName();
    strcpy (cstr, surfadj_filename.c_str());

    /*--- Write file name with extension if unsteady or steady ---*/
    if (config->GetTime_Marching() == HARMONIC_BALANCE)
    	SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iTimeIter));

    if ((config->GetTime_Marching() && config->GetTime_Domain()) ||
        (config->GetTime_Marching() == HARMONIC_BALANCE)) {
      if ((SU2_TYPE::Int(iTimeIter) >= 0)    && (SU2_TYPE::Int(iTimeIter) < 10))    SPRINTF (buffer, "_0000%d.csv", SU2_TYPE::Int(iTimeIter));
      if ((SU2_TYPE::Int(iTimeIter) >= 10)   && (SU2_TYPE::Int(iTimeIter) < 100))   SPRINTF (buffer, "_000%d.csv",  SU2_TYPE::Int(iTimeIter));
      if ((SU2_TYPE::Int(iTimeIter) >= 100)  && (SU2_TYPE::Int(iTimeIter) < 1000))  SPRINTF (buffer, "_00%d.csv",   SU2_TYPE::Int(iTimeIter));
      if ((SU2_TYPE::Int(iTimeIter) >= 1000) && (SU2_TYPE::Int(iTimeIter) < 10000)) SPRINTF (buffer, "_0%d.csv",    SU2_TYPE::Int(iTimeIter));
      if (SU2_TYPE::Int(iTimeIter) >= 10000) SPRINTF (buffer, "_%d.csv", SU2_TYPE::Int(iTimeIter));
    }
    else
      SPRINTF (buffer, ".csv");

    strcat (cstr, buffer);

    /*--- Read the sensitivity file ---*/

    string::size_type position;

    Surface_file.open(cstr, ios::in);

    /*--- File header ---*/

    getline(Surface_file, text_line);

    vector<string> split_line;

    char delimiter = ',';
    split_line = PrintingToolbox::split(text_line, delimiter);

    for (unsigned short iField = 0; iField < split_line.size(); iField++){
      PrintingToolbox::trim(split_line[iField]);
    }

    std::vector<string>::iterator it = std::find(split_line.begin(), split_line.end(), "\"Surface_Sensitivity\"");

    if (it == split_line.end()){
      SU2_MPI::Error("Surface sensitivity not found in file.", CURRENT_FUNCTION);
    }

    unsigned short sens_index = std::distance(split_line.begin(), it);

    while (getline(Surface_file, text_line)) {
      for (icommas = 0; icommas < 50; icommas++) {
        position = text_line.find( ",", 0 );
        if (position!=string::npos) text_line.erase (position,1);
      }
      stringstream  point_line(text_line);
      point_line >> iPoint;

      for (int i = 1; i <= sens_index; i++)
        point_line >> Sensitivity;

      if (PointInDomain[iPoint]) {

        /*--- Find the vertex for the Point and Marker ---*/

        iMarker = Point2Vertex[iPoint][0];
        iVertex = Point2Vertex[iPoint][1];

        /*--- Increment the auxiliary variable with the contribution of
         this unsteady timestep. For steady problems, this reduces to
         a single sensitivity value multiplied by 1.0. ---*/

        vertex[iMarker][iVertex]->AddAuxVar(Sensitivity*(delta_T/total_T));
      }

    }
    Surface_file.close();
  }

  delete[] Point2Vertex;
  delete[] PointInDomain;

}

void CPhysicalGeometry::SetSensitivity(CConfig *config) {

  ifstream restart_file;
  string filename = config->GetSolution_AdjFileName();

  su2double AoASens;
  unsigned short nTimeIter;
  unsigned long index;
  string::size_type position;
  int counter = 0;

  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  if (config->GetTime_Domain()) {
    nTimeIter = config->GetnTime_Iter();
  }else {
    nTimeIter = 1;
  }

  if (rank == MASTER_NODE)
    cout << "Reading in sensitivity at iteration " << nTimeIter-1 << "."<< endl;

  /*--- Read all lines in the restart file ---*/
  long iPoint_Local; unsigned long iPoint_Global = 0; string text_line;

  iPoint_Global = 0;

  filename = config->GetSolution_AdjFileName();

  filename = config->GetObjFunc_Extension(filename);


  if (config->GetRead_Binary_Restart()) {

    filename = config->GetFilename(filename, ".dat", nTimeIter-1);

    char str_buf[CGNS_STRING_SIZE], fname[100];
    unsigned short iVar;
    strcpy(fname, filename.c_str());
    int nRestart_Vars = 5, nFields;
    int *Restart_Vars = new int[5];
    passivedouble *Restart_Data = NULL;
    int Restart_Iter = 0;
    passivedouble Restart_Meta_Passive[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    su2double Restart_Meta[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + fname, CURRENT_FUNCTION);
    }

    /*--- First, read the number of variables and points. ---*/

    ret = fread(Restart_Vars, sizeof(int), nRestart_Vars, fhw);
    if (ret != (unsigned long)nRestart_Vars) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (Restart_Vars[0] != 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }

    /*--- Store the number of fields for simplicity. ---*/

    nFields = Restart_Vars[1];

    /*--- Read the variable names from the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. We pad the beginning of the
     variable string vector with the Point_ID tag that wasn't written. ---*/

    config->fields.push_back("Point_ID");
    for (iVar = 0; iVar < nFields; iVar++) {
      ret = fread(str_buf, sizeof(char), CGNS_STRING_SIZE, fhw);
      if (ret != (unsigned long)CGNS_STRING_SIZE) {
        SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
      }
      config->fields.push_back(str_buf);
    }

    /*--- For now, create a temp 1D buffer to read the data from file. ---*/

    Restart_Data = new passivedouble[nFields*GetnPointDomain()];

    /*--- Read in the data for the restart at all local points. ---*/

    ret = fread(Restart_Data, sizeof(passivedouble), nFields*GetnPointDomain(), fhw);
    if (ret != (unsigned long)nFields*GetnPointDomain()) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Compute (negative) displacements and grab the metadata. ---*/

    ret = sizeof(int) + 8*sizeof(passivedouble);
    fseek(fhw,-ret, SEEK_END);

    /*--- Read the external iteration. ---*/

    ret = fread(&Restart_Iter, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Read the metadata. ---*/

    ret = fread(Restart_Meta_Passive, sizeof(passivedouble), 8, fhw);
    if (ret != 8) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Close the file. ---*/

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    SU2_MPI::Status status;
    MPI_Datatype etype, filetype;
    MPI_Offset disp;
    unsigned long iPoint_Global, iChar;
    string field_buf;

    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- First, read the number of variables and points (i.e., cols and rows),
     which we will need in order to read the file later. Also, read the
     variable string names here. Only the master rank reads the header. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, Restart_Vars, nRestart_Vars, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(Restart_Vars, nRestart_Vars, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (Restart_Vars[0] != 535532) {

      SU2_MPI::Error(string("File ") + string(fname) + string(" is not a binary SU2 restart file.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }

    /*--- Store the number of fields for simplicity. ---*/

    nFields = Restart_Vars[1];

    /*--- Read the variable names from the file. Note that we are adopting a
     fixed length of 33 for the string length to match with CGNS. This is
     needed for when we read the strings later. ---*/

    char *mpi_str_buf = new char[nFields*CGNS_STRING_SIZE];
    if (rank == MASTER_NODE) {
      disp = nRestart_Vars*sizeof(int);
      MPI_File_read_at(fhw, disp, mpi_str_buf, nFields*CGNS_STRING_SIZE,
                       MPI_CHAR, MPI_STATUS_IGNORE);
    }

    /*--- Broadcast the string names of the variables. ---*/

    SU2_MPI::Bcast(mpi_str_buf, nFields*CGNS_STRING_SIZE, MPI_CHAR,
                   MASTER_NODE, MPI_COMM_WORLD);

    /*--- Now parse the string names and load into the config class in case
     we need them for writing visualization files (SU2_SOL). ---*/

    config->fields.push_back("Point_ID");
    for (iVar = 0; iVar < nFields; iVar++) {
      index = iVar*CGNS_STRING_SIZE;
      for (iChar = 0; iChar < (unsigned long)CGNS_STRING_SIZE; iChar++) {
        str_buf[iChar] = mpi_str_buf[index + iChar];
      }
      field_buf.append(str_buf);
      config->fields.push_back(field_buf.c_str());
      field_buf.clear();
    }

    /*--- Free string buffer memory. ---*/

    delete [] mpi_str_buf;

    /*--- We're writing only su2doubles in the data portion of the file. ---*/

    etype = MPI_DOUBLE;

    /*--- We need to ignore the 4 ints describing the nVar_Restart and nPoints,
     along with the string names of the variables. ---*/

    disp = nRestart_Vars*sizeof(int) + CGNS_STRING_SIZE*nFields*sizeof(char);

    /*--- Define a derived datatype for this rank's set of non-contiguous data
     that will be placed in the restart. Here, we are collecting each one of the
     points which are distributed throughout the file in blocks of nVar_Restart data. ---*/

    int *blocklen = new int[GetnPointDomain()];
    int *displace = new int[GetnPointDomain()];

    counter = 0;
    for (iPoint_Global = 0; iPoint_Global < GetGlobal_nPointDomain(); iPoint_Global++ ) {
      if (GetGlobal_to_Local_Point(iPoint_Global) > -1) {
        blocklen[counter] = nFields;
        displace[counter] = iPoint_Global*nFields;
        counter++;
      }
    }
    MPI_Type_indexed(GetnPointDomain(), blocklen, displace, MPI_DOUBLE, &filetype);
    MPI_Type_commit(&filetype);

    /*--- Set the view for the MPI file write, i.e., describe the location in
     the file that this rank "sees" for writing its piece of the restart file. ---*/

    MPI_File_set_view(fhw, disp, etype, filetype, (char*)"native", MPI_INFO_NULL);

    /*--- For now, create a temp 1D buffer to read the data from file. ---*/

    Restart_Data = new passivedouble[nFields*GetnPointDomain()];

    /*--- Collective call for all ranks to read from their view simultaneously. ---*/

    MPI_File_read_all(fhw, Restart_Data, nFields*GetnPointDomain(), MPI_DOUBLE, &status);

    /*--- Free the derived datatype. ---*/

    MPI_Type_free(&filetype);

    /*--- Reset the file view before writing the metadata. ---*/

    MPI_File_set_view(fhw, 0, MPI_BYTE, MPI_BYTE, (char*)"native", MPI_INFO_NULL);

    /*--- Access the metadata. ---*/

    if (rank == MASTER_NODE) {

      /*--- External iteration. ---*/
      disp = (nRestart_Vars*sizeof(int) + nFields*CGNS_STRING_SIZE*sizeof(char) +
              nFields*Restart_Vars[2]*sizeof(passivedouble));
      MPI_File_read_at(fhw, disp, &Restart_Iter, 1, MPI_INT, MPI_STATUS_IGNORE);

      /*--- Additional doubles for AoA, AoS, etc. ---*/

      disp = (nRestart_Vars*sizeof(int) + nFields*CGNS_STRING_SIZE*sizeof(char) +
              nFields*Restart_Vars[2]*sizeof(passivedouble) + 1*sizeof(int));
      MPI_File_read_at(fhw, disp, Restart_Meta_Passive, 8, MPI_DOUBLE, MPI_STATUS_IGNORE);

    }

    /*--- Communicate metadata. ---*/

    SU2_MPI::Bcast(&Restart_Iter, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Copy to a su2double structure (because of the SU2_MPI::Bcast
              doesn't work with passive data)---*/

    for (unsigned short iVar = 0; iVar < 8; iVar++)
      Restart_Meta[iVar] = Restart_Meta_Passive[iVar];

    SU2_MPI::Bcast(Restart_Meta, 8, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);

    /*--- All ranks close the file after writing. ---*/

    MPI_File_close(&fhw);

    delete [] blocklen;
    delete [] displace;

#endif

    std::vector<string>::iterator itx = std::find(config->fields.begin(), config->fields.end(), "Sensitivity_x");
    std::vector<string>::iterator ity = std::find(config->fields.begin(), config->fields.end(), "Sensitivity_y");
    std::vector<string>::iterator itz = std::find(config->fields.begin(), config->fields.end(), "Sensitivity_z");

    if (itx == config->fields.end()){
      SU2_MPI::Error("Sensitivity x not found in file.", CURRENT_FUNCTION);
    }
    if (ity == config->fields.end()){
      SU2_MPI::Error("Sensitivity y not found in file.", CURRENT_FUNCTION);
    }
    if (nDim == 3){
      if (itz == config->fields.end()){
        SU2_MPI::Error("Sensitivity z not found in file.", CURRENT_FUNCTION);
      }
    }

    unsigned short sens_x_idx = std::distance(config->fields.begin(), itx);
    unsigned short sens_y_idx = std::distance(config->fields.begin(), ity);
    unsigned short sens_z_idx = 0;
    if (nDim == 3)
      sens_z_idx = std::distance(config->fields.begin(), itz);

    /*--- Load the data from the binary restart. ---*/

    counter = 0;
    for (iPoint_Global = 0; iPoint_Global < GetGlobal_nPointDomain(); iPoint_Global++ ) {

      /*--- Retrieve local index. If this node from the restart file lives
       on the current processor, we will load and instantiate the vars. ---*/

      iPoint_Local = GetGlobal_to_Local_Point(iPoint_Global);

      if (iPoint_Local > -1) {

        /*--- We need to store this point's data, so jump to the correct
         offset in the buffer of data from the restart file and load it. ---*/

        index = counter*nFields + sens_x_idx - 1;
        Sensitivity(iPoint_Local,0) = Restart_Data[index];
        index = counter*nFields + sens_y_idx - 1;
        Sensitivity(iPoint_Local,1) = Restart_Data[index];

        if (nDim == 3){
          index = counter*nFields + sens_z_idx - 1;
          Sensitivity(iPoint_Local,2) = Restart_Data[index];
        }
        /*--- Increment the overall counter for how many points have been loaded. ---*/
        counter++;
      }
    }

    /*--- Lastly, load the AoA sensitivity from the binary metadata. ---*/

    config->SetAoA_Sens(Restart_Meta[4]);

  } else {

    filename = config->GetFilename(filename, ".csv", nTimeIter-1);

    /*--- First, check that this is not a binary restart file. ---*/

    char fname[100];
    strcpy(fname, filename.c_str());
    int magic_number;

#ifndef HAVE_MPI

    /*--- Serial binary input. ---*/

    FILE *fhw;
    fhw = fopen(fname,"rb");
    size_t ret;

    /*--- Error check for opening the file. ---*/

    if (!fhw) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- Attempt to read the first int, which should be our magic number. ---*/

    ret = fread(&magic_number, sizeof(int), 1, fhw);
    if (ret != 1) {
      SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
    }

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {
      SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }

    fclose(fhw);

#else

    /*--- Parallel binary input using MPI I/O. ---*/

    MPI_File fhw;
    int ierr;

    /*--- All ranks open the file using MPI. ---*/

    ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

    /*--- Error check opening the file. ---*/

    if (ierr) {
      SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
    }

    /*--- Have the master attempt to read the magic number. ---*/

    if (rank == MASTER_NODE)
      MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

    /*--- Broadcast the number of variables to all procs and store clearly. ---*/

    SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    /*--- Check that this is an SU2 binary file. SU2 binary files
     have the hex representation of "SU2" as the first int in the file. ---*/

    if (magic_number == 535532) {

      SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
                     string("SU2 reads/writes binary restart files by default.\n") +
                     string("Note that backward compatibility for ASCII restart files is\n") +
                     string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
    }

    MPI_File_close(&fhw);

#endif

  restart_file.open(filename.data(), ios::in);
  if (restart_file.fail()) {
    SU2_MPI::Error(string("There is no adjoint restart file ") + filename, CURRENT_FUNCTION);
  }

  /*--- The first line is the header ---*/

  getline (restart_file, text_line);

  vector<string> fields = PrintingToolbox::split(text_line, ',');

  for (unsigned short iField = 0; iField < fields.size(); iField++){
    PrintingToolbox::trim(fields[iField]);
  }

  std::vector<string>::iterator itx = std::find(fields.begin(), fields.end(), "\"Sensitivity_x\"");
  std::vector<string>::iterator ity = std::find(fields.begin(), fields.end(), "\"Sensitivity_y\"");
  std::vector<string>::iterator itz = std::find(fields.begin(), fields.end(), "\"Sensitivity_z\"");

  if (itx == fields.end()){
    SU2_MPI::Error("Sensitivity x not found in file.", CURRENT_FUNCTION);
  }
  if (ity ==fields.end()){
    SU2_MPI::Error("Sensitivity y not found in file.", CURRENT_FUNCTION);
  }
  if (nDim == 3){
    if (itz == fields.end()){
      SU2_MPI::Error("Sensitivity z not found in file.", CURRENT_FUNCTION);
    }
  }

  unsigned short sens_x_idx = std::distance(fields.begin(), itx);
  unsigned short sens_y_idx = std::distance(fields.begin(), ity);
  unsigned short sens_z_idx = 0;
  if (nDim == 3)
    sens_z_idx = std::distance(fields.begin(), itz);


  for (iPoint_Global = 0; iPoint_Global < GetGlobal_nPointDomain(); iPoint_Global++ ) {

    getline (restart_file, text_line);

    vector<string> point_line = PrintingToolbox::split(text_line, ',');

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    iPoint_Local = GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {
      Sensitivity(iPoint_Local,0) = PrintingToolbox::stod(point_line[sens_x_idx]);
      Sensitivity(iPoint_Local,1) = PrintingToolbox::stod(point_line[sens_y_idx]);
      if (nDim == 3)
        Sensitivity(iPoint_Local,2) = PrintingToolbox::stod(point_line[sens_z_idx]);
    }

  }

  /*--- Read AoA sensitivity ---*/

  while (getline (restart_file, text_line)) {
    position = text_line.find ("SENS_AOA=",0);
    if (position != string::npos) {
      text_line.erase (0,9); AoASens = atof(text_line.c_str());
      config->SetAoA_Sens(AoASens);
    }
  }

  restart_file.close();

  }

}

void CPhysicalGeometry::ReadUnorderedSensitivity(CConfig *config) {

  /*--- This routine makes SU2_DOT more interoperable with other
   packages so that folks can customize their workflows. For example, one
   may want to compute flow and adjoint with package A, deform the mesh
   and project the sensitivities with SU2, and control the actual shape
   parameterization with package C. This routine allows SU2_DOT to read
   in an additional format for volume sensitivities that looks like:

    x0, y0, z0, dj/dx, dj/dy, dj/dz
    x1, y1, z1, dj/dx, dj/dy, dj/dz
    ...
    xN, yN, zN, dj/dx, dj/dy, dj/dz

   with N being the number of grid points. This is a format already used
   in other packages. Note that the nodes can be in any order in the file. ---*/

  unsigned short iDim;
  unsigned long iPoint, pointID;
  unsigned long unmatched = 0, iPoint_Found = 0, iPoint_Ext = 0;

  su2double Coor_External[3] = {0.0,0.0,0.0}, Sens_External[3] = {0.0,0.0,0.0};
  su2double dist;
  int rankID;

  string filename, text_line;
  ifstream external_file;
  ofstream sens_file;

  if (rank == MASTER_NODE)
    cout << "Parsing unordered ASCII volume sensitivity file."<< endl;

  /*--- Allocate space for the sensitivity and initialize. ---*/

  Sensitivity.resize(nPoint,nDim) = su2double(0.0);

  /*--- Get the filename for the unordered ASCII sensitivity file input. ---*/

  filename = config->GetDV_Unordered_Sens_Filename();
  external_file.open(filename.data(), ios::in);
  if (external_file.fail()) {
    SU2_MPI::Error(string("There is no unordered ASCII sensitivity file ") +
                   filename, CURRENT_FUNCTION);
  }

  /*--- Allocate the vectors to hold boundary node coordinates
   and its local ID. ---*/

  vector<su2double>     Coords(nDim*nPointDomain);
  vector<unsigned long> PointIDs(nPointDomain);

  /*--- Retrieve and store the coordinates of owned interior nodes
   and their local point IDs. ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
    PointIDs[iPoint] = iPoint;
    for (iDim = 0; iDim < nDim; iDim++)
      Coords[iPoint*nDim + iDim] = node[iPoint]->GetCoord(iDim);
  }

  /*--- Build the ADT of all interior nodes. ---*/

  CADTPointsOnlyClass VertexADT(nDim, nPointDomain,
                                Coords.data(), PointIDs.data(), true);

  /*--- Loop over all interior mesh nodes owned by this rank and find the
   matching point with minimum distance. Once we have the match, store the
   sensitivities from the file for that node. ---*/

  if (VertexADT.IsEmpty()) {

    SU2_MPI::Error("No external points given to ADT.", CURRENT_FUNCTION);

  } else {

    /*--- Read the input sensitivity file and locate the point matches
     using the ADT search, on a line-by-line basis. ---*/

    iPoint_Found = 0; iPoint_Ext  = 0;
    while (getline (external_file, text_line)) {

      /*--- First, check that the line has 6 entries, otherwise throw out. ---*/

      istringstream point_line(text_line);
      vector<string> tokens((istream_iterator<string>(point_line)),
                             istream_iterator<string>());

      if (tokens.size() == 6) {

        istringstream point_line(text_line);

        /*--- Get the coordinates and sensitivity for this line. ---*/

        for (iDim = 0; iDim < nDim; iDim++) point_line >> Coor_External[iDim];
        for (iDim = 0; iDim < nDim; iDim++) point_line >> Sens_External[iDim];

        /*--- Locate the nearest node to this external point. If it is on
         our rank, then store the sensitivity value. ---*/

        VertexADT.DetermineNearestNode(&Coor_External[0], dist,
                                       pointID, rankID);

        if (rankID == rank) {

          /*--- Store the sensitivities at the matched local node. ---*/

          for (iDim = 0; iDim < nDim; iDim++)
            Sensitivity(pointID,iDim) = Sens_External[iDim];

          /*--- Keep track of how many points we match. ---*/

          iPoint_Found++;

          /*--- Keep track of points with poor matches for reporting. ---*/

          if (dist > 1e-10) unmatched++;

        }

        /*--- Increment counter for total points in the external file. ---*/

        iPoint_Ext++;

      }
    }

    /*--- Close the external file. ---*/

    external_file.close();

    /*--- We have not received all nodes in the input file. Throw an error. ---*/

    if ((iPoint_Ext < GetGlobal_nPointDomain()) && (rank == MASTER_NODE)) {
      sens_file.open(config->GetDV_Unordered_Sens_Filename().data(), ios::out);
      sens_file.close();
      SU2_MPI::Error("Not enough points in the input sensitivity file.",
                     CURRENT_FUNCTION);
    }

    /*--- Check for points with a poor match and report the count. ---*/

    unsigned long myUnmatched = unmatched; unmatched = 0;
    SU2_MPI::Allreduce(&myUnmatched, &unmatched, 1,
                       MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    if ((unmatched > 0) && (rank == MASTER_NODE)) {
      cout << " Warning: there are " << unmatched;
      cout << " points with a match distance > 1e-10." << endl;
    }

  }

}

void CPhysicalGeometry::Check_Periodicity(CConfig *config) {

  /*--- Check for the presence of any periodic BCs and disable multigrid
   for now if found. ---*/

  if ((config->GetnMarker_Periodic() != 0) && (config->GetnMGLevels() > 0)) {
    if (rank == MASTER_NODE)
      cout << "WARNING: Periodicity has been detected. Disabling multigrid. "<< endl;
    config->SetMGLevels(0);
  }

}

su2double CPhysicalGeometry::Compute_MaxThickness(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config, vector<su2double> &Xcoord_Airfoil,
                                                  vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, jVertex, n, Trailing_Point, Leading_Point;
  su2double Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, zp1, zpn, MaxThickness_Value = 0, Thickness, Length, Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, MaxDistance, Distance, AoA;
  vector<su2double> Xcoord, Ycoord, Zcoord, Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;

  /*--- Translate to the origin ---*/

  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];

  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }

  /*--- Rotate the airfoil ---*/

  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);

  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];
    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
  }

  /*--- Identify upper and lower side, and store the value of the normal --*/

  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));

    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;

    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;

    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];

    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);

    unsigned short index = 2;

    /*--- Removing the trailing edge from list of points that we are going to use in the interpolation,
			to be sure that a blunt trailing edge do not affect the interpolation ---*/

    if ((Normal[index] >= 0.0) && (fabs(Xcoord_Airfoil_[iVertex]) > MaxDistance*0.01)) {
      Xcoord.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord.push_back(Zcoord_Airfoil_[iVertex]);
    }

  }

  /*--- Order the arrays using the X component ---*/

  for (iVertex = 0; iVertex < Xcoord.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord.size() - 1 - iVertex; jVertex++) {
      if (Xcoord[jVertex] > Xcoord[jVertex+1]) {
        auxXCoord = Xcoord[jVertex]; Xcoord[jVertex] = Xcoord[jVertex+1]; Xcoord[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord[jVertex]; Ycoord[jVertex] = Ycoord[jVertex+1]; Ycoord[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord[jVertex]; Zcoord[jVertex] = Zcoord[jVertex+1]; Zcoord[jVertex+1] = auxZCoord;
      }
    }
  }

  n = Xcoord.size();
  if (n > 1) {
    zp1 = (Zcoord[1]-Zcoord[0])/(Xcoord[1]-Xcoord[0]);
    zpn = (Zcoord[n-1]-Zcoord[n-2])/(Xcoord[n-1]-Xcoord[n-2]);
    Z2coord.resize(n+1);
    SetSpline(Xcoord, Zcoord, n, zp1, zpn, Z2coord);

    /*--- Compute the thickness (we add a fabs because we can not guarantee the
     right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/

    MaxThickness_Value = 0.0;
    for (iVertex = 0; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
      if (Zcoord_Normal[iVertex] < 0.0) {
        Thickness = fabs(Zcoord_Airfoil_[iVertex] - GetSpline(Xcoord, Zcoord, Z2coord, n, Xcoord_Airfoil_[iVertex]));
        if (Thickness > MaxThickness_Value) { MaxThickness_Value = Thickness; }
      }
    }
  }
  else { MaxThickness_Value = 0.0; }

  return MaxThickness_Value;

}

su2double CPhysicalGeometry::Compute_Dihedral(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                              su2double *LeadingEdge_i, su2double *TrailingEdge_i) {

  // su2double Dihedral_Leading = atan((LeadingEdge_i[2] - LeadingEdge_im1[2]) / (LeadingEdge_i[1] - LeadingEdge_im1[1]))*180/PI_NUMBER;
  su2double Dihedral_Trailing = atan((TrailingEdge_i[2] - TrailingEdge_im1[2]) / (TrailingEdge_i[1] - TrailingEdge_im1[1]))*180/PI_NUMBER;

  // su2double Dihedral = 0.5*(Dihedral_Leading + Dihedral_Trailing);

  return Dihedral_Trailing;

}

su2double CPhysicalGeometry::Compute_Curvature(su2double *LeadingEdge_im1, su2double *TrailingEdge_im1,
                                               su2double *LeadingEdge_i, su2double *TrailingEdge_i,
                                               su2double *LeadingEdge_ip1, su2double *TrailingEdge_ip1) {

  su2double A[2], B[2], C[2], BC[2], AB[2], AC[2], BC_MOD, AB_MOD,  AC_MOD, AB_CROSS_AC;

  // A[0] = LeadingEdge_im1[1];     A[1] = LeadingEdge_im1[2];
  // B[0] = LeadingEdge_i[1];           B[1] = LeadingEdge_i[2];
  // C[0] = LeadingEdge_ip1[1];      C[1] = LeadingEdge_ip1[2];

  // BC[0] = C[0] - B[0]; BC[1] = C[1] - B[1];
  // AB[0] = B[0] - A[0]; AB[1] = B[1] - A[1];
  // AC[0] = C[0] - A[0]; AC[1] = C[1] - A[1];
  // BC_MOD = sqrt(BC[0]*BC[0] + BC[1]*BC[1] );
  // AB_MOD = sqrt(AB[0]*AB[0] + AB[1]*AB[1] );
  // AC_MOD = sqrt(AC[0]*AC[0] + AC[1]*AC[1] );
  // AB_CROSS_AC = AB[0]* AC[1] - AB[1]* AC[0];

  // su2double Curvature_Leading = fabs(1.0/(0.5*BC_MOD*AB_MOD*AC_MOD/AB_CROSS_AC));

  A[0] = TrailingEdge_im1[1];      A[1] = TrailingEdge_im1[2];
  B[0] = TrailingEdge_i[1];                B[1] = TrailingEdge_i[2];
  C[0] = TrailingEdge_ip1[1];      C[1] = TrailingEdge_ip1[2];

  BC[0] = C[0] - B[0]; BC[1] = C[1] - B[1];
  AB[0] = B[0] - A[0]; AB[1] = B[1] - A[1];
  AC[0] = C[0] - A[0]; AC[1] = C[1] - A[1];
  BC_MOD = sqrt(BC[0]*BC[0] + BC[1]*BC[1] );
  AB_MOD = sqrt(AB[0]*AB[0] + AB[1]*AB[1] );
  AC_MOD = sqrt(AC[0]*AC[0] + AC[1]*AC[1] );
  AB_CROSS_AC = AB[0]* AC[1] - AB[1]* AC[0];

  su2double Curvature_Trailing = fabs(1.0/(0.5*BC_MOD*AB_MOD*AC_MOD/AB_CROSS_AC));

  // su2double Curvature = 0.5*(Curvature_Leading + Curvature_Trailing);

  return Curvature_Trailing;

}

su2double CPhysicalGeometry::Compute_Twist(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                           vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance, Twist = 0.0;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  Twist = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;

  return Twist;

}

void CPhysicalGeometry::Compute_Wing_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal,
                                                vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {

    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  LeadingEdge[0] = Xcoord_Airfoil[Leading_Point];
  LeadingEdge[1] = Ycoord_Airfoil[Leading_Point];
  LeadingEdge[2] = Zcoord_Airfoil[Leading_Point];

  TrailingEdge[0] = Xcoord_Airfoil[Trailing_Point];
  TrailingEdge[1] = Ycoord_Airfoil[Trailing_Point];
  TrailingEdge[2] = Zcoord_Airfoil[Trailing_Point];

}

void CPhysicalGeometry::Compute_Fuselage_LeadingTrailing(su2double *LeadingEdge, su2double *TrailingEdge, su2double *Plane_P0, su2double *Plane_Normal,
                                                vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance;

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0));
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  LeadingEdge[0] = Xcoord_Airfoil[Leading_Point];
  LeadingEdge[1] = Ycoord_Airfoil[Leading_Point];
  LeadingEdge[2] = Zcoord_Airfoil[Leading_Point];

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Zcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  TrailingEdge[0] = 0.5*(Xcoord_Airfoil[Trailing_Point]+Xcoord_Airfoil[Leading_Point]);
  TrailingEdge[1] = 0.5*(Ycoord_Airfoil[Trailing_Point]+Ycoord_Airfoil[Leading_Point]);
  TrailingEdge[2] = 0.5*(Zcoord_Airfoil[Trailing_Point]+Zcoord_Airfoil[Leading_Point]);

}

su2double CPhysicalGeometry::Compute_Chord(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                           vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex, Trailing_Point;
  su2double MaxDistance, Distance, Chord = 0.0;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/
  MaxDistance = 0.0; Trailing_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {

    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; }
  }

  Chord = MaxDistance;

  return Chord;

}

su2double CPhysicalGeometry::Compute_Width(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                           vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point;
  su2double MaxDistance, Distance, Width = 0.0;

  MaxDistance = 0.0; Trailing_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = fabs(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point]);
    if (MaxDistance < Distance) { MaxDistance = Distance; }
  }

  Width = MaxDistance;
  return Width;

}

su2double CPhysicalGeometry::Compute_WaterLineWidth(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config,
                                                    vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point;
  su2double MinDistance, Distance, WaterLineWidth = 0.0;
  su2double WaterLine = config->GetGeo_Waterline_Location();

  MinDistance = 1E10; WaterLineWidth = 0; Trailing_Point = 0;
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = fabs(Zcoord_Airfoil[iVertex] - WaterLine);
    if (Distance < MinDistance) {
    	MinDistance = Distance;
    	WaterLineWidth = fabs(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point]);
    }
  }

  return WaterLineWidth;

}

su2double CPhysicalGeometry::Compute_Height(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                            vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point;
  su2double MaxDistance, Distance, Height = 0.0;

  MaxDistance = 0.0; Trailing_Point = 0;
  for (iVertex = 1; iVertex < Zcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));
    if (MaxDistance < Distance) { MaxDistance = Distance; }
  }

  Height = MaxDistance;

  return Height;

}

su2double CPhysicalGeometry::Compute_LERadius(su2double *Plane_P0, su2double *Plane_Normal, vector<su2double> &Xcoord_Airfoil,
                                              vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {

  unsigned long iVertex, Trailing_Point, Leading_Point;
  su2double MaxDistance, Distance, LERadius = 0.0, X1, X2, X3, Y1, Y2, Y3, Ma, Mb, Xc, Yc, Radius;

  /*--- Find the leading and trailing edges and compute the radius of curvature ---*/

  MaxDistance = 0.0; Trailing_Point = 0;  Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {

    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  X1 = Xcoord_Airfoil[Leading_Point-3];
  Y1 = Zcoord_Airfoil[Leading_Point-3];

  X2 = Xcoord_Airfoil[Leading_Point];
  Y2 = Zcoord_Airfoil[Leading_Point];

  X3 = Xcoord_Airfoil[Leading_Point+3];
  Y3 = Zcoord_Airfoil[Leading_Point+3];

  if (X2 != X1) Ma = (Y2-Y1) / (X2-X1); else Ma = 0.0;
  if (X3 != X2) Mb = (Y3-Y2) / (X3-X2); else Mb = 0.0;

  if (Mb != Ma) Xc = (Ma*Mb*(Y1-Y3)+Mb*(X1+X2)-Ma*(X2+X3))/(2.0*(Mb-Ma)); else Xc = 0.0;
  if (Ma != 0.0) Yc = -(1.0/Ma)*(Xc-0.5*(X1+X2))+0.5*(Y1+Y2); else Yc = 0.0;

  Radius = sqrt((Xc-X1)*(Xc-X1)+(Yc-Y1)*(Yc-Y1));
  if (Radius != 0.0) LERadius = 1.0/Radius; else LERadius = 0.0;

  return LERadius;

}

su2double CPhysicalGeometry::Compute_Thickness(su2double *Plane_P0, su2double *Plane_Normal, su2double Location, CConfig *config,
                                               vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil, su2double &ZLoc) {

  unsigned long iVertex, jVertex, n_Upper, n_Lower, Trailing_Point, Leading_Point;
  su2double Thickness_Location, Normal[3], Tangent[3], BiNormal[3], auxXCoord, auxYCoord, auxZCoord, Thickness_Value = 0.0, Length,
            Xcoord_Trailing, Ycoord_Trailing, Zcoord_Trailing, ValCos, ValSin, XValue, ZValue, zp1, zpn, Chord, MaxDistance, Distance, AoA;

  vector<su2double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Z2coord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord_Lower,
                    Z2coord, Xcoord_Normal, Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;

  su2double Zcoord_Up, Zcoord_Down, ZLoc_, YLoc_;

  /*--- Find the leading and trailing edges and compute the angle of attack ---*/

  MaxDistance = 0.0; Trailing_Point = 0; Leading_Point = 0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Distance = sqrt(pow(Xcoord_Airfoil[iVertex] - Xcoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Ycoord_Airfoil[iVertex] - Ycoord_Airfoil[Trailing_Point], 2.0) +
                    pow(Zcoord_Airfoil[iVertex] - Zcoord_Airfoil[Trailing_Point], 2.0));

    if (MaxDistance < Distance) { MaxDistance = Distance; Leading_Point = iVertex; }
  }

  AoA = atan((Zcoord_Airfoil[Leading_Point] - Zcoord_Airfoil[Trailing_Point]) / (Xcoord_Airfoil[Trailing_Point] - Xcoord_Airfoil[Leading_Point]))*180/PI_NUMBER;
  Chord = MaxDistance;

  /*--- Translate to the origin ---*/

  Xcoord_Trailing = Xcoord_Airfoil[0];
  Ycoord_Trailing = Ycoord_Airfoil[0];
  Zcoord_Trailing = Zcoord_Airfoil[0];

  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    Xcoord_Airfoil_.push_back(Xcoord_Airfoil[iVertex] - Xcoord_Trailing);
    Ycoord_Airfoil_.push_back(Ycoord_Airfoil[iVertex] - Ycoord_Trailing);
    Zcoord_Airfoil_.push_back(Zcoord_Airfoil[iVertex] - Zcoord_Trailing);
  }

  /*--- Rotate the airfoil ---*/

  ValCos = cos(AoA*PI_NUMBER/180.0);
  ValSin = sin(AoA*PI_NUMBER/180.0);

  for (iVertex = 0; iVertex < Xcoord_Airfoil.size(); iVertex++) {
    XValue = Xcoord_Airfoil_[iVertex];
    ZValue = Zcoord_Airfoil_[iVertex];

    Xcoord_Airfoil_[iVertex] = XValue*ValCos - ZValue*ValSin;
    Zcoord_Airfoil_[iVertex] = ZValue*ValCos + XValue*ValSin;
  }

  /*--- Identify upper and lower side, and store the value of the normal --*/

  for (iVertex = 1; iVertex < Xcoord_Airfoil_.size(); iVertex++) {
    Tangent[0] = Xcoord_Airfoil_[iVertex] - Xcoord_Airfoil_[iVertex-1];
    Tangent[1] = Ycoord_Airfoil_[iVertex] - Ycoord_Airfoil_[iVertex-1];
    Tangent[2] = Zcoord_Airfoil_[iVertex] - Zcoord_Airfoil_[iVertex-1];
    Length = sqrt(pow(Tangent[0], 2.0) + pow(Tangent[1], 2.0) + pow(Tangent[2], 2.0));
    Tangent[0] /= Length; Tangent[1] /= Length; Tangent[2] /= Length;

    BiNormal[0] = Plane_Normal[0];
    BiNormal[1] = Plane_Normal[1];
    BiNormal[2] = Plane_Normal[2];
    Length = sqrt(pow(BiNormal[0], 2.0) + pow(BiNormal[1], 2.0) + pow(BiNormal[2], 2.0));
    BiNormal[0] /= Length; BiNormal[1] /= Length; BiNormal[2] /= Length;

    Normal[0] = Tangent[1]*BiNormal[2] - Tangent[2]*BiNormal[1];
    Normal[1] = Tangent[2]*BiNormal[0] - Tangent[0]*BiNormal[2];
    Normal[2] = Tangent[0]*BiNormal[1] - Tangent[1]*BiNormal[0];

    Xcoord_Normal.push_back(Normal[0]); Ycoord_Normal.push_back(Normal[1]); Zcoord_Normal.push_back(Normal[2]);

    unsigned short index = 2;

    if (Normal[index] >= 0.0) {
      Xcoord_Upper.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Upper.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Upper.push_back(Zcoord_Airfoil_[iVertex]);
    }
    else {
      Xcoord_Lower.push_back(Xcoord_Airfoil_[iVertex]);
      Ycoord_Lower.push_back(Ycoord_Airfoil_[iVertex]);
      Zcoord_Lower.push_back(Zcoord_Airfoil_[iVertex]);
    }

  }

  /*--- Order the arrays using the X component ---*/

  for (iVertex = 0; iVertex < Xcoord_Upper.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Upper.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Upper[jVertex] > Xcoord_Upper[jVertex+1]) {
        auxXCoord = Xcoord_Upper[jVertex]; Xcoord_Upper[jVertex] = Xcoord_Upper[jVertex+1]; Xcoord_Upper[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Upper[jVertex]; Ycoord_Upper[jVertex] = Ycoord_Upper[jVertex+1]; Ycoord_Upper[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Upper[jVertex]; Zcoord_Upper[jVertex] = Zcoord_Upper[jVertex+1]; Zcoord_Upper[jVertex+1] = auxZCoord;
      }
    }
  }

  /*--- Order the arrays using the X component ---*/

  for (iVertex = 0; iVertex < Xcoord_Lower.size(); iVertex++) {
    for (jVertex = 0; jVertex < Xcoord_Lower.size() - 1 - iVertex; jVertex++) {
      if (Xcoord_Lower[jVertex] > Xcoord_Lower[jVertex+1]) {
        auxXCoord = Xcoord_Lower[jVertex]; Xcoord_Lower[jVertex] = Xcoord_Lower[jVertex+1]; Xcoord_Lower[jVertex+1] = auxXCoord;
        auxYCoord = Ycoord_Lower[jVertex]; Ycoord_Lower[jVertex] = Ycoord_Lower[jVertex+1]; Ycoord_Lower[jVertex+1] = auxYCoord;
        auxZCoord = Zcoord_Lower[jVertex]; Zcoord_Lower[jVertex] = Zcoord_Lower[jVertex+1]; Zcoord_Lower[jVertex+1] = auxZCoord;
      }
    }
  }

  n_Upper = Xcoord_Upper.size();
  if (n_Upper > 1) {
    zp1 = (Zcoord_Upper[1]-Zcoord_Upper[0])/(Xcoord_Upper[1]-Xcoord_Upper[0]);
    zpn = (Zcoord_Upper[n_Upper-1]-Zcoord_Upper[n_Upper-2])/(Xcoord_Upper[n_Upper-1]-Xcoord_Upper[n_Upper-2]);
    Z2coord_Upper.resize(n_Upper+1);
    SetSpline(Xcoord_Upper, Zcoord_Upper, n_Upper, zp1, zpn, Z2coord_Upper);
  }

  n_Lower = Xcoord_Lower.size();
  if (n_Lower > 1) {
    zp1 = (Zcoord_Lower[1]-Zcoord_Lower[0])/(Xcoord_Lower[1]-Xcoord_Lower[0]);
    zpn = (Zcoord_Lower[n_Lower-1]-Zcoord_Lower[n_Lower-2])/(Xcoord_Lower[n_Lower-1]-Xcoord_Lower[n_Lower-2]);
    Z2coord_Lower.resize(n_Lower+1);
    SetSpline(Xcoord_Lower, Zcoord_Lower, n_Lower, zp1, zpn, Z2coord_Lower);
  }

  if ((n_Upper > 1) && (n_Lower > 1)) {

    Thickness_Location = - Chord*(1.0-Location);

    Zcoord_Up = GetSpline(Xcoord_Upper, Zcoord_Upper, Z2coord_Upper, n_Upper, Thickness_Location);
    Zcoord_Down = GetSpline(Xcoord_Lower, Zcoord_Lower, Z2coord_Lower, n_Lower, Thickness_Location);

    YLoc_ = Thickness_Location;
    ZLoc_ = 0.5*(Zcoord_Up + Zcoord_Down);

    ZLoc = sin(-AoA*PI_NUMBER/180.0)*YLoc_ + cos(-AoA*PI_NUMBER/180.0)*ZLoc_ + Zcoord_Trailing;

    /*--- Compute the thickness (we add a fabs because we can not guarantee the
     right sorting of the points and the upper and/or lower part of the airfoil is not well defined) ---*/

    Thickness_Value = fabs(Zcoord_Up - Zcoord_Down);

  }
  else { Thickness_Value = 0.0; }

  return Thickness_Value;

}

su2double CPhysicalGeometry::Compute_Area(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config,
                                          vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex;
  su2double Area_Value = 0.0;
  vector<su2double> Xcoord_Upper, Ycoord_Upper, Zcoord_Upper, Xcoord_Lower, Ycoord_Lower, Zcoord_Lower, Z2coord, Xcoord_Normal,
                    Ycoord_Normal, Zcoord_Normal, Xcoord_Airfoil_, Ycoord_Airfoil_, Zcoord_Airfoil_;
  su2double DeltaZ, DeltaX, X, Z;

  /*--- Use the Green theorem to evaluate the area (the points have been sortered),
   we assume that the airfoil is in the X-Z plane  ---*/

  Area_Value = 0.0;

  for (iVertex = 0; iVertex < Xcoord_Airfoil.size()-1; iVertex++) {
    X = 0.5*(Xcoord_Airfoil[iVertex]+Xcoord_Airfoil[iVertex+1]);
    Z = 0.5*(Zcoord_Airfoil[iVertex]+Zcoord_Airfoil[iVertex+1]);
    DeltaX = Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex];
    DeltaZ = Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex];
    Area_Value += 0.5*( X*DeltaZ-Z*DeltaX);
  }

  X = 0.5*(Xcoord_Airfoil[Xcoord_Airfoil.size()-1]+Xcoord_Airfoil[0]);
  Z = 0.5*(Zcoord_Airfoil[Xcoord_Airfoil.size()-1]+Zcoord_Airfoil[0]);
  DeltaX = Xcoord_Airfoil[0] - Xcoord_Airfoil[Xcoord_Airfoil.size()-1];
  DeltaZ = Zcoord_Airfoil[0] - Zcoord_Airfoil[Xcoord_Airfoil.size()-1];
  Area_Value += 0.5 * (X*DeltaZ-Z*DeltaX);

  Area_Value = fabs(Area_Value);

  return Area_Value;

}

su2double CPhysicalGeometry::Compute_Length(su2double *Plane_P0, su2double *Plane_Normal, CConfig *config,
                                          vector<su2double> &Xcoord_Airfoil, vector<su2double> &Ycoord_Airfoil, vector<su2double> &Zcoord_Airfoil) {
  unsigned long iVertex;
  su2double Length_Value = 0.0, Length_Value_ = 0.0;
  su2double DeltaZ, DeltaX;

  /*--- Not that in a symmetry plane configuration there is an extra edge that connects
   the two extremes, and we really don't now the curve orientation. We will evaluate
   both distance and picked the smallest one ---*/

  Length_Value = 0.0;
  for (iVertex = 0; iVertex < Xcoord_Airfoil.size()-2; iVertex++) {
    DeltaX = Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex];
    DeltaZ = Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex];
    Length_Value += sqrt(DeltaX*DeltaX + DeltaZ*DeltaZ);
  }

  Length_Value_ = 0.0;
  for (iVertex = 1; iVertex < Xcoord_Airfoil.size()-1; iVertex++) {
    DeltaX = Xcoord_Airfoil[iVertex+1] - Xcoord_Airfoil[iVertex];
    DeltaZ = Zcoord_Airfoil[iVertex+1] - Zcoord_Airfoil[iVertex];
    Length_Value_ += sqrt(DeltaX*DeltaX + DeltaZ*DeltaZ);
  }

  Length_Value = min(Length_Value, Length_Value_);

  return Length_Value;

}

void CPhysicalGeometry::Compute_Wing(CConfig *config, bool original_surface,
                                     su2double &Wing_Volume, su2double &Wing_MinMaxThickness, su2double &Wing_MaxMaxThickness, su2double &Wing_MinChord, su2double &Wing_MaxChord,
                                     su2double &Wing_MinLERadius, su2double &Wing_MaxLERadius,
                                     su2double &Wing_MinToC, su2double &Wing_MaxToC, su2double &Wing_ObjFun_MinToC, su2double &Wing_MaxTwist, su2double &Wing_MaxCurvature,
                                     su2double &Wing_MaxDihedral) {

  unsigned short iPlane, iDim, nPlane = 0;
  unsigned long iVertex;
  su2double MinPlane, MaxPlane, dPlane, *Area, *MaxThickness, *ToC, *Chord, *LERadius, *Twist, *Curvature, *Dihedral, SemiSpan;
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  ofstream Wing_File, Section_File;

  /*--- Make a large number of section cuts for approximating volume ---*/

  nPlane = config->GetnWingStations();
  SemiSpan = config->GetSemiSpan();

  /*--- Allocate memory for the section cutting ---*/

  Area = new su2double [nPlane];
  MaxThickness = new su2double [nPlane];
  Chord = new su2double [nPlane];
  LERadius = new su2double [nPlane];
  ToC = new su2double [nPlane];
  Twist = new su2double [nPlane];
  Curvature = new su2double [nPlane];
  Dihedral = new su2double [nPlane];

  su2double **LeadingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    LeadingEdge[iPlane] = new su2double[nDim];

  su2double **TrailingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    TrailingEdge[iPlane] = new su2double[nDim];

  su2double **Plane_P0 = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_P0[iPlane] = new su2double[nDim];

  su2double **Plane_Normal = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_Normal[iPlane] = new su2double[nDim];

  MinPlane = config->GetStations_Bounds(0); MaxPlane = config->GetStations_Bounds(1);
  dPlane = fabs((MaxPlane - MinPlane)/su2double(nPlane-1));

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
    Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;

    if (config->GetGeo_Description() == WING) {
      Plane_Normal[iPlane][1] = 1.0;
      Plane_P0[iPlane][1] = MinPlane + iPlane*dPlane;
    }

    if (config->GetGeo_Description() == TWOD_AIRFOIL) {
      Plane_Normal[iPlane][2] = 1.0;
      Plane_P0[iPlane][2] = MinPlane + iPlane*dPlane;
    }

  }

  /*--- Allocate some vectors for storing airfoil coordinates ---*/

  Xcoord_Airfoil   = new vector<su2double>[nPlane];
  Ycoord_Airfoil   = new vector<su2double>[nPlane];
  Zcoord_Airfoil   = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  /*--- Create the section slices through the geometry ---*/

  for (iPlane = 0; iPlane < nPlane; iPlane++) {

    ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane],
                           -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, NULL, Xcoord_Airfoil[iPlane],
                           Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                           Variable_Airfoil[iPlane], original_surface, config);

   }

  /*--- Compute airfoil characteristic only in the master node ---*/

  if (rank == MASTER_NODE) {

    /*--- Write an output file---*/

    if (config->GetTabular_FileFormat() == TAB_CSV) {
      Wing_File.open("wing_description.csv", ios::out);
      if (config->GetSystemMeasurements() == US)
        Wing_File << "\"yCoord/SemiSpan\",\"Area (in^2)\",\"Max. Thickness (in)\",\"Chord (in)\",\"Leading Edge Radius (1/in)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge XLoc/SemiSpan\",\"Leading Edge ZLoc/SemiSpan\",\"Trailing Edge XLoc/SemiSpan\",\"Trailing Edge ZLoc/SemiSpan\"" << endl;
      else
        Wing_File << "\"yCoord/SemiSpan\",\"Area (m^2)\",\"Max. Thickness (m)\",\"Chord (m)\",\"Leading Edge Radius (1/m)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge XLoc/SemiSpan\",\"Leading Edge ZLoc/SemiSpan\",\"Trailing Edge XLoc/SemiSpan\",\"Trailing Edge ZLoc/SemiSpan\"" << endl;
    }
    else {
      Wing_File.open("wing_description.dat", ios::out);
      Wing_File << "TITLE = \"Wing description\"" << endl;
      if (config->GetSystemMeasurements() == US)
        Wing_File << "VARIABLES = \"<greek>h</greek>\",\"Area (in<sup>2</sup>)\",\"Max. Thickness (in)\",\"Chord (in)\",\"Leading Edge Radius (1/in)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge XLoc/SemiSpan\",\"Leading Edge ZLoc/SemiSpan\",\"Trailing Edge XLoc/SemiSpan\",\"Trailing Edge ZLoc/SemiSpan\"" << endl;
      else
        Wing_File << "VARIABLES = \"<greek>h</greek>\",\"Area (m<sup>2</sup>)\",\"Max. Thickness (m)\",\"Chord (m)\",\"Leading Edge Radius (1/m)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Curvature (1/m)\",\"Dihedral (deg)\",\"Leading Edge XLoc/SemiSpan\",\"Leading Edge ZLoc/SemiSpan\",\"Trailing Edge XLoc/SemiSpan\",\"Trailing Edge ZLoc/SemiSpan\"" << endl;
      Wing_File << "ZONE T= \"Baseline wing\"" << endl;
    }


    /*--- Evaluate  geometrical quatities that do not require any kind of filter, local to each point ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        LeadingEdge[iPlane][iDim]  = 0.0;
        TrailingEdge[iPlane][iDim] = 0.0;
      }

      Area[iPlane]                = 0.0;
      MaxThickness[iPlane]        = 0.0;
      Chord[iPlane]               = 0.0;
      LERadius[iPlane]            = 0.0;
      ToC[iPlane]                 = 0.0;
      Twist[iPlane]               = 0.0;

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        Compute_Wing_LeadingTrailing(LeadingEdge[iPlane], TrailingEdge[iPlane], Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        MaxThickness[iPlane] = Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Chord[iPlane] = Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Twist[iPlane] = Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        LERadius[iPlane] = Compute_LERadius(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        ToC[iPlane] = MaxThickness[iPlane] / Chord[iPlane];


      }

    }

    /*--- Evaluate  geometrical quatities that have been computed using a filtered value (they depend on more than one point) ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      Curvature[iPlane]      = 0.0;
      Dihedral[iPlane]       = 0.0;

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        if ((iPlane == 0) || (iPlane == nPlane-1)) Curvature[iPlane] = 0.0;
        else Curvature[iPlane] = Compute_Curvature(LeadingEdge[iPlane-1], TrailingEdge[iPlane-1],
                                                   LeadingEdge[iPlane], TrailingEdge[iPlane],
                                                   LeadingEdge[iPlane+1], TrailingEdge[iPlane+1]);

        if (iPlane == 0) Dihedral[iPlane] = 0.0;
        else Dihedral[iPlane] = Compute_Dihedral(LeadingEdge[iPlane-1], TrailingEdge[iPlane-1],
                                                 LeadingEdge[iPlane], TrailingEdge[iPlane]);

      }

    }

    /*--- Set the curvature and dihedral angles at the extremes ---*/

    if (nPlane > 1) {
      if ((Xcoord_Airfoil[0].size() != 0) && (Xcoord_Airfoil[1].size() != 0)) {
        Curvature[0] = Curvature[1]; Dihedral[0] = Dihedral[1];
      }
      if ((Xcoord_Airfoil[nPlane-1].size() != 0) && (Xcoord_Airfoil[nPlane-2].size() != 0)) {
        Curvature[nPlane-1] = Curvature[nPlane-2];
      }
    }

    /*--- Plot the geometrical quatities ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
    	if (Xcoord_Airfoil[iPlane].size() > 1) {
    		if (config->GetTabular_FileFormat() == TAB_CSV) {
    			Wing_File  << Ycoord_Airfoil[iPlane][0]/SemiSpan <<", "<< Area[iPlane] <<", "<< MaxThickness[iPlane] <<", "<< Chord[iPlane] <<", "<< LERadius[iPlane] <<", "<< ToC[iPlane]
    			           <<", "<< Twist[iPlane] <<", "<< Curvature[iPlane] <<", "<< Dihedral[iPlane]
    			           <<", "<< LeadingEdge[iPlane][0]/SemiSpan <<", "<< LeadingEdge[iPlane][2]/SemiSpan
    			           <<", "<< TrailingEdge[iPlane][0]/SemiSpan <<", "<< TrailingEdge[iPlane][2]/SemiSpan << endl;
    		}
    		else  {
    			Wing_File  << Ycoord_Airfoil[iPlane][0]/SemiSpan <<" "<< Area[iPlane] <<" "<< MaxThickness[iPlane] <<" "<< Chord[iPlane] <<" "<< LERadius[iPlane] <<" "<< ToC[iPlane]
    			           <<" "<< Twist[iPlane] <<" "<< Curvature[iPlane]  <<" "<< Dihedral[iPlane]
    			           <<" "<< LeadingEdge[iPlane][0]/SemiSpan <<" "<< LeadingEdge[iPlane][2]/SemiSpan
    			           <<" "<< TrailingEdge[iPlane][0]/SemiSpan <<" "<< TrailingEdge[iPlane][2]/SemiSpan << endl;

    		}
    	}
    }

    Wing_File.close();

    Section_File.open("wing_slices.dat", ios::out);

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (iPlane == 0) {
        Section_File << "TITLE = \"Aircraft Slices\"" << endl;
        if (config->GetSystemMeasurements() == US)
          Section_File << "VARIABLES = \"x (in)\", \"y (in)\", \"z (in)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
        else Section_File << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
      }

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        Section_File << "ZONE T=\"<greek>h</greek> = " << Ycoord_Airfoil[iPlane][0]/SemiSpan << " \", I= " << Xcoord_Airfoil[iPlane].size() << ", F=POINT" << endl;

        for (iVertex = 0; iVertex < Xcoord_Airfoil[iPlane].size(); iVertex++) {

          /*--- Move to the origin  ---*/

          su2double XValue_ = Xcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][0];
          su2double ZValue_ = Zcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][2];

          /*--- Rotate the airfoil and divide by the chord ---*/

          su2double ValCos = cos(Twist[iPlane]*PI_NUMBER/180.0);
          su2double ValSin = sin(Twist[iPlane]*PI_NUMBER/180.0);

          su2double XValue = (XValue_*ValCos - ZValue_*ValSin) / Chord[iPlane];
          su2double ZValue = (ZValue_*ValCos + XValue_*ValSin) / Chord[iPlane];

          /*--- Write the file ---*/

          Section_File  << Xcoord_Airfoil[iPlane][iVertex] << " " << Ycoord_Airfoil[iPlane][iVertex] << " " << Zcoord_Airfoil[iPlane][iVertex] << " " << XValue  << " " << ZValue << endl;
        }
      }

    }

    Section_File.close();


    /*--- Compute the wing volume using a composite Simpson's rule ---*/

    Wing_Volume = 0.0;
    for (iPlane = 0; iPlane < nPlane-2; iPlane+=2) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
        Wing_Volume += (1.0/3.0)*dPlane*(Area[iPlane] + 4.0*Area[iPlane+1] + Area[iPlane+2]);
      }
    }

    /*--- Evaluate Max and Min quantities ---*/

    Wing_MaxMaxThickness = -1E6; Wing_MinMaxThickness = 1E6; Wing_MinChord = 1E6; Wing_MaxChord = -1E6;
    Wing_MinLERadius = 1E6; Wing_MaxLERadius = -1E6; Wing_MinToC = 1E6; Wing_MaxToC = -1E6;
    Wing_MaxTwist = -1E6; Wing_MaxCurvature = -1E6; Wing_MaxDihedral = -1E6;

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      if (MaxThickness[iPlane] != 0.0) Wing_MinMaxThickness = min(Wing_MinMaxThickness, MaxThickness[iPlane]);
      Wing_MaxMaxThickness = max(Wing_MaxMaxThickness, MaxThickness[iPlane]);
      if (Chord[iPlane] != 0.0) Wing_MinChord = min(Wing_MinChord, Chord[iPlane]);
      Wing_MaxChord = max(Wing_MaxChord, Chord[iPlane]);
      if (LERadius[iPlane] != 0.0) Wing_MinLERadius = min(Wing_MinLERadius, LERadius[iPlane]);
      Wing_MaxLERadius = max(Wing_MaxLERadius, LERadius[iPlane]);
      if (ToC[iPlane] != 0.0) Wing_MinToC = min(Wing_MinToC, ToC[iPlane]);
      Wing_MaxToC = max(Wing_MaxToC, ToC[iPlane]);
      Wing_ObjFun_MinToC = sqrt((Wing_MinToC - 0.07)*(Wing_MinToC - 0.07));
      Wing_MaxTwist = max(Wing_MaxTwist, fabs(Twist[iPlane]));
      Wing_MaxCurvature = max(Wing_MaxCurvature, Curvature[iPlane]);
      Wing_MaxDihedral = max(Wing_MaxDihedral, fabs(Dihedral[iPlane]));
    }

  }

  /*--- Free memory for the section cuts ---*/

  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  delete [] Variable_Airfoil;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] LeadingEdge[iPlane];
  delete [] LeadingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] TrailingEdge[iPlane];
  delete [] TrailingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_P0[iPlane];
  delete [] Plane_P0;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_Normal[iPlane];
  delete [] Plane_Normal;

  delete [] Area;
  delete [] MaxThickness;
  delete [] Chord;
  delete [] LERadius;
  delete [] ToC;
  delete [] Twist;
  delete [] Curvature;
  delete [] Dihedral;

}

void CPhysicalGeometry::Compute_Fuselage(CConfig *config, bool original_surface,
                                         su2double &Fuselage_Volume, su2double &Fuselage_WettedArea,
                                         su2double &Fuselage_MinWidth, su2double &Fuselage_MaxWidth,
                                         su2double &Fuselage_MinWaterLineWidth, su2double &Fuselage_MaxWaterLineWidth,
                                         su2double &Fuselage_MinHeight, su2double &Fuselage_MaxHeight,
                                         su2double &Fuselage_MaxCurvature) {

  unsigned short iPlane, iDim, nPlane = 0;
  unsigned long iVertex;
  su2double MinPlane, MaxPlane, dPlane, *Area, *Length, *Width, *WaterLineWidth, *Height, *Curvature;
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  ofstream Fuselage_File, Section_File;

  /*--- Make a large number of section cuts for approximating volume ---*/

  nPlane = config->GetnWingStations();

  /*--- Allocate memory for the section cutting ---*/

  Area = new su2double [nPlane];
  Length = new su2double [nPlane];
  Width = new su2double [nPlane];
  WaterLineWidth = new su2double [nPlane];
  Height = new su2double [nPlane];
  Curvature = new su2double [nPlane];

  su2double **LeadingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    LeadingEdge[iPlane] = new su2double[nDim];

  su2double **TrailingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    TrailingEdge[iPlane] = new su2double[nDim];

  su2double **Plane_P0 = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_P0[iPlane] = new su2double[nDim];

  su2double **Plane_Normal = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_Normal[iPlane] = new su2double[nDim];

  MinPlane = config->GetStations_Bounds(0); MaxPlane = config->GetStations_Bounds(1);
  dPlane = fabs((MaxPlane - MinPlane)/su2double(nPlane-1));

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
    Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;

    Plane_Normal[iPlane][0] = 1.0;
    Plane_P0[iPlane][0] = MinPlane + iPlane*dPlane;

  }

  /*--- Allocate some vectors for storing airfoil coordinates ---*/

  Xcoord_Airfoil   = new vector<su2double>[nPlane];
  Ycoord_Airfoil   = new vector<su2double>[nPlane];
  Zcoord_Airfoil   = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  /*--- Create the section slices through the geometry ---*/

  for (iPlane = 0; iPlane < nPlane; iPlane++) {

    ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane],
                           -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, NULL, Xcoord_Airfoil[iPlane],
                           Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                           Variable_Airfoil[iPlane], original_surface, config);
  }

  /*--- Compute the area at each section ---*/

  if (rank == MASTER_NODE) {

    /*--- Write an output file---*/

    if (config->GetTabular_FileFormat() == TAB_CSV) {
      Fuselage_File.open("fuselage_description.csv", ios::out);
      if (config->GetSystemMeasurements() == US)
        Fuselage_File << "\"x (in)\",\"Area (in^2)\",\"Length (in)\",\"Width (in)\",\"Waterline width (in)\",\"Height (in)\",\"Curvature (1/in)\",\"Generatrix Curve X (in)\",\"Generatrix Curve Y (in)\",\"Generatrix Curve Z (in)\",\"Axis Curve X (in)\",\"Axis Curve Y (in)\",\"Axis Curve Z (in)\"" << endl;
      else
        Fuselage_File << "\"x (m)\",\"Area (m^2)\",\"Length (m)\",\"Width (m)\",\"Waterline width (m)\",\"Height (m)\",\"Curvature (1/in)\",\"Generatrix Curve X (m)\",\"Generatrix Curve Y (m)\",\"Generatrix Curve Z (m)\",\"Axis Curve X (m)\",\"Axis Curve Y (m)\",\"Axis Curve Z (m)\"" << endl;
    }
    else {
      Fuselage_File.open("fuselage_description.dat", ios::out);
      Fuselage_File << "TITLE = \"Fuselage description\"" << endl;
      if (config->GetSystemMeasurements() == US)
        Fuselage_File << "VARIABLES = \"x (in)\",\"Area (in<sup>2</sup>)\",\"Length (in)\",\"Width (in)\",\"Waterline width (in)\",\"Height (in)\",\"Curvature (1/in)\",\"Generatrix Curve X (in)\",\"Generatrix Curve Y (in)\",\"Generatrix Curve Z (in)\",\"Axis Curve X (in)\",\"Axis Curve Y (in)\",\"Axis Curve Z (in)\"" << endl;
      else
        Fuselage_File << "VARIABLES = \"x (m)\",\"Area (m<sup>2</sup>)\",\"Length (m)\",\"Width (m)\",\"Waterline width (m)\",\"Height (m)\",\"Curvature (1/m)\",\"Generatrix Curve X (m)\",\"Generatrix Curve Y (m)\",\"Generatrix Curve Z (m)\",\"Axis Curve X (m)\",\"Axis Curve Y (m)\",\"Axis Curve Z (m)\"" << endl;
      Fuselage_File << "ZONE T= \"Baseline fuselage\"" << endl;
    }


    /*--- Evaluate  geometrical quatities that do not require any kind of filter, local to each point ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        LeadingEdge[iPlane][iDim]  = 0.0;
        TrailingEdge[iPlane][iDim] = 0.0;
      }

      Area[iPlane]   = 0.0;
      Length[iPlane]   = 0.0;
      Width[iPlane]  = 0.0;
      WaterLineWidth[iPlane]  = 0.0;
      Height[iPlane] = 0.0;

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        Compute_Fuselage_LeadingTrailing(LeadingEdge[iPlane], TrailingEdge[iPlane], Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Length[iPlane] = Compute_Length(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Width[iPlane] = Compute_Width(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        WaterLineWidth[iPlane] = Compute_WaterLineWidth(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Height[iPlane] = Compute_Height(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

      }

    }

    /*--- Evaluate  geometrical quatities that have been computed using a filtered value (they depend on more than one point) ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      Curvature[iPlane] = 0.0;

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        if ((iPlane == 0) || (iPlane == nPlane-1)) Curvature[iPlane] = 0.0;
        else Curvature[iPlane] = Compute_Curvature(LeadingEdge[iPlane-1], TrailingEdge[iPlane-1],
                                                   LeadingEdge[iPlane], TrailingEdge[iPlane],
                                                   LeadingEdge[iPlane+1], TrailingEdge[iPlane+1]);

      }

    }

    /*--- Set the curvature and dihedral angles at the extremes ---*/

    if (nPlane > 1) {
      if ((Xcoord_Airfoil[0].size() != 0) && (Xcoord_Airfoil[1].size() != 0)) {
        Curvature[0] = Curvature[1];
      }
      if ((Xcoord_Airfoil[nPlane-1].size() != 0) && (Xcoord_Airfoil[nPlane-2].size() != 0)) {
        Curvature[nPlane-1] = Curvature[nPlane-2];
      }
    }

    /*--- Plot the geometrical quatities ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
    	if (Xcoord_Airfoil[iPlane].size() > 1) {
    		if (config->GetTabular_FileFormat() == TAB_CSV) {
    			Fuselage_File  << -Ycoord_Airfoil[iPlane][0] <<", "<< Area[iPlane] <<", "<< Length[iPlane] <<", "<< Width[iPlane] <<", "<< WaterLineWidth[iPlane] <<", "<< Height[iPlane] <<", "<< Curvature[iPlane]
    			           <<", "<< -LeadingEdge[iPlane][1] <<", "<< LeadingEdge[iPlane][0]  <<", "<< LeadingEdge[iPlane][2]
    			           <<", "<< -TrailingEdge[iPlane][1] <<", "<< TrailingEdge[iPlane][0]  <<", "<< TrailingEdge[iPlane][2]  << endl;
    		}
    		else  {
    			Fuselage_File  << -Ycoord_Airfoil[iPlane][0] <<" "<< Area[iPlane] <<" "<< Length[iPlane] <<" "<< Width[iPlane] <<" "<< WaterLineWidth[iPlane] <<" "<< Height[iPlane] <<" "<< Curvature[iPlane]
    			           <<" "<< -LeadingEdge[iPlane][1] <<" "<< LeadingEdge[iPlane][0]  <<" "<< LeadingEdge[iPlane][2]
    			           <<" "<< -TrailingEdge[iPlane][1] <<" "<< TrailingEdge[iPlane][0]  <<" "<< TrailingEdge[iPlane][2] << endl;
    		}
    	}
    }

    Fuselage_File.close();

    Section_File.open("fuselage_slices.dat", ios::out);

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (iPlane == 0) {
        Section_File << "TITLE = \"Aircraft Slices\"" << endl;
        if (config->GetSystemMeasurements() == US)
          Section_File << "VARIABLES = \"x (in)\", \"y (in)\", \"z (in)\"" << endl;
        else Section_File << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\"" << endl;
      }

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        Section_File << "ZONE T=\"X = " << -Ycoord_Airfoil[iPlane][0] << " \", I= " << Ycoord_Airfoil[iPlane].size() << ", F=POINT" << endl;

        for (iVertex = 0; iVertex < Xcoord_Airfoil[iPlane].size(); iVertex++) {

          /*--- Write the file ---*/

          Section_File  << -Ycoord_Airfoil[iPlane][iVertex] << " " << Xcoord_Airfoil[iPlane][iVertex] << " " << Zcoord_Airfoil[iPlane][iVertex] << endl;
        }
      }

    }

    Section_File.close();


    /*--- Compute the fuselage volume using a composite Simpson's rule ---*/

    Fuselage_Volume = 0.0;
    for (iPlane = 0; iPlane < nPlane-2; iPlane+=2) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
      	Fuselage_Volume += (1.0/3.0)*dPlane*(Area[iPlane] + 4.0*Area[iPlane+1] + Area[iPlane+2]);
      }
    }

    /*--- Compute the fuselage wetted area ---*/

    Fuselage_WettedArea = 0.0;
    if (Xcoord_Airfoil[0].size() > 1) Fuselage_WettedArea += (1.0/2.0)*dPlane*Length[0];
    for (iPlane = 1; iPlane < nPlane-1; iPlane++) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
      	Fuselage_WettedArea += dPlane*Length[iPlane];
      }
    }
    if (Xcoord_Airfoil[nPlane-1].size() > 1) Fuselage_WettedArea += (1.0/2.0)*dPlane*Length[nPlane-1];

    /*--- Evaluate Max and Min quantities ---*/

    Fuselage_MaxWidth = -1E6; Fuselage_MinWidth = 1E6;
    Fuselage_MaxWaterLineWidth = -1E6; Fuselage_MinWaterLineWidth = 1E6;
    Fuselage_MaxHeight = -1E6; Fuselage_MinHeight = 1E6;
    Fuselage_MaxCurvature = -1E6;

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      if (Width[iPlane] != 0.0) Fuselage_MinWidth = min(Fuselage_MinWidth, Width[iPlane]);
      Fuselage_MaxWidth = max(Fuselage_MaxWidth, Width[iPlane]);
      if (WaterLineWidth[iPlane] != 0.0) Fuselage_MinWaterLineWidth = min(Fuselage_MinWaterLineWidth, WaterLineWidth[iPlane]);
      Fuselage_MaxWaterLineWidth = max(Fuselage_MaxWaterLineWidth, WaterLineWidth[iPlane]);
      if (Height[iPlane] != 0.0) Fuselage_MinHeight = min(Fuselage_MinHeight, Height[iPlane]);
      Fuselage_MaxHeight = max(Fuselage_MaxHeight, Height[iPlane]);
      Fuselage_MaxCurvature = max(Fuselage_MaxCurvature, Curvature[iPlane]);
    }

  }

  /*--- Free memory for the section cuts ---*/

  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  delete [] Variable_Airfoil;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] LeadingEdge[iPlane];
  delete [] LeadingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] TrailingEdge[iPlane];
  delete [] TrailingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_P0[iPlane];
  delete [] Plane_P0;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_Normal[iPlane];
  delete [] Plane_Normal;

  delete [] Area;
  delete [] Length;
  delete [] Width;
  delete [] WaterLineWidth;
  delete [] Height;
  delete [] Curvature;

}

void CPhysicalGeometry::Compute_Nacelle(CConfig *config, bool original_surface,
                                        su2double &Nacelle_Volume, su2double &Nacelle_MinMaxThickness,
                                        su2double &Nacelle_MinChord, su2double &Nacelle_MaxChord,
                                        su2double &Nacelle_MaxMaxThickness, su2double &Nacelle_MinLERadius,
                                        su2double &Nacelle_MaxLERadius, su2double &Nacelle_MinToC, su2double &Nacelle_MaxToC,
                                        su2double &Nacelle_ObjFun_MinToC, su2double &Nacelle_MaxTwist) {

  unsigned short iPlane, iDim, nPlane = 0;
  unsigned long iVertex;
  su2double Angle, MinAngle, MaxAngle, dAngle, *Area, *MaxThickness, *ToC, *Chord, *LERadius, *Twist;
  vector<su2double> *Xcoord_Airfoil, *Ycoord_Airfoil, *Zcoord_Airfoil, *Variable_Airfoil;
  ofstream Nacelle_File, Section_File;


  /*--- Make a large number of section cuts for approximating volume ---*/

  nPlane = config->GetnWingStations();

  /*--- Allocate memory for the section cutting ---*/

  Area = new su2double [nPlane];
  MaxThickness = new su2double [nPlane];
  Chord = new su2double [nPlane];
  LERadius = new su2double [nPlane];
  ToC = new su2double [nPlane];
  Twist = new su2double [nPlane];

  su2double **LeadingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    LeadingEdge[iPlane] = new su2double[nDim];

  su2double **TrailingEdge = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    TrailingEdge[iPlane] = new su2double[nDim];

  su2double **Plane_P0 = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_P0[iPlane] = new su2double[nDim];

  su2double **Plane_Normal = new su2double*[nPlane];
  for (iPlane = 0; iPlane < nPlane; iPlane++ )
    Plane_Normal[iPlane] = new su2double[nDim];

  MinAngle = config->GetStations_Bounds(0); MaxAngle = config->GetStations_Bounds(1);
  dAngle = fabs((MaxAngle - MinAngle)/su2double(nPlane-1));

  for (iPlane = 0; iPlane < nPlane; iPlane++) {
    Plane_Normal[iPlane][0] = 0.0;    Plane_P0[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = 0.0;    Plane_P0[iPlane][1] = 0.0;
    Plane_Normal[iPlane][2] = 0.0;    Plane_P0[iPlane][2] = 0.0;

    /*--- Apply roll to cut the nacelle ---*/

    Angle = MinAngle + iPlane*dAngle*PI_NUMBER/180.0;

    if (Angle <= 0) Angle = 1E-6;
    if (Angle >= 360) Angle = 359.999999;

    Plane_Normal[iPlane][0] = 0.0;
    Plane_Normal[iPlane][1] = -sin(Angle);
    Plane_Normal[iPlane][2] = cos(Angle);

    /*--- Apply tilt angle to the plane ---*/

    su2double Tilt_Angle = config->GetNacelleLocation(3)*PI_NUMBER/180;
    su2double Plane_NormalX_Tilt = Plane_Normal[iPlane][0]*cos(Tilt_Angle) + Plane_Normal[iPlane][2]*sin(Tilt_Angle);
    su2double Plane_NormalY_Tilt = Plane_Normal[iPlane][1];
    su2double Plane_NormalZ_Tilt = Plane_Normal[iPlane][2]*cos(Tilt_Angle) - Plane_Normal[iPlane][0]*sin(Tilt_Angle);

    /*--- Apply toe angle to the plane ---*/

    su2double Toe_Angle = config->GetNacelleLocation(4)*PI_NUMBER/180;
    su2double Plane_NormalX_Tilt_Toe = Plane_NormalX_Tilt*cos(Toe_Angle) - Plane_NormalY_Tilt*sin(Toe_Angle);
    su2double Plane_NormalY_Tilt_Toe = Plane_NormalX_Tilt*sin(Toe_Angle) + Plane_NormalY_Tilt*cos(Toe_Angle);
    su2double Plane_NormalZ_Tilt_Toe = Plane_NormalZ_Tilt;

    /*--- Update normal vector ---*/

    Plane_Normal[iPlane][0] = Plane_NormalX_Tilt_Toe;
    Plane_Normal[iPlane][1] = Plane_NormalY_Tilt_Toe;
    Plane_Normal[iPlane][2] = Plane_NormalZ_Tilt_Toe;

    /*--- Point in the plane ---*/

    Plane_P0[iPlane][0] = config->GetNacelleLocation(0);
    Plane_P0[iPlane][1] = config->GetNacelleLocation(1);
    Plane_P0[iPlane][2] = config->GetNacelleLocation(2);

  }

  /*--- Allocate some vectors for storing airfoil coordinates ---*/

  Xcoord_Airfoil   = new vector<su2double>[nPlane];
  Ycoord_Airfoil   = new vector<su2double>[nPlane];
  Zcoord_Airfoil   = new vector<su2double>[nPlane];
  Variable_Airfoil = new vector<su2double>[nPlane];

  /*--- Create the section slices through the geometry ---*/

  for (iPlane = 0; iPlane < nPlane; iPlane++) {

    ComputeAirfoil_Section(Plane_P0[iPlane], Plane_Normal[iPlane],
                           -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, NULL, Xcoord_Airfoil[iPlane],
                           Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                           Variable_Airfoil[iPlane], original_surface, config);

  }

  /*--- Compute airfoil characteristic only in the master node ---*/

  if (rank == MASTER_NODE) {

    /*--- Write an output file---*/

    if (config->GetTabular_FileFormat() == TAB_CSV) {
      Nacelle_File.open("nacelle_description.csv", ios::out);
      if (config->GetSystemMeasurements() == US)
        Nacelle_File << "\"Theta (deg)\",\"Area (in^2)\",\"Max. Thickness (in)\",\"Chord (in)\",\"Leading Edge Radius (1/in)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Leading Edge XLoc\",\"Leading Edge ZLoc\",\"Trailing Edge XLoc\",\"Trailing Edge ZLoc\"" << endl;
      else
        Nacelle_File << "\"Theta (deg)\",\"Area (m^2)\",\"Max. Thickness (m)\",\"Chord (m)\",\"Leading Edge Radius (1/m)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Curvature (1/in)\",\"Dihedral (deg)\",\"Leading Edge XLoc\",\"Leading Edge ZLoc\",\"Trailing Edge XLoc\",\"Trailing Edge ZLoc\"" << endl;
    }
    else {
      Nacelle_File.open("nacelle_description.dat", ios::out);
      Nacelle_File << "TITLE = \"Nacelle description\"" << endl;
      if (config->GetSystemMeasurements() == US)
        Nacelle_File << "VARIABLES = \"<greek>q</greek> (deg)\",\"Area (in<sup>2</sup>)\",\"Max. Thickness (in)\",\"Chord (in)\",\"Leading Edge Radius (1/in)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Leading Edge XLoc\",\"Leading Edge ZLoc\",\"Trailing Edge XLoc\",\"Trailing Edge ZLoc\"" << endl;
      else
        Nacelle_File << "VARIABLES = \"<greek>q</greek> (deg)\",\"Area (m<sup>2</sup>)\",\"Max. Thickness (m)\",\"Chord (m)\",\"Leading Edge Radius (1/m)\",\"Max. Thickness/Chord\",\"Twist (deg)\",\"Leading Edge XLoc\",\"Leading Edge ZLoc\",\"Trailing Edge XLoc\",\"Trailing Edge ZLoc\"" << endl;
      Nacelle_File << "ZONE T= \"Baseline nacelle\"" << endl;
    }


    /*--- Evaluate  geometrical quatities that do not require any kind of filter, local to each point ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      for (iDim = 0; iDim < nDim; iDim++) {
        LeadingEdge[iPlane][iDim]  = 0.0;
        TrailingEdge[iPlane][iDim] = 0.0;
      }

      Area[iPlane]                = 0.0;
      MaxThickness[iPlane]        = 0.0;
      Chord[iPlane]               = 0.0;
      LERadius[iPlane]            = 0.0;
      ToC[iPlane]                 = 0.0;
      Twist[iPlane]               = 0.0;

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        Compute_Wing_LeadingTrailing(LeadingEdge[iPlane], TrailingEdge[iPlane], Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Area[iPlane] = Compute_Area(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        MaxThickness[iPlane] = Compute_MaxThickness(Plane_P0[iPlane], Plane_Normal[iPlane], config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Chord[iPlane] = Compute_Chord(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        Twist[iPlane] = Compute_Twist(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        LERadius[iPlane] = Compute_LERadius(Plane_P0[iPlane], Plane_Normal[iPlane], Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        ToC[iPlane] = MaxThickness[iPlane] / Chord[iPlane];


      }

    }

    /*--- Plot the geometrical quatities ---*/

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      su2double theta_deg = atan2(Plane_Normal[iPlane][1], -Plane_Normal[iPlane][2])/PI_NUMBER*180 + 180;

      if (Xcoord_Airfoil[iPlane].size() > 1) {
        if (config->GetTabular_FileFormat() == TAB_CSV) {
          Nacelle_File  << theta_deg <<", "<< Area[iPlane] <<", "<< MaxThickness[iPlane] <<", "<< Chord[iPlane] <<", "<< LERadius[iPlane] <<", "<< ToC[iPlane]
          <<", "<< Twist[iPlane] <<", "<< LeadingEdge[iPlane][0] <<", "<< LeadingEdge[iPlane][2]
          <<", "<< TrailingEdge[iPlane][0] <<", "<< TrailingEdge[iPlane][2] << endl;
        }
        else  {
          Nacelle_File  << theta_deg <<" "<< Area[iPlane] <<" "<< MaxThickness[iPlane] <<" "<< Chord[iPlane] <<" "<< LERadius[iPlane] <<" "<< ToC[iPlane]
          <<" "<< Twist[iPlane] <<" "<< LeadingEdge[iPlane][0] <<" "<< LeadingEdge[iPlane][2]
          <<" "<< TrailingEdge[iPlane][0] <<" "<< TrailingEdge[iPlane][2] << endl;

        }
      }

    }

    Nacelle_File.close();

    Section_File.open("nacelle_slices.dat", ios::out);

    for (iPlane = 0; iPlane < nPlane; iPlane++) {

      if (iPlane == 0) {
        Section_File << "TITLE = \"Nacelle Slices\"" << endl;
        if (config->GetSystemMeasurements() == US)
          Section_File << "VARIABLES = \"x (in)\", \"y (in)\", \"z (in)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
        else Section_File << "VARIABLES = \"x (m)\", \"y (m)\", \"z (m)\", \"x<sub>2D</sub>/c\", \"y<sub>2D</sub>/c\"" << endl;
      }

      if (Xcoord_Airfoil[iPlane].size() > 1) {

        su2double theta_deg = atan2(Plane_Normal[iPlane][1], -Plane_Normal[iPlane][2])/PI_NUMBER*180 + 180;
        su2double Angle = theta_deg*PI_NUMBER/180 - 0.5*PI_NUMBER;

        Section_File << "ZONE T=\"<greek>q</greek> = " << theta_deg << " deg\", I= " << Xcoord_Airfoil[iPlane].size() << ", F=POINT" << endl;

        for (iVertex = 0; iVertex < Xcoord_Airfoil[iPlane].size(); iVertex++) {

          /*--- Move to the origin  ---*/

          su2double XValue_ = Xcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][0];
          su2double ZValue_ = Zcoord_Airfoil[iPlane][iVertex] - LeadingEdge[iPlane][2];

          /*--- Rotate the airfoil and divide by the chord ---*/

          su2double ValCos = cos(Twist[iPlane]*PI_NUMBER/180.0);
          su2double ValSin = sin(Twist[iPlane]*PI_NUMBER/180.0);

          su2double XValue = (XValue_*ValCos - ZValue_*ValSin) / Chord[iPlane];
          su2double ZValue = (ZValue_*ValCos + XValue_*ValSin) / Chord[iPlane];

          su2double XCoord = Xcoord_Airfoil[iPlane][iVertex] + config->GetNacelleLocation(0);
          su2double YCoord = (Ycoord_Airfoil[iPlane][iVertex]*cos(Angle) - Zcoord_Airfoil[iPlane][iVertex]*sin(Angle)) + config->GetNacelleLocation(1);
          su2double ZCoord = (Zcoord_Airfoil[iPlane][iVertex]*cos(Angle) + Ycoord_Airfoil[iPlane][iVertex]*sin(Angle)) + config->GetNacelleLocation(2);

          /*--- Write the file ---*/

          Section_File  << XCoord << " " << YCoord << " " << ZCoord << " " << XValue  << " " << ZValue << endl;
        }
      }

    }

    Section_File.close();


    /*--- Compute the wing volume using a composite Simpson's rule ---*/

    Nacelle_Volume = 0.0;
    for (iPlane = 0; iPlane < nPlane-2; iPlane+=2) {
      if (Xcoord_Airfoil[iPlane].size() > 1) {
        Nacelle_Volume += (1.0/3.0)*dAngle*(Area[iPlane] + 4.0*Area[iPlane+1] + Area[iPlane+2]);
      }
    }

    /*--- Evaluate Max and Min quantities ---*/

    Nacelle_MaxMaxThickness = -1E6; Nacelle_MinMaxThickness = 1E6; Nacelle_MinChord = 1E6; Nacelle_MaxChord = -1E6;
    Nacelle_MinLERadius = 1E6; Nacelle_MaxLERadius = -1E6; Nacelle_MinToC = 1E6; Nacelle_MaxToC = -1E6;
    Nacelle_MaxTwist = -1E6;

    for (iPlane = 0; iPlane < nPlane; iPlane++) {
      if (MaxThickness[iPlane] != 0.0) Nacelle_MinMaxThickness = min(Nacelle_MinMaxThickness, MaxThickness[iPlane]);
      Nacelle_MaxMaxThickness = max(Nacelle_MaxMaxThickness, MaxThickness[iPlane]);
      if (Chord[iPlane] != 0.0) Nacelle_MinChord = min(Nacelle_MinChord, Chord[iPlane]);
      Nacelle_MaxChord = max(Nacelle_MaxChord, Chord[iPlane]);
      if (LERadius[iPlane] != 0.0) Nacelle_MinLERadius = min(Nacelle_MinLERadius, LERadius[iPlane]);
      Nacelle_MaxLERadius = max(Nacelle_MaxLERadius, LERadius[iPlane]);
      if (ToC[iPlane] != 0.0) Nacelle_MinToC = min(Nacelle_MinToC, ToC[iPlane]);
      Nacelle_MaxToC = max(Nacelle_MaxToC, ToC[iPlane]);
      Nacelle_ObjFun_MinToC = sqrt((Nacelle_MinToC - 0.07)*(Nacelle_MinToC - 0.07));
      Nacelle_MaxTwist = max(Nacelle_MaxTwist, fabs(Twist[iPlane]));
    }

  }

  /*--- Free memory for the section cuts ---*/

  delete [] Xcoord_Airfoil;
  delete [] Ycoord_Airfoil;
  delete [] Zcoord_Airfoil;
  delete [] Variable_Airfoil;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] LeadingEdge[iPlane];
  delete [] LeadingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] TrailingEdge[iPlane];
  delete [] TrailingEdge;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_P0[iPlane];
  delete [] Plane_P0;

  for (iPlane = 0; iPlane < nPlane; iPlane++)
    delete [] Plane_Normal[iPlane];
  delete [] Plane_Normal;

  delete [] Area;
  delete [] MaxThickness;
  delete [] Chord;
  delete [] LERadius;
  delete [] ToC;
  delete [] Twist;

}

std::unique_ptr<CADTElemClass> CPhysicalGeometry::ComputeViscousWallADT(const CConfig *config) const{

  /*--------------------------------------------------------------------------*/
  /*--- Step 1: Create the coordinates and connectivity of the linear      ---*/
  /*---         subelements of the local boundaries that must be taken     ---*/
  /*---         into account in the wall distance computation.             ---*/
  /*--------------------------------------------------------------------------*/

  /* Initialize an array for the mesh points, which eventually contains the
     mapping from the local nodes to the number used in the connectivity of the
     local boundary faces. However, in a first pass it is an indicator whether
     or not a mesh point is on a local wall boundary. */
  vector<unsigned long> meshToSurface(nPoint, 0);

  /* Define the vectors for the connectivity of the local linear subelements,
     the element ID's, the element type and marker ID's. */
  vector<unsigned long> surfaceConn;
  vector<unsigned long> elemIDs;
  vector<unsigned short> VTK_TypeElem;
  vector<unsigned short> markerIDs;

  /* Loop over the boundary markers. */

  for(unsigned short iMarker=0; iMarker<config->GetnMarker_All(); ++iMarker) {


    /* Check for a viscous wall. */
    if( config->GetViscous_Wall(iMarker)) {

      /* Loop over the surface elements of this marker. */
      for(unsigned long iElem=0; iElem < nElem_Bound[iMarker]; iElem++) {

        /* Set the flag of the mesh points on this surface to true. */
        for (unsigned short iNode = 0; iNode < bound[iMarker][iElem]->GetnNodes(); iNode++) {
          unsigned long iPoint = bound[iMarker][iElem]->GetNode(iNode);
          meshToSurface[iPoint] = 1;
        }
        /* Determine the necessary data from the corresponding standard face,
          such as the number of linear subfaces, the number of DOFs per
          linear subface and the corresponding local connectivity. */
        const unsigned short VTK_Type      = bound[iMarker][iElem]->GetVTK_Type();
        const unsigned short nDOFsPerElem  = bound[iMarker][iElem]->GetnNodes();

          /* Loop over the nodes of element and store the required data. */

        markerIDs.push_back(iMarker);
        VTK_TypeElem.push_back(VTK_Type);
        elemIDs.push_back(iElem);

        for (unsigned short iNode = 0; iNode < nDOFsPerElem; iNode++)
          surfaceConn.push_back(bound[iMarker][iElem]->GetNode(iNode));
      }
    }
  }

  /*--- Create the coordinates of the local points on the viscous surfaces and
        create the final version of the mapping from all volume points to the
        points on the viscous surfaces. ---*/
  vector<su2double> surfaceCoor;
  unsigned long nVertex_SolidWall = 0;

  for(unsigned long i=0; i<nPoint; ++i) {
    if( meshToSurface[i] ) {
      meshToSurface[i] = nVertex_SolidWall++;

      for(unsigned short k=0; k<nDim; ++k)
        surfaceCoor.push_back(node[i]->GetCoord(k));
    }
  }

  /*--- Change the surface connectivity, such that it corresponds to
        the entries in surfaceCoor rather than in meshPoints. ---*/
  for(unsigned long i=0; i<surfaceConn.size(); ++i)
    surfaceConn[i] = meshToSurface[surfaceConn[i]];

  /*--------------------------------------------------------------------------*/
  /*--- Step 2: Build the ADT, which is an ADT of bounding boxes of the    ---*/
  /*---         surface elements. A nearest point search does not give     ---*/
  /*---         accurate results, especially not for the integration       ---*/
  /*---         points of the elements close to a wall boundary.           ---*/
  /*--------------------------------------------------------------------------*/

  std::unique_ptr<CADTElemClass> WallADT(new CADTElemClass(nDim, surfaceCoor, surfaceConn, VTK_TypeElem,
                                                           markerIDs, elemIDs, true));

  return WallADT;

}

void CPhysicalGeometry::SetWallDistance(const CConfig *config, CADTElemClass *WallADT){


  /*--------------------------------------------------------------------------*/
  /*--- Step 3: Loop over all interior mesh nodes and compute minimum      ---*/
  /*---        distance to a solid wall element                           ---*/
  /*--------------------------------------------------------------------------*/

  if (!WallADT->IsEmpty()){
    /*--- Solid wall boundary nodes are present. Compute the wall
     distance for all nodes. ---*/

    for (unsigned long iPoint=0; iPoint<GetnPoint(); ++iPoint) {
      unsigned short markerID;
      unsigned long  elemID;
      int            rankID;
      su2double      dist;

      WallADT->DetermineNearestElement(node[iPoint]->GetCoord(), dist, markerID,
                                       elemID, rankID);
      if (dist < node[iPoint]->GetWall_Distance())
        node[iPoint]->SetWall_Distance(dist);
    }
  }
}
