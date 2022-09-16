/*!
 * \file CCGNSMeshReaderBase.cpp
 * \brief Helper class for the reading of a CGNS grid file.
 * \author T. Economon
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

#include "../../../include/geometry/meshreader/CCGNSMeshReaderBase.hpp"
#include "../../../include/toolboxes/CLinearPartitioner.hpp"

CCGNSMeshReaderBase::CCGNSMeshReaderBase(CConfig        *val_config,
                                         unsigned short val_iZone,
                                         unsigned short val_nZone)
: CMeshReaderFVM(val_config, val_iZone, val_nZone) {

#ifdef HAVE_CGNS
  /*--- We use val_iZone with +1 for the 1-based indexing in CGNS. ---*/
  cgnsZone = val_iZone + 1;
  nZones   = val_nZone;

#else
  SU2_MPI::Error(string(" SU2 built without CGNS support. \n") +
                 string(" To use CGNS, build SU2 accordingly."),
                 CURRENT_FUNCTION);
#endif
}

CCGNSMeshReaderBase::~CCGNSMeshReaderBase(void) { }

#ifdef HAVE_CGNS
void CCGNSMeshReaderBase::OpenCGNSFile(string val_filename) {
  
  /*--- For proper support of the high order elements, at least version 3.3
        of CGNS must be used. Check this. ---*/
#if CGNS_VERSION >= 3300

  /*--- Check whether the supplied file is truly a CGNS file. ---*/
  
  int file_type;
  if (cg_is_cgns(val_filename.c_str(), &file_type) != CG_OK) {
    SU2_MPI::Error(val_filename +
                   string(" was not found or is not a properly formatted") +
                   string(" CGNS file.\nNote that SU2 expects unstructured") +
                   string(" CGNS files in ADF data format."),
                   CURRENT_FUNCTION);
  }
  
  /*--- Open the CGNS file for reading. The value of cgnsFileID returned
   is the specific index number for this file and will be
   repeatedly used in the function calls. ---*/
  
  if (cg_open(val_filename.c_str(), CG_MODE_READ, &cgnsFileID))
    cg_error_exit();
  if (rank == MASTER_NODE) {
    cout << "Reading the CGNS file: ";
    cout << val_filename.c_str() << "." << endl;
  }

#else

  SU2_MPI::Error("CGNS version 3.3 or higher is necessary",
                 CURRENT_FUNCTION);

#endif
}

void CCGNSMeshReaderBase::ReadCGNSDatabaseMetadata() {
  
  /*--- Get the number of databases. This is the highest node
   in the CGNS heirarchy. ---*/
  
  int nbases;
  if (cg_nbases(cgnsFileID, &nbases)) cg_error_exit();
  if (rank == MASTER_NODE)
    cout << "CGNS file contains " << nbases << " database(s)." << endl;
  
  /*--- Check if there is more than one database. Throw an
   error if there is because this reader can currently
   only handle one database. ---*/
  
  if ( nbases > 1 ) {
    SU2_MPI::Error("CGNS reader currently can only handle 1 database.",
                   CURRENT_FUNCTION);
  }
  
  /*--- Read the database. Note that the CGNS indexing starts at 1. ---*/
  
  int cell_dim, phys_dim;
  char basename[CGNS_STRING_SIZE];
  if (cg_base_read(cgnsFileID, cgnsBase, basename, &cell_dim, &phys_dim))
    cg_error_exit();
  if (rank == MASTER_NODE) {
    cout << "Database " << cgnsBase << ", " << basename << ": ";
    cout << " cell dimension of " << cell_dim << ", physical ";
    cout << "dimension of " << phys_dim << "." << endl;
  }
  
  /*--- Set the number of dimensions baed on cell_dim. ---*/
  
  dimension = (unsigned short)cell_dim;
  
}

void CCGNSMeshReaderBase::ReadCGNSZoneMetadata() {
  
  /*--- First, check all sections to find the element types and to
   classify them as either surface or volume elements. We will also
   perform some error checks here to avoid partitioning issues. ---*/
  
  /*--- Get the number of zones for this base. ---*/
  
  int nzones;
  if (cg_nzones(cgnsFileID, cgnsBase, &nzones)) cg_error_exit();
  if (rank == MASTER_NODE) {
    cout <<  nzones << " total zone(s)." << endl;
  }
  
  /*--- Check if there is more than one zone. Until we enable it, we
   will require a single zone CGNS file. Multizone problems can still
   be run with CGNS by using separate CGNS files for each zone. ---*/
  
  if ( nzones > 1 ) {
    SU2_MPI::Error(string("CGNS reader currently expects only 1 zone per CGNS file.") +
                   string("Multizone problems can be run with separate CGNS files for each zone."), CURRENT_FUNCTION);
  }
  
  /*--- Read the basic information for this zone, including
   the name and the number of vertices, cells, and
   boundary cells which are stored in the cgsize variable. ---*/
  
  vector<cgsize_t> cgsize(3);
  ZoneType_t zonetype;
  char zonename[CGNS_STRING_SIZE];
  if (cg_zone_read(cgnsFileID, cgnsBase, cgnsZone, zonename, cgsize.data()))
    cg_error_exit();
  
  /*--- Rename the zone size information for clarity.
   NOTE: The number of cells here may be only the number of
   interior elements or it may be the total. This needs to
   be counted explicitly later. ---*/
  
  numberOfGlobalPoints = cgsize[0];
  int nElemCGNS        = cgsize[1];
  
  /*--- Get some additional information about the current zone. ---*/
  
  if (cg_zone_type(cgnsFileID, cgnsBase, cgnsZone, &zonetype)) cg_error_exit();
  
  /*--- Check for an unstructured mesh. Throw an error if not found. ---*/
  
  if (zonetype != Unstructured)
    SU2_MPI::Error("Structured CGNS zone found while unstructured expected.",
                   CURRENT_FUNCTION);
  
  /*--- Print current zone info to the console. ---*/
  
  if (rank == MASTER_NODE) {
    cout << "Zone " << cgnsZone << ", " << zonename << ": ";
    cout << numberOfGlobalPoints << " total vertices, ";
    cout << nElemCGNS << " total elements." << endl;
  }
  
  /*--- Retrieve the number of grids in this zone. For now, we know
   this is one, but to be more general, this will need to check and
   allow for a loop over all grids. ---*/
  
  int ngrids;
  if (cg_ngrids(cgnsFileID, cgnsBase, cgnsZone, &ngrids)) cg_error_exit();
  if (ngrids > 1) {
    SU2_MPI::Error("CGNS reader currently handles only 1 grid per zone.",
                   CURRENT_FUNCTION);
  }
  
}

void CCGNSMeshReaderBase::ReadCGNSSectionMetadata() {
  
  /*--- Begin section for retrieving the connectivity info. ---*/
  
  if ((rank == MASTER_NODE) && (size > SINGLE_NODE))
    cout << "Distributing connectivity across all ranks." << endl;
  
  /*--- First check the number of sections. ---*/
  
  if (cg_nsections(cgnsFileID, cgnsBase, cgnsZone, &nSections)) cg_error_exit();
  if (rank == MASTER_NODE) {
    cout << "Number of connectivity sections is ";
    cout << nSections << "." << endl;
  }
  
  /*--- Prepare several data structures to hold the various
   pieces of information describing each section. ---*/
  
  isInterior.resize(nSections);
  nElems.resize(nSections,0);
  elemOffset.resize(nSections+1, 0); elemOffset[0] = 0;
  connElems.resize(nSections);
  sectionNames.resize(nSections, vector<char>(CGNS_STRING_SIZE));
  numberOfGlobalElements = 0;

  for (int s = 0; s < nSections; s++) {
    
    /*--- Read the connectivity details for this section. ---*/
    
    int nbndry, parent_flag, vtk_type;
    cgsize_t startE, endE, sizeNeeded;
    ElementType_t elemType;
    if (cg_section_read(cgnsFileID, cgnsBase, cgnsZone, s+1,
                        sectionNames[s].data(), &elemType, &startE, &endE,
                        &nbndry, &parent_flag)) cg_error_exit();
    
    /*--- Compute the total element count in this section (global). ---*/
    
    unsigned long element_count = (endE-startE+1);
    
    /* Get the details for the CGNS element type in this section. */
    
    string elem_name = GetCGNSElementType(elemType, vtk_type);
    
    /* We assume that each section contains interior elements by default.
     If we find 1D elements in a 2D problem or 2D elements in a 3D
     problem, then we know the section must contain boundary elements.
     We assume that each section is composed of either entirely interior
     or entirely boundary elements. */
    
    isInterior[s] = true;
    
    if (elemType == MIXED) {
      
      /* For a mixed section, we check the type of the first element
       so that we can correctly label this section as an interior or
       boundary element section. Here, we also assume that a section
       can not hold both interior and boundary elements. First, get
       the size required to read a single element from the section. */
      
      if (cg_ElementPartialSize(cgnsFileID, cgnsBase, cgnsZone, s+1, startE,
                                startE, &sizeNeeded) != CG_OK)
        cg_error_exit();
      
      /* A couple of auxiliary vectors for mixed element sections. */
      
      vector<cgsize_t> connElemCGNS(sizeNeeded);
      vector<cgsize_t> connOffsetCGNS(2,0);
      
      /* Retrieve the connectivity information for the first element. */
      
      if (cg_poly_elements_partial_read(cgnsFileID, cgnsBase, cgnsZone, s+1,
                                        startE, startE, connElemCGNS.data(),
                                        connOffsetCGNS.data(), NULL) != CG_OK)
        cg_error_exit();
      
      /* The element type is in the first position of the connectivity
       information that we retrieved from the CGNS file. */
      
      elemType = ElementType_t(connElemCGNS[0]);
      
    }
    
    /* Check for 1D elements in 2D problems, or for 2D elements in
     3D problems. If found, mark the section as a boundary section. */

    if(dimension == 2) {
      switch( elemType ) {
        case BAR_2: case BAR_3: case BAR_4: case BAR_5:
          isInterior[s] = false;
          break;
        default:  // To avoid a compiler warning.
          break;
      }
    }
    else {
      switch( elemType ) {
        case TRI_3:  case TRI_6:   case TRI_9:   case TRI_10:
        case TRI_12: case TRI_15:
        case QUAD_4:  case QUAD_8:     case QUAD_9:  case QUAD_12:
        case QUAD_16: case QUAD_P4_16: case QUAD_25:
          isInterior[s] = false;
          break;
        default:  // To avoid a compiler warning.
          break;
      }
    }  

    /*--- Increment the global element offset for each section
     based on whether or not this is a surface or volume section.
     We also keep a running count of the total elements globally. ---*/
    
    elemOffset[s+1] = elemOffset[s];
    if (!isInterior[s]) elemOffset[s+1] += element_count;
    else numberOfGlobalElements += element_count;
    
    /*--- Print some information to the console. ---*/
    
    if (rank == MASTER_NODE) {
      cout << "Section " << string(sectionNames[s].data());
      cout << " contains " << element_count << " elements";
      cout << " of type " << elem_name << "." <<endl;
    }
    
  }
  
}

void CCGNSMeshReaderBase::ReadCGNSPointCoordinates() {
  
  /*--- Compute the number of points that will be on each processor.
   This is a linear partitioning with the addition of a simple load
   balancing for any remainder points. ---*/
  
  CLinearPartitioner pointPartitioner(numberOfGlobalPoints,0);
  
  /*--- Store the local number of nodes for this rank. ---*/
  
  numberOfLocalPoints = pointPartitioner.GetSizeOnRank(rank);
  
  /*--- Create buffer to hold the grid coordinates for our rank. ---*/
  
  localPointCoordinates.resize(dimension);
  for (int k = 0; k < dimension; k++)
    localPointCoordinates[k].resize(numberOfLocalPoints, 0.0);
  
  /*--- Set the value of range_max to the total number of nodes in
   the unstructured mesh. Also allocate memory for the temporary array
   that will hold the grid coordinates as they are extracted. Note the
   +1 for CGNS convention. ---*/
  
  cgsize_t range_min = (cgsize_t)pointPartitioner.GetFirstIndexOnRank(rank)+1;
  cgsize_t range_max = (cgsize_t)pointPartitioner.GetLastIndexOnRank(rank);
  
  /*--- Loop over each set of coordinates. ---*/
  
  for (int k = 0; k < dimension; k++) {
    
    /*--- Read the coordinate info. This will retrieve the
     data type (either RealSingle or RealDouble) as
     well as the coordname which will specify the
     type of data that it is based in the SIDS convention.
     This might be "CoordinateX," for instance. ---*/
    
    char coordname[CGNS_STRING_SIZE];
    DataType_t datatype;
    if (cg_coord_info(cgnsFileID, cgnsBase, cgnsZone, k+1,
                      &datatype, coordname)) cg_error_exit();
    if (rank == MASTER_NODE) {
      cout << "Loading " << coordname;
      if (size > SINGLE_NODE) {
        cout << " values into linear partitions." << endl;
      } else {
        cout << " values." << endl;
      }
    }
    
    /*--- Check the coordinate name to decide the index for storage. ---*/
    
    unsigned short indC = 0;
    if      (string(coordname) == "CoordinateX") indC = 0;
    else if (string(coordname) == "CoordinateY") indC = 1;
    else if (string(coordname) == "CoordinateZ") indC = 2;
    else
      SU2_MPI::Error(string("Unknown coordinate name, ") + coordname +
                     string(", in the CGNS file."), CURRENT_FUNCTION);
    
    /*--- Now read our rank's chunk of coordinates from the file.
     Ask for datatype RealDouble and let CGNS library do the translation
     when RealSingle is found. ---*/
    
    if (cg_coord_read(cgnsFileID, cgnsBase, cgnsZone, coordname, RealDouble,
                      &range_min, &range_max, localPointCoordinates[indC].data()))
      cg_error_exit();
  }
  
}

string CCGNSMeshReaderBase::GetCGNSElementType(ElementType_t val_elem_type,
                                               int           &val_vtk_type) {
  
  /* Check the CGNS element type and return the string name
   for the element and the associated VTK type index. */
  
  string elem_name;
  switch (val_elem_type) {
    case NODE:
      elem_name      = "Vertex";
      val_vtk_type   = 1;
      SU2_MPI::Error("Vertex elements detected. Please remove.",
                     CURRENT_FUNCTION);
      break;
    case BAR_2: case BAR_3: case BAR_4: case BAR_5:
      elem_name     = "Line";
      val_vtk_type  = 3;
      if (dimension == 3)
        SU2_MPI::Error("Line elements detected in a 3D mesh. Please remove.",
                       CURRENT_FUNCTION);
      break;
    case TRI_3:  case TRI_6:   case TRI_9:   case TRI_10:
    case TRI_12: case TRI_15:
      elem_name     = "Triangle";
      val_vtk_type  = 5;
      break;
    case QUAD_4:  case QUAD_8:     case QUAD_9:  case QUAD_12:
    case QUAD_16: case QUAD_P4_16: case QUAD_25:
      elem_name     = "Quadrilateral";
      val_vtk_type  = 9;
      break;
    case TETRA_4:  case TETRA_10: case TETRA_16: case TETRA_20:
    case TETRA_22: case TETRA_34: case TETRA_35:
      elem_name     = "Tetrahedron";
      val_vtk_type  = 10;
      break;
    case HEXA_8:  case HEXA_20:  case HEXA_27: case HEXA_32:
    case HEXA_56: case HEXA_64: case HEXA_44:  case HEXA_98:
    case HEXA_125:
      elem_name     = "Hexahedron";
      val_vtk_type  = 12;
      break;
    case PENTA_6:  case PENTA_15: case PENTA_18: case PENTA_24:
    case PENTA_38: case PENTA_40: case PENTA_33: case PENTA_66:
    case PENTA_75:
      elem_name     = "Prism";
      val_vtk_type  = 13;
      break;
    case PYRA_5:  case PYRA_14: case PYRA_13:    case PYRA_21:
    case PYRA_29: case PYRA_30: case PYRA_P4_29: case PYRA_50:
    case PYRA_55:
      elem_name     = "Pyramid";
      val_vtk_type  = 14;
      break;
    case MIXED:
      elem_name     = "Mixed";
      val_vtk_type  = -1;
      break;
    default:
      char buf[100];
      SPRINTF(buf, "Unsupported or unknown CGNS element type: (type %d)\n",
              val_elem_type);
      SU2_MPI::Error(string(buf), CURRENT_FUNCTION);
      break;
  }
  
  return elem_name;
  
}
#endif
