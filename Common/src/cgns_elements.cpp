/*!
 * \file cgns_elements.cpp
 * \brief CGNS element definitions and conversions to the SU2 standard.
 * \author E. van der Weide
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/geometry_structure.hpp"

#ifdef HAVE_CGNS
#if CGNS_VERSION >= 3300

void CGNSElementTypeClass::DetermineMetaData(const unsigned short nDim,
                                             const int            fn,
                                             const int            iBase,
                                             const int            iZone,
                                             const int            iConn) {

  /* Store the connectivity ID. */
  connID = iConn;

  /* Read the element type and range from the CGNS file. */
  char cgnsname[CGNS_STRING_SIZE];
  int nBndry, parentFlag;
  if(cg_section_read(fn, iBase, iZone, iConn, cgnsname, &elemType,
                     &indBeg, &indEnd, &nBndry, &parentFlag) != CG_OK) cg_error_exit();

  /* Determine the number of elements present in this connectivity section. */
  nElem = indEnd - indBeg + 1;

  /* Determine the dimension of this element type, i.e. the number of
     parametric coordinates and conclude whether or not this is a volume or
     surface connectivity. Note that it is possible that it is neither of these,
     e.g. line elements for 3D. This information is not needed for the DG
     flow solver and is therefore ignored. */
  const unsigned short nDimElem = DetermineElementDimension(fn, iBase, iZone);
  volumeConn  = (nDimElem == nDim);
  surfaceConn = (nDimElem == nDim-1);
}

void CGNSElementTypeClass::ReadConnectivityRange(const int           fn,
                                                 const int           iBase,
                                                 const int           iZone,
                                                 const unsigned long offsetRank,
                                                 const unsigned long nElemRank,
                                                 const unsigned long startingElemIDRank,
                                                 CPrimalGrid         **&elem,
                                                 unsigned long       &locElemCount,
                                                 unsigned long       &nDOFsLoc) {

  /* Determine the index range to be read for this rank. */
  const cgsize_t iBeg = indBeg + offsetRank;
  const cgsize_t iEnd = iBeg + nElemRank -1;

  /* Determine the size of the vector needed to read the connectivity
     data from the CGNS file. */
  cgsize_t sizeNeeded;
  if(cg_ElementPartialSize(fn, iBase, iZone, connID, iBeg, iEnd,
                           &sizeNeeded) != CG_OK) cg_error_exit();

  /* Allocate the memory for the connectivity and read the data. */
  vector<cgsize_t> connCGNSVec(sizeNeeded);
  if(cg_elements_partial_read(fn, iBase, iZone, connID, iBeg, iEnd,
                              connCGNSVec.data(), NULL) != CG_OK) cg_error_exit();

  /* Define the variables needed to convert the connectivities from CGNS to
     SU2 format. Note that the vectors are needed to support a connectivity
     section with mixed element types. */
  vector<ElementType_t>  CGNS_Type;
  vector<unsigned short> VTK_Type;
  vector<unsigned short> nPoly;
  vector<unsigned short> nDOFs;

  vector<vector<unsigned short> > SU2ToCGNS;

  /* Definition of variables used in the loop below. */
  cgsize_t *connCGNS = connCGNSVec.data();
  ElementType_t typeElem = elemType;
  vector<unsigned long> connSU2;

  /* Loop over the elements just read. */
  for(unsigned long i=0; i<nElemRank; ++i, ++locElemCount) {

    /* Determine the element type for this element if this is a mixed
       connectivity and set the pointer to the actual connectivity data. */
    if(elemType == MIXED) {
      typeElem = (ElementType_t) connCGNS[0];
      ++connCGNS;
    }

    /* Determine the index in the stored vectors (CGNS_Type, VTK_Type, etc),
       which corresponds to this element. If the type is not stored yet,
       a new entry in these vectors will be created. */
    const unsigned short ind = IndexInStoredTypes(typeElem, CGNS_Type, VTK_Type,
                                                  nPoly, nDOFs, SU2ToCGNS);

    /* Resize the connSU2 vector to the appropriate size. */
    connSU2.resize(nDOFs[ind]);

    /* Create the connectivity used in SU2 by carrying out the renumbering.
       Note that in CGNS the numbering starts at 1, while in SU2 it starts
       at 0. This explains the addition of -1. */
    for(unsigned short j=0; j<nDOFs[ind]; ++j)
      connSU2[j] = connCGNS[SU2ToCGNS[ind][j]] - 1;

    /* Set the pointer for connCGNS for the data of the next element. */
    connCGNS += nDOFs[ind];

    /* Determine the global element ID of this element. */
    const unsigned long globElemID = startingElemIDRank + locElemCount;

    /* Create a new FEM primal grid element, which stores the required
       information. Note that the second to last argument contains the local
       offset of the DOFs. This will be corrected to the global offset
       afterwards in CPhysicalGeometry::Read_CGNS_Format_Parallel_FEM.
       Also note that in CGNS it is not possible (at least at the moment)
       to specify a different polynomial degree for the grid and solution. */
    elem[locElemCount] = new CPrimalGridFEM(globElemID, VTK_Type[ind], nPoly[ind],
                                            nPoly[ind], nDOFs[ind], nDOFs[ind],
                                            nDOFsLoc, connSU2.data());

    /* Update the counter nDOFsLoc. */
    nDOFsLoc += nDOFs[ind];
  }
}

unsigned short CGNSElementTypeClass::DetermineElementDimension(const int fn,
                                                               const int iBase,
                                                               const int iZone) {

  /*--- Determine the element type and return the appropriate element
        dimension, i.e. the number of parametric coordinates needed to
        describe the element. ---*/
  switch( elemType ) {
    case NODE:
      return 0;

    case BAR_2: case BAR_3: case BAR_4: case BAR_5:
      return 1;

    case TRI_3:  case TRI_6:   case TRI_9:   case TRI_10:
    case TRI_12: case TRI_15:  case QUAD_4:  case QUAD_8:
    case QUAD_9: case QUAD_12: case QUAD_16: case QUAD_P4_16:
    case QUAD_25:
      return 2;

    case TETRA_4:  case TETRA_10:   case TETRA_16: case TETRA_20:
    case TETRA_22: case TETRA_34:   case TETRA_35: case PYRA_5:
    case PYRA_14:  case PYRA_13:    case PYRA_21:  case PYRA_29:
    case PYRA_30:  case PYRA_P4_29: case PYRA_50:  case PYRA_55:
    case PENTA_6:  case PENTA_15:   case PENTA_18: case PENTA_24:
    case PENTA_38: case PENTA_40:   case PENTA_33: case PENTA_66:
    case PENTA_75: case HEXA_8:     case HEXA_20:  case HEXA_27:
    case HEXA_32:  case HEXA_56:    case HEXA_64:  case HEXA_44:
    case HEXA_98:  case HEXA_125:
      return 3;

    case MIXED:
      return DetermineElementDimensionMixed(fn, iBase, iZone);

    default:
      int rank = MASTER_NODE;
#ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
      if(rank == MASTER_NODE) {
        cout << endl << endl << "   !!! Error !!!" << endl;
        cout << "Unsupported CGNS element type, " << elemType
             << ", encountered. " << endl;
        cout << " Now exiting..." << endl << endl;
      }

#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }

  /* Just return a value to avoid a compiler warning. The program should
     never reach this position. */
  return 0;
}

unsigned short CGNSElementTypeClass::DetermineElementDimensionMixed(const int fn,
                                                                    const int iBase,
                                                                    const int iZone) {
  /* Determine the required size of the vector to read the
     data of the first element of this connectivity. */
  cgsize_t sizeNeeded;
  if(cg_ElementPartialSize(fn, iBase, iZone, connID, indBeg, indBeg,
                           &sizeNeeded) != CG_OK) cg_error_exit();

  /* Read the data of the first element in this section. */
  vector<cgsize_t> buf(sizeNeeded);
  if(cg_elements_partial_read(fn, iBase, iZone, connID, indBeg, indBeg,
                              buf.data(), NULL) != CG_OK) cg_error_exit();

  /* The first entry of buf contains the element type. Copy this value
     temporarily into the member variable elemType and determine the
     element dimension of this element. */
  elemType = (ElementType_t) buf[0];

  unsigned short nDimElem = DetermineElementDimension(fn, iBase, iZone);

  /* Reset the element type to MIXED and return the element dimension. */
  elemType = MIXED;
  return nDimElem;
}

unsigned short CGNSElementTypeClass::IndexInStoredTypes(
                                  const ElementType_t             typeElem,
                                  vector<ElementType_t>           &CGNS_Type,
                                  vector<unsigned short>          &VTK_Type,
                                  vector<unsigned short>          &nPoly,
                                  vector<unsigned short>          &nDOFs,
                                  vector<vector<unsigned short> > &SU2ToCGNS) {

  /* Loop over the available types and check if the current type is present.
     If so, break the loop, such that the correct index is stored. */
  unsigned short ind;
  for(ind=0; ind<CGNS_Type.size(); ++ind)
    if(typeElem == CGNS_Type[ind]) break;

  /* If the current element type is not present yet, the data for this
     type must be created. */
  if(ind == CGNS_Type.size()) {

    unsigned short VTKElem, nPolyElem, nDOFsElem;
    vector<unsigned short> SU2ToCGNSElem;
    CreateDataElement(typeElem, VTKElem, nPolyElem, nDOFsElem, SU2ToCGNSElem);

    CGNS_Type.push_back(typeElem);
    VTK_Type.push_back(VTKElem);
    nPoly.push_back(nPolyElem);
    nDOFs.push_back(nDOFsElem);
    SU2ToCGNS.push_back(SU2ToCGNSElem);
  }

  /* Return the index. */
  return ind;
}

#endif
#endif
