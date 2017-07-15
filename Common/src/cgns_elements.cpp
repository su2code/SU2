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
}

unsigned short CGNSElementTypeClass::DetermineElementDimensionMixed(const int fn,
                                                                    const int iBase,
                                                                    const int iZone) {
  /* Read the data of the first element in this section.
     Reserve enough memory for the worst case scenario. */
  cgsize_t buf[200];
  if(cg_elements_partial_read(fn, iBase, iZone, connID, indBeg, indBeg,
                              buf, NULL) != CG_OK) cg_error_exit();

  /* The first entry of buf contains the element type. Copy this value
     temporarily into the member variable elemType and determine the
     element dimension of this element. */
  elemType = (ElementType_t) buf[0];

  unsigned short nDimElem = DetermineElementDimension(fn, iBase, iZone);

  /* Reset the element type to MIXED and return the element dimension. */
  elemType = MIXED;
  return nDimElem;
}

#endif
#endif
