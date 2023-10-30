/*!
 * \file fem_cgns_elements.cpp
 * \brief CGNS element definitions and conversions to the SU2 standard.
 * \author E. van der Weide
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

#ifdef HAVE_CGNS
#include "../../include/fem/fem_cgns_elements.hpp"
#include "../../include/fem/geometry_structure_fem_part.hpp"
#include "../../include/parallelization/mpi_structure.hpp"

#include <cmath>
#include <climits>
#include <algorithm>

using namespace std;

#if CGNS_VERSION >= 3300

void CCGNSElementType::DetermineMetaData(const unsigned short nDim, const int fn, const int iBase, const int iZone,
                                         const int iConn) {
  /* Store the connectivity ID. */
  connID = iConn;

  /* Read the element type and range from the CGNS file. */
  char cgnsname[CGNS_STRING_SIZE];
  int nBndry, parentFlag;
  if (cg_section_read(fn, iBase, iZone, iConn, cgnsname, &elemType, &indBeg, &indEnd, &nBndry, &parentFlag) != CG_OK)
    cg_error_exit();

  /* Store the name of this connectivity. */
  connName = cgnsname;

  /* Determine the number of elements present in this connectivity section. */
  nElem = indEnd - indBeg + 1;

  /* Determine the dimension of this element type, i.e. the number of
     parametric coordinates and conclude whether or not this is a volume or
     surface connectivity. Note that it is possible that it is neither of these,
     e.g. line elements for 3D. This information is not needed for the DG
     flow solver and is therefore ignored. */
  const unsigned short nDimElem = DetermineElementDimension(fn, iBase, iZone);
  volumeConn = (nDimElem == nDim);
  surfaceConn = (nDimElem == nDim - 1);
}

void CCGNSElementType::ReadBoundaryConnectivityRange(const int fn, const int iBase, const int iZone,
                                                     const unsigned long offsetRank, const unsigned long nBoundElemRank,
                                                     const unsigned long startingBoundElemIDRank,
                                                     unsigned long& locBoundElemCount,
                                                     vector<CBoundaryFace>& boundElems) {
  /* Determine the index range to be read for this rank. */
  const cgsize_t iBeg = indBeg + offsetRank;
  const cgsize_t iEnd = iBeg + nBoundElemRank - 1;

  /* Determine the size of the vector needed to read the connectivity
     data from the CGNS file. */
  cgsize_t sizeNeeded;
  if (cg_ElementPartialSize(fn, iBase, iZone, connID, iBeg, iEnd, &sizeNeeded) != CG_OK) cg_error_exit();

  /* Allocate the memory for the connectivity and read the data. */
  vector<cgsize_t> connCGNSVec(sizeNeeded);
  if (elemType == MIXED) {
    vector<cgsize_t> connCGNSOffsetVec(iEnd - iBeg + 2);
    if (cg_poly_elements_partial_read(fn, iBase, iZone, connID, iBeg, iEnd, connCGNSVec.data(),
                                      connCGNSOffsetVec.data(), nullptr) != CG_OK)
      cg_error_exit();
  } else {
    if (cg_elements_partial_read(fn, iBase, iZone, connID, iBeg, iEnd, connCGNSVec.data(), nullptr) != CG_OK)
      cg_error_exit();
  }

  /* Define the variables needed to convert the connectivities from CGNS to
     SU2 format. Note that the vectors are needed to support a connectivity
     section with mixed element types. */
  vector<ElementType_t> CGNS_Type;
  vector<unsigned short> VTK_Type;
  vector<unsigned short> nPoly;
  vector<unsigned short> nDOFs;

  vector<vector<unsigned short> > SU2ToCGNS;

  /* Definition of variables used in the loop below. */
  cgsize_t* connCGNS = connCGNSVec.data();
  ElementType_t typeElem = elemType;
  vector<unsigned long> connSU2;

  /* Loop over the elements just read. */
  for (unsigned long i = 0; i < nBoundElemRank; ++i, ++locBoundElemCount) {
    /* Determine the element type for this element if this is a mixed
       connectivity and set the pointer to the actual connectivity data. */
    if (elemType == MIXED) {
      typeElem = (ElementType_t)connCGNS[0];
      ++connCGNS;
    }

    /* Determine the index in the stored vectors (CGNS_Type, VTK_Type, etc),
       which corresponds to this element. If the type is not stored yet,
       a new entry in these vectors will be created. */
    const unsigned short ind = IndexInStoredTypes(typeElem, CGNS_Type, VTK_Type, nPoly, nDOFs, SU2ToCGNS);

    /* Resize the connSU2 vector to the appropriate size. */
    connSU2.resize(nDOFs[ind]);

    /* Create the connectivity used in SU2 by carrying out the renumbering.
       Note that in CGNS the numbering starts at 1, while in SU2 it starts
       at 0. This explains the addition of -1. */
    for (unsigned short j = 0; j < nDOFs[ind]; ++j) connSU2[j] = connCGNS[SU2ToCGNS[ind][j]] - 1;

    /* Set the pointer for connCGNS for the data of the next element. */
    connCGNS += nDOFs[ind];

    /* Determine the global boundary element ID of this element. */
    const unsigned long globBoundElemID = startingBoundElemIDRank + locBoundElemCount;

    /* Create an object of CBoundaryFace and store it in boundElems.
       Note that the corresponding domain element is not known yet and
       is therefore set to ULONG_MAX. */
    CBoundaryFace thisBoundFace;
    thisBoundFace.VTK_Type = VTK_Type[ind];
    thisBoundFace.nPolyGrid = nPoly[ind];
    thisBoundFace.nDOFsGrid = nDOFs[ind];
    thisBoundFace.globalBoundElemID = globBoundElemID;
    thisBoundFace.domainElementID = ULONG_MAX;
    thisBoundFace.Nodes = connSU2;

    boundElems.push_back(thisBoundFace);
  }
}

void CCGNSElementType::ReadConnectivityRange(const int fn, const int iBase, const int iZone,
                                             const unsigned long offsetRank, const unsigned long nElemRank,
                                             const unsigned long startingElemIDRank, CPrimalGrid**& elem,
                                             unsigned long& locElemCount, unsigned long& nDOFsLoc) {
  /* Determine the index range to be read for this rank. */
  const cgsize_t iBeg = indBeg + offsetRank;
  const cgsize_t iEnd = iBeg + nElemRank - 1;

  /* Determine the size of the vector needed to read the connectivity
     data from the CGNS file. */
  cgsize_t sizeNeeded;
  if (cg_ElementPartialSize(fn, iBase, iZone, connID, iBeg, iEnd, &sizeNeeded) != CG_OK) cg_error_exit();

  /* Allocate the memory for the connectivity and read the data. */
  vector<cgsize_t> connCGNSVec(sizeNeeded);
  if (elemType == MIXED) {
    vector<cgsize_t> connCGNSOffsetVec(iEnd - iBeg + 2);
    if (cg_poly_elements_partial_read(fn, iBase, iZone, connID, iBeg, iEnd, connCGNSVec.data(),
                                      connCGNSOffsetVec.data(), nullptr) != CG_OK)
      cg_error_exit();

  } else {
    if (cg_elements_partial_read(fn, iBase, iZone, connID, iBeg, iEnd, connCGNSVec.data(), nullptr) != CG_OK)
      cg_error_exit();
  }

  /* Define the variables needed to convert the connectivities from CGNS to
     SU2 format. Note that the vectors are needed to support a connectivity
     section with mixed element types. */
  vector<ElementType_t> CGNS_Type;
  vector<unsigned short> VTK_Type;
  vector<unsigned short> nPoly;
  vector<unsigned short> nDOFs;

  vector<vector<unsigned short> > SU2ToCGNS;

  /* Definition of variables used in the loop below. */
  cgsize_t* connCGNS = connCGNSVec.data();
  ElementType_t typeElem = elemType;
  vector<unsigned long> connSU2;

  /* Loop over the elements just read. */
  for (unsigned long i = 0; i < nElemRank; ++i, ++locElemCount) {
    /* Determine the element type for this element if this is a mixed
       connectivity and set the pointer to the actual connectivity data. */
    if (elemType == MIXED) {
      typeElem = (ElementType_t)connCGNS[0];
      ++connCGNS;
    }

    /* Determine the index in the stored vectors (CGNS_Type, VTK_Type, etc),
       which corresponds to this element. If the type is not stored yet,
       a new entry in these vectors will be created. */
    const unsigned short ind = IndexInStoredTypes(typeElem, CGNS_Type, VTK_Type, nPoly, nDOFs, SU2ToCGNS);

    /* Resize the connSU2 vector to the appropriate size. */
    connSU2.resize(nDOFs[ind]);

    /* Create the connectivity used in SU2 by carrying out the renumbering.
       Note that in CGNS the numbering starts at 1, while in SU2 it starts
       at 0. This explains the addition of -1. */
    for (unsigned short j = 0; j < nDOFs[ind]; ++j) connSU2[j] = connCGNS[SU2ToCGNS[ind][j]] - 1;

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
    elem[locElemCount] = new CPrimalGridFEM(globElemID, VTK_Type[ind], nPoly[ind], nPoly[ind], nDOFs[ind], nDOFs[ind],
                                            nDOFsLoc, connSU2.data());

    /* Update the counter nDOFsLoc. */
    nDOFsLoc += nDOFs[ind];
  }
}

unsigned short CCGNSElementType::DetermineElementDimension(const int fn, const int iBase, const int iZone) {
  /*--- Determine the element type and return the appropriate element
        dimension, i.e. the number of parametric coordinates needed to
        describe the element. ---*/
  switch (elemType) {
    case NODE:
      return 0;

    case BAR_2:
    case BAR_3:
    case BAR_4:
    case BAR_5:
      return 1;

    case TRI_3:
    case TRI_6:
    case TRI_9:
    case TRI_10:
    case TRI_12:
    case TRI_15:
    case QUAD_4:
    case QUAD_8:
    case QUAD_9:
    case QUAD_12:
    case QUAD_16:
    case QUAD_P4_16:
    case QUAD_25:
      return 2;

    case TETRA_4:
    case TETRA_10:
    case TETRA_16:
    case TETRA_20:
    case TETRA_22:
    case TETRA_34:
    case TETRA_35:
    case PYRA_5:
    case PYRA_14:
    case PYRA_13:
    case PYRA_21:
    case PYRA_29:
    case PYRA_30:
    case PYRA_P4_29:
    case PYRA_50:
    case PYRA_55:
    case PENTA_6:
    case PENTA_15:
    case PENTA_18:
    case PENTA_24:
    case PENTA_38:
    case PENTA_40:
    case PENTA_33:
    case PENTA_66:
    case PENTA_75:
    case HEXA_8:
    case HEXA_20:
    case HEXA_27:
    case HEXA_32:
    case HEXA_56:
    case HEXA_64:
    case HEXA_44:
    case HEXA_98:
    case HEXA_125:
      return 3;

    case MIXED:
      return DetermineElementDimensionMixed(fn, iBase, iZone);

    default:
      ostringstream message;
      message << "Unsupported CGNS element type, " << elemType << ", encountered.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }

  /* Just return a value to avoid a compiler warning. The program should
     never reach this position. */
  return 0;
}

unsigned short CCGNSElementType::DetermineElementDimensionMixed(const int fn, const int iBase, const int iZone) {
  /* Determine the required size of the vector to read the
     data of the first element of this connectivity. */
  cgsize_t sizeNeeded;
  if (cg_ElementPartialSize(fn, iBase, iZone, connID, indBeg, indBeg, &sizeNeeded) != CG_OK) cg_error_exit();

  /* Read the data of the first element in this section. */
  vector<cgsize_t> buf(sizeNeeded);
  vector<cgsize_t> buf_offset(2, 0);
  if (cg_poly_elements_partial_read(fn, iBase, iZone, connID, indBeg, indBeg, buf.data(), buf_offset.data(), nullptr) !=
      CG_OK)
    cg_error_exit();

  /* The first entry of buf contains the element type. Copy this value
     temporarily into the member variable elemType and determine the
     element dimension of this element. */
  elemType = (ElementType_t)buf[0];

  unsigned short nDimElem = DetermineElementDimension(fn, iBase, iZone);

  /* Reset the element type to MIXED and return the element dimension. */
  elemType = MIXED;
  return nDimElem;
}

unsigned short CCGNSElementType::IndexInStoredTypes(const ElementType_t typeElem, vector<ElementType_t>& CGNS_Type,
                                                    vector<unsigned short>& VTK_Type, vector<unsigned short>& nPoly,
                                                    vector<unsigned short>& nDOFs,
                                                    vector<vector<unsigned short> >& SU2ToCGNS) {
  /* Loop over the available types and check if the current type is present.
     If so, break the loop, such that the correct index is stored. */
  unsigned short ind;
  for (ind = 0; ind < CGNS_Type.size(); ++ind)
    if (typeElem == CGNS_Type[ind]) break;

  /* If the current element type is not present yet, the data for this
     type must be created. */
  if (ind == CGNS_Type.size()) {
    unsigned short VTKElem, nPolyElem, nDOFsElem;
    vector<unsigned short> SU2ToCGNSElem;
    CreateDataElementType(typeElem, VTKElem, nPolyElem, nDOFsElem, SU2ToCGNSElem);

    CGNS_Type.push_back(typeElem);
    VTK_Type.push_back(VTKElem);
    nPoly.push_back(nPolyElem);
    nDOFs.push_back(nDOFsElem);
    SU2ToCGNS.push_back(SU2ToCGNSElem);
  }

  /* Return the index. */
  return ind;
}

void CCGNSElementType::CreateDataElementType(const ElementType_t typeElem, unsigned short& VTK_Type,
                                             unsigned short& nPoly, unsigned short& nDOFs,
                                             vector<unsigned short>& SU2ToCGNS) {
  /*--- Determine the element type and call the corresponding function
        to determine the actual data. ---*/
  switch (typeElem) {
    case NODE:
      CreateDataNODE(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case BAR_2:
      CreateDataBAR_2(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case BAR_3:
      CreateDataBAR_3(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case BAR_4:
      CreateDataBAR_4(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case BAR_5:
      CreateDataBAR_5(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TRI_3:
      CreateDataTRI_3(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TRI_6:
      CreateDataTRI_6(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TRI_10:
      CreateDataTRI_10(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TRI_15:
      CreateDataTRI_15(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case QUAD_4:
      CreateDataQUAD_4(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case QUAD_9:
      CreateDataQUAD_9(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case QUAD_16:
      CreateDataQUAD_16(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case QUAD_25:
      CreateDataQUAD_25(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TETRA_4:
      CreateDataTETRA_4(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TETRA_10:
      CreateDataTETRA_10(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TETRA_20:
      CreateDataTETRA_20(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case TETRA_35:
      CreateDataTETRA_35(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PYRA_5:
      CreateDataPYRA_5(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PYRA_14:
      CreateDataPYRA_14(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PYRA_30:
      CreateDataPYRA_30(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PYRA_55:
      CreateDataPYRA_55(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PENTA_6:
      CreateDataPENTA_6(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PENTA_18:
      CreateDataPENTA_18(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PENTA_40:
      CreateDataPENTA_40(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case PENTA_75:
      CreateDataPENTA_75(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case HEXA_8:
      CreateDataHEXA_8(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case HEXA_27:
      CreateDataHEXA_27(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case HEXA_64:
      CreateDataHEXA_64(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;
    case HEXA_125:
      CreateDataHEXA_125(VTK_Type, nPoly, nDOFs, SU2ToCGNS);
      break;

    default:
      /* Print an error message that this element type is not supported and exit. */
      ostringstream message;
      message << "CGNS element type " << typeElem << " not supported.";
      SU2_MPI::Error(message.str(), CURRENT_FUNCTION);
  }
}

/*------------------------------------------------------------------------*/
/*--- Below this line there are only the conversion routines between   ---*/
/*--- CGNS and SU2 format for the specific elements.                   ---*/
/*------------------------------------------------------------------------*/

void CCGNSElementType::CreateDataNODE(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                      vector<unsigned short>& SU2ToCGNS) {
  /* Set the required data for a NODE. */
  VTK_Type = VERTEX;
  nPoly = 0;
  nDOFs = 1;
  SU2ToCGNS.push_back(0);
}

void CCGNSElementType::CreateDataBAR_2(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The BAR_2 element is a linear element. The numbering of the nodes is
     the same for CGNS and SU2. */
  VTK_Type = LINE;
  nPoly = 1;
  nDOFs = 2;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
}

void CCGNSElementType::CreateDataBAR_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The BAR_3 element is a quadratic element. SU2 numbers to nodes with
     increasing parametric value, while in CGNS the internal node is
     numbered last. */
  VTK_Type = LINE;
  nPoly = 2;
  nDOFs = 3;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 2;
  SU2ToCGNS[2] = 1;
}

void CCGNSElementType::CreateDataBAR_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The BAR_4 element is a cubic element. SU2 numbers to nodes with
     increasing parametric value, while in CGNS the internal nodes are
     numbered last. */
  VTK_Type = LINE;
  nPoly = 3;
  nDOFs = 4;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 2;
  SU2ToCGNS[2] = 3;
  SU2ToCGNS[3] = 1;
}

void CCGNSElementType::CreateDataBAR_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The BAR_5 element is a quartic element. SU2 numbers to nodes with
     increasing parametric value, while in CGNS the internal nodes are
     numbered last. */
  VTK_Type = LINE;
  nPoly = 4;
  nDOFs = 5;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 2;
  SU2ToCGNS[2] = 3;
  SU2ToCGNS[3] = 4;
  SU2ToCGNS[4] = 1;
}

void CCGNSElementType::CreateDataTRI_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The TRI_3 element is a linear triangle. The node numbering is the same
     in SU2 and CGNS. */
  VTK_Type = TRIANGLE;
  nPoly = 1;
  nDOFs = 3;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 2;
}

void CCGNSElementType::CreateDataTRI_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                       vector<unsigned short>& SU2ToCGNS) {
  /* The TRI_6 element is a quadratic triangle. In CGNS the nodes are numbered
     as follows: - First the vertex nodes.
                 - Second the edge nodes.
                 - Third the face nodes (not present for TRI_6).
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = TRIANGLE;
  nPoly = 2;
  nDOFs = 6;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 3;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 5;
  SU2ToCGNS[4] = 4;
  SU2ToCGNS[5] = 2;
}

void CCGNSElementType::CreateDataTRI_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The TRI_10 element is a cubic triangle. In CGNS the nodes are numbered
     as follows: - First the vertex nodes.
                 - Second the edge nodes.
                 - Third the face nodes.
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = TRIANGLE;
  nPoly = 3;
  nDOFs = 10;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 3;
  SU2ToCGNS[2] = 4;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 8;
  SU2ToCGNS[5] = 9;
  SU2ToCGNS[6] = 5;
  SU2ToCGNS[7] = 7;
  SU2ToCGNS[8] = 6;
  SU2ToCGNS[9] = 2;
}

void CCGNSElementType::CreateDataTRI_15(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The TRI_15 element is a quartic triangle. In CGNS the nodes are numbered
     as follows: - First the vertex nodes.
                 - Second the edge nodes.
                 - Third the face nodes.
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge.
     Also note that in the CGNS standard the face nodes are not placed at
     the location for uniform spacing. This effect is currently ignored,
     i.e. it is assumed that the spacing in parameter space is uniform. */
  VTK_Type = TRIANGLE;
  nPoly = 4;
  nDOFs = 15;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 3;
  SU2ToCGNS[2] = 4;
  SU2ToCGNS[3] = 5;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 11;
  SU2ToCGNS[6] = 12;
  SU2ToCGNS[7] = 13;
  SU2ToCGNS[8] = 6;
  SU2ToCGNS[9] = 10;
  SU2ToCGNS[10] = 14;
  SU2ToCGNS[11] = 7;
  SU2ToCGNS[12] = 9;
  SU2ToCGNS[13] = 8;
  SU2ToCGNS[14] = 2;
}

void CCGNSElementType::CreateDataQUAD_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The QUAD_4 element is a linear quadrilateral. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes (not present for QUAD_4).
                          - Third the face nodes (not present for QUAD_4).
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = QUADRILATERAL;
  nPoly = 1;
  nDOFs = 4;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 3;
  SU2ToCGNS[3] = 2;
}

void CCGNSElementType::CreateDataQUAD_9(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The QUAD_9 element is a quadratic quadrilateral. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = QUADRILATERAL;
  nPoly = 2;
  nDOFs = 9;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 7;
  SU2ToCGNS[4] = 8;
  SU2ToCGNS[5] = 5;
  SU2ToCGNS[6] = 3;
  SU2ToCGNS[7] = 6;
  SU2ToCGNS[8] = 2;
}

void CCGNSElementType::CreateDataQUAD_16(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The QUAD_16 element is a cubic quadrilateral. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = QUADRILATERAL;
  nPoly = 3;
  nDOFs = 16;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 5;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 11;
  SU2ToCGNS[5] = 12;
  SU2ToCGNS[6] = 13;
  SU2ToCGNS[7] = 6;
  SU2ToCGNS[8] = 10;
  SU2ToCGNS[9] = 15;
  SU2ToCGNS[10] = 14;
  SU2ToCGNS[11] = 7;
  SU2ToCGNS[12] = 3;
  SU2ToCGNS[13] = 9;
  SU2ToCGNS[14] = 8;
  SU2ToCGNS[15] = 2;
}

void CCGNSElementType::CreateDataQUAD_25(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The QUAD_25 element is a quartic quadrilateral. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
     In SU2 the node numbering takes place along increasing i- and j-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the other edge. */
  VTK_Type = QUADRILATERAL;
  nPoly = 4;
  nDOFs = 25;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 5;
  SU2ToCGNS[3] = 6;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 15;
  SU2ToCGNS[6] = 16;
  SU2ToCGNS[7] = 17;
  SU2ToCGNS[8] = 18;
  SU2ToCGNS[9] = 7;
  SU2ToCGNS[10] = 14;
  SU2ToCGNS[11] = 23;
  SU2ToCGNS[12] = 24;
  SU2ToCGNS[13] = 19;
  SU2ToCGNS[14] = 8;
  SU2ToCGNS[15] = 13;
  SU2ToCGNS[16] = 22;
  SU2ToCGNS[17] = 21;
  SU2ToCGNS[18] = 20;
  SU2ToCGNS[19] = 9;
  SU2ToCGNS[20] = 3;
  SU2ToCGNS[21] = 12;
  SU2ToCGNS[22] = 11;
  SU2ToCGNS[23] = 10;
  SU2ToCGNS[24] = 2;
}

void CCGNSElementType::CreateDataTETRA_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The TETRA_4 element is a linear tetrahedron. The node numbering is
     the same in SU2 and CGNS. */
  VTK_Type = TETRAHEDRON;
  nPoly = 1;
  nDOFs = 4;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 2;
  SU2ToCGNS[3] = 3;
}

void CCGNSElementType::CreateDataTETRA_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The TETRA_10 element is a quadratic tetrahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes (not present in TETRA_10).
                          - Fourth the volume nodes (not present in TETRA_10).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = TETRAHEDRON;
  nPoly = 2;
  nDOFs = 10;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 6;
  SU2ToCGNS[4] = 5;
  SU2ToCGNS[5] = 2;
  SU2ToCGNS[6] = 7;
  SU2ToCGNS[7] = 8;
  SU2ToCGNS[8] = 9;
  SU2ToCGNS[9] = 3;
}

void CCGNSElementType::CreateDataTETRA_20(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The TETRA_20 element is a cubic tetrahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes (not present in TETRA_20).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = TETRAHEDRON;
  nPoly = 3;
  nDOFs = 20;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 5;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 9;
  SU2ToCGNS[5] = 16;
  SU2ToCGNS[6] = 6;
  SU2ToCGNS[7] = 8;
  SU2ToCGNS[8] = 7;
  SU2ToCGNS[9] = 2;
  SU2ToCGNS[10] = 10;
  SU2ToCGNS[11] = 17;
  SU2ToCGNS[12] = 12;
  SU2ToCGNS[13] = 19;
  SU2ToCGNS[14] = 18;
  SU2ToCGNS[15] = 14;
  SU2ToCGNS[16] = 11;
  SU2ToCGNS[17] = 13;
  SU2ToCGNS[18] = 15;
  SU2ToCGNS[19] = 3;
}

void CCGNSElementType::CreateDataTETRA_35(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The TETRA_35 element is a quartic tetrahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. Also note that in the CGNS standard the face and volume
     nodes are not placed at the location for uniform spacing. This effect is
     currently ignored, i.e. it is assumed that the spacing in parameter space
     is uniform. */
  VTK_Type = TETRAHEDRON;
  nPoly = 4;
  nDOFs = 35;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 4;
  SU2ToCGNS[2] = 5;
  SU2ToCGNS[3] = 6;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 12;
  SU2ToCGNS[6] = 22;
  SU2ToCGNS[7] = 23;
  SU2ToCGNS[8] = 7;
  SU2ToCGNS[9] = 11;
  SU2ToCGNS[10] = 24;
  SU2ToCGNS[11] = 8;
  SU2ToCGNS[12] = 10;
  SU2ToCGNS[13] = 9;
  SU2ToCGNS[14] = 2;
  SU2ToCGNS[15] = 13;
  SU2ToCGNS[16] = 25;
  SU2ToCGNS[17] = 26;
  SU2ToCGNS[18] = 16;
  SU2ToCGNS[19] = 32;
  SU2ToCGNS[20] = 34;
  SU2ToCGNS[21] = 28;
  SU2ToCGNS[22] = 31;
  SU2ToCGNS[23] = 29;
  SU2ToCGNS[24] = 19;
  SU2ToCGNS[25] = 14;
  SU2ToCGNS[26] = 27;
  SU2ToCGNS[27] = 17;
  SU2ToCGNS[28] = 33;
  SU2ToCGNS[29] = 30;
  SU2ToCGNS[30] = 20;
  SU2ToCGNS[31] = 15;
  SU2ToCGNS[32] = 18;
  SU2ToCGNS[33] = 21;
  SU2ToCGNS[34] = 3;
}

void CCGNSElementType::CreateDataPYRA_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The PYRA_5 element is a linear pyramid. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes (not present in PYRA_5).
                          - Third the face nodes (not present in PYRA_5).
                          - Fourth the volume nodes (not present in PYRA_5).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = PYRAMID;
  nPoly = 1;
  nDOFs = 5;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 3;
  SU2ToCGNS[3] = 2;
  SU2ToCGNS[4] = 4;
}

void CCGNSElementType::CreateDataPYRA_14(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The PYRA_14 element is a quadratic pyramid. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes (not present in PYRA_14).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = PYRAMID;
  nPoly = 2;
  nDOFs = 14;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 5;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 8;
  SU2ToCGNS[4] = 13;
  SU2ToCGNS[5] = 6;
  SU2ToCGNS[6] = 3;
  SU2ToCGNS[7] = 7;
  SU2ToCGNS[8] = 2;
  SU2ToCGNS[9] = 9;
  SU2ToCGNS[10] = 10;
  SU2ToCGNS[11] = 12;
  SU2ToCGNS[12] = 11;
  SU2ToCGNS[13] = 4;
}

void CCGNSElementType::CreateDataPYRA_30(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The PYRA_30 element is a cubic pyramid. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = PYRAMID;
  nPoly = 3;
  nDOFs = 30;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 5;
  SU2ToCGNS[2] = 6;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 12;
  SU2ToCGNS[5] = 21;
  SU2ToCGNS[6] = 22;
  SU2ToCGNS[7] = 7;
  SU2ToCGNS[8] = 11;
  SU2ToCGNS[9] = 24;
  SU2ToCGNS[10] = 23;
  SU2ToCGNS[11] = 8;
  SU2ToCGNS[12] = 3;
  SU2ToCGNS[13] = 20;
  SU2ToCGNS[14] = 9;
  SU2ToCGNS[15] = 3;
  SU2ToCGNS[16] = 13;
  SU2ToCGNS[17] = 25;
  SU2ToCGNS[18] = 15;
  SU2ToCGNS[19] = 28;
  SU2ToCGNS[20] = 29;
  SU2ToCGNS[21] = 26;
  SU2ToCGNS[22] = 19;
  SU2ToCGNS[23] = 27;
  SU2ToCGNS[24] = 17;
  SU2ToCGNS[25] = 14;
  SU2ToCGNS[26] = 16;
  SU2ToCGNS[27] = 20;
  SU2ToCGNS[28] = 18;
  SU2ToCGNS[29] = 4;
}

void CCGNSElementType::CreateDataPYRA_55(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The PYRA_55 element is a quartic pyramid. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. Also note that in the CGNS standard the triangular face
     and volume nodes are not placed at the location for uniform spacing. This
     effect is currently ignored, i.e. it is assumed that the spacing in
     parameter space is uniform.*/
  VTK_Type = PYRAMID;
  nPoly = 4;
  nDOFs = 55;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 5;
  SU2ToCGNS[2] = 6;
  SU2ToCGNS[3] = 7;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 16;
  SU2ToCGNS[6] = 29;
  SU2ToCGNS[7] = 30;
  SU2ToCGNS[8] = 31;
  SU2ToCGNS[9] = 8;
  SU2ToCGNS[10] = 15;
  SU2ToCGNS[11] = 36;
  SU2ToCGNS[12] = 37;
  SU2ToCGNS[13] = 32;
  SU2ToCGNS[14] = 9;
  SU2ToCGNS[15] = 14;
  SU2ToCGNS[16] = 35;
  SU2ToCGNS[17] = 34;
  SU2ToCGNS[18] = 33;
  SU2ToCGNS[19] = 10;
  SU2ToCGNS[20] = 3;
  SU2ToCGNS[21] = 13;
  SU2ToCGNS[22] = 12;
  SU2ToCGNS[23] = 11;
  SU2ToCGNS[24] = 2;
  SU2ToCGNS[25] = 17;
  SU2ToCGNS[26] = 38;
  SU2ToCGNS[27] = 39;
  SU2ToCGNS[28] = 20;
  SU2ToCGNS[29] = 48;
  SU2ToCGNS[30] = 50;
  SU2ToCGNS[31] = 51;
  SU2ToCGNS[32] = 41;
  SU2ToCGNS[33] = 47;
  SU2ToCGNS[34] = 53;
  SU2ToCGNS[35] = 52;
  SU2ToCGNS[36] = 42;
  SU2ToCGNS[37] = 26;
  SU2ToCGNS[38] = 45;
  SU2ToCGNS[39] = 44;
  SU2ToCGNS[40] = 23;
  SU2ToCGNS[41] = 18;
  SU2ToCGNS[42] = 40;
  SU2ToCGNS[43] = 21;
  SU2ToCGNS[44] = 49;
  SU2ToCGNS[45] = 54;
  SU2ToCGNS[46] = 43;
  SU2ToCGNS[47] = 27;
  SU2ToCGNS[48] = 46;
  SU2ToCGNS[49] = 24;
  SU2ToCGNS[50] = 19;
  SU2ToCGNS[51] = 22;
  SU2ToCGNS[52] = 28;
  SU2ToCGNS[53] = 25;
  SU2ToCGNS[54] = 4;
}

void CCGNSElementType::CreateDataPENTA_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The PENTA_6 element is a linear prism. The node numbering is
     the same in SU2 and CGNS. */
  VTK_Type = PRISM;
  nPoly = 1;
  nDOFs = 6;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 2;
  SU2ToCGNS[3] = 3;
  SU2ToCGNS[4] = 4;
  SU2ToCGNS[5] = 5;
}

void CCGNSElementType::CreateDataPENTA_18(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The PENTA_18 element is a quadratic prism. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes (not present in PENTA_18).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = PRISM;
  nPoly = 2;
  nDOFs = 18;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 6;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 8;
  SU2ToCGNS[4] = 7;
  SU2ToCGNS[5] = 2;
  SU2ToCGNS[6] = 9;
  SU2ToCGNS[7] = 15;
  SU2ToCGNS[8] = 10;
  SU2ToCGNS[9] = 17;
  SU2ToCGNS[10] = 16;
  SU2ToCGNS[11] = 11;
  SU2ToCGNS[12] = 3;
  SU2ToCGNS[13] = 12;
  SU2ToCGNS[14] = 4;
  SU2ToCGNS[15] = 14;
  SU2ToCGNS[16] = 13;
  SU2ToCGNS[17] = 5;
}

void CCGNSElementType::CreateDataPENTA_40(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The PENTA_40 element is a cubic prism. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = PRISM;
  nPoly = 3;
  nDOFs = 40;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 6;
  SU2ToCGNS[2] = 7;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 11;
  SU2ToCGNS[5] = 24;
  SU2ToCGNS[6] = 8;
  SU2ToCGNS[7] = 10;
  SU2ToCGNS[8] = 9;
  SU2ToCGNS[9] = 2;
  SU2ToCGNS[10] = 12;
  SU2ToCGNS[11] = 25;
  SU2ToCGNS[12] = 26;
  SU2ToCGNS[13] = 14;
  SU2ToCGNS[14] = 34;
  SU2ToCGNS[15] = 38;
  SU2ToCGNS[16] = 29;
  SU2ToCGNS[17] = 33;
  SU2ToCGNS[18] = 30;
  SU2ToCGNS[19] = 16;
  SU2ToCGNS[20] = 13;
  SU2ToCGNS[21] = 28;
  SU2ToCGNS[22] = 27;
  SU2ToCGNS[23] = 15;
  SU2ToCGNS[24] = 35;
  SU2ToCGNS[25] = 39;
  SU2ToCGNS[26] = 32;
  SU2ToCGNS[27] = 36;
  SU2ToCGNS[28] = 31;
  SU2ToCGNS[29] = 17;
  SU2ToCGNS[30] = 3;
  SU2ToCGNS[31] = 18;
  SU2ToCGNS[32] = 19;
  SU2ToCGNS[33] = 4;
  SU2ToCGNS[34] = 23;
  SU2ToCGNS[35] = 37;
  SU2ToCGNS[36] = 20;
  SU2ToCGNS[37] = 22;
  SU2ToCGNS[38] = 21;
  SU2ToCGNS[39] = 5;
}

void CCGNSElementType::CreateDataPENTA_75(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The PENTA_75 element is a quartic prism. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. Also note that in the CGNS standard the triangular face
     and volume nodes are not placed at the location for uniform spacing. This
     effect is currently ignored, i.e. it is assumed that the spacing in
     parameter space is uniform. */
  VTK_Type = PRISM;
  nPoly = 4;
  nDOFs = 75;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 6;
  SU2ToCGNS[2] = 7;
  SU2ToCGNS[3] = 8;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 14;
  SU2ToCGNS[6] = 33;
  SU2ToCGNS[7] = 34;
  SU2ToCGNS[8] = 9;
  SU2ToCGNS[9] = 13;
  SU2ToCGNS[10] = 35;
  SU2ToCGNS[11] = 10;
  SU2ToCGNS[12] = 12;
  SU2ToCGNS[13] = 11;
  SU2ToCGNS[14] = 2;
  SU2ToCGNS[15] = 15;
  SU2ToCGNS[16] = 36;
  SU2ToCGNS[17] = 37;
  SU2ToCGNS[18] = 38;
  SU2ToCGNS[19] = 18;
  SU2ToCGNS[20] = 56;
  SU2ToCGNS[21] = 66;
  SU2ToCGNS[22] = 67;
  SU2ToCGNS[23] = 45;
  SU2ToCGNS[24] = 55;
  SU2ToCGNS[25] = 68;
  SU2ToCGNS[26] = 46;
  SU2ToCGNS[27] = 54;
  SU2ToCGNS[28] = 47;
  SU2ToCGNS[29] = 21;
  SU2ToCGNS[30] = 16;
  SU2ToCGNS[31] = 43;
  SU2ToCGNS[32] = 44;
  SU2ToCGNS[33] = 39;
  SU2ToCGNS[34] = 19;
  SU2ToCGNS[35] = 57;
  SU2ToCGNS[36] = 69;
  SU2ToCGNS[37] = 70;
  SU2ToCGNS[38] = 52;
  SU2ToCGNS[39] = 62;
  SU2ToCGNS[40] = 71;
  SU2ToCGNS[41] = 53;
  SU2ToCGNS[42] = 61;
  SU2ToCGNS[43] = 48;
  SU2ToCGNS[44] = 22;
  SU2ToCGNS[45] = 17;
  SU2ToCGNS[46] = 42;
  SU2ToCGNS[47] = 41;
  SU2ToCGNS[48] = 40;
  SU2ToCGNS[49] = 20;
  SU2ToCGNS[50] = 58;
  SU2ToCGNS[51] = 72;
  SU2ToCGNS[52] = 73;
  SU2ToCGNS[53] = 51;
  SU2ToCGNS[54] = 59;
  SU2ToCGNS[55] = 74;
  SU2ToCGNS[56] = 50;
  SU2ToCGNS[57] = 60;
  SU2ToCGNS[58] = 49;
  SU2ToCGNS[59] = 23;
  SU2ToCGNS[60] = 3;
  SU2ToCGNS[61] = 24;
  SU2ToCGNS[62] = 25;
  SU2ToCGNS[63] = 26;
  SU2ToCGNS[64] = 4;
  SU2ToCGNS[65] = 32;
  SU2ToCGNS[66] = 63;
  SU2ToCGNS[67] = 64;
  SU2ToCGNS[68] = 27;
  SU2ToCGNS[69] = 31;
  SU2ToCGNS[70] = 65;
  SU2ToCGNS[71] = 28;
  SU2ToCGNS[72] = 30;
  SU2ToCGNS[73] = 29;
  SU2ToCGNS[74] = 5;
}

void CCGNSElementType::CreateDataHEXA_8(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                        vector<unsigned short>& SU2ToCGNS) {
  /* The HEXA_8 element is a linear hexahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes (not present in HEXA_8).
                          - Third the face nodes (not present in HEXA_8).
                          - Fourth the volume nodes (not present in HEXA_8).
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = HEXAHEDRON;
  nPoly = 1;
  nDOFs = 8;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 1;
  SU2ToCGNS[2] = 3;
  SU2ToCGNS[3] = 2;
  SU2ToCGNS[4] = 4;
  SU2ToCGNS[5] = 5;
  SU2ToCGNS[6] = 7;
  SU2ToCGNS[7] = 6;
}

void CCGNSElementType::CreateDataHEXA_27(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The HEXA_27 element is a quadratic hexahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = HEXAHEDRON;
  nPoly = 2;
  nDOFs = 27;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 8;
  SU2ToCGNS[2] = 1;
  SU2ToCGNS[3] = 11;
  SU2ToCGNS[4] = 20;
  SU2ToCGNS[5] = 9;
  SU2ToCGNS[6] = 3;
  SU2ToCGNS[7] = 10;
  SU2ToCGNS[8] = 2;
  SU2ToCGNS[9] = 12;
  SU2ToCGNS[10] = 21;
  SU2ToCGNS[11] = 13;
  SU2ToCGNS[12] = 24;
  SU2ToCGNS[13] = 26;
  SU2ToCGNS[14] = 22;
  SU2ToCGNS[15] = 15;
  SU2ToCGNS[16] = 23;
  SU2ToCGNS[17] = 14;
  SU2ToCGNS[18] = 4;
  SU2ToCGNS[19] = 16;
  SU2ToCGNS[20] = 5;
  SU2ToCGNS[21] = 19;
  SU2ToCGNS[22] = 25;
  SU2ToCGNS[23] = 17;
  SU2ToCGNS[24] = 7;
  SU2ToCGNS[25] = 18;
  SU2ToCGNS[26] = 6;
}

void CCGNSElementType::CreateDataHEXA_64(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                         vector<unsigned short>& SU2ToCGNS) {
  /* The HEXA_64 element is a cubic hexahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = HEXAHEDRON;
  nPoly = 3;
  nDOFs = 64;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 8;
  SU2ToCGNS[2] = 9;
  SU2ToCGNS[3] = 1;
  SU2ToCGNS[4] = 15;
  SU2ToCGNS[5] = 32;
  SU2ToCGNS[6] = 33;
  SU2ToCGNS[7] = 10;
  SU2ToCGNS[8] = 14;
  SU2ToCGNS[9] = 35;
  SU2ToCGNS[10] = 34;
  SU2ToCGNS[11] = 11;
  SU2ToCGNS[12] = 3;
  SU2ToCGNS[13] = 13;
  SU2ToCGNS[14] = 12;
  SU2ToCGNS[15] = 2;
  SU2ToCGNS[16] = 16;
  SU2ToCGNS[17] = 36;
  SU2ToCGNS[18] = 37;
  SU2ToCGNS[19] = 18;
  SU2ToCGNS[20] = 49;
  SU2ToCGNS[21] = 56;
  SU2ToCGNS[22] = 57;
  SU2ToCGNS[23] = 40;
  SU2ToCGNS[24] = 48;
  SU2ToCGNS[25] = 59;
  SU2ToCGNS[26] = 58;
  SU2ToCGNS[27] = 41;
  SU2ToCGNS[28] = 22;
  SU2ToCGNS[29] = 45;
  SU2ToCGNS[30] = 44;
  SU2ToCGNS[31] = 20;
  SU2ToCGNS[32] = 17;
  SU2ToCGNS[33] = 39;
  SU2ToCGNS[34] = 38;
  SU2ToCGNS[35] = 19;
  SU2ToCGNS[36] = 50;
  SU2ToCGNS[37] = 60;
  SU2ToCGNS[38] = 61;
  SU2ToCGNS[39] = 43;
  SU2ToCGNS[40] = 51;
  SU2ToCGNS[41] = 63;
  SU2ToCGNS[42] = 62;
  SU2ToCGNS[43] = 42;
  SU2ToCGNS[44] = 23;
  SU2ToCGNS[45] = 46;
  SU2ToCGNS[46] = 47;
  SU2ToCGNS[47] = 21;
  SU2ToCGNS[48] = 4;
  SU2ToCGNS[49] = 24;
  SU2ToCGNS[50] = 25;
  SU2ToCGNS[51] = 5;
  SU2ToCGNS[52] = 31;
  SU2ToCGNS[53] = 52;
  SU2ToCGNS[54] = 53;
  SU2ToCGNS[55] = 26;
  SU2ToCGNS[56] = 30;
  SU2ToCGNS[57] = 55;
  SU2ToCGNS[58] = 54;
  SU2ToCGNS[59] = 27;
  SU2ToCGNS[60] = 7;
  SU2ToCGNS[61] = 29;
  SU2ToCGNS[62] = 28;
  SU2ToCGNS[63] = 6;
}

void CCGNSElementType::CreateDataHEXA_125(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                                          vector<unsigned short>& SU2ToCGNS) {
  /* The HEXA_125 element is a quartic hexahedron. In CGNS the nodes are
     numbered as follows: - First the vertex nodes.
                          - Second the edge nodes.
                          - Third the face nodes.
                          - Fourth the volume nodes.
     In SU2 the node numbering takes place along increasing i-, j- and k-lines,
     where the i-direction is defined from node 0 to 1 and the j-direction
     from node 0 along the second edge and the k-direction from node 0 along
     the third edge. */
  VTK_Type = HEXAHEDRON;
  nPoly = 4;
  nDOFs = 125;

  SU2ToCGNS.resize(nDOFs);
  SU2ToCGNS[0] = 0;
  SU2ToCGNS[1] = 8;
  SU2ToCGNS[2] = 9;
  SU2ToCGNS[3] = 10;
  SU2ToCGNS[4] = 1;
  SU2ToCGNS[5] = 19;
  SU2ToCGNS[6] = 44;
  SU2ToCGNS[7] = 45;
  SU2ToCGNS[8] = 46;
  SU2ToCGNS[9] = 11;
  SU2ToCGNS[10] = 18;
  SU2ToCGNS[11] = 51;
  SU2ToCGNS[12] = 52;
  SU2ToCGNS[13] = 47;
  SU2ToCGNS[14] = 12;
  SU2ToCGNS[15] = 17;
  SU2ToCGNS[16] = 50;
  SU2ToCGNS[17] = 49;
  SU2ToCGNS[18] = 48;
  SU2ToCGNS[19] = 13;
  SU2ToCGNS[20] = 3;
  SU2ToCGNS[21] = 16;
  SU2ToCGNS[22] = 15;
  SU2ToCGNS[23] = 14;
  SU2ToCGNS[24] = 2;
  SU2ToCGNS[25] = 20;
  SU2ToCGNS[26] = 53;
  SU2ToCGNS[27] = 54;
  SU2ToCGNS[28] = 55;
  SU2ToCGNS[29] = 23;
  SU2ToCGNS[30] = 82;
  SU2ToCGNS[31] = 98;
  SU2ToCGNS[32] = 99;
  SU2ToCGNS[33] = 100;
  SU2ToCGNS[34] = 62;
  SU2ToCGNS[35] = 81;
  SU2ToCGNS[36] = 105;
  SU2ToCGNS[37] = 106;
  SU2ToCGNS[38] = 101;
  SU2ToCGNS[39] = 63;
  SU2ToCGNS[40] = 80;
  SU2ToCGNS[41] = 104;
  SU2ToCGNS[42] = 103;
  SU2ToCGNS[43] = 102;
  SU2ToCGNS[44] = 64;
  SU2ToCGNS[45] = 29;
  SU2ToCGNS[46] = 73;
  SU2ToCGNS[47] = 72;
  SU2ToCGNS[48] = 71;
  SU2ToCGNS[49] = 26;
  SU2ToCGNS[50] = 21;
  SU2ToCGNS[51] = 60;
  SU2ToCGNS[52] = 61;
  SU2ToCGNS[53] = 56;
  SU2ToCGNS[54] = 24;
  SU2ToCGNS[55] = 83;
  SU2ToCGNS[56] = 107;
  SU2ToCGNS[57] = 108;
  SU2ToCGNS[58] = 109;
  SU2ToCGNS[59] = 69;
  SU2ToCGNS[60] = 88;
  SU2ToCGNS[61] = 114;
  SU2ToCGNS[62] = 115;
  SU2ToCGNS[63] = 110;
  SU2ToCGNS[64] = 70;
  SU2ToCGNS[65] = 87;
  SU2ToCGNS[66] = 113;
  SU2ToCGNS[67] = 112;
  SU2ToCGNS[68] = 111;
  SU2ToCGNS[69] = 65;
  SU2ToCGNS[70] = 30;
  SU2ToCGNS[71] = 74;
  SU2ToCGNS[72] = 79;
  SU2ToCGNS[73] = 78;
  SU2ToCGNS[74] = 27;
  SU2ToCGNS[75] = 22;
  SU2ToCGNS[76] = 59;
  SU2ToCGNS[77] = 58;
  SU2ToCGNS[78] = 57;
  SU2ToCGNS[79] = 25;
  SU2ToCGNS[80] = 84;
  SU2ToCGNS[81] = 116;
  SU2ToCGNS[82] = 117;
  SU2ToCGNS[83] = 118;
  SU2ToCGNS[84] = 68;
  SU2ToCGNS[85] = 85;
  SU2ToCGNS[86] = 123;
  SU2ToCGNS[87] = 124;
  SU2ToCGNS[88] = 119;
  SU2ToCGNS[89] = 67;
  SU2ToCGNS[90] = 86;
  SU2ToCGNS[91] = 122;
  SU2ToCGNS[92] = 121;
  SU2ToCGNS[93] = 120;
  SU2ToCGNS[94] = 66;
  SU2ToCGNS[95] = 31;
  SU2ToCGNS[96] = 75;
  SU2ToCGNS[97] = 76;
  SU2ToCGNS[98] = 77;
  SU2ToCGNS[99] = 28;
  SU2ToCGNS[100] = 4;
  SU2ToCGNS[101] = 32;
  SU2ToCGNS[102] = 33;
  SU2ToCGNS[103] = 34;
  SU2ToCGNS[104] = 5;
  SU2ToCGNS[105] = 43;
  SU2ToCGNS[106] = 89;
  SU2ToCGNS[107] = 90;
  SU2ToCGNS[108] = 91;
  SU2ToCGNS[109] = 35;
  SU2ToCGNS[110] = 42;
  SU2ToCGNS[111] = 96;
  SU2ToCGNS[112] = 97;
  SU2ToCGNS[113] = 92;
  SU2ToCGNS[114] = 36;
  SU2ToCGNS[115] = 41;
  SU2ToCGNS[116] = 95;
  SU2ToCGNS[117] = 94;
  SU2ToCGNS[118] = 93;
  SU2ToCGNS[119] = 37;
  SU2ToCGNS[120] = 7;
  SU2ToCGNS[121] = 40;
  SU2ToCGNS[122] = 39;
  SU2ToCGNS[123] = 38;
  SU2ToCGNS[124] = 6;
}

#endif
#endif
