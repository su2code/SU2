/*!
 * \file fem_cgns_elements.hpp
 * \brief Headers of the classes and functions for reading CGNS files
 *        with high order elements.
 *        The functions are in the <i>cgns_elements.cpp</i> file.
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

#pragma once

#include "../parallelization/mpi_structure.hpp"

#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif

#include "../geometry/primal_grid/CPrimalGridFEM.hpp"

/* Forward declaration of CBoundaryFace to avoid problems. */
class CBoundaryFace;

#ifdef HAVE_CGNS
#if CGNS_VERSION >= 3300

/*!
 * \class CCGNSElementType
 * \brief Class which stores the CGNS element type info for a connectivity section.
 */

class CCGNSElementType {
 public:
  int connID;             /*!< \brief CGNS connectivity ID of this connectivity. */
  ElementType_t elemType; /*!< \brief Element type according to the CGNS convention,
                                      possibly MIXED. */
  cgsize_t indBeg;        /*!< \brief Index of the first element in the CGNS connectivity. */
  cgsize_t indEnd;        /*!< \brief Index of the last element in the CGNS connectivity. */
  cgsize_t nElem;         /*!< \brief Number of elements present for this element type. */

  std::string connName; /*!< \brief Name of this connectivity. */

  bool volumeConn;  /*!< \brief Whether or not this is a volume connectivity. */
  bool surfaceConn; /*!< \brief Whether or not this is a surface connectivity. */

  /* Standard constructor, nothing to be done. */
  CCGNSElementType() {}

  /* Destructor, nothing to be done. */
  ~CCGNSElementType() {}

  /*--- Member function, which determines the meta data for this element type. ---*/
  void DetermineMetaData(const unsigned short nDim, const int fn, const int iBase, const int iZone, const int iConn);

  /*--- Member function, which reads the required boundary connectivity range. ---*/
  void ReadBoundaryConnectivityRange(const int fn, const int iBase, const int iZone, const unsigned long offsetRank,
                                     const unsigned long nBoundElemRank, const unsigned long startingBoundElemIDRank,
                                     unsigned long& locBoundElemCount, std::vector<CBoundaryFace>& boundElems);

  /*--- Member function, which reads the required connectivity range. ---*/
  void ReadConnectivityRange(const int fn, const int iBase, const int iZone, const unsigned long offsetRank,
                             const unsigned long nElemRank, const unsigned long startingElemIDRank, CPrimalGrid**& elem,
                             unsigned long& locElemCount, unsigned long& nDOFsLoc);

 private:
  /*--- Member function, which creates the required data for the given
        element type. ---*/
  void CreateDataElementType(const ElementType_t typeElem, unsigned short& VTK_Type, unsigned short& nPoly,
                             unsigned short& nDOFs, std::vector<unsigned short>& SU2ToCGNS);

  /*--- Member function, which determines the element dimension, i.e. the
        number of parametric coordinates. ---*/
  unsigned short DetermineElementDimension(const int fn, const int iBase, const int iZone);

  /*--- Member function, which determines the element dimension when the
        connectivity is mixed. ---*/
  unsigned short DetermineElementDimensionMixed(const int fn, const int iBase, const int iZone);

  /*--- Member function, which determines the corresponding index of the
        given element in the stored types. If not present, a new index
        is created. ---*/
  unsigned short IndexInStoredTypes(const ElementType_t typeElem, std::vector<ElementType_t>& CGNS_Type,
                                    std::vector<unsigned short>& VTK_Type, std::vector<unsigned short>& nPoly,
                                    std::vector<unsigned short>& nDOFs,
                                    std::vector<std::vector<unsigned short> >& SU2ToCGNS);

  /*--- Functions to create the conversion data from CGNS format to SU2 format
        for all the supported CGNS elements. ---*/
  void CreateDataNODE(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                      std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataBAR_2(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataBAR_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataBAR_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataBAR_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTRI_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTRI_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTRI_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTRI_15(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataQUAD_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataQUAD_9(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataQUAD_16(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataQUAD_25(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTETRA_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTETRA_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTETRA_20(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataTETRA_35(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPYRA_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPYRA_14(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPYRA_30(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPYRA_55(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPENTA_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPENTA_18(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPENTA_40(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataPENTA_75(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataHEXA_8(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataHEXA_27(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataHEXA_64(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         std::vector<unsigned short>& SU2ToCGNS);

  void CreateDataHEXA_125(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          std::vector<unsigned short>& SU2ToCGNS);
};
#endif
#endif

using namespace std;
