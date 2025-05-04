/*!
 * \file CCGNSElementType.hpp
 * \brief Header file for the class CCGNSElementType.
 *        The implementations are in the <i>CCGNSElementType.cpp</i> file.
 * \author E. van der Weide
 * \version 8.2.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include <vector>

#ifdef HAVE_CGNS
#include "cgnslib.h"
#endif

#ifdef HAVE_CGNS

using namespace std;

/*!
 * \class CCGNSElementType
 * \brief Class used to convert the CGNS format to SU2 format for high order elements.
 * \author: E. van der Weide
 */
class CCGNSElementType {
 private:
  vector<ElementType_t> CGNSTypeStored;            /*!< \brief CGNS element types for which the data is stored. */
  vector<unsigned short> VTKTypeStored;            /*!< \brief VTK type of the element. */
  vector<unsigned short> nPolyStored;              /*!< \brief Polynomial degree of the element. */
  vector<unsigned short> nDOFsStored;              /*!< \brief Number of DOFs of the element. */
  vector<vector<unsigned short> > SU2ToCGNSStored; /*!< \brief Double vector, which stores the conversion
                                                               from SU2 to CGNS for the type in local numbering. */
 public:
  /*--- Standard constructor, nothing to be done. ---*/
  CCGNSElementType() = default;

  /*--- Destructor, nothing to be done. ---*/
  ~CCGNSElementType() = default;

  /*!
   * \brief Converts the connectivity information from CGNS to SU2 format.
   * \param[in]  val_elemType  - CGNS elements type to be converted.
   * \param[in]  val_globalID  - Global ID of this element.
   * \param[in]  connCGNS      - Array with the connectivity of the element in CGNS format.
   * \param[out] connSU2       - Vector with the connectivity and meta data in SU2 format.
   */
  void CGNSToSU2(const ElementType_t val_elemType, const unsigned long val_globalID, const cgsize_t* connCGNS,
                 std::vector<unsigned long>& connSU2);

 private:
  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a node.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataNODE(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                      vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Bar2 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataBAR_2(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Bar3 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataBAR_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Bar4 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataBAR_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Bar5 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataBAR_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tri3 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTRI_3(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tri6 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTRI_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                       vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tri10 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTRI_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tri15 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTRI_15(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Quad4 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataQUAD_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Quad9 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataQUAD_9(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Quad16 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataQUAD_16(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Quad25 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataQUAD_25(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tetra4 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTETRA_4(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tetra10 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTETRA_10(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tetra20 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTETRA_20(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Tetra35 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataTETRA_35(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Pyra5 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPYRA_5(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Pyra14 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPYRA_14(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Pyra30 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPYRA_30(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Pyra55 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPYRA_55(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Penta6 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPENTA_6(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Penta18 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPENTA_18(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Penta40 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPENTA_40(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Penta75 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataPENTA_75(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Hexa8 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataHEXA_8(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                        vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Hexa27 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataHEXA_27(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Hexa64 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataHEXA_64(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                         vector<unsigned short>& SU2ToCGNS);

  /*!
   * \brief Converts the connectivity from CGNS to SU2 for a Hexa125 element.
   * \param[out] VTK_Type  - Corresponding VTK type
   * \param[out] nPoly     - Polynomial degree
   * \param[out] nDOFs     - Number of DOFs of the element.
   * \param[out] SU2ToCGNS - Vector containing the mapping from SU2 to CGNS.
   */
  void CreateDataHEXA_125(unsigned short& VTK_Type, unsigned short& nPoly, unsigned short& nDOFs,
                          vector<unsigned short>& SU2ToCGNS);
};
#endif
