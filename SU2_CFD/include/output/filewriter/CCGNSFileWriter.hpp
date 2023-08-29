/*!
 * \file CCGNSFileWriter.hpp
 * \brief Headers for CGNS file writer class.
 * \author G. Baldan
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

#ifdef HAVE_CGNS
#ifdef __MINGW32__
#include <basetsd.h>
#endif
#include "cgnslib.h"
#endif

#include "CFileWriter.hpp"

class CCGNSFileWriter final : public CFileWriter {
 private:
  const bool isSurface; /*!< \brief True if surface file. */

#ifdef HAVE_CGNS
  int cgnsFileID; /*!< \brief CGNS file identifier. */
  int cgnsBase;   /*!< \brief CGNS database index. */
  int cgnsZone;   /*!< \brief CGNS zone index. */
  int cgnsFields; /*!< \brief CGNS flow solution index. */

  int nZones;    /*!< \brief Total number of zones in the CGNS file. */
  int nSections; /*!< \brief Total number of sections in the CGNS file. */

  unsigned short nDim;        /*!< \brief Problem dimension. */
  unsigned long nLocalPoints; /*!< \brief Local number of points. */
  cgsize_t GlobalPoint;       /*!< \brief Total number of points. */
  cgsize_t GlobalElem;        /*!< \brief Total number of elements. */

  typedef float dataPrecision;            /*!< \brief Define data precision of output (float or double). */
  const DataType_t dataType = RealSingle; /*!< \brief Datatype of fields can be RealSingle or RealDouble. */

  vector<cgsize_t> sendBufferConnectivity; /*!< \brief Send buffer for connectivity data. */
  vector<cgsize_t> recvBufferConnectivity; /*!< \brief Receive buffer for connectivity data. */
  vector<dataPrecision> recvBufferField;   /*!< \brief Send buffer for field data. */
  vector<dataPrecision> sendBufferField;   /*!< \brief Receive buffer for field data. */

  cgsize_t cumulative; /*!< \brief Cumulative number of elements written. */
#endif
 public:
  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valDataSorter - The parallel sorted data to write.
   * \param[in] isSurf - True if it is a surface file.
   */
  CCGNSFileWriter(CParallelDataSorter* valDataSorter, bool isSurf = false);

  /*!
   * \brief Write sorted data to file in CGNS file format.
   * \param[in] val_filename - The name of the file.
   */
  void WriteData(string val_filename) override ;

 private:
#ifdef HAVE_CGNS
  /*!
   * \brief Initialize CGNS mesh file.
   */
  void InitializeMeshFile(const string& val_filename);

  /*!
   * \brief Write i-th coordinate to file in CGNS file format.
   * \param[in] iField - the output field ID.
   * \param[in] FieldName - Field name in the CGNS.
   */
  void WriteField(int iField, const string& FieldName);

  /*!
   * \brief Write connectivity to file for GEO_TYPE in CGNS file format.
   * \param[in] type - GEO_TYPE.
   * \param[in] SectionName - Section name in the CGNS file.
   */
  void WriteConnectivity(GEO_TYPE type, const string& SectionName);

  /*!
   * \brief Initialize flow solution in the CGNS file.
   */
  void InitializeFields();

  /*!
   * \brief Call a generic CGNS function.
   * \param[in] ier - error value.
   */
  static inline void CallCGNS(const int& ier) {
    if (ier) cg_error_exit();
  }

  /*!
   * \brief Return the CGNS element type (ElementType_t).
   * \param[in] elementType - GEO_TYPE.
   */
  static inline ElementType_t GetCGNSType(unsigned short elementType) {
    switch (elementType) {
      case LINE:
        return BAR_2;
      case TRIANGLE:
        return TRI_3;
      case QUADRILATERAL:
        return QUAD_4;
      case TETRAHEDRON:
        return TETRA_4;
      case HEXAHEDRON:
        return HEXA_8;
      case PYRAMID:
        return PYRA_5;
      case PRISM:
        return PENTA_6;
      default:
        assert(false && "Invalid element type.");
        return ElementTypeNull;
    }
  }
#endif
};
