/*!
 * \file CCGNSFileWriter.hpp
 * \brief Headers fo paraview binary file writer class.
 * \author G. Baldan
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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
#include "cgnslib.h"
#endif

#include "CFileWriter.hpp"

class CCGNSFileWriter final : public CFileWriter {
 private:
  bool isSurface; /*!< \brief True if surface file. */

#ifdef HAVE_CGNS
  int cgnsFileID; /*!< \brief CGNS file identifier. */
  int cgnsBase;   /*!< \brief CGNS database index. */
  int cgnsZone;   /*!< \brief CGNS zone index. */
  int cgnsFields; /*!< \brief CGNS flow solution index. */

  int nZones;    /*!< \brief Total number of zones in the CGNS file. */
  int nSections; /*!< \brief Total number of sections in the CGNS file. */

  unsigned short nDim; /*!< \brief Problem dimension. */
  int nLocalPoints;    /*!< \brief Local number of points. */
  int GlobalPoint;     /*!< \brief Total number of points. */
  int GlobalElem;      /*!< \brief Total number of elements. */

  cgsize_t cumulative; /*!< \brief Cumulative number of elements written. */
#endif

 public:
  /*!
   * \brief File extension
   */
  const static string fileExt;

  /*!
   * \brief Construct a file writer using field names and the data sorter.
   * \param[in] valFileName - The name of the file.
   * \param[in] valDataSorter - The parallel sorted data to write.
   * \param[in] isSurf - True if it is a surface file.
   */
  CCGNSFileWriter(string valFileName, CParallelDataSorter* valDataSorter, bool isSurf = false);

  /*!
   * \brief Destructor
   */
  ~CCGNSFileWriter() override;

  /*!
   * \brief Write sorted data to file in CGNS file format.
   */
  void Write_Data() override;
#ifdef HAVE_CGNS
  /*!
   * \brief Initialize CGNS mesh file.
   */
  void initializeMeshFile();

  /*!
   * \brief Write i-th coordinate to file in CGNS file format.
   * \param[in] CoordinateNumber - Coordinate number according to CGNS.
   * \param[in] CoordinateName - Coordinate name according to CGNS.
   */
  void WriteCoordinate(int CoordinateNumber, const string& CoordinateName);

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
   * \brief Write i-th coordinate to file in CGNS file format.
   * \param[in] iField - the output field ID.
   * \param[in] FieldName - Field name in the CGNS.
   */
  void WriteField(int iField, const string& FieldName);

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
