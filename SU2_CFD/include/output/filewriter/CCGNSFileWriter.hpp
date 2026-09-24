/*!
 * \file CCGNSFileWriter.hpp
 * \brief Headers for CGNS file writer class.
 * \author G. Baldan
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
#ifdef HAVE_MPI
#include "pcgnslib.h"
#endif
#endif

#include <cstdint>
#include <limits>

#include "CFileWriter.hpp"

class CConfig;
class CGeometry;
class CFVMDataSorter;

class CCGNSFileWriter final : public CFileWriter {
 private:
  const bool isSurface;             /*!< \brief True if surface file. */
  const bool doublePrecisionFields; /*!< \brief True to write the fields in double precision instead of single. */

  /*!
   * \brief Boundary elements of one marker owned by this rank, written as a boundary section of a volume file.
   */
  struct BoundaryMarker {
    string name;                /*!< \brief Marker tag. */
    unsigned short kindBC;      /*!< \brief SU2 boundary condition kind. */
    vector<unsigned long> conn; /*!< \brief VTK type followed by the (1-based) output ids of the nodes, per element. */
  };
  vector<BoundaryMarker> boundaryMarkers; /*!< \brief Markers written as boundaries of a volume file. */

  vector<string> surfaceMarkers; /*!< \brief Markers written as one zone each in a surface file. */
  CConfig* config = nullptr;     /*!< \brief Config, to sort the surface data of each marker. */
  CGeometry* geometry = nullptr; /*!< \brief Geometry, to sort the surface data of each marker. */

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

  cgsize_t cumulative; /*!< \brief Cumulative number of elements written. */

  /*--- Max connectivity entries per section, so that readers using 32-bit sizes can read it. ---*/
  static constexpr cgsize_t maxSectionEntries = std::numeric_limits<int32_t>::max();
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
  CCGNSFileWriter(CParallelDataSorter* valDataSorter, bool isSurf = false, bool doublePrecision = false);

  /*!
   * \brief Write sorted data to file in CGNS file format.
   * \param[in] val_filename - The name of the file.
   */
  void WriteData(string val_filename) override ;

  /*!
   * \brief Add the boundaries to a volume file: one boundary section, BC and family per marker, named as the marker.
   * \param[in] valConfig - Definition of the problem.
   * \param[in] valGeometry - Geometrical definition of the problem.
   * \param[in] volumeSorter - The volume data sorter, to find the boundary elements owned by this rank.
   */
  void SetBoundaryMarkers(CConfig* valConfig, CGeometry* valGeometry, const CFVMDataSorter* volumeSorter);

  /*!
   * \brief Write a surface file with one zone per plotted marker, named as the marker. The data of the surface
   *        sorter is sorted again for each marker when the file is written.
   * \param[in] valConfig - Definition of the problem.
   * \param[in] valGeometry - Geometrical definition of the problem.
   */
  void SetSurfaceMarkers(CConfig* valConfig, CGeometry* valGeometry);

 private:
#ifdef HAVE_CGNS
  /*!
   * \brief Create the CGNS file and its base.
   * \param[in] val_filename - The name of the file.
   */
  void InitializeMeshFile(const string& val_filename);

  /*!
   * \brief Write a zone with the data currently held by the data sorter.
   * \param[in] zoneName - Name of the zone.
   */
  void WriteZone(const string& zoneName);

  /*!
   * \brief Create a zone for the data currently held by the data sorter.
   * \param[in] zoneName - Name of the zone.
   */
  void InitializeZone(const string& zoneName);

  /*!
   * \brief Write the boundary sections, BCs and families of the markers set with SetBoundaryMarkers.
   */
  void WriteBoundaries();

  /*!
   * \brief Write i-th coordinate to file in CGNS file format. Coordinates are always written in double precision,
   *        the fields follow the precision requested by the user.
   * \param[in] iField - the output field ID.
   * \param[in] FieldName - Field name in the CGNS.
   */
  void WriteField(int iField, const string& FieldName);

  /*!
   * \brief Write i-th coordinate or field to file with the given data type.
   * \param[in] iField - the output field ID.
   * \param[in] FieldName - Field name in the CGNS.
   * \param[in] dataType - CGNS data type matching T, i.e. RealSingle for float and RealDouble for double.
   */
  template <class T>
  void WriteFieldOfType(int iField, const string& FieldName, DataType_t dataType);

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
   * \brief Create an element section, collectively when built with MPI.
   * \param[in] name - Name of the section.
   * \param[in] type - CGNS element type of the section.
   * \param[in] start - First element of the section.
   * \param[in] end - Last element of the section.
   * \param[out] section - Index of the section.
   * \returns CGNS error code.
   */
  int SectionWrite(const string& name, ElementType_t type, cgsize_t start, cgsize_t end, int* section);

  /*!
   * \brief Write the connectivity of the elements [start, end] of a section, the elements of this rank.
   * \param[in] section - Index of the section.
   * \param[in] start - First element written by this rank.
   * \param[in] end - Last element written by this rank.
   * \param[in] elements - Connectivity, nullptr if this rank writes no element of the section.
   * \returns CGNS error code.
   */
  int ElementsWriteData(int section, cgsize_t start, cgsize_t end, const cgsize_t* elements);

  /*!
   * \brief Call a generic CGNS function.
   * \param[in] ier - error value.
   */
  static inline void CallCGNS(const int& ier) {
    if (ier) cg_error_exit();
  }

  /*!
   * \brief Return the CGNS boundary condition type of an SU2 boundary condition kind.
   * \param[in] kindBC - SU2 boundary condition kind.
   */
  static BCType_t GetCGNSBCType(unsigned short kindBC);

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
