/*!
 * \file CCGNSMeshReaderFEM.hpp
 * \brief Header file for the class CCGNSMeshReaderFEM.
 *        The implementations are in the <i>CCGNSMeshReaderFEM.cpp</i> file.
 * \author T. Economon
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

#include "CCGNSMeshReaderBase.hpp"

/*!
 * \class CCGNSMeshReaderFEM
 * \brief Reads a CGNS zone into linear partitions for the finite element solver (FEM).
 * \author: T. Economon
 */
class CCGNSMeshReaderFEM final : public CCGNSMeshReaderBase {
 private:
  /*!
   * \brief Communicates the grid points to the MPI rank where they are needed.
   */
  void CommPointCoordinates();

#ifdef HAVE_CGNS

  /*!
   * \brief Reads the connectivity range from a CGNS section and convert it to the internal format.
   * \param[in]    val_section    - CGNS section index.
   * \param[in]    val_firstIndex - Global index of the first element to be stored on this rank.
   * \param[in]    val_lastIndex  - Global index of the last element (not included) to be stored on this rank.
   * \param[inout] elemCount      - Counter, which keeps track how many global elements are stored.
   * \param[inout] localElemCount - Counter, which keeps track how many local elements are stored.
   * \param[inout] localConn      - Vector where the connectivity must be stored.
   */
  void ReadCGNSConnectivityRangeSection(const int val_section, const unsigned long val_firstIndex,
                                        const unsigned long val_lastIndex, unsigned long& elemCount,
                                        unsigned long& localElemCount, vector<unsigned long>& localConn);

  /*!
   * \brief Reads the interior volume elements from one section of a CGNS zone into linear partitions across all ranks.
   */
  void ReadCGNSVolumeElementConnectivity();

  /*!
   * \brief Reads the surface (boundary) elements from one section of a CGNS zone into linear partitions across all
   * ranks.
   */
  void ReadCGNSSurfaceElementConnectivity();

  /*!
   * \brief Reads the connectivity from a CGNS surface section and select the relevant faces.
   * \param[in]  val_section    - CGNS section index.
   * \param[in]  localFaces     - The faces of the locally stored volume elements.
   * \param[out] nSurfElem      - Number of local surface elements stored for this surface section.
   * \param[out] surfConn       - Vector to store the connectivity of the surface elements to be stored.
   */
  void ReadCGNSSurfaceSection(const int val_section, const vector<CFaceOfElement>& localFaces, unsigned long& nSurfElem,
                              vector<unsigned long>& surfConn);
#endif

 public:
  /*!
   * \brief Constructor of the CCGNSMeshReaderFEM class.
   */
  CCGNSMeshReaderFEM(const CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Destructor of the CCGNSMeshReaderFEM class.
   */
  ~CCGNSMeshReaderFEM(void) override;
};
