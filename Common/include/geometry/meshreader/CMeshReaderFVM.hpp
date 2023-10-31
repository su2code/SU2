
/*!
 * \file CMeshReaderFVM.hpp
 * \brief Header file for the class CMeshReaderFVM.
 *        The implementations are in the <i>CMeshReaderFVM.cpp</i> file.
 * \author T. Economon
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

#include <string.h>

#include "../../parallelization/mpi_structure.hpp"
#include "../../CConfig.hpp"

/*!
 * \class CMeshReaderFVM
 * \brief Base class for the mesh zone readers of the finite volume solver (FVM).
 * \author T. Economon
 */
class CMeshReaderFVM {
 protected:
  const int rank; /*!< \brief MPI Rank. */
  const int size; /*!< \brief MPI Size. */

  const CConfig* config = nullptr; /*!< \brief Local pointer to the config parameter object. */

  unsigned short dimension = 0; /*!< \brief Dimension of the problem (2 or 3). */

  unsigned long numberOfLocalPoints =
      0; /*!< \brief Number of local grid points within the linear partition on this rank. */
  unsigned long numberOfGlobalPoints = 0; /*!< \brief Number of global grid points within the mesh file. */
  vector<vector<passivedouble> >
      localPointCoordinates; /*!< \brief Vector holding the coordinates from the mesh file for the local grid points.
                                First index is dimension, second is point index. */

  unsigned long numberOfLocalElements =
      0; /*!< \brief Number of local elements within the linear partition on this rank. */
  unsigned long numberOfGlobalElements = 0;             /*!< \brief Number of global elements within the mesh file. */
  vector<unsigned long> localVolumeElementConnectivity; /*!< \brief Vector containing the element connectivity from the
                                                           mesh file for the local elements. */

  unsigned long numberOfMarkers = 0; /*!< \brief Total number of markers contained within the mesh file. */
  vector<string> markerNames;        /*!< \brief String names for all markers in the mesh file. */
  vector<vector<unsigned long> >
      surfaceElementConnectivity; /*!< \brief Vector containing the surface element connectivity from the mesh file on a
                                     per-marker basis. Only the master node reads and stores this connectivity. */

 public:
  /*!
   * \brief Constructor of the CMeshReaderFVM class.
   * \param[in] val_config - config object for the current zone.
   * \param[in] val_iZone  - Current zone index.
   * \param[in] val_nZone  - Total number of zones.
   */
  CMeshReaderFVM(const CConfig* val_config, unsigned short val_iZone, unsigned short val_nZone);

  /*!
   * \brief Get the physical dimension of the problem (2 or 3).
   * \returns Physical dimension of the problem.
   */
  inline unsigned short GetDimension() const { return dimension; }

  /*!
   * \brief Get the local point coordinates (linearly partitioned).
   * \returns Local point coordinates (linear partitioned).
   */
  inline const vector<vector<passivedouble> >& GetLocalPointCoordinates() const { return localPointCoordinates; }

  /*!
   * \brief Get the surface element connectivity for the specified marker. Only the master node owns the surface
   * connectivity. \param[in] val_iMarker - current marker index. \returns Surface element connecitivity for a marker
   * from the master rank.
   */
  inline const vector<unsigned long>& GetSurfaceElementConnectivityForMarker(int val_iMarker) const {
    return surfaceElementConnectivity[val_iMarker];
  }

  /*!
   * \brief Get the number surface elements for the specified marker.
   * \param[in] val_iMarker - current marker index.
   * \returns Number of surface elements for a marker.
   */
  inline unsigned long GetNumberOfSurfaceElementsForMarker(int val_iMarker) const {
    return (unsigned long)surfaceElementConnectivity[val_iMarker].size() / SU2_CONN_SIZE;
  }

  /*!
   * \brief Get the local volume element connectivity (linearly partitioned).
   * \returns Local volume element connectivity (linearly partitioned).
   */
  inline const vector<unsigned long>& GetLocalVolumeElementConnectivity() const {
    return localVolumeElementConnectivity;
  }

  /*!
   * \brief Get the total number of markers in the mesh zone.
   * \returns Total number of markers in the mesh zone.
   */
  inline unsigned long GetNumberOfMarkers() const { return numberOfMarkers; }

  /*!
   * \brief Get the vector of string names for all markers in the mesh zone.
   * \returns Vector of string names for all markers in the mesh zone.
   */
  inline const vector<string>& GetMarkerNames() const { return markerNames; }

  /*!
   * \brief Get the number of local grid points within the linear partition on this rank.
   * \returns Number of local grid points within the linear partition on this rank.
   */
  inline unsigned long GetNumberOfLocalPoints() const { return numberOfLocalPoints; }

  /*!
   * \brief Get the number of global grid points within the mesh file.
   * \returns Number of global grid points within the mesh file.
   */
  inline unsigned long GetNumberOfGlobalPoints() const { return numberOfGlobalPoints; }

  /*!
   * \brief Get the number of local elements within the linear partition on this rank.
   * \returns Number of local elements within the linear partition on this rank.
   */
  inline unsigned long GetNumberOfLocalElements() const { return numberOfLocalElements; }

  /*!
   * \brief Get the number of global elements within the mesh file.
   * \returns Number of global elements within the mesh file.
   */
  inline unsigned long GetNumberOfGlobalElements() const { return numberOfGlobalElements; }
};
