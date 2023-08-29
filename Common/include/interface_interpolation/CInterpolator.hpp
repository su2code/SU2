/*!
 * \file CInterpolator.hpp
 * \brief Base class for multiphysics interpolation.
 * \author H. Kline
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

#include "../../include/basic_types/datatype_structure.hpp"
#include "../../include/containers/C2DContainer.hpp"
#include "../../include/containers/container_decorators.hpp"
#include <vector>
#include <algorithm>

class CConfig;
class CGeometry;

using namespace std;

/*!
 * \class CInterpolator
 * \ingroup Interfaces
 * \brief Main class for defining the interpolator, it requires
 *        a child class for each particular interpolation method.
 * \author H. Kline
 */
class CInterpolator {
 protected:
  const int rank;            /*!< \brief MPI Rank. */
  const int size;            /*!< \brief MPI Size. */
  const unsigned donorZone;  /*!< \brief Index of donor zone. */
  const unsigned targetZone; /*!< \brief Index of target zone. */

  unsigned long MaxLocalVertex_Donor, /*!< \brief Maximum vertices per processor. */
      Buffer_Send_nVertex_Donor[1],   /*!< \brief Buffer to send number of vertices on the local processor. */
      *Buffer_Receive_nVertex_Donor;  /*!< \brief Buffer to store the number of vertices per processor on the Donor
                                         domain. */

  su2vector<unsigned long> Buffer_Send_GlobalPoint;    /*!< \brief Buffer to send global point indices. */
  su2vector<unsigned long> Buffer_Receive_GlobalPoint; /*!< \brief Buffer to receive global point indices. */

  su2activematrix Buffer_Send_Coord;    /*!< \brief Buffer to send coordinate values. */
  su2activematrix Buffer_Receive_Coord; /*!< \brief Buffer to receive coordinate values. */

  /*! \brief Buffer to receive the number of surface-connected edges, for each vertex. */
  su2vector<unsigned long> Buffer_Receive_nLinkedNodes;
  /*! \brief Buffer to receive the index of the Receive_LinkedNodes buffer where corresponding list of linked nodes
   * begins. */
  su2vector<unsigned long> Buffer_Receive_StartLinkedNodes;
  /*! \brief Buffer to receive the list of surface-connected nodes, for each vertex.
   * \details The vertices are ordered as in Buffer_Receive_nLinkedNodes and Buffer_Receive_StartLinkedNodes, but for
   * each*/
  su2vector<unsigned long> Buffer_Receive_LinkedNodes;
  /*! \brief Buffer to receive the rank that owns the vertex. */
  su2vector<unsigned long> Buffer_Receive_Proc;

  unsigned long nGlobalVertex_Target, /*!< \brief Global number of vertex of the target boundary. */
      nLocalVertex_Target,            /*!< \brief Number of vertex of the target boundary owned by the thread. */
      nGlobalVertex_Donor,            /*!< \brief Global number of vertex of the donor boundary. */
      nLocalVertex_Donor,             /*!< \brief Number of vertex of the donor boundary owned by the thread. */
      nGlobalVertex,     /*!< \brief Dummy variable to temporarily store the global number of vertex of a boundary. */
      nLocalLinkedNodes; /*!< \brief Dummy variable to temporarily store the number of vertex of a boundary. */

  CGeometry**** const Geometry;     /*! \brief Vector which stores n zones of geometry. */
  CGeometry* const donor_geometry;  /*! \brief Donor geometry. */
  CGeometry* const target_geometry; /*! \brief Target geometry. */

 public:
  struct CDonorInfo {
    vector<int> processor;
    vector<unsigned long> globalPoint;
    vector<su2double> coefficient;

    unsigned long nDonor() const { return processor.size(); }

    void resize(size_t nDonor) {
      processor.resize(nDonor);
      globalPoint.resize(nDonor);
      coefficient.resize(nDonor);
    }
  };
  vector<vector<CDonorInfo> > targetVertices; /*! \brief Donor information per marker per vertex of the target. */

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry_container - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone.
   * \param[in] jZone - index of the target zone.
   */
  CInterpolator(CGeometry**** geometry_container, const CConfig* const* config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief No default construction allowed to force zones and geometry to always be set.
   */
  CInterpolator(void) = delete;

  /*!
   * \brief Destructor of the class, nothing is deleted, derived classes need to manage the MPI buffers.
   */
  virtual ~CInterpolator(void) = default;

  /*!
   * \brief Set up transfer matrix defining relation between two meshes
   * \note Main method that derived classes must implement.
   * \param[in] config - Definition of the particular problem.
   */
  virtual void SetTransferCoeff(const CConfig* const* config) = 0;

  /*!
   * \brief Print information about the interpolation.
   */
  virtual void PrintStatistics(void) const {}

  /*!
   * \brief Check whether an interface should be processed or not, i.e. if it is part of the zones.
   * \param[in] val_markDonor  - Marker tag from donor zone.
   * \param[in] val_markTarget - Marker tag from target zone.
   */
  static bool CheckInterfaceBoundary(int val_markDonor, int val_markTarget);

  /*!
   * \brief Check whether two zones have a common interface.
   * \param[in] donor - Configuration of the donor zone.
   * \param[in] target - Configuration of the target zone.
   */
  static bool CheckZonesInterface(const CConfig* donor, const CConfig* target);

 protected:
  /*!
   * \brief Reconstruct the boundary connectivity from parallel partitioning and broadcasts it to all threads.
   * \param[in] val_zone   - index of the zone
   * \param[in] val_marker - index of the marker
   */
  void ReconstructBoundary(unsigned long val_zone, int val_marker);

  /*!
   * \brief Determine array sizes used to collect and send coordinate and global point information.
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Determine_ArraySize(int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

  /*!
   * \brief Collect and communicate vertex info: coord, global point.
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Collect_VertexInfo(int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

  /*!
   * \brief Collect all donor elements in an interface pair.
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] nDim - number of physical dimensions.
   * \param[in] compress - Squeeze the information (Allgatherv instead of Allgather).
   * \param[out] allNumElem - Number of donor element per rank.
   * \param[out] numNodes - Number of nodes for each element.
   * \param[out] idxNodes - Index (global) of those nodes.
   * \return Number of collected donor elements.
   * \note The last two outputs are always sized for Allgather.
   */
  unsigned long Collect_ElementInfo(int markDonor, unsigned short nDim, bool compress,
                                    vector<unsigned long>& allNumElem, vector<unsigned short>& numNodes,
                                    su2matrix<long>& idxNodes) const;
};
