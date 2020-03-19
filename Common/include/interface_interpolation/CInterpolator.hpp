/*!
 * \file CInterpolator.hpp
 * \brief Base class for multiphysics interpolation.
 * \author H. Kline
 * \version 7.0.2 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/datatype_structure.hpp"

class CConfig;
class CGeometry;

/*!
 * \class CInterpolator
 * \brief Main class for defining the interpolator, it requires
 * a child class for each particular interpolation method
 * \author H. Kline
 */
class CInterpolator {
protected:
  const int rank; 	         /*!< \brief MPI Rank. */
  const int size;       	   /*!< \brief MPI Size. */
  const unsigned donorZone;  /*!< \brief Index of donor zone. */
  const unsigned targetZone; /*!< \brief Index of target zone. */

  unsigned long
  MaxLocalVertex_Donor,      /*!< \brief Maximum vertices per processor*/
  nGlobalFace_Donor,         /*!< \brief Number of global donor faces*/
  nGlobalFaceNodes_Donor,    /*!< \brief Number of global donor face nodes*/
  MaxFace_Donor,             /*!< \brief Maximum faces per processor*/
  MaxFaceNodes_Donor;        /*!< \brief Maximum nodes associated with faces per processor*/

  unsigned long
  *Buffer_Receive_nVertex_Donor,     /*!< \brief Buffer to store the number of vertices per processor on the Donor domain */
  *Buffer_Receive_nFace_Donor,       /*!< \brief Buffer to store the number of faces per processor*/
  *Buffer_Receive_nFaceNodes_Donor,  /*!< \brief Buffer to store the number of nodes associated with faces per processor*/
  Buffer_Send_nVertex_Donor[1],      /*!< \brief Buffer to send number of vertices on the local processor*/
  Buffer_Send_nFace_Donor[1],        /*!< \brief Buffer to send number of faces on the local processor*/
  Buffer_Send_nFaceNodes_Donor[1],   /*!< \brief Buffer to send the number of nodes assocated with faces per processor*/
  *Buffer_Send_FaceIndex,            /*!< \brief Buffer to send indices pointing to the node indices that define the faces*/
  *Buffer_Receive_FaceIndex,         /*!< \brief Buffer to receive indices pointing to the node indices that define the faces*/
  *Buffer_Send_FaceNodes,            /*!< \brief Buffer to send indices pointing to the location of node information in other buffers, defining faces*/
  *Buffer_Receive_FaceNodes,         /*!< \brief Buffer to receive indices pointing to the location of node information in other buffers, defining faces*/
  *Buffer_Send_FaceProc,             /*!< \brief Buffer to send processor which stores the node indicated in Buffer_Receive_FaceNodes*/
  *Buffer_Receive_FaceProc;          /*!< \brief Buffer to receive processor which stores the node indicated in Buffer_Receive_FaceNodes*/

  long *Buffer_Send_GlobalPoint,     /*!< \brief Buffer to send global point indices*/
  *Buffer_Receive_GlobalPoint;       /*!< \brief Buffer to receive global point indices*/

  su2double *Buffer_Send_Coord,      /*!< \brief Buffer to send coordinate values*/
  *Buffer_Send_Normal,               /*!< \brief Buffer to send normal vector values */
  *Buffer_Receive_Coord,             /*!< \brief Buffer to receive coordinate values*/
  *Buffer_Receive_Normal;            /*!< \brief Buffer to receive normal vector values*/

  unsigned long
  *Receive_GlobalPoint,              /*!< \brief Buffer to receive Global point indexes*/
  *Buffer_Receive_nLinkedNodes,      /*!< \brief Buffer to receive the number of edges connected to each node*/
  *Buffer_Receive_LinkedNodes,       /*!< \brief Buffer to receive the list of notes connected to the nodes through an edge*/
  *Buffer_Receive_StartLinkedNodes,  /*!< \brief Buffer to receive the index of the Receive_LinkedNodes buffer where corresponding list of linked nodes begins */
  *Buffer_Receive_Proc;              /*!< \brief Buffer to receive the thread that owns the node*/

  unsigned long
  nGlobalVertex_Target,              /*!< \brief Global number of vertex of the target boundary*/
  nLocalVertex_Target,               /*!< \brief Number of vertex of the target boundary owned by the thread*/
  nGlobalVertex_Donor,               /*!< \brief Global number of vertex of the donor boundary*/
  nLocalVertex_Donor,                /*!< \brief Number of vertex of the donor boundary owned by the thread*/
  nGlobalVertex,                     /*!< \brief Dummy variable to temporarily store the global number of vertex of a boundary*/
  nLocalLinkedNodes;                 /*!< \brief Dummy variable to temporarily store the number of vertex of a boundary*/

  CGeometry**** const Geometry;      /*! \brief Vector which stores n zones of geometry. */
  CGeometry* const donor_geometry;   /*! \brief Donor geometry. */
  CGeometry* const target_geometry;  /*! \brief Target geometry. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - index of the donor zone
   * \param[in] jZone - index of the target zone
   */
  CInterpolator(CGeometry ****geometry_container, const CConfig* const* config, unsigned int iZone, unsigned int jZone);

  /*!
   * \brief No default construction allowed.
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
  virtual void Set_TransferCoeff(const CConfig* const* config) = 0;

  /*!
   * \brief Find the index of the interface marker shared by that zone
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker_interface - Interface tag.
   */
  static int Find_InterfaceMarker(const CConfig *config, unsigned short val_marker_interface);

  /*!
   * \brief Check whether an interface should be processed or not, i.e. if it is part of the zones.
   * \param[in] val_markDonor  - Marker tag from donor zone.
   * \param[in] val_markTarget - Marker tag from target zone.
   */
  static bool CheckInterfaceBoundary(int val_markDonor, int val_markTarget);

protected:
  /*!
   * \brief Recontstruct the boundary connectivity from parallel partitioning and broadcasts it to all threads
   * \param[in] val_zone   - index of the zone
   * \param[in] val_marker - index of the marker
   */
  void ReconstructBoundary(unsigned long val_zone, int val_marker);

  /*!
   * \brief compute squared distance between 2 points
   * \param[in] nDim - number of dimensions
   * \param[in] point_i - coordinates of point i
   * \param[in] point_j - coordinates of point j
   */
  static inline su2double PointsSquareDistance(unsigned short nDim, const su2double *point_i, const su2double *point_j) {
    su2double d = 0.0;
    for(unsigned short iDim = 0; iDim < nDim; iDim++)
      d += pow(point_j[iDim] - point_i[iDim], 2);
    return d;
  }

  /*!
   * \brief compute distance between 2 points
   * \param[in] nDim - number of dimensions
   * \param[in] point_i - coordinates of point i
   * \param[in] point_j - coordinates of point j
   */
  static inline su2double PointsDistance(unsigned short nDim, const su2double *point_i, const su2double *point_j) {
    return sqrt(PointsSquareDistance(nDim, point_i, point_j));
  }

  /*!
   * \brief Determine array sizes used to collect and send coordinate and global point
   * information.
   * \param[in] faces - boolean that determines whether or not to set face information as well
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

  /*!
   * \brief Collect and communicate vertex info: coord, global point, and if faces=true the normal vector
   * \param[in] faces - boolean that determines whether or not to set face information as well
   * \param[in] markDonor - Index of the boundary on the donor domain.
   * \param[in] markTarget - Index of the boundary on the target domain.
   * \param[in] nVertexDonor - Number of vertices on the donor boundary.
   * \param[in] nDim - number of physical dimensions.
   */
  void Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim);

};
