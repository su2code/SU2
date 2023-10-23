/*!
 * \file CADTPointsOnlyClass.hpp
 * \brief Class for storing an ADT of only points in an arbitrary number of dimensions.
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
#include "./CADTBaseClass.hpp"

/*!
 * \class CADTPointsOnlyClass
 * \ingroup ADT
 * \brief  Class for storing an ADT of only points in an arbitrary number of dimensions.
 * \author E. van der Weide
 */
class CADTPointsOnlyClass : public CADTBaseClass {
 private:
  vector<su2double> coorPoints;        /*!< \brief Vector, which contains the coordinates
                                                   of the points in the ADT. */
  vector<unsigned long> localPointIDs; /*!< \brief Vector, which contains the local point ID's
                                                   of the points in the ADT. */
  vector<int> ranksOfPoints;           /*!< \brief Vector, which contains the ranks
                                                   of the points in the ADT. */
 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] nDim       Number of spatial dimensions of the problem.
   * \param[in] nPoints    Number of local points to be stored in the ADT.
   * \param[in] coor       Coordinates of the local points.
   * \param[in] pointID    Local point IDs of the local points.
   * \param[in] globalTree Whether or not a global tree must be built. If false
                           a local ADT is built.
   */
  CADTPointsOnlyClass(unsigned short nDim, unsigned long nPoints, const su2double* coor, const unsigned long* pointID,
                      const bool globalTree);

  /*!
   * \brief Function, which determines the nearest node in the ADT for the given coordinate.
   * \note This simply forwards the call to the implementation function selecting the right
   *       working variables for the current thread.
   * \param[in]  coor    Coordinate for which the nearest node in the ADT must be determined.
   * \param[out] dist    Distance to the nearest node in the ADT.
   * \param[out] pointID Local point ID of the nearest node in the ADT.
   * \param[out] rankID  Rank on which the nearest node in the ADT is stored.
   */
  inline void DetermineNearestNode(const su2double* coor, su2double& dist, unsigned long& pointID, int& rankID) {
    const auto iThread = omp_get_thread_num();
    DetermineNearestNode_impl(FrontLeaves[iThread], FrontLeavesNew[iThread], coor, dist, pointID, rankID);
  }

  /*!
   * \brief Default constructor of the class, disabled.
   */
  CADTPointsOnlyClass() = delete;

 private:
  /*!
   * \brief Implementation of DetermineNearestNode.
   * \note Working variables (first two) passed explicitly for thread safety.
   */
  void DetermineNearestNode_impl(vector<unsigned long>& frontLeaves, vector<unsigned long>& frontLeavesNew,
                                 const su2double* coor, su2double& dist, unsigned long& pointID, int& rankID) const;
};
