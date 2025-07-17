/*!
 * \file CGraphPartitioning.hpp
 * \brief Headers for the classes realted to the algorithms that are used
                to divide the matrix acyclic graph into parallel partitions.
 * \author A. Raj
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

#include "../geometry/CGeometry.hpp"

/*!
 * \class CGraphPartitioning
 * \brief Abstract base class for defining graph partitioning algorithms.
 * \author A. Raj
 *
 * In order to use certain parallel algorithms in the solution process -
 * whether with linear solvers or preconditioners - we require the matrix
 * to be partitioned into certain parallel divisions. These maybe in the form
 * of levels, blocks, colors and so on. Since a number of different algorithms
 * can be used to split the graph, we've introduced a base class containing the
 * "Partition" member function from which child classes of the specific
 * algorithm can be derived. Currently, we are only using direct declarations
 * of the derived classes in the code. However, this method was chosen as it
 * allows us to pass different child class algorithms to a single implementation
 * of the function that requires it - similar to the CMatrixVectorProduct class.
 */

template <class ScalarType>

class CGraphPartitioning {
 public:
  virtual ~CGraphPartitioning() = 0;
  virtual void Partition(vector<ScalarType>& pointList, vector<ScalarType>& partitionOffsets,
                         vector<ScalarType>& chainPtr, unsigned short chainLimit) = 0;
};
template <class ScalarType>
CGraphPartitioning<ScalarType>::~CGraphPartitioning() {}

template <class ScalarType>

class CLevelScheduling final : public CGraphPartitioning<ScalarType> {
 private:
  ScalarType nPointDomain;
  CPoint* nodes;

 public:
  ScalarType nLevels;
  ScalarType maxLevelWidth;
  vector<ScalarType> levels;

  /*!
   * \brief constructor of the class
   * \param[in] nPointDomain_ref - number of points associated with the problem
   * \param[in] nodes_ref - represents the relationships between the points
   */
  inline CLevelScheduling<ScalarType>(ScalarType nPointDomain_ref, CPoint* nodes_ref)
      : nPointDomain(nPointDomain_ref), nodes(nodes_ref) {
    nLevels = 0ul;
    maxLevelWidth = 0ul;
  }

  CLevelScheduling() = delete;  // Removing default constructor

  /*!
   * \brief Divides the levels into groups of chains depending on the preset GPU block and warp size.
   * \param[in] levelOffsets - Represents the vector array containing the ordered list of starting rows of each level.
   * \param[in] chainPtr - Represents the vector array containing the ordered list of starting levels of each chain.
   * \param[in] rowsPerBlock - Represents the maximum number of rows that can be accomodated per CUDA block.
   */
  void CalculateChain(vector<ScalarType> levelOffsets, vector<ScalarType>& chainPtr, unsigned short rowsPerBlock) {
    ScalarType levelWidth = 0;

    /*This is not a magic number. We are simply initializing
    the point array with its first element that is always zero.*/
    chainPtr.push_back(0);

    for (ScalarType iLevel = 0ul; iLevel < nLevels; iLevel++) {
      levelWidth = levelOffsets[iLevel + 1] - levelOffsets[iLevel];
      maxLevelWidth = std::max(levelWidth, maxLevelWidth);

      if (levelWidth > rowsPerBlock) {
        if (chainPtr.back() != iLevel) {
          chainPtr.push_back(iLevel);
        }

        chainPtr.push_back(iLevel + 1);
      }
    }

    chainPtr.push_back(nLevels);
  }

  /*!
   * \brief Reorders the points according to the levels
   * \param[in] pointList - Ordered array that contains the list of all mesh points.
   * \param[in] inversePointList - Array utilized to access the index of each point in pointList.
   * \param[in] levelOffsets - Vector array containing the ordered list of starting rows of each level.
   */
  void Reorder(vector<ScalarType>& pointList, vector<ScalarType>& inversePointList, vector<ScalarType> levelOffsets) {
    for (auto localPoint = 0ul; localPoint < nPointDomain; ++localPoint) {
      const auto globalPoint = pointList[localPoint];
      inversePointList[levelOffsets[levels[localPoint]]++] = globalPoint;
    }

    pointList = std::move(inversePointList);
  }

  /*!
   * \brief Reorders the points according to the levels
   * \param[in] pointList - Ordered array that contains the list  of all mesh points.
   * \param[in] levelOffsets - Vector array containing the ordered list of starting rows of each level.
   * \param[in] chainPtr - Represents the vector array containing the ordered list of starting levels of each chain.
   * \param[in] rowsPerBlock - Represents the maximum number of rows that can be accomodated per CUDA block.
   */
  void Partition(vector<ScalarType>& pointList, vector<ScalarType>& levelOffsets, vector<ScalarType>& chainPtr,
                 unsigned short rowsPerBlock) override {
    vector<ScalarType> inversePointList;
    inversePointList.reserve(nPointDomain);
    levels.reserve(nPointDomain);

    for (auto point = 0ul; point < nPointDomain; point++) {
      inversePointList[pointList[point]] = point;
      levels[point] = 0;
    }

    //  Local Point - Ordering of the points post the RCM ordering
    //  Global Point - Original order of the points before the RCM ordering

    for (auto localPoint = 0ul; localPoint < nPointDomain; ++localPoint) {
      const auto globalPoint = pointList[localPoint];

      for (auto adjPoints = 0u; adjPoints < nodes->GetnPoint(globalPoint); adjPoints++) {
        const auto adjGlobalPoint = nodes->GetPoint(globalPoint, adjPoints);

        if (adjGlobalPoint < nPointDomain) {
          const auto adjLocalPoint = inversePointList[adjGlobalPoint];

          if (adjLocalPoint < localPoint) {
            levels[localPoint] = std::max(levels[localPoint], levels[adjLocalPoint] + 1);
          }
        }
      }

      nLevels = std::max(nLevels, levels[localPoint] + 1);
    }

    levelOffsets.resize(nLevels + 1);
    for (auto iPoint = 0ul; iPoint < nPointDomain; iPoint++) {
      ++levelOffsets[levels[iPoint] + 1];
    }

    for (auto iLevel = 2ul; iLevel <= nLevels; ++iLevel) {
      levelOffsets[iLevel] += levelOffsets[iLevel - 1];
    }

    Reorder(pointList, inversePointList, levelOffsets);

    CalculateChain(levelOffsets, chainPtr, rowsPerBlock);
  }
};
