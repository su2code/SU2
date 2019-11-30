/*!
 * \file graph_toolbox.hpp
 * \brief Functions and classes to build/represent sparse graphs or sparse patterns.
 * \author P. Gomes, F. Palacios, A. Bueno, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "C2DContainer.hpp"

#include <set>
#include <vector>
#include <cassert>

/*!
 * \enum ConnectivityType
 * \brief In FVM points are connected by the edges (faces) of the grid.
 *        In FEM, two points are connected if they have an element in common.
 */
enum class ConnectivityType {FiniteVolume=0, FiniteElement=1};


/*!
 * \class CCompressedSparsePattern
 * \brief A simple class to store adjacency information in a
 * compressed format suitable for sparse matrix operations.
 * If built for row-major storage the inner indices are column indices
 * and the pattern should be used as (row,icol), otherwise as (col,irow).
 */
template<typename Index_t>
class CCompressedSparsePattern {
  static_assert(std::is_integral<Index_t>::value,"");

private:
  su2vector<Index_t> m_outerPtr; /*!< \brief Start positions of the inner indices for each outer index. */
  su2vector<Index_t> m_innerIdx; /*!< \brief Inner indices of the non zero entries. */
  su2vector<Index_t> m_diagPtr;  /*!< \brief Position of the diagonal entry. */

public:
  using CEdgeToNonZeroMap =
    C2DContainer<unsigned long, Index_t, StorageType::RowMajor, 64, DynamicSize, 2>;

  /*!
   * \brief Construct from vector-like objects of any type with
   *        methods "size()" and "data()" (returning a pointer).
   * \param[in] outerPtr - Outer index pointers.
   * \param[in] innerIdx - Inner indices.
   */
  template<class T>
  CCompressedSparsePattern(const T& outerPtr, const T& innerIdx)
  {
    m_outerPtr.resize(outerPtr.size());
    for(Index_t i=0; i<outerPtr.size(); ++i)
      m_outerPtr(i) = outerPtr.data()[i];

    m_innerIdx.resize(innerIdx.size());
    for(Index_t i=0; i<innerIdx.size(); ++i)
      m_innerIdx(i) = innerIdx.data()[i];

    /*--- perform a basic sanity check ---*/
    assert(m_innerIdx.size() == m_outerPtr(m_outerPtr.size()-1));
  }

  /*!
   * \brief Build a list of pointers to the diagonal entries of the pattern.
   */
  void buildDiagPtr() {
    if(!m_diagPtr.empty()) return;

    m_diagPtr.resize(getOuterSize());
    for(Index_t k = 0; k < getOuterSize(); ++k)
      m_diagPtr(k) = findInnerIdx(k,k);
  }

  /*!
   * \return True if the pattern is empty, i.e. has not been built yet.
   */
  inline bool empty() const {
    return m_outerPtr.empty() || m_innerIdx.empty();
  }

  /*!
   * \return Number of rows/columns.
   */
  inline Index_t getOuterSize() const {
    return m_outerPtr.size()-1;
  }

  /*!
   * \return Number of non zero entries.
   */
  inline Index_t getNumNonZeros() const {
    return m_innerIdx.size();
  }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \return Number of inner indices associated with the outer index.
   */
  inline Index_t getNumNeighbors(Index_t iOuterIdx) const {
    return m_outerPtr(iOuterIdx+1) - m_outerPtr(iOuterIdx);
  }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \param[in] iNeigh - Relative position of the inner index.
   * \return The index of the iNeigh'th inner index associated with the outer index.
   */
  inline Index_t getInnerIdx(Index_t iOuterIdx, Index_t iNeigh) const {
    assert(iNeigh >= 0 && iNeigh < getNumNeighbors(iOuterIdx));
    return m_innerIdx(m_outerPtr(iOuterIdx) + iNeigh);
  }

  /*!
   * \param[in] iOuterIdx - Outer index (row/col).
   * \param[in] iInnerIdx - Inner index (col/row).
   * \return True if (iOuterIdx,iInnerIdx) exists, i.e. is non zero.
   */
  inline bool isNonZero(Index_t iOuterIdx, Index_t iInnerIdx) const {
    for(Index_t k = m_outerPtr(iOuterIdx); k < m_outerPtr(iOuterIdx+1); ++k)
      if(m_innerIdx(k) == iInnerIdx) return true;
    return false;
  }

  /*!
   * \param[in] iOuterIdx - Outer index (row/col).
   * \param[in] iInnerIdx - Inner index (col/row).
   * \return Absolute position of non zero entry (iOuterIdx,iInnerIdx).
   * \note Method segfaults if (iOuterIdx,iInnerIdx) is not a non zero.
   *       On DEBUG builds this requirement is asserted.
   */
  inline Index_t findInnerIdx(Index_t iOuterIdx, Index_t iInnerIdx) const {
    assert(isNonZero(iOuterIdx, iInnerIdx) && "Error, j does not belong to NZ(i).");
    Index_t k = m_outerPtr(iOuterIdx);
    while(m_innerIdx(k) != iInnerIdx) ++k;
    return k;
  }

  /*!
   * \param[in] iDiagIdx - Diagonal index (row == col).
   * \return Absolute position of the diagonal entry.
   */
  inline Index_t getDiagPtr(Index_t iDiagIdx) const {
    return m_diagPtr(iDiagIdx);
  }

  /*!
   * \return Raw pointer to the outer pointer vector.
   */
  inline const Index_t* outerPtr() const {
    assert(!empty() && "Sparse pattern has not been built.");
    return m_outerPtr.data();
  }

  /*!
   * \return Raw pointer to the inner index vector.
   */
  inline const Index_t* innerIdx() const {
    assert(!empty() && "Sparse pattern has not been built.");
    return m_innerIdx.data();
  }

  /*!
   * \return Raw pointer to the diagonal pointer vector.
   */
  inline const Index_t* diagPtr() const {
    assert(!m_diagPtr.empty() && "Diagonal map has not been built.");
    return m_diagPtr.data();
  }
};


/*!
 * \brief Build a sparse pattern from geometry information, of type FVM or FEM,
 *        for a given fill-level. At fill-level N, the immediate neighbors of the
 *        points in level N-1 are also considered neighbors of the base point.
 *        The resulting pattern is that of A^{N+1} where A is the sparse matrix
 *        of immediate neighbors.
 * \param[in] geometry - Definition of the grid.
 * \param[in] type - Of connectivity.
 * \param[in] fillLevel - Target degree of neighborhood (immediate neighbors always added).
 * \return Compressed-Storage-Row sparse pattern.
 */
template<class Geometry_t, typename Index_t>
CCompressedSparsePattern<Index_t> buildCSRPattern(Geometry_t& geometry,
                                                  ConnectivityType type,
                                                  Index_t fillLevel)
{
  Index_t nPoint = geometry.GetnPoint();

  std::vector<Index_t> outerPtr(nPoint+1);
  std::vector<Index_t> innerIdx;
  innerIdx.reserve(nPoint); // at least this much space is needed

  for(Index_t iPoint = 0; iPoint < nPoint; ++iPoint)
  {
    /*--- Inner indices for iPoint start here. ---*/
    outerPtr[iPoint] = innerIdx.size();

    /*--- Use a set to avoid duplication and keep ascending order. ---*/
    std::set<Index_t> neighbors;

    /*--- Insert base point. ---*/
    neighbors.insert(iPoint);

    /*--- Neighbors added in previous level. ---*/
    std::set<Index_t> addedNeighbors(neighbors);

    for(Index_t iLevel = 0; ; ++iLevel)
    {
      /*--- New points added in this level. ---*/
      std::set<Index_t> newNeighbors;

      /*--- For each point previously added, add its level 0
       *    neighbors, not duplicating any existing neighbor. ---*/
      for(auto jPoint : addedNeighbors)
      {
        auto point = geometry.node[jPoint];

        if(type == ConnectivityType::FiniteVolume)
        {
          /*--- For FVM we know the neighbors of point j directly. ---*/
          for(unsigned short iNeigh = 0; iNeigh < point->GetnPoint(); ++iNeigh)
          {
            Index_t kPoint = point->GetPoint(iNeigh);

            if(neighbors.count(kPoint) == 0) // no duplication
              newNeighbors.insert(kPoint);
          }
        }
        else // FiniteElement
        {
          /*--- For FEM we need the nodes of all elements that contain point j. ---*/
          for(unsigned short iNeigh = 0; iNeigh < point->GetnElem(); ++iNeigh)
          {
            auto elem = geometry.elem[point->GetElem(iNeigh)];

            for(unsigned short iNode = 0; iNode < elem->GetnNode(); ++iNode)
            {
              Index_t kPoint = elem->GetNode(iNode);

              if(neighbors.count(kPoint) == 0) // no duplication
                newNeighbors.insert(kPoint);
            }
          }
        }
      }

      neighbors.insert(newNeighbors.begin(), newNeighbors.end());

      if(iLevel >= fillLevel) break;

      /*--- For the next level we get the neighbours of the new points. ---*/
      addedNeighbors = newNeighbors;
    }

    /*--- Store final sparse pattern for iPoint. ---*/
    innerIdx.insert(innerIdx.end(), neighbors.begin(), neighbors.end());
  }
  outerPtr.back() = innerIdx.size();

  /*--- Return pattern as CCompressedSparsePattern object. ---*/
  return CCompressedSparsePattern<Index_t>(outerPtr, innerIdx);
}


/*!
 * \brief Build a lookup table of the absolute positions of the non zero entries
 *        of a compressed sparse pattern, accessed when visiting the FVM edges
 *        of a grid. The table can then be used for fast access (avoids searches)
 *        to the non zero entries of a sparse matrix associated with the pattern.
 * \param[in] geometry - Definition of the grid.
 * \param[in] pattern - Sparse pattern.
 * \return nEdge by 2 matrix.
 */
template<class Geometry_t, typename Index_t>
CCompressedSparsePattern<Index_t>::CEdgeToNonZeroMap mapEdgesToSparsePattern(
  Geometry_t& geometry, const CCompressedSparsePattern<Index_t>& pattern)
{
  assert(!pattern.empty());

  pattern::CEdgeToNonZeroMap edgeMap(geometry.GetnEdge(),2);

  for(Index_t iEdge = 0; iEdge < geometry.GetnEdge(); ++iEdge)
  {
    Index_t iPoint = geometry.edge[iEdge]->GetNode(0);
    Index_t jPoint = geometry.edge[iEdge]->GetNode(1);

    edgeMap(iEdge,0) = pattern.findInnerIdx(iPoint,jPoint);
    edgeMap(iEdge,1) = pattern.findInnerIdx(jPoint,iPoint);
  }

  return edgeMap;
}
