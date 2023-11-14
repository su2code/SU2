/*!
 * \file graph_toolbox.hpp
 * \brief Functions and classes to build/represent sparse graphs or sparse patterns.
 * \author P. Gomes
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

#include "../containers/C2DContainer.hpp"
#include "../parallelization/omp_structure.hpp"

#include <set>
#include <vector>
#include <limits>
#include <cassert>
#include <algorithm>
#include <numeric>

/// \addtogroup Graph
/// @{

/*!
 * \enum ConnectivityType
 * \brief In FVM points are connected by the edges (faces) of the grid.
 *        In FEM, two points are connected if they have an element in common.
 */
enum class ConnectivityType { FiniteVolume = 0, FiniteElement = 1 };

/*!
 * \class CCompressedSparsePattern
 * \brief A simple class to store adjacency information in a
 * compressed format suitable for sparse matrix operations.
 * If built for row-major storage the inner indices are column indices
 * and the pattern should be used as (row,icol), otherwise as (col,irow).
 */
template <typename Index_t>
class CCompressedSparsePattern {
  static_assert(std::is_integral<Index_t>::value, "");

 private:
  su2vector<Index_t> m_outerPtr;       /*!< \brief Start positions of the inner indices for each outer index. */
  su2vector<Index_t> m_innerIdx;       /*!< \brief Inner indices of the non zero entries. */
  su2vector<Index_t> m_diagPtr;        /*!< \brief Position of the diagonal entry. */
  su2vector<Index_t> m_innerIdxTransp; /*!< \brief Position of the transpose non zero entries, requires symmetry. */

 public:
  using IndexType = Index_t;

  /*!
   * \brief Type to allow range for loops over inner indices.
   */
  struct CInnerIter {
    const IndexType* const m_first = nullptr;
    const IndexType* const m_last = nullptr;
    CInnerIter(const IndexType* first, const IndexType* last) : m_first(first), m_last(last) {}
    const IndexType* begin() const { return m_first; }
    const IndexType* end() const { return m_last; }
  };

  /*!
   * \brief Default construction.
   */
  CCompressedSparsePattern() = default;

  /*!
   * \brief Construct with default inner indices.
   * \param[in] outerPtrBegin - Start of outer pointers.
   * \param[in] outerPtrEnd - End of outer pointers.
   * \param[in] defaultInnerIdx - Default value for inner indices.
   */
  template <class Iterator>
  CCompressedSparsePattern(Iterator outerPtrBegin, Iterator outerPtrEnd, Index_t defaultInnerIdx) {
    const auto size = outerPtrEnd - outerPtrBegin;
    m_outerPtr.resize(size);
    Index_t k = 0;
    for (auto it = outerPtrBegin; it != outerPtrEnd; ++it) m_outerPtr(k++) = *it;

    m_innerIdx.resize(m_outerPtr(size - 1)) = defaultInnerIdx;
  }

  /*!
   * \brief Construct from rvalue refs.
   * \note This is the most efficient constructor as no data copy occurs.
   * \param[in] outerPtr - Outer index pointers.
   * \param[in] innerIdx - Inner indices.
   */
  CCompressedSparsePattern(su2vector<Index_t>&& outerPtr, su2vector<Index_t>&& innerIdx)
      : m_outerPtr(outerPtr), m_innerIdx(innerIdx) {
    /*--- perform a basic sanity check ---*/
    assert(m_innerIdx.size() == m_outerPtr(m_outerPtr.size() - 1));
  }

  /*!
   * \brief Construct from vector-like objects of any type with
   *        methods "size()" and "data()" (returning a pointer).
   * \param[in] outerPtr - Outer index pointers.
   * \param[in] innerIdx - Inner indices.
   */
  template <class T>
  CCompressedSparsePattern(const T& outerPtr, const T& innerIdx) {
    m_outerPtr.resize(outerPtr.size());
    for (Index_t i = 0; i < outerPtr.size(); ++i) m_outerPtr(i) = outerPtr.data()[i];

    m_innerIdx.resize(innerIdx.size());
    for (Index_t i = 0; i < innerIdx.size(); ++i) m_innerIdx(i) = innerIdx.data()[i];

    /*--- perform a basic sanity check ---*/
    assert(m_innerIdx.size() == m_outerPtr(m_outerPtr.size() - 1));
  }

  /*!
   * \brief Build from a "list of lists" type object
   * \param[in] lil - An object with operator [][] and method "size" e.g. vector<vector<?> >.
   */
  template <class T>
  CCompressedSparsePattern(const T& lil) {
    m_outerPtr.resize(lil.size() + 1);
    m_outerPtr(0) = 0;
    for (Index_t i = 1; i < Index_t(m_outerPtr.size()); ++i) m_outerPtr(i) = m_outerPtr(i - 1) + lil[i - 1].size();

    m_innerIdx.resize(m_outerPtr(lil.size()));
    Index_t k = 0;
    for (Index_t i = 0; i < Index_t(lil.size()); ++i)
      for (Index_t j = 0; j < Index_t(lil[i].size()); ++j) m_innerIdx(k++) = lil[i][j];
  }

  /*!
   * \brief Build a list of pointers to the diagonal entries of the pattern.
   */
  void buildDiagPtr() {
    if (!m_diagPtr.empty()) return;

    m_diagPtr.resize(getOuterSize());

    SU2_OMP_PARALLEL_(for schedule(static,roundUpDiv(getOuterSize(),omp_get_max_threads())))
    for (Index_t k = 0; k < getOuterSize(); ++k) m_diagPtr(k) = findInnerIdx(k, k);
    END_SU2_OMP_PARALLEL
  }

  /*!
   * \brief Build a list of pointers to the transpose entries of the pattern, requires symmetry.
   */
  void buildTransposePtr() {
    if (!m_innerIdxTransp.empty()) return;

    m_innerIdxTransp.resize(getNumNonZeros());

    SU2_OMP_PARALLEL_(for schedule(static,roundUpDiv(getOuterSize(),omp_get_max_threads())))
    for (Index_t i = 0; i < getOuterSize(); ++i) {
      for (Index_t k = m_outerPtr(i); k < m_outerPtr(i + 1); ++k) {
        auto j = m_innerIdx(k);
        m_innerIdxTransp(k) = findInnerIdx(j, i);
        assert(m_innerIdxTransp(k) != m_innerIdx.size() && "The pattern is not symmetric.");
      }
    }
    END_SU2_OMP_PARALLEL
  }

  /*!
   * \return True if the pattern is empty, i.e. has not been built yet.
   */
  inline bool empty() const { return m_outerPtr.empty() || m_innerIdx.empty(); }

  /*!
   * \return Number of rows/columns.
   */
  inline Index_t getOuterSize() const { return m_outerPtr.size() - 1; }

  /*!
   * \return Number of non zero entries.
   */
  inline Index_t getNumNonZeros() const { return m_innerIdx.size(); }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \return Number of inner indices associated with the outer index.
   */
  inline Index_t getNumNonZeros(Index_t iOuterIdx) const { return m_outerPtr(iOuterIdx + 1) - m_outerPtr(iOuterIdx); }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \param[in] iNonZero - Relative position of the inner index.
   * \return The index of the i'th inner index associated with the outer index.
   */
  inline Index_t getInnerIdx(Index_t iOuterIdx, Index_t iNonZero) const {
    assert(iNonZero >= 0 && iNonZero < getNumNonZeros(iOuterIdx));
    return m_innerIdx(m_outerPtr(iOuterIdx) + iNonZero);
  }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \param[in] iNonZero - Relative position of the inner index.
   * \return The index of the i'th inner index associated with the outer index.
   */
  inline Index_t& getInnerIdx(Index_t iOuterIdx, Index_t iNonZero) {
    assert(iNonZero >= 0 && iNonZero < getNumNonZeros(iOuterIdx));
    return m_innerIdx(m_outerPtr(iOuterIdx) + iNonZero);
  }

  /*!
   * \param[in] iOuterIdx - Outer index.
   * \return Iterator to inner dimension to use in range for loops.
   */
  inline CInnerIter getInnerIter(Index_t iOuterIdx) const {
    return CInnerIter(m_innerIdx.data() + m_outerPtr(iOuterIdx), m_innerIdx.data() + m_outerPtr(iOuterIdx + 1));
  }

  /*!
   * \param[in] iOuterIdx - Outer index (row/col).
   * \param[in] iInnerIdx - Inner index (col/row).
   * \return Absolute position of non zero entry (iOuterIdx,iInnerIdx),
   *         or NNZ if position does not belong to the pattern.
   */
  inline Index_t findInnerIdx(Index_t iOuterIdx, Index_t iInnerIdx) const {
    for (Index_t k = m_outerPtr(iOuterIdx); k < m_outerPtr(iOuterIdx + 1); ++k)
      if (m_innerIdx(k) == iInnerIdx) return k;
    return m_innerIdx.size();
  }

  /*!
   * \param[in] iOuterIdx - Outer index (row/col).
   * \param[in] iInnerIdx - Inner index (col/row).
   * \return True if (iOuterIdx,iInnerIdx) exists, i.e. is non zero.
   */
  inline bool isNonZero(Index_t iOuterIdx, Index_t iInnerIdx) const {
    return findInnerIdx(iOuterIdx, iInnerIdx) < m_innerIdx.size();
  }

  /*!
   * \param[in] iOuterIdx - Outer index (row/col).
   * \param[in] iInnerIdx - Inner index (col/row).
   * \return Absolute position of non zero entry (iOuterIdx,iInnerIdx).
   * \note This method is only safe if the entry exists.
   */
  inline Index_t quickFindInnerIdx(Index_t iOuterIdx, Index_t iInnerIdx) const {
    assert(isNonZero(iOuterIdx, iInnerIdx) && "Error, j does not belong to NZ(i).");
    Index_t k = m_outerPtr(iOuterIdx);
    while (m_innerIdx(k) != iInnerIdx) ++k;
    return k;
  }

  /*!
   * \param[in] iDiagIdx - Diagonal index (row == col).
   * \return Absolute position of the diagonal entry.
   */
  inline Index_t getDiagPtr(Index_t iDiagIdx) const { return m_diagPtr(iDiagIdx); }

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
   * \return Raw pointer to the inner index vector, offset for a given outer index.
   */
  inline const Index_t* innerIdx(Index_t iOuterIdx) const {
    assert(!empty() && "Sparse pattern has not been built.");
    return m_innerIdx.data() + m_outerPtr(iOuterIdx);
  }

  /*!
   * \return Raw pointer to the diagonal pointer vector.
   */
  inline const Index_t* diagPtr() const {
    assert(!m_diagPtr.empty() && "Diagonal map has not been built.");
    return m_diagPtr.data();
  }

  /*!
   * \return Raw pointer to the transpose pointer vector.
   */
  inline const su2vector<Index_t>& transposePtr() const {
    assert(!m_innerIdxTransp.empty() && "Transpose map has not been built.");
    return m_innerIdxTransp;
  }

  /*!
   * \return The minimum inner index.
   */
  Index_t getMinInnerIdx() const {
    Index_t idx = std::numeric_limits<Index_t>::max();
    for (Index_t k = 0; k < m_innerIdx.size(); ++k) idx = std::min(idx, m_innerIdx(k));
    return idx;
  }

  /*!
   * \return The maximum inner index.
   */
  Index_t getMaxInnerIdx() const {
    Index_t idx = std::numeric_limits<Index_t>::min();
    for (Index_t k = 0; k < m_innerIdx.size(); ++k) idx = std::max(idx, m_innerIdx(k));
    return idx;
  }
};

/*!
 * \brief Alias a type of container as the edge map class, this is a N by 2 container
 *        that maps the two sparse entries referenced by an edge (ij and ji) to two
 *        non zero entries of a sparse pattern.
 */
template <typename Index_t>
using CEdgeToNonZeroMap = C2DContainer<unsigned long, Index_t, StorageType::RowMajor, 64, DynamicSize, 2>;

using CCompressedSparsePatternUL = CCompressedSparsePattern<unsigned long>;
using CCompressedSparsePatternL = CCompressedSparsePattern<long>;
using CEdgeToNonZeroMapUL = CEdgeToNonZeroMap<unsigned long>;

/*!
 * \brief Build a sparse pattern from geometry information, of type FVM or FEM,
 *        for a given fill-level. At fill-level N, the immediate neighbors of the
 *        points in level N-1 are also considered neighbors of the base point.
 *        The resulting pattern is that of A^{N+1} where A is the sparse matrix
 *        of immediate neighbors.
 * \note Algorithm is equivalent to the implementation by F. Palacios,
 *       A. Bueno, and T. Economon from CSysMatrix.
 * \param[in] geometry - Definition of the grid.
 * \param[in] type - Of connectivity.
 * \param[in] fillLvl - Target degree of neighborhood (immediate neighbors always added).
 * \return Compressed-Storage-Row sparse pattern.
 */
template <class Geometry_t, typename Index_t>
CCompressedSparsePattern<Index_t> buildCSRPattern(Geometry_t& geometry, ConnectivityType type, Index_t fillLvl) {
  Index_t nPoint = geometry.GetnPoint();

  std::vector<Index_t> outerPtr(nPoint + 1);
  std::vector<Index_t> innerIdx;
  innerIdx.reserve(nPoint);  // at least this much space is needed

  for (Index_t iPoint = 0; iPoint < nPoint; ++iPoint) {
    /*--- Inner indices for iPoint start here. ---*/
    outerPtr[iPoint] = innerIdx.size();

    /*--- Use a set to avoid duplication and keep ascending order. ---*/
    std::set<Index_t> neighbors;

    /*--- Insert base point. ---*/
    neighbors.insert(iPoint);

    /*--- Neighbors added in previous level. ---*/
    std::set<Index_t> addedNeighbors(neighbors);

    for (Index_t iLevel = 0;; ++iLevel) {
      /*--- New points added in this level. ---*/
      std::set<Index_t> newNeighbors;

      /*--- For each point previously added, add its level 0
       *    neighbors, not duplicating any existing neighbor. ---*/
      for (auto jPoint : addedNeighbors) {
        if (type == ConnectivityType::FiniteVolume) {
          /*--- For FVM we know the neighbors of point j directly. ---*/
          for (Index_t kPoint : geometry.nodes->GetPoints(jPoint))
            if (neighbors.count(kPoint) == 0)  // no duplication
              newNeighbors.insert(kPoint);
        } else  // FiniteElement
        {
          /*--- For FEM we need the nodes of all elements that contain point j. ---*/
          for (auto iElem : geometry.nodes->GetElems(jPoint)) {
            auto elem = geometry.elem[iElem];

            for (unsigned short iNode = 0; iNode < elem->GetnNodes(); ++iNode) {
              Index_t kPoint = elem->GetNode(iNode);

              if (neighbors.count(kPoint) == 0)  // no duplication
                newNeighbors.insert(kPoint);
            }
          }
        }
      }

      neighbors.insert(newNeighbors.begin(), newNeighbors.end());

      if (iLevel >= fillLvl) break;

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
template <class Geometry_t, typename Index_t>
CEdgeToNonZeroMap<Index_t> mapEdgesToSparsePattern(Geometry_t& geometry,
                                                   const CCompressedSparsePattern<Index_t>& pattern) {
  assert(!pattern.empty());

  CEdgeToNonZeroMap<Index_t> edgeMap(geometry.GetnEdge(), 2);

  for (Index_t iEdge = 0; iEdge < geometry.GetnEdge(); ++iEdge) {
    Index_t iPoint = geometry.edges->GetNode(iEdge, 0);
    Index_t jPoint = geometry.edges->GetNode(iEdge, 1);

    edgeMap(iEdge, 0) = pattern.quickFindInnerIdx(iPoint, jPoint);
    edgeMap(iEdge, 1) = pattern.quickFindInnerIdx(jPoint, iPoint);
  }

  return edgeMap;
}

/*!
 * \brief Create the natural coloring (equivalent to the normal sequential loop
 *        order) for a given number of inner indexes.
 * \note This is to reduce overhead in "OpenMP-ready" code when only 1 thread is used.
 * \param[in] numInnerIndexes - Number of indexes that are to be colored.
 * \return Natural (sequential) coloring of the inner indices.
 */
template <class T = CCompressedSparsePatternUL, class Index_t = typename T::IndexType>
T createNaturalColoring(Index_t numInnerIndexes) {
  /*--- One color. ---*/
  su2vector<Index_t> outerPtr(2);
  outerPtr(0) = 0;
  outerPtr(1) = numInnerIndexes;

  /*--- Containing all indexes in ascending order. ---*/
  su2vector<Index_t> innerIdx(numInnerIndexes);
  std::iota(innerIdx.data(), innerIdx.data() + numInnerIndexes, 0);

  return T(std::move(outerPtr), std::move(innerIdx));
}

/*!
 * \brief Color contiguous groups of outer indices of a sparse pattern such that
 *        within each color, any two groups do not have inner indices in common.
 * \note  Within a group, two outer indices will generally have common inner indices.
 *        The coloring is returned as a compressed sparse pattern where the colors
 *        are outer indices, and the outer indices of the input pattern are the
 *        inner indices of the coloring. A simple greedy algorithm is used.
 *        Using a sparse pattern as input allows "anything" to be colored e.g.
 *        FVM edges, FEM elements, the rows/columns of a sparse matrix, etc.
 * \note  The worst that can happen in this method is needing an unreasonable number
 *        of colors, or too much memory due to a large range of the inner indices.
 *        The last two template parameters limit both, in case of failure an empty
 *        pattern is returned.
 * \param[in] pattern - Sparse pattern to be colored.
 * \param[in] groupSize - Size of the outer index groups, default 1.
 * \param[in] balanceColors - Try to balance number of indexes per color,
 *            tends to result in worse locality (thus false by default).
 * \param[out] indexColor - Optional, vector with colors given to the outer indices.
 * \return Coloring in the same type of the input pattern.
 */
template <typename Color_t = char, size_t MaxColors = 64, size_t MaxMB = 128, class T>
T colorSparsePattern(const T& pattern, size_t groupSize = 1, bool balanceColors = false,
                     std::vector<Color_t>* indexColor = nullptr) {
  static_assert(std::is_integral<Color_t>::value, "");
  static_assert(std::numeric_limits<Color_t>::max() >= MaxColors, "");

  using Index_t = typename T::IndexType;

  const Index_t grpSz = groupSize;
  const Index_t nOuter = pattern.getOuterSize();

  /*--- Trivial case. ---*/
  if (groupSize >= nOuter) return createNaturalColoring(nOuter);

  const Index_t minIdx = pattern.getMinInnerIdx();
  const Index_t nInner = pattern.getMaxInnerIdx() + 1 - minIdx;

  /*--- Check the max memory condition (<< 23 is to count bits). ---*/
  if (size_t(nInner) > (MaxMB << 23)) return T();

  /*--- Vector with the color given to each outer index. ---*/
  std::vector<Color_t> idxColor(nOuter);

  /*--- Start with one color, with no indices assigned. ---*/
  std::vector<Index_t> colorSize(1, 0);
  Color_t nColor = 1;

  {
    /*--- For each color keep track of the inner indices that are in it. ---*/
    std::vector<std::vector<bool> > innerInColor;
    innerInColor.emplace_back(nInner, false);

    /*--- Order in which we look for space in the colors to insert a new group. ---*/
    std::vector<Color_t> searchOrder(MaxColors);

    auto outerPtr = pattern.outerPtr();
    auto innerIdx = pattern.innerIdx();

    for (Index_t iOuter = 0; iOuter < nOuter; iOuter += grpSz) {
      Index_t grpEnd = std::min(iOuter + grpSz, nOuter);

      searchOrder.resize(nColor);
      std::iota(searchOrder.begin(), searchOrder.end(), 0);

      /*--- Balance sizes by looking for space in smaller colors first. ---*/
      if (balanceColors) {
        std::sort(searchOrder.begin(), searchOrder.end(),
                  [&colorSize](Color_t a, Color_t b) { return colorSize[a] < colorSize[b]; });
      }

      auto it = searchOrder.begin();

      for (; it != searchOrder.end(); ++it) {
        bool free = true;
        /*--- Traverse entire group as a large outer index. ---*/
        for (Index_t k = outerPtr[iOuter]; k < outerPtr[grpEnd] && free; ++k) {
          free = !innerInColor[*it][innerIdx[k] - minIdx];
        }
        /*--- If none of the inner indices in the group appears in
         *    this color yet, it is assigned to the group. ---*/
        if (free) break;
      }

      Color_t color;

      if (it != searchOrder.end()) {
        /*--- Found a free color. ---*/
        color = *it;
      } else {
        /*--- No color was free, make space for a new one. ---*/
        color = nColor++;
        if (nColor == MaxColors) return T();
        colorSize.push_back(0);
        innerInColor.emplace_back(nInner, false);
      }

      /*--- Assign color to group. ---*/
      for (Index_t k = iOuter; k < grpEnd; ++k) idxColor[k] = color;

      /*--- Mark the inner indices of the group as belonging to the color. ---*/
      for (Index_t k = outerPtr[iOuter]; k < outerPtr[grpEnd]; ++k) {
        innerInColor[color][innerIdx[k] - minIdx] = true;
      }

      /*--- Update count for the assigned color. ---*/
      colorSize[color] += grpEnd - iOuter;
    }
  }  // matrix of bools goes out of scope

  /*--- Compress the coloring information. ---*/

  su2vector<Index_t> colorPtr(nColor + 1);
  colorPtr(0) = 0;
  su2vector<Index_t> outerIdx(nOuter);

  Index_t k = 0;
  for (Color_t color = 0; color < nColor; ++color) {
    colorPtr(color + 1) = colorPtr(color) + colorSize[color];

    for (Index_t iOuter = 0; iOuter < nOuter; ++iOuter)
      if (idxColor[iOuter] == color) outerIdx(k++) = iOuter;
  }

  /*--- Optional return of the direct color information. ---*/
  if (indexColor) *indexColor = std::move(idxColor);

  /*--- Move compressed coloring into result pattern instance. ---*/
  return T(std::move(colorPtr), std::move(outerIdx));
}

/*!
 * \brief A way to represent one grid color that allows range-for syntax.
 */
template <typename T = unsigned long>
struct GridColor {
  static_assert(std::is_integral<T>::value, "");

  const T size;
  T groupSize;
  const T* const indices;

  GridColor(const T* idx = nullptr, T sz = 0, T grp = 0) : size(sz), groupSize(grp), indices(idx) {}

  inline const T* begin() const { return indices; }
  inline const T* end() const { return indices + size; }
};

/*!
 * \brief A way to represent natural coloring {0,1,2,...,size-1} with zero
 * overhead (behaves like looping with an integer index, after optimization...).
 */
template <typename T = unsigned long>
struct DummyGridColor {
  static_assert(std::is_integral<T>::value, "");

  T size;
  struct {
    inline T operator[](T i) const { return i; }
  } indices;

  DummyGridColor(T sz = 0) : size(sz) {}

  struct IteratorLikeInt {
    T i;
    inline IteratorLikeInt(T pos = 0) : i(pos) {}
    inline IteratorLikeInt& operator++() {
      ++i;
      return *this;
    }
    inline IteratorLikeInt operator++(int) {
      auto j = i++;
      return IteratorLikeInt(j);
    }
    inline T operator*() const { return i; }
    inline T operator->() const { return i; }
    inline bool operator==(const IteratorLikeInt& other) const { return i == other.i; }
    inline bool operator!=(const IteratorLikeInt& other) const { return i != other.i; }
  };
  inline IteratorLikeInt begin() const { return IteratorLikeInt(0); }
  inline IteratorLikeInt end() const { return IteratorLikeInt(size); }
};

/*!
 * \brief Computes the efficiency of a grid coloring for given number of threads and chunk size.
 */
template <class SparsePattern>
su2double coloringEfficiency(const SparsePattern& coloring, int numThreads, int chunkSize) {
  using Index_t = typename SparsePattern::IndexType;

  /*--- Ideally compute time is proportional to total work over number of threads. ---*/
  su2double ideal = coloring.getNumNonZeros() / su2double(numThreads);

  /*--- In practice the total work is quantized first by colors and then by chunks. ---*/
  Index_t real = 0;
  for (Index_t color = 0; color < coloring.getOuterSize(); ++color)
    real += chunkSize * roundUpDiv(roundUpDiv(coloring.getNumNonZeros(color), chunkSize), numThreads);

  return ideal / real;
}

/// @}
