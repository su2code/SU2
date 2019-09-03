/*!
 * \file CVertexMap.hpp
 * \brief An index to index lookup vector.
 * \author P. Gomes
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include <type_traits>
#include <cassert>

/*!
 * \class CVertexMap
 * \brief A lookup type map, maps indices in a large range to indices in a smaller one.
 *
 * The usage case is:
 *   1 - Initialize via Reset for the large range size.
 *   2 - Use SetIsVertex to define what large range indices have small range correspondence.
 *   3 - Call Build, vertices are given a sequential index in increasing order of iPoint,
 *       not in order of calls to SetIsVertex. Keep this in mind when allocating vertex data.
 *   4 - Use GetVertexIndex to convert a point index to vertex index.
 *
 * Only consider using this class if you have no reasonable way to keep track of vertices.
 * For example if you need to assume every point might be a vertex, but you do not want to
 * allocate data for the entire set of points.
 *
 * \note For efficiency use the smallest type that can fit the maximum number of vertices.
 */
template<typename T = unsigned long>
class CVertexMap {
  static_assert(std::is_unsigned<T>::value && std::is_integral<T>::value,
                "Vertex map requires an unsigned integral type (e.g. unsigned).");

 private:
  su2vector<T> Map;      /*!< \brief Map from range 0-(nPoint-1) to 1-nVertex. */
  bool isValid = false;  /*!< \brief Set to true when it is safe to use the accessors. */
  T nVertex = 0;         /*!< \brief Number of vertices. */

 public:
  /*!
   * \brief Check if the current mapping is valid.
   */
  inline bool GetIsValid() const { return isValid; }

  /*!
   * \brief Get the number of vertices currently mapped.
   */
  inline T GetnVertex() const { return nVertex; }

  /*!
   * \brief Reset the map for size nPoint, marks every point as not-a-vertex.
   */
  void Reset(unsigned long nPoint) {
    Map.resize(nPoint) = 0;
    nVertex = 0;
    isValid = true;
  }

  /*!
   * \brief Set the vertex status of a point.
   */
  inline bool SetIsVertex(unsigned long iPoint, bool isVertex) {
    /*--- Invalidate map if change is requested as that destroys it. ---*/
    if (isVertex != bool(Map(iPoint))) {
      isValid = false;
      Map(iPoint) = T(isVertex);
    }
    return isValid;
  }

  /*!
   * \brief Get wheter a point is marked as vertex.
   */
  inline bool GetIsVertex(unsigned long iPoint) const {
    return (Map(iPoint) != 0);
  }

  /*!
   * \brief Build the point to vertex map.
   */
  T Build() {
    if (!isValid) {
      /*--- The map is 1 based, the accessors correct
       accordingly when returning vertex indices. ---*/
      nVertex = 0;

      for (unsigned long iPoint = 0; iPoint < Map.size(); ++iPoint)
        if (Map(iPoint)!=0)
          Map(iPoint) = ++nVertex;

      isValid = true;
    }
    return nVertex;
  }

  /*!
   * \brief Convert a point index to vertex index.
   * \param[in,out] iVertex - On entry point index, on exit vertex index.
   * \return True if conversion is successful (i.e. point is vertex).
   */
  inline bool GetVertexIndex(unsigned long &iVertex) const {
    assert(isValid && "Vertex map is not in valid state.");
    iVertex = Map(iVertex);
    if(iVertex==0) return false; // not a vertex
    iVertex--;    // decrement for 0 based
    return true;  // is a vertex
  }

};
