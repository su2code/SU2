/*!
 * \file CMeshFEM_Base.hpp
 * \brief Base class definition for a mesh object for the FEM solver.
 *        The implementations are in the <i>CMeshFEM_Base.cpp</i> file.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
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

#include "../../../include/toolboxes/classes_multiple_integers.hpp"
#include "../../fem/CFEMStandardElementBase.hpp"
#include "../../fem/CGemmBase.hpp"
#include "../CGeometry.hpp"
#include "CPointFEM.hpp"
#include "CBoundaryFEM.hpp"

using namespace std;

/*!
 * \class CMeshFEM_Base
 * \brief Base class for the FEM solver.
 * \author E. van der Weide
 * \version 7.0.8 "Blackbird"
 */
class CMeshFEM_Base: public CGeometry {
protected:
  unsigned long nVolElemOwned{0};  /*!< \brief Number of owned local volume elements. */

  vector<CPointFEM>    meshPoints;   /*!< \brief Vector of the points of the FEM mesh. */
  vector<CBoundaryFEM> boundaries;   /*!< \brief Vector of the boundaries of the FEM mesh. */

  vector<CFEMStandardElementBase *> standardVolumeElementsGrid;  /*!< \brief Vector of standard volume
                                                                             elements for the grid. */
  vector<CFEMStandardElementBase *> standardSurfaceElementsGrid; /*!< \brief Vector of standard surface
                                                                             elements for the grid. */

  vector<CGemmBase *> gemmTypesFaces;   /*!< \brief Vector of gemm types that occur for the faces. */

  vector<int> ranksRecv;     /*!< \brief Vector of ranks, from which this rank will receive halo
                                         information. Self communication is included. */
  vector<int> ranksSend;     /*!< \brief Vector of ranks, to which this rank will send halo
                                         information. Self communication is included. */

  vector<vector<unsigned long> > entitiesSend; /*!< \brief Vector of vector, which contains the entities that
                                                           must be sent. Self communication is included. For DG
                                                           an entitity is an element, for regular FEM an entity
                                                           is a DOF. */
  vector<vector<unsigned long> > entitiesRecv; /*!< \brief Vector of vector, which contains the entities that
                                                           must be received. Self communication is included. For DG
                                                           an entity is an element, for regular FEM an entity
                                                           is a DOF. */

  vector<unsigned short> rotPerMarkers; /*!< \brief Vector, which contains the indices of the rotational
                                                    periodic markers. */
  vector<vector<unsigned long> > rotPerHalos; /*!< \brief Vector of vector, which contains the indices of
                                                          the halo elements for which a rotationally periodic
                                                          correction must be applied. */
public:
  /*!
   * \brief Constructor of the class.
   */
 CMeshFEM_Base(void) : CGeometry() { }

  /*!
   * \overload
   * \brief Redistributes the grid over the ranks and creates the halo layer.
   * \param[in] geometry - The linear distributed grid that must be redistributed.
   * \param[in] config   - Definition of the particular problem.
   */
  CMeshFEM_Base(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CMeshFEM_Base(void);

protected:

  /*!
   * \brief Function, which creates the standard elements for the grid.
   * \param[in] elemTypes   - Information about the element types to be created.
   * \param[in] locGridDOFs - Location of the grid DOFs, either LGL or equidistant.
   */
  void CreateStandardVolumeElementsGrid(const vector<CUnsignedShort4T> &elemTypes,
                                        const unsigned short           locGridDOFs);
};
