/*!
 * \file CMeshFEM_Base.hpp
 * \brief Base class definition for a mesh object for the FEM solver.
 *        The implementations are in the <i>CMeshFEM_Base.cpp</i> file.
 * \author E. van der Weide
 * \version 7.1.0 "Blackbird"
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
 * \version 7.1.0 "Blackbird"
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

  /*-----------------------------------------------------------------------------------*/
  /*---                     Constructors and destructors.                           ---*/
  /*-----------------------------------------------------------------------------------*/

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

  /*-----------------------------------------------------------------------------------*/
  /*---                  Inline public member functions.                            ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which makes available the boundaries of the local FEM mesh.
   * \return  Pointer to the boundaries of the local FEM mesh.
   */
  inline CBoundaryFEM *GetBoundaries(void) {return boundaries.data();}

  /*!
   * \brief Function, which makes available the mesh points of the local FEM mesh.
   * \return  Pointer to the mesh points of the local FEM mesh.
   */
  inline CPointFEM *GetMeshPoints(void) {return meshPoints.data();}

  /*!
   * \brief Function, which makes available the number of mesh points of the local FEM mesh.
   * \return  Number of mesh points of the local FEM mesh.
   */
  inline unsigned long GetNMeshPoints(void) {return meshPoints.size();}

  /*!
   * \brief Function, which makes available the number of owned volume elements in the local FEM mesh.
   * \return  Number of owned volume elements of the local FEM mesh.
   */
  inline unsigned long GetNVolElemOwned(void) const {return nVolElemOwned;}

  /*!
   * \brief Function, which makes available the vector of receive ranks as
   *        a const reference.
   * \return  Const reference to the vector of ranks.
   */
  inline const vector<int>& GetRanksRecv(void) const {return ranksRecv;}

  /*!
   * \brief Function, which makes available the vector of send ranks as
   *        a const reference.
   * \return  Const reference to the vector of ranks.
   */
  inline const vector<int>& GetRanksSend(void) const {return ranksSend;}

  /*!
   * \brief Function, which makes available the vector of vectors containing the receive
   *        entities as a const reference.
   * \return  Const reference to the vector of vectors of receive entities.
   */
  inline const vector<vector<unsigned long> >& GetEntitiesRecv(void) const {return entitiesRecv;}

  /*!
   * \brief Function, which makes available the vector of vectors containing the send
   *        entities as a const reference.
   * \return  Const reference to the vector of vectors of send entities.
   */
  inline const vector<vector<unsigned long> >& GetEntitiesSend(void) const {return entitiesSend;}

protected:

  /*-----------------------------------------------------------------------------------*/
  /*---                     Proteced member functions.                              ---*/
  /*-----------------------------------------------------------------------------------*/

  /*!
   * \brief Function, which creates the standard elements for the grid.
   * \param[in] elemTypes   - Information about the element types to be created.
   * \param[in] locGridDOFs - Location of the grid DOFs, either LGL or equidistant.
   */
  void CreateStandardVolumeElementsGrid(const vector<CUnsignedShort4T> &elemTypes,
                                        const unsigned short           locGridDOFs);

  /*!
   * \brief Compute an ADT including the coordinates of all viscous markers
   * \param[in] config - Definition of the particular problem.
   * \return pointer to the ADT
   */
  std::unique_ptr<CADTElemClass> ComputeViscousWallADT(const CConfig *config) const override;
};
