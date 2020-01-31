/*!
 * \file CMultiGridGeometry.hpp
 * \brief Headers of the multigrid geometry class.
 * \author F. Palacios, T. Economon
 * \version 7.0.1 "Blackbird"
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

#include "CGeometry.hpp"

/*!
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delicated part is the
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios
 */
class CMultiGridGeometry final : public CGeometry {

public:
  /*--- This is to suppress Woverloaded-virtual, omitting it has no negative impact. ---*/
  using CGeometry::SetVertex;
  using CGeometry::SetMeshFile;
  using CGeometry::SetControlVolume;
  using CGeometry::SetBoundControlVolume;
  using CGeometry::SetPoint_Connectivity;

  /*!
   * \brief Constructor of the class.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Level of the multigrid.
   * \param[in] iZone - Current zone in the mesh.
   */
  CMultiGridGeometry(CGeometry **geometry, CConfig *config_container, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CMultiGridGeometry(void);

  /*!
   * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
   * \param[in] CVPoint - Control volume to be agglomerated.
   * \param[in] marker_seed - Marker of the seed.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return <code>TRUE</code> or <code>FALSE</code> depending if the control volume can be agglomerated.
   */
  bool SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, CGeometry *fine_grid, CConfig *config);

  /*!
   * \brief Determine if a can be agglomerated using geometrical criteria.
   * \param[in] iPoint - Seed point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  bool GeometricalCheck(unsigned long iPoint, CGeometry *fine_grid, CConfig *config);

  /*!
   * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
   * \param[in] Suitable_Indirect_Neighbors - List of Indirect Neighbours that can be agglomerated.
   * \param[in] iPoint - Seed point.
   * \param[in] Index_CoarseCV - Index of agglomerated point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   */
  void SetSuitableNeighbors(vector<unsigned long> *Suitable_Indirect_Neighbors, unsigned long iPoint,
                            unsigned long Index_CoarseCV, CGeometry *fine_grid);

  /*!
   * \brief Set boundary vertex.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetVertex(CGeometry *geometry, CConfig *config) override;

  /*!
   * \brief Set points which surround a point.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SetPoint_Connectivity(CGeometry *geometry) override;

  /*!
   * \brief Set the edge structure of the agglomerated control volume.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] action - Allocate or not the new elements.
   */
  void SetControlVolume(CConfig *config, CGeometry *geometry, unsigned short action) override;

  /*!
   * \brief Mach the near field boundary condition.
   * \param[in] config - Definition of the particular problem.
   */
  void MatchNearField(CConfig *config) override;

  /*!
   * \brief Mach the near field boundary condition.
   * \param[in] config - Definition of the particular problem.
   */
  void MatchActuator_Disk(CConfig *config) override;

  /*!
   * \brief Mach the periodic boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_periodic - Index of the first periodic face in a pair.
   */
  void MatchPeriodic(CConfig *config, unsigned short val_periodic) override;

  /*!
   * \brief Set boundary vertex structure of the agglomerated control volume.
   * \param[in] config - Definition of the particular problem.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] action - Allocate or not the new elements.
   */
  void SetBoundControlVolume(CConfig *config, CGeometry *geometry, unsigned short action) override;

  /*!
   * \brief Set a representative coordinates of the agglomerated control volume.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void SetCoord(CGeometry *geometry) override;

  /*!
   * \brief Set a representative wall normal heat flux of the agglomerated control volume on a particular boundary marker.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   */
  void SetMultiGridWallHeatFlux(CGeometry *geometry, unsigned short val_marker) override;

  /*!
   * \brief Set a representative wall temperature of the agglomerated control volume on a particular boundary marker.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   */
  void SetMultiGridWallTemperature(CGeometry *geometry, unsigned short val_marker) override;

  /*!
   * \brief Set the grid velocity at each node in the coarse mesh level based
   *        on a restriction from a finer mesh.
   * \param[in] fine_mesh - Geometry container for the finer mesh level.
   * \param[in] config - Definition of the particular problem.
   */
  void SetRestricted_GridVelocity(CGeometry *fine_mesh, CConfig *config) override;

  /*!
   * \brief Find and store the closest neighbor to a vertex.
   * \param[in] config - Definition of the particular problem.
   */
  void FindNormal_Neighbor(CConfig *config) override;

};

