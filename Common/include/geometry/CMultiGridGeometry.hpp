/*!
 * \file CMultiGridGeometry.hpp
 * \brief Headers of the multigrid geometry class.
 * \author F. Palacios, T. Economon
 * \version 8.4.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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
class CMultiGridQueue;

/*!
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delegated part is the
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios
 */
class CMultiGridGeometry final : public CGeometry {
 private:
  /*!
   * \brief Determine if a CVPoint can be agglomerated, if it has the same marker point as the seed.
   * \param[in] CVPoint - Control volume to be agglomerated.
   * \param[in] marker_seed - Marker of the seed.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return <code>TRUE</code> or <code>FALSE</code> depending if the control volume can be agglomerated.
   */
  bool SetBoundAgglomeration(unsigned long CVPoint, vector<short> marker_seed, const CGeometry* fine_grid,
                             const CConfig* config) const;

  /*!
   * \brief Determine if a Point can be agglomerated using geometrical criteria.
   * \param[in] iPoint - Seed point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  bool GeometricalCheck(unsigned long iPoint, const CGeometry* fine_grid, const CConfig* config) const;

  /*!
   * \brief Determine if a CVPoint can be agglomerated, if it has the same marker point as the seed.
   * \param[out] Suitable_Indirect_Neighbors - List of Indirect Neighbours that can be agglomerated.
   * \param[in] iPoint - Seed point.
   * \param[in] Index_CoarseCV - Index of agglomerated point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   */
  void SetSuitableNeighbors(vector<unsigned long>& Suitable_Indirect_Neighbors, unsigned long iPoint,
                            unsigned long Index_CoarseCV, const CGeometry* fine_grid) const;

  /*!
   * \brief Compute surface straightness for multigrid geometry.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeSurfStraightness(CConfig* config);

  /*!
   * \brief Compute local curvature at a boundary vertex on Euler wall.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] iPoint - Point index.
   * \param[in] iMarker - Marker index.
   * \return Maximum angle (in degrees) between this vertex normal and adjacent vertex normals.
   */
  su2double ComputeLocalCurvature(const CGeometry* fine_grid, unsigned long iPoint, unsigned short iMarker) const;

 public:
  /*--- This is to suppress Woverloaded-virtual, omitting it has no negative impact. ---*/
  using CGeometry::SetBoundControlVolume;
  using CGeometry::SetControlVolume;
  using CGeometry::SetPoint_Connectivity;
  using CGeometry::SetVertex;

  /*!
   * \brief Constructor of the class.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Level of the multigrid.
   */
  CMultiGridGeometry(CGeometry* fine_grid, CConfig* config, unsigned short iMesh);

  /*!
   * \brief Set boundary vertex.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetVertex(const CGeometry* fine_grid, const CConfig* config) override;

  /*!
   * \brief Set points which surround a point.
   * \param[in] fine_grid - Geometrical definition of the child grid.
   */
  void SetPoint_Connectivity(const CGeometry* fine_grid) override;

  /*!
   * \brief Set the edge structure of the agglomerated control volume.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] action - Allocate or not the new elements.
   */
  void SetControlVolume(const CGeometry* fine_grid, unsigned short action) override;

  /*!
   * \brief Set boundary vertex structure of the agglomerated control volume.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] action - Allocate or not the new elements.
   */
  void SetBoundControlVolume(const CGeometry* fine_grid, const CConfig* config, unsigned short action) override;

  /*!
   * \brief Set a representative coordinates of the agglomerated control volume.
   * \param[in] fine_grid - Geometrical definition of the problem.
   */
  void SetCoord(const CGeometry* fine_grid) override;

  /*!
   * \brief Set the grid velocity at each node in the coarse mesh level based
   *        on a restriction from a finer mesh.
   * \param[in] fine_grid - Geometry container for the finer mesh level.
   */
  void SetRestricted_GridVelocity(const CGeometry* fine_grid) override;

  /*!
   * \brief Find and store the closest, most normal, neighbor to a vertex.
   * \param[in] config - Definition of the particular problem.
   */
  void FindNormal_Neighbor(const CConfig* config) override;

  /*!
   * \brief Mach the near field boundary condition.
   * \param[in] config - Definition of the particular problem.
   */
  void MatchActuator_Disk(const CConfig* config) override;

  /*!
   * \brief Mach the periodic boundary conditions.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_periodic - Index of the first periodic face in a pair.
   */
  void MatchPeriodic(const CConfig* config, unsigned short val_periodic) override;

  /*!
   * \brief Set a representative wall normal heat flux of the agglomerated control volume on a particular boundary
   * marker. \param[in] fine_grid - Geometrical definition of the problem. \param[in] val_marker - Index of the boundary
   * marker.
   */
  void SetMultiGridWallHeatFlux(const CGeometry* fine_grid, unsigned short val_marker) override;

  /*!
   * \brief Set a representative wall temperature of the agglomerated control volume on a particular boundary marker.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   */
  void SetMultiGridWallTemperature(const CGeometry* fine_grid, unsigned short val_marker) override;

  /*!
   * \brief Validate that halo CV coordinates match corresponding domain CVs on remote ranks (debug feature).
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Multigrid level for reporting.
   */
  void ValidateHaloCoordinates(const CConfig* config, unsigned short iMesh) const;
};
