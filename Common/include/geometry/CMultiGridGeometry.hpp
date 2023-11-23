/*!
 * \file CMultiGridGeometry.hpp
 * \brief Headers of the multigrid geometry class.
 * \author F. Palacios, T. Economon
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

#include "CGeometry.hpp"

/*!
 * \class CMultiGridGeometry
 * \brief Class for defining the multigrid geometry, the main delicated part is the
 *        agglomeration stage, which is done in the declaration.
 * \author F. Palacios
 */
class CMultiGridGeometry final : public CGeometry {
 private:
  /*!
   * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
   * \param[in] CVPoint - Control volume to be agglomerated.
   * \param[in] marker_seed - Marker of the seed.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return <code>TRUE</code> or <code>FALSE</code> depending if the control volume can be agglomerated.
   */
  bool SetBoundAgglomeration(unsigned long CVPoint, short marker_seed, const CGeometry* fine_grid,
                             const CConfig* config) const;

  /*!
   * \brief Determine if a can be agglomerated using geometrical criteria.
   * \param[in] iPoint - Seed point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  bool GeometricalCheck(unsigned long iPoint, const CGeometry* fine_grid, const CConfig* config) const;

  /*!
   * \brief Determine if a CVPoint van be agglomerated, if it have the same marker point as the seed.
   * \param[out] Suitable_Indirect_Neighbors - List of Indirect Neighbours that can be agglomerated.
   * \param[in] iPoint - Seed point.
   * \param[in] Index_CoarseCV - Index of agglomerated point.
   * \param[in] fine_grid - Geometrical definition of the problem.
   */
  void SetSuitableNeighbors(vector<unsigned long>& Suitable_Indirect_Neighbors, unsigned long iPoint,
                            unsigned long Index_CoarseCV, const CGeometry* fine_grid) const;

  /*!
   * \brief Set a representative wall value of the agglomerated control volumes on a particular boundary marker.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] val_marker - Index of the boundary marker.
   * \param[in] wall_quantity - Object with methods Get(iVertex_fine) and Set(iVertex_coarse, val).
   */
  template <class T>
  void SetMultiGridWallQuantity(const CGeometry* fine_grid, unsigned short val_marker, T& wall_quantity) {
    for (auto iVertex = 0ul; iVertex < nVertex[val_marker]; iVertex++) {
      const auto Point_Coarse = vertex[val_marker][iVertex]->GetNode();

      if (!nodes->GetDomain(Point_Coarse)) continue;

      su2double Area_Parent = 0.0;

      /*--- Compute area parent by taking into account only volumes that are on the marker. ---*/
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
        const auto Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        const auto isVertex =
            fine_grid->nodes->GetDomain(Point_Fine) && (fine_grid->nodes->GetVertex(Point_Fine, val_marker) != -1);
        if (isVertex) {
          Area_Parent += fine_grid->nodes->GetVolume(Point_Fine);
        }
      }

      su2double Quantity_Coarse = 0.0;

      /*--- Loop again to average coarser value. ---*/
      for (auto iChildren = 0u; iChildren < nodes->GetnChildren_CV(Point_Coarse); iChildren++) {
        const auto Point_Fine = nodes->GetChildren_CV(Point_Coarse, iChildren);
        const auto isVertex =
            fine_grid->nodes->GetDomain(Point_Fine) && (fine_grid->nodes->GetVertex(Point_Fine, val_marker) != -1);
        if (isVertex) {
          const auto Vertex_Fine = fine_grid->nodes->GetVertex(Point_Fine, val_marker);
          const auto Area_Children = fine_grid->nodes->GetVolume(Point_Fine);
          Quantity_Coarse += wall_quantity.Get(Vertex_Fine) * Area_Children / Area_Parent;
        }
      }

      /*--- Set the value at the coarse level. ---*/
      wall_quantity.Set(iVertex, Quantity_Coarse);
    }
  }

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
   * \param[in] action - Allocate or not the new elements.
   */
  void SetBoundControlVolume(const CGeometry* fine_grid, unsigned short action) override;

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
   * \brief Find and store the closest neighbor to a vertex.
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
};
