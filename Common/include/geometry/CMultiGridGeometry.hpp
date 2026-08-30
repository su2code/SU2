/*!
 * \file CMultiGridGeometry.hpp
 * \brief Headers of the multigrid geometry class.
 * \author F. Palacios, T. Economon
 * \version 8.5.0 "Harrier"
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
 * \brief Class for defining the multigrid geometry, the main dedicated part is the
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
   * \brief Compute local curvature at a boundary vertex on Euler wall.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] iPoint - Point index.
   * \param[in] iMarker - Marker index.
   * \return Maximum angle (in degrees) between this vertex normal and adjacent vertex normals.
   */
  su2double ComputeLocalCurvature(const CGeometry* fine_grid, unsigned long iPoint, unsigned short iMarker) const;

  /*!
   * \brief Agglomerate high-aspect-ratio interior cells along implicit lines from wall vertices.
   * \param[in,out] Index_CoarseCV - Current coarse CV index, incremented as new coarse CVs are created.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Configuration.
   */
  void AgglomerateImplicitLines(unsigned long& Index_CoarseCV, const CGeometry* fine_grid, const CConfig* config);

  /*!
   * \brief Per-node dual-grid stiffness data used by the implicit-line agglomeration: the weakest and
   *        strongest coupling at each node, and the neighbour the strongest one leads to. Their ratio is
   *        the local cell aspect ratio, available on every multigrid level unlike CGeometry::Aspect_Ratio.
   */
  struct CNodeStiffness {
    vector<su2double> wMin, wMax;   /*!< \brief Weakest and strongest edge coupling at each node. */
    vector<unsigned long> jStiffest; /*!< \brief Neighbour across the strongest edge. */

    /*!< \brief Local aspect ratio at a node, 1.0 where it could not be measured. */
    su2double AspectRatio(unsigned long iPoint) const {
      return (wMin[iPoint] > 0.0) ? wMax[iPoint] / wMin[iPoint] : su2double(1.0);
    }
  };

  /*!
   * \brief Measure the dual-grid coupling at every node of a grid.
   * \param[in] fine_grid - Grid to measure.
   * \return Weakest/strongest coupling per node, see CNodeStiffness.
   */
  CNodeStiffness ComputeNodeStiffness(const CGeometry* fine_grid) const;

  /*!
   * \brief PHASE A of the implicit-line agglomeration: grow wall-normal lines of nodes from the
   *        boundaries that carry a stretched layer. Lines are node-disjoint.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Configuration.
   * \param[in] stiff - Node coupling from ComputeNodeStiffness.
   * \return One vector per line, each starting at its boundary node.
   */
  vector<vector<unsigned long>> BuildImplicitLines(const CGeometry* fine_grid, const CConfig* config,
                                                  const CNodeStiffness& stiff) const;

  /*!
   * \brief PHASE B of the implicit-line agglomeration: partition the lines into compact bundles that
   *        share a wall footprint, by repeated pairwise matching.
   * \param[in] lines - Lines from BuildImplicitLines.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Configuration.
   * \param[out] adj - Line adjacency inherited from the boundary nodes, reused by PHASE C.
   * \return One vector of line indices per bundle.
   */
  vector<vector<unsigned long>> BundleImplicitLines(const vector<vector<unsigned long>>& lines,
                                                    const CGeometry* fine_grid, const CConfig* config,
                                                    vector<vector<unsigned long>>& adj) const;

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
