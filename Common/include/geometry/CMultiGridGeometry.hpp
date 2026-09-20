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
                             const CConfig* config, const vector<char>& mixedBC) const;

  /*!
   * \brief Find nodes where two boundary conditions of different type meet. These are never
   *        agglomerated, since a coarse CV holding one would average both conditions.
   * \param[in] fine_grid - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \return One flag per fine grid point, set where that point must stay on its own.
   */
  vector<char> FindMixedBoundaryNodes(const CGeometry* fine_grid, const CConfig* config) const;

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
   * \brief Pave the domain with columns grown from the boundary patches by one breadth first walk,
   *        so that every node goes to the patch that reaches it in the fewest steps.
   * \param[in,out] Index_CoarseCV - Current coarse CV index, incremented as new coarse CVs are created.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Configuration.
   * \param[in] iMesh - Multigrid level being built, used to label the summary.
   * \param[in] mixedBC - Nodes that must stay on their own, from FindMixedBoundaryNodes.
   * \param[in] onPhysBoundary - Nodes carrying a physical boundary condition, excluding SEND_RECEIVE.
   * \param[in] onPeriodic - Nodes on a periodic marker, which a column never claims.
   * \param[out] neverGrewCV - Coarse CV index of every column that never grew past its seed layer.
   * \return Summary of the paving, empty except on the master rank.
   */
  string PaveAdvancingFronts(unsigned long& Index_CoarseCV, const CGeometry* fine_grid, const CConfig* config,
                             unsigned short iMesh, const vector<char>& mixedBC, const vector<char>& onPhysBoundary,
                             const vector<char>& onPeriodic, vector<unsigned long>& neverGrewCV);

  /*!
   * \brief Boundary nodes that seed a column, ordered with the most anisotropic first.
   */
  struct CFrontSeeds {
    vector<unsigned long> node;                    /*!< \brief Seed node on the boundary. */
    vector<std::array<su2double, MAXNDIM>> normal; /*!< \brief Unit normal there, pointing into the domain. */
    vector<su2double> strength;                    /*!< \brief Local anisotropy, the ordering key. */
    vector<char> tier;                             /*!< \brief Paving order: 0 viscous wall, 1 Euler
                                                        wall, 2 anything else. */
    unsigned long nRefusedCurvature = 0;           /*!< \brief Euler wall nodes the curvature limit kept out. */
  };

  /*!
   * \brief Collect the boundary nodes that seed a column: those where the stiffest edge is also the
   *        edge most nearly along the boundary normal, so the mesh is layered against this boundary.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Definition of the particular problem.
   * \return Seed nodes and their inward boundary normals, most anisotropic first.
   */
  CFrontSeeds SeedFrontNodes(const CGeometry* fine_grid, const CConfig* config) const;

  /*!
   * \brief Element strip a column follows while the mesh is extruded along it.
   */
  struct CWalkState {
    unsigned long elem;  /*!< \brief Element the column stands in, or NO_ELEM once it is lost. */
    unsigned short face; /*!< \brief Local face it entered that element through. */
  };

  /*!
   * \brief Partition the seed nodes into compact surface patches. A patch is a boundary face where
   *        the primal grid gives one, otherwise it is built by repeated pairwise matching.
   * \param[in] seeds - Seed nodes from SeedFrontNodes.
   * \param[in] fine_grid - Fine grid geometry.
   * \param[in] config - Definition of the particular problem.
   * \param[in] mixedBC - Nodes that must stay on their own, from FindMixedBoundaryNodes.
   * \param[out] walk - Element and face each patch starts from, NO_ELEM where it has none.
   * \return One vector of indices into seeds.node per patch, at most two entries in 2D, four in 3D.
   */
  vector<vector<unsigned long>> BuildFrontPatches(const CFrontSeeds& seeds, const CGeometry* fine_grid,
                                                  const CConfig* config, const vector<char>& mixedBC,
                                                  vector<CWalkState>& walk) const;

  string levelReport; /*!< \brief Console summary for this level. */

 public:
  /*!
   * \brief Get the console summary for this level, held back so it does not interleave with the
   *        multigrid table.
   * \return Summary text, empty except on the master rank.
   */
  const string& GetLevelReport() const { return levelReport; }

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
