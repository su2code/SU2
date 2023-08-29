/*!
 * \file CMeshSolver.hpp
 * \brief Declaration and inlines of the class to compute the deformation of
 *        the volumetric numerical grid using the linear elasticity solver.
 * \author Ruben Sanchez, based on CVolumetricMovement developments (F. Palacios, A. Bueno, T. Economon, S. Padron)
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

#include "CFEASolver.hpp"
#include "../variables/CMeshBoundVariable.hpp"
#include "../variables/CMeshElement.hpp"

/*!
 * \brief Mesh deformation solver (pseudo elasticity).
 * \ingroup Elasticity_Equations
 */
class CMeshSolver final : public CFEASolver {
protected:

  bool time_domain;
  bool multizone;

  bool stiffness_set;          /*!< \brief Element-based stiffness is set. */

  unsigned long ElemCounter;   /*!< \brief Error (negative volume) counter. */

  /*!
   * \brief Minimum/Maximum distance and volume (in reference and current (deformed) coords).
   */
  su2double MinVolume, MinVolume_Ref, MinVolume_Curr;
  su2double MaxVolume, MaxVolume_Ref, MaxVolume_Curr;
  su2double MinDistance, MaxDistance;

  vector<CMeshElement> element; /*!< \brief Vector which stores element information for each problem. */

  /*!
   * \brief Compute the min and max volume of the elements in the domain.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] updated - Boolean, computes the volumes with the updated coordinates.
   * \return Value of the length of the smallest edge of the grid.
   */
  void SetMinMaxVolume(CGeometry *geometry, CConfig *config, bool updated);

  /*!
   * \brief Compute the min and max volume of the elements in the domain.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetWallDistance(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the value of the coordinates after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateGridCoord(CGeometry *geometry, const CConfig *config);

  /*!
   * \brief Compute the grid velocity form the displacements of the mesh.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeGridVelocity(CGeometry **geometry, const CConfig *config) const;

  /*!
   * \brief Compute the grid velocity form the velocity at deformable boundary.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeGridVelocity_FromBoundary(CGeometry **geometry, CNumerics **numerics, CConfig *config);

  /*!
   * \brief Check the boundary vertex that are going to be moved.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] velocity_transfer - Boolean for deforming displacement or velocity
   */
  void SetBoundaryDisplacements(CGeometry *geometry, CConfig *config, bool velocity_transfer);

  /*!
   * \brief Apply forced displacement boundary conditions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Index of the marker.
   * \param[in] velocity - Boolean for deforming displacement or velocity.
   */
  void BC_Deforming(CGeometry *geometry, const CConfig *config, unsigned short val_marker, bool velocity);

  /*!
   * \brief Load the geometries at the previous time states n and nM1.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void RestartOldGeometry(CGeometry *geometry, const CConfig *config);

public:
  /*!
   * \brief Constructor of the class.
   */
  CMeshSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Grid deformation using the linear elasticity equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] numerics - Numerics used in the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void DeformMesh(CGeometry **geometry,
                  CNumerics **numerics,
                  CConfig *config) override;

  /*!
   * \brief Set the stiffness of the mesh.
   * \param[in] numerics - Numerics used in the solution.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMesh_Stiffness(CNumerics **numerics,
                         CConfig *config) override;

  /*!
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  inline su2double Get_ValCoord(const CGeometry*,
                                unsigned long indexNode,
                                unsigned short iDim) const override {
    return static_cast<const CMeshBoundVariable*>(nodes)->GetMesh_Coord(indexNode,iDim);
  }

  /*!
   * \brief Move the mesh in time.
   */
  void SetDualTime_Mesh(void) override;

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry,
                   CSolver ***solver,
                   CConfig *config,
                   int val_iter,
                   bool val_update_geo) override;

  /*!
   * \brief Get minimun volume in the mesh
   * \return
   */
  inline su2double GetMinimum_Volume() const override {return MinVolume_Curr;}

  /*!
   * \brief Get maximum volume in the mesh
   * \return
   */
  inline su2double GetMaximum_Volume() const override {return MaxVolume_Curr;}

  /*!
   * \brief Pitching definition for deforming mesh
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current time iteration number
   */
  void Surface_Pitching(CGeometry *geometry, CConfig *config, unsigned long iter);

  /*!
   * \brief Rotating definition for deforming mesh
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current time iteration number
   */
  void Surface_Rotating(CGeometry *geometry, CConfig *config, unsigned long iter);

  /*!
   * \brief Plunging definition for deforming mesh
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current time iteration number
   */
  void Surface_Plunging(CGeometry *geometry, CConfig *config, unsigned long iter);

  /*!
   * \brief Translating definition for deforming mesh
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iter - Current time iteration number
   */
  void Surface_Translating(CGeometry *geometry, CConfig *config, unsigned long iter);

};
