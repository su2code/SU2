/*!
 * \file CMeshSolver.hpp
 * \brief Declaration and inlines of the class to compute the deformation of
 *        the volumetric numerical grid using the linear elasticity solver.
 * \author Ruben Sanchez, based on CVolumetricMovement developments (F. Palacios, A. Bueno, T. Economon, S. Padron)
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

#include "../solver_structure.hpp"

class CMeshSolver : public CFEASolver {
protected:

  bool time_domain;        /*!< \brief Number of dimensions. */
  bool multizone;

  bool stiffness_set;          /*!< \brief Element-based stiffness is set. */

  su2double *Coordinate;       /*!< \brief Auxiliary nDim vector. */

  su2double MinVolume_Ref,     /*!< \brief Minimum volume in  to make zeros and impose boundary conditions. */
            MinVolume_Curr;

  su2double MaxVolume_Ref,
            MaxVolume_Curr;

  su2double MinDistance;
  su2double MaxDistance;

  su2double E;                  /*!< \brief Young's modulus of elasticity. */
  su2double Nu;                 /*!< \brief Poisson's ratio. */

  su2double Mu;                 /*!< \brief Lame's coeficient. */
  su2double Lambda;             /*!< \brief Lame's coeficient. */

public:

  CMeshElement* element;         /*!< \brief Vector which stores element information for each problem. */

  /*!
   * \brief Constructor of the class.
   */
  CMeshSolver(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMeshSolver(void);

  /*!
   * \brief Grid deformation using the linear elasticity equations.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void DeformMesh(CGeometry **geometry, CNumerics **numerics, CConfig *config);

  /*!
   * \brief Set the stiffness of the mesh.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetMesh_Stiffness(CGeometry **geometry, CNumerics **numerics, CConfig *config);

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
   * \brief Get the value of the reference coordinate to set on the element structure.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] indexNode - Index of the node.
   * \param[in] iDim - Dimension required.
   */
  inline su2double Get_ValCoord(CGeometry *geometry, unsigned long indexNode, unsigned short iDim) {
    return nodes->GetMesh_Coord(indexNode,iDim);
  }

  /*!
   * \brief Update the value of the coordinates after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateGridCoord(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the dual grid after the grid movement (edges and control volumes).
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateDualGrid(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Compute the grid velocity form the displacements of the mesh.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeGridVelocity(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Update the coarse multigrid levels after the grid movement.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void UpdateMultiGrid(CGeometry **geometry, CConfig *config);

  /*!
   * \brief Check the boundary vertex that are going to be moved.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  void SetBoundaryDisplacements(CGeometry *geometry, CNumerics *numerics, CConfig *config);

  /*!
   * \brief Move the mesh in time.
   */
  void SetDualTime_Mesh(void);

  /*!
   * \brief Load a solution from a restart file.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all of the solvers.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_iter - Current external iteration number.
   * \param[in] val_update_geo - Flag for updating coords and grid velocity.
   */
  void LoadRestart(CGeometry **geometry, CSolver ***solver, CConfig *config, int val_iter, bool val_update_geo);

  /*!
   * \brief Load the geometries at the previous time states n and nM1.
   * \param[in] geometry - Geometrical definition of the problem.
   */
  void Restart_OldGeometry(CGeometry *geometry, CConfig *config);

  /*!
   * \brief Get minimun volume in the mesh
   * \return 
   */
  su2double GetMinimum_Volume(){return MinVolume_Curr;}
  
  /*!
   * \brief Get maximum volume in the mesh
   * \return 
   */
  su2double GetMaximum_Volume(){return MaxVolume_Curr;}
  
};
