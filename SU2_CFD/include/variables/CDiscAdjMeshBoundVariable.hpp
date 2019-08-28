/*!
 * \file CDiscAdjMeshVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the adjoint variables of the mesh movement.
 * \author Ruben Sanchez
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

#include "CDiscAdjMeshVariable.hpp"

class CDiscAdjMeshBoundVariable : public CDiscAdjMeshVariable {
protected:

  su2double* Bound_Disp_Sens;     /*!< \brief Store the reference coordinates of the mesh. */
  su2double* Bound_Disp_Direct;   /*!< \brief Store the reference boundary displacements of the mesh. */

  su2double* Solution_BGS_k;      /*!< \brief BGS solution to compute overall convergence. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_coor - Values of the coordinates (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjMeshBoundVariable(su2double *val_coor, unsigned short val_nDim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMeshBoundVariable(void);

  /*!
   * \brief Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline su2double* GetBoundDisp_Direct(void) { return Bound_Disp_Direct; }

  /*!
   * \brief Set the solution for the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline void SetBoundDisp_Direct(const su2double *val_BoundDisp) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Bound_Disp_Direct[iDim] = val_BoundDisp[iDim];
  }

  /*!
   * \brief Set the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] val_sens - Pointer to the sensitivities of the boundary displacements.
   */
  inline void SetBoundDisp_Sens(const su2double *val_sens) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Bound_Disp_Sens[iDim] = val_sens[iDim];
  }

  /*!
   * \brief Get the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord_Sens[nDim]
   * \return Value of the original Mesh_Coord_Sens iDim.
   */
  inline su2double GetBoundDisp_Sens(unsigned short iDim) const final { return Bound_Disp_Sens[iDim]; }

  /*!
   * \brief Determine whether the node is a moving vertex.
   * \return True. The node is at the boundary.
   */
  inline bool Get_isVertex(void) const final { return true; }

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  inline void Set_BGSSolution_k(void) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++)
      Solution_BGS_k[iVar] = Bound_Disp_Sens[iVar];
  }

  /*!
   * \brief Get the value of the solution in the previous BGS subiteration.
   * \param[out] val_solution - solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution_k(unsigned short iDim) { return Solution_BGS_k[iDim];}

};
