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

#include "CVariable.hpp"
#include "../../../Common/include/toolboxes/CVertexMap.hpp"

class CDiscAdjMeshBoundVariable final : public CVariable {
private:

  MatrixType Bound_Disp_Sens;     /*!< \brief Store the reference coordinates of the mesh. */
  MatrixType Bound_Disp_Direct;   /*!< \brief Store the reference boundary displacements of the mesh. */

  MatrixType Solution_BGS_k;      /*!< \brief BGS solution to compute overall convergence. */

  CVertexMap<unsigned> VertexMap; /*!< \brief Object that controls accesses to the variables of this class. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Values of the coordinates (initialization value).
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjMeshBoundVariable(Idx_t npoint, Idx_t ndim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjMeshBoundVariable() = default;

  /*!
   * \brief Allocate member variables for points marked as vertex (via "Set_isVertex").
   * \param[in] config - Definition of the particular problem.
   */
  void AllocateBoundaryVariables(CConfig *config);

  /*!
   * \brief Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline const su2double* GetBoundDisp_Direct(Idx_t iPoint) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return nullptr;
    return Bound_Disp_Direct[iPoint];
  }

  /*!
   * \brief Set the solution for the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline void SetBoundDisp_Direct(Idx_t iPoint, const su2double *val_BoundDisp) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (Idx_t iDim = 0; iDim < nDim; iDim++)
      Bound_Disp_Direct(iPoint,iDim) = val_BoundDisp[iDim];
  }

  /*!
   * \brief Set the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] val_sens - Pointer to the sensitivities of the boundary displacements.
   */
  inline void SetBoundDisp_Sens(Idx_t iPoint, const su2double *val_sens) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (Idx_t iDim = 0; iDim < nDim; iDim++)
      Bound_Disp_Sens(iPoint,iDim) = val_sens[iDim];
  }

  /*!
   * \brief Get the value of the sensitivity with respect to the undeformed coordinates.
   * \param[in] iDim - Index of Mesh_Coord_Sens[nDim]
   * \return Value of the original Mesh_Coord_Sens iDim.
   */
  inline su2double GetBoundDisp_Sens(Idx_t iPoint, Idx_t iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Bound_Disp_Sens(iPoint,iDim);
  }

  /*!
   * \brief Get whether a node is on the boundary
   */
  inline bool Get_isVertex(Idx_t iPoint) const override {
    return VertexMap.GetIsVertex(iPoint);
  }

  /*!
   * \brief Set whether a node is on the boundary
   */
  inline void Set_isVertex(Idx_t iPoint, bool isVertex) override {
    VertexMap.SetIsVertex(iPoint,isVertex);
  }

  /*!
   * \brief Set the value of the solution in the previous BGS subiteration.
   */
  void Set_BGSSolution_k() override;

  /*!
   * \brief Get the value of the solution in the previous BGS subiteration.
   * \param[out] val_solution - solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution_k(Idx_t iPoint, Idx_t iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Solution_BGS_k(iPoint,iDim);
  }

};
