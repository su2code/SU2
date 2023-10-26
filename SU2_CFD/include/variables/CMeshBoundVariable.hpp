/*!
 * \file CMeshBoundVariable.hpp
 * \brief Declaration and inlines of the class
 *        to define the variables of the mesh movement at the moving boundaries.
 * \author Ruben Sanchez
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

#include "CMeshVariable.hpp"
#include "../../../Common/include/containers/CVertexMap.hpp"

class CMeshBoundVariable final : public CMeshVariable {
private:

  MatrixType Boundary_Displacement;  /*!< \brief Store the reference coordinates of the mesh. */
  MatrixType Boundary_Velocity;      /*!< \brief Store the boundary velocities of the mesh. */
  CVertexMap<unsigned> VertexMap;    /*!< \brief Object that controls accesses to the variables of this class. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_coor - Values of the coordinates (initialization value).
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CMeshBoundVariable(unsigned long npoint, unsigned long ndim, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CMeshBoundVariable() override = default;

  /*!
   * \brief Allocate member variables for points marked as vertex (via "Set_isVertex").
   * \param[in] config - Definition of the particular problem.
   */
  void AllocateBoundaryVariables(CConfig *config);

  /*!
   * \brief Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary displacement.
   */
  inline su2double GetBound_Disp(unsigned long iPoint, unsigned long iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Boundary_Displacement(iPoint,iDim);
  }

  /*!
   * \brief Set the boundary displacements.
   * \param[in] val_BoundDisp - Pointer to the boundary displacements.
   */
  inline void SetBound_Disp(unsigned long iPoint, const su2double *val_BoundDisp) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Boundary_Displacement(iPoint,iDim) = val_BoundDisp[iDim];
  }

  /*!
   * \brief Set the boundary displacement.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundDisp - Value of the boundary displacements.
   */
  inline void SetBound_Disp(unsigned long iPoint, unsigned long iDim, su2double val_BoundDisp) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    Boundary_Displacement(iPoint,iDim) = val_BoundDisp;
  }

  /*!
   * \brief Get the value of the displacement imposed at the boundary.
   * \return Value of the boundary velocity.
   */
  inline su2double GetBound_Vel(unsigned long iPoint, unsigned long iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Boundary_Velocity(iPoint,iDim);
  }

    /*!
   * \brief Set the boundary displacements.
   * \param[in] val_BoundVel - Pointer to the boundary velocities.
   */
  inline void SetBound_Vel(unsigned long iPoint, const su2double *val_BoundVel) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Boundary_Velocity(iPoint,iDim) = val_BoundVel[iDim];
  }

  /*!
   * \brief Set the boundary velocity.
   * \param[in] iDim - Index of the dimension of interest.
   * \param[in] val_BoundVel - Value of the boundary velocities.
   */
  inline void SetBound_Vel(unsigned long iPoint, unsigned long iDim, su2double val_BoundVel) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    Boundary_Velocity(iPoint,iDim) = val_BoundVel;
  }

  /*!
   * \brief Register the boundary displacements of the mesh.
   */
  void Register_BoundDisp() override;

  /*!
   * \brief Recover the value of the adjoint of the boundary displacements.
   */
  inline void GetAdjoint_BoundDisp(unsigned long iPoint, su2double *adj_disp) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iVar = 0; iVar < nVar; iVar++) {
        adj_disp[iVar] = SU2_TYPE::GetDerivative(Boundary_Displacement(iPoint,iVar));
    }
  }

  /*!
   * \brief Get whether a node is on the boundary
   */
  inline bool Get_isVertex(unsigned long iPoint) const override {
    return VertexMap.GetIsVertex(iPoint);
  }

  /*!
   * \brief Set whether a node is on the boundary
   */
  inline void Set_isVertex(unsigned long iPoint, bool isVertex) override {
    VertexMap.SetIsVertex(iPoint,isVertex);
  }

};
