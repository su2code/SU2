/*!
 * \file CFEABoundVariable.hpp
 * \brief Class for defining the variables on the FEA boundaries for FSI applications.
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

#include "CFEAVariable.hpp"
#include "../../../Common/include/containers/CVertexMap.hpp"

/*!
 * \class CFEABoundVariable
 * \brief Class that adds storage of boundary variables (tractions) to CFEAVariable.
 * \note Member variables are allocated only for points marked as "vertex" i.e. on a boundary.
 * A map is constructed so that variables can be referenced by iPoint instead of iVertex.
 * \ingroup Structural Finite Element Analysis Variables
 * \author R. Sanchez.
 * \version 8.0.0 "Harrier"
 */
class CFEABoundVariable final : public CFEAVariable {
protected:

  MatrixType FlowTraction;         /*!< \brief Traction from the fluid field. */
  MatrixType FlowTraction_n;       /*!< \brief Traction from the fluid field at time n. */

  MatrixType Residual_Ext_Surf;    /*!< \brief Term of the residual due to external forces. */
  MatrixType Residual_Ext_Surf_n;  /*!< \brief Term of the residual due to external forces at time n. */

  CVertexMap<unsigned> VertexMap;  /*!< \brief Object that controls accesses to the variables of this class. */

  bool fsi_analysis = false;       /*!< \brief If flow tractions are available. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEABoundVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CFEABoundVariable() override = default;

  /*!
   * \brief Allocate member variables for points marked as vertex (via "Set_isVertex").
   * \param[in] config - Definition of the particular problem.
   */
  void AllocateBoundaryVariables(CConfig *config);

  /*!
   * \brief Add surface load to the residual term
   */
  inline void Add_SurfaceLoad_Res(unsigned long iPoint, const su2double *val_surfForce) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Residual_Ext_Surf(iPoint,iVar) += val_surfForce[iVar];
  }

  /*!
   * \brief Get the residual term due to surface load
   */
  inline su2double Get_SurfaceLoad_Res(unsigned long iPoint, unsigned long iVar) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Residual_Ext_Surf(iPoint,iVar);
  }

  /*!
   * \brief Clear the surface load residual
   */
  inline void Clear_SurfaceLoad_Res() override;

  /*!
   * \brief Store the surface load as the load for the previous time step.
   */
  void Set_SurfaceLoad_Res_n() override;

  /*!
   * \brief Get the surface load from the previous time step.
   */
  inline su2double Get_SurfaceLoad_Res_n(unsigned long iPoint, unsigned long iVar) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return Residual_Ext_Surf_n(iPoint,iVar);
  }

  /*!
   * \brief Set the flow traction at a node on the structural side
   */
  inline void Set_FlowTraction(unsigned long iPoint, const su2double *val_flowTraction) override {
    if (!fsi_analysis) return;
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iVar = 0; iVar < nVar; iVar++) FlowTraction(iPoint,iVar) = val_flowTraction[iVar];
  }

  /*!
   * \brief Add a value to the flow traction at a node on the structural side
   */
  inline void Add_FlowTraction(unsigned long iPoint, const su2double *val_flowTraction) override {
    if (!fsi_analysis) return;
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    for (unsigned long iVar = 0; iVar < nVar; iVar++) FlowTraction(iPoint,iVar) += val_flowTraction[iVar];
  }

  /*!
   * \brief Get the residual term due to the flow traction
   */
  inline su2double Get_FlowTraction(unsigned long iPoint, unsigned long iVar) const override {
    if (!fsi_analysis) return 0.0;
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return FlowTraction(iPoint,iVar);
  }

  /*!
   * \brief Set the value of the flow traction at the previous time step.
   */
  void Set_FlowTraction_n() override;

  /*!
   * \brief Retrieve the value of the flow traction from the previous time step.
   */
  inline su2double Get_FlowTraction_n(unsigned long iPoint, unsigned long iVar) const override {
    if (!fsi_analysis) return 0.0;
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return FlowTraction_n(iPoint,iVar);
  }

  /*!
   * \brief Register the flow tractions as input variable.
   */
  void RegisterFlowTraction(bool reset) override;

  /*!
   * \brief Extract the flow traction derivatives.
   */
  inline su2double ExtractFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim) const override {
    if (!fsi_analysis) return 0.0;
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return SU2_TYPE::GetDerivative(FlowTraction(iPoint,iDim));
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
