/*!
 * \file CDiscAdjFEABoundVariable.hpp
 * \brief Main class for defining the variables of the adjoint FEA solver at the boundary.
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

#include "CDiscAdjVariable.hpp"
#include "../../../Common/include/containers/CVertexMap.hpp"

/*!
 * \class CDiscAdjFEABoundVariable
 * \ingroup DiscAdj
 * \brief Main class for defining the variables on the FEA boundaries for adjoint applications.
 * \author R. Sanchez.
 * \version 8.0.0 "Harrier"
 */
class CDiscAdjFEABoundVariable final : public CDiscAdjVariable {
private:

  MatrixType FlowTraction_Sens;        /*!< \brief Adjoint of the flow tractions. */
  MatrixType SourceTerm_DispAdjoint;   /*!< \brief Source term applied into the displacement
                                                   adjoint coming from external solvers. */
  MatrixType SourceTerm_VelAdjoint;

  CVertexMap<unsigned> VertexMap;      /*!< \brief Object that controls accesses to the variables of this class. */

public:
  /*!
   * \overload
   * \param[in] sol - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEABoundVariable(const su2double *sol, unsigned long npoint, unsigned long ndim,
                           unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEABoundVariable() override = default;

  /*!
   * \brief Allocate member variables for points marked as vertex (via "Set_isVertex").
   * \param[in] config - Definition of the particular problem.
   */
  void AllocateBoundaryVariables(CConfig *config);

  /*!
   * \brief Set the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    FlowTraction_Sens(iPoint,iDim) = val;
  }

  /*!
   * \brief Get the FSI force sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetFlowTractionSensitivity(unsigned long iPoint, unsigned long iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return FlowTraction_Sens(iPoint,iDim);
  }

  /*!
   * \brief Set the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \param[in] val - value of the source term
   */
  inline void SetSourceTerm_DispAdjoint(unsigned long iPoint, unsigned long iDim, su2double val) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    SourceTerm_DispAdjoint(iPoint,iDim) = val;
  }
  inline void SetSourceTerm_VelAdjoint(unsigned long iPoint, unsigned long iDim, su2double val) override {
    if (!VertexMap.GetVertexIndex(iPoint)) return;
    SourceTerm_VelAdjoint(iPoint,iDim) = val;
  }

  /*!
   * \brief Get the source term applied into the displacement adjoint coming from external solvers
   * \param[in] iDim - spacial component
   * \return value of the source term
   */
  inline su2double GetSourceTerm_DispAdjoint(unsigned long iPoint, unsigned long iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return SourceTerm_DispAdjoint(iPoint,iDim);
  }
  inline su2double GetSourceTerm_VelAdjoint(unsigned long iPoint, unsigned long iDim) const override {
    if (!VertexMap.GetVertexIndex(iPoint)) return 0.0;
    return SourceTerm_VelAdjoint(iPoint,iDim);
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
