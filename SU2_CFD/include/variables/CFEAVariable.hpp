/*!
 * \file CFEAVariable.hpp
 * \brief Class for defining the variables of the FEM structural problem.
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

/*!
 * \class CFEAVariable
 * \brief Class for defining the variables of the FEM structural problem.
 * \ingroup Structural Finite Element Analysis Variables
 * \author F. Palacios, R. Sanchez.
 * \version 7.1.1 "Blackbird"
 */
class CFEAVariable : public CVariable {
protected:

  MatrixType Stress;                /*!< \brief Stress tensor. */

  MatrixType Residual_Ext_Body;     /*!< \brief Term of the residual due to body forces */

  VectorType VonMises_Stress;       /*!< \brief Von Mises stress. */

  MatrixType Solution_Pred;         /*!< \brief Predictor of the solution for FSI purposes */
  MatrixType Solution_Pred_Old;     /*!< \brief Predictor of the solution at time n for FSI purposes */

  MatrixType Reference_Geometry;    /*!< \brief Reference solution for optimization problems */

  MatrixType Prestretch;            /*!< \brief Prestretch geometry */

  /*!
   * \brief Constructor of the class.
   * \note This class is not supposed to be instantiated, it is only a building block for CFEABoundVariable
   * \param[in] val_fea - Values of the fea solution (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CFEAVariable(const su2double *val_fea, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CFEAVariable() override = default;

  /*!
   * \brief Get the value of the stress.
   * \return Value of the stress.
   */
  inline const su2double *GetStress_FEM(unsigned long iPoint) const final { return Stress[iPoint]; }

  /*!
   * \brief Set the value of the stress at the node
   * \param[in] iVar - index of the stress term
   * \param[in] val_stress - value of the stress
   */
  inline void SetStress_FEM(unsigned long iPoint, unsigned long iVar, su2double val_stress) final { Stress(iPoint,iVar) = val_stress; }

  /*!
   * \brief Add a certain value to the value of the stress at the node
   * \param[in] iVar - index of the stress term
   * \param[in] val_stress - value of the stress
   */
  inline void AddStress_FEM(unsigned long iPoint, unsigned long iVar, su2double val_stress) final { Stress(iPoint,iVar) += val_stress; }

  /*!
   * \brief Add body forces to the residual term.
   */
  inline void Add_BodyForces_Res(unsigned long iPoint, const su2double *val_bodyForce) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Residual_Ext_Body(iPoint,iVar) += val_bodyForce[iVar];
  }

  /*!
   * \brief Clear the surface load residual
   */
  inline void Clear_BodyForces_Res(unsigned long iPoint) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Residual_Ext_Body(iPoint,iVar) = 0.0;
  }

  /*!
   * \brief Get the body forces.
   */
  inline su2double Get_BodyForces_Res(unsigned long iPoint, unsigned long iVar) const final { return Residual_Ext_Body(iPoint,iVar); }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_pred - Pointer to the residual vector.
   */
  inline void SetSolution_Pred(unsigned long iPoint, const su2double *val_solution_pred) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Pred(iPoint,iVar) = val_solution_pred[iVar];
  }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline const su2double *GetSolution_Pred(unsigned long iPoint) const final { return Solution_Pred[iPoint]; }

  /*!
   * \brief Set the value of the old solution.
   * \param[in] val_solution_pred_old - Pointer to the residual vector.
   */
  inline void SetSolution_Pred_Old(unsigned long iPoint, const su2double *val_solution_pred_old) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Pred_Old(iPoint,iVar) = val_solution_pred_old[iVar];
  }

  /*!
   * \brief Get the solution at time n.
   * \return Pointer to the solution (at time n) vector.
   */
  inline const su2double *GetSolution_Pred_Old(unsigned long iPoint) const final { return Solution_Pred_Old[iPoint]; }

  /*!
   * \brief A virtual member.
   */
  inline void SetPrestretch(unsigned long iPoint, unsigned long iVar, su2double val_prestretch) final {
    Prestretch(iPoint,iVar) = val_prestretch;
  }

  /*!
   * \brief A virtual member.
   */
  inline const su2double *GetPrestretch(unsigned long iPoint) const final { return Prestretch[iPoint]; }

  /*!
   * \brief A virtual member.
   */
  inline su2double GetPrestretch(unsigned long iPoint, unsigned long iVar) const final { return Prestretch(iPoint,iVar); }

  /*!
   * \brief Set the value of the Von Mises stress.
   * \param[in] val_stress - Value of the Von Mises stress.
   */
  inline void SetVonMises_Stress(unsigned long iPoint, su2double val_stress) final { VonMises_Stress(iPoint) = val_stress; }

  /*!
   * \brief Get the value of the Von Mises stress.
   * \return Value of the Von Mises stress.
   */
  inline su2double GetVonMises_Stress(unsigned long iPoint) const final { return VonMises_Stress(iPoint); }

  /*!
   * \brief Set the reference geometry.
   * \return Pointer to the solution (at time n) vector.
   */
  inline void SetReference_Geometry(unsigned long iPoint, unsigned long iVar, su2double ref_geometry) final {
    Reference_Geometry(iPoint,iVar) = ref_geometry;
  }

  /*!
   * \brief Get the pointer to the reference geometry
   */
  inline const su2double* GetReference_Geometry(unsigned long iPoint) const final { return Reference_Geometry[iPoint]; }

  /*!
   * \brief Register the variables in the solution time_n array as input/output variable.
   * \param[in] input - input or output variables.
   */
  void Register_femSolution_time_n(bool input, bool push_index) final;

};
