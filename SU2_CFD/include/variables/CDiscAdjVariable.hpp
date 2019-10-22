/*!
 * \file CDiscAdjVariable.hpp
 * \brief Main class for defining the variables of the adjoint solver.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CDiscAdjVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring.
 */
class CDiscAdjVariable final : public CVariable {
private:
  MatrixType Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  MatrixType Solution_Direct;
  MatrixType DualTime_Derivative;
  MatrixType DualTime_Derivative_n;

  MatrixType Cross_Term_Derivative;
  MatrixType Geometry_CrossTerm_Derivative;
  MatrixType Geometry_CrossTerm_Derivative_Flow;

  MatrixType Solution_Geometry;
  MatrixType Solution_Geometry_Old;
  MatrixType Geometry_Direct;

  MatrixType Solution_BGS;
  MatrixType Solution_Geometry_BGS_k;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] sol - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjVariable(const su2double* sol, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjVariable() = default;

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) override { Sensitivity(iPoint,iDim) = val;}

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const override { return Sensitivity(iPoint,iDim); }

  inline void SetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) override { DualTime_Derivative(iPoint,iVar) = der; }

  inline void SetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar, su2double der) override { DualTime_Derivative_n(iPoint,iVar) = der; }

  inline su2double GetDual_Time_Derivative(unsigned long iPoint, unsigned long iVar) const override { return DualTime_Derivative(iPoint,iVar); }

  inline su2double GetDual_Time_Derivative_n(unsigned long iPoint, unsigned long iVar) const override { return DualTime_Derivative_n(iPoint,iVar); }

  inline void SetSolution_Direct(unsigned long iPoint, const su2double *val_solution_direct) override {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_Direct(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct(unsigned long iPoint) override { return Solution_Direct[iPoint]; }

  /*!
   * \brief Set the restart geometry (coordinate of the converged solution)
   * \param[in] val_geometry_direct - Value of the restart coordinate.
   */
  inline void SetGeometry_Direct(unsigned long iPoint, const su2double *val_geometry_direct) override {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Geometry_Direct(iPoint,iDim) = val_geometry_direct[iDim];
  }

  /*!
   * \brief Get the restart geometry (coordinate of the converged solution).
   * \return Pointer to the restart coordinate vector.
   */
  inline su2double *GetGeometry_Direct(unsigned long iPoint) override { return Geometry_Direct[iPoint]; }

  /*!
   * \brief Get the restart geometry (coordinate of the converged solution).
   * \return Coordinate iDim of the geometry_direct vector.
   */
  inline su2double GetGeometry_Direct(unsigned long iPoint, unsigned long iDim) const override { return Geometry_Direct(iPoint,iDim); }

  /*!
   * \brief Get the geometry solution.
   * \param[in] iDim - Index of the coordinate.
   * \return Value of the solution for the index <i>iDim</i>.
   */
  inline su2double GetSolution_Geometry(unsigned long iPoint, unsigned long iDim) const override { return Solution_Geometry(iPoint,iDim); }

  /*!
   * \brief Set the value of the mesh solution (adjoint).
   * \param[in] val_solution_geometry - Solution of the problem (acceleration).
   */
  inline void SetSolution_Geometry(unsigned long iPoint, const su2double *val_solution_geometry) override {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Geometry(iPoint,iDim) = val_solution_geometry[iDim];
  }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] val_solution_geometry - Solution of the problem (acceleration).
   */
  inline void SetSolution_Geometry(unsigned long iPoint, unsigned long iVar, su2double val_solution_geometry) override {
    Solution_Geometry(iPoint,iVar) = val_solution_geometry;
  }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetGeometry_CrossTerm_Derivative(unsigned long iPoint, unsigned long iVar) const override {
    return Geometry_CrossTerm_Derivative(iPoint,iVar);
  }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] der - cross term derivative.
   */
  inline void SetGeometry_CrossTerm_Derivative(unsigned long iPoint, unsigned long iDim, su2double der) override {
    Geometry_CrossTerm_Derivative(iPoint,iDim) = der;
  }

  /*!
   * \brief Get the mesh cross term derivative from the flow solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetGeometry_CrossTerm_Derivative_Flow(unsigned long iPoint, unsigned long iVar) const override {
    return Geometry_CrossTerm_Derivative_Flow(iPoint,iVar);
  }

  /*!
   * \brief Set the value of the mesh cross term derivative from the flow solution (adjoint).
   * \param[in] der - cross term derivative.
   */
  inline void SetGeometry_CrossTerm_Derivative_Flow(unsigned long iPoint, unsigned long iDim, su2double der) override {
    Geometry_CrossTerm_Derivative_Flow(iPoint,iDim) = der;
  }

  /*!
   * \brief Set the value of the mesh solution (adjoint).
   */
  void Set_OldSolution_Geometry() override;

  /*!
   * \brief Get the value of the old geometry solution (adjoint).
   * \param[out] val_solution - old adjoint solution for coordinate iDim
   */
  inline su2double Get_OldSolution_Geometry(unsigned long iPoint, unsigned long iDim) const override {
    return Solution_Geometry_Old(iPoint,iDim);
  }

  /*!
   * \brief Set the value of the adjoint solution in the current BGS subiteration.
   */
  inline void Set_BGSSolution(unsigned long iPoint, unsigned long iDim, su2double val_solution) override {
    Solution_BGS(iPoint,iDim) = val_solution;
  }

  /*!
   * \brief Get the value of the adjoint solution in the previous BGS subiteration.
   * \param[out] val_solution - adjoint solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution(unsigned long iPoint, unsigned long iDim) const override { return Solution_BGS(iPoint,iDim);}

  /*!
   * \brief Set the contribution of crossed terms into the derivative.
   */
  inline void SetCross_Term_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) override {
    Cross_Term_Derivative(iPoint,iVar) = der;
  }

  /*!
   * \brief Get the contribution of crossed terms into the derivative.
   */
  inline su2double GetCross_Term_Derivative(unsigned long iPoint, unsigned long iVar) const override {
    return Cross_Term_Derivative(iPoint,iVar);
  }

};
