/*!
 * \file CDiscAdjFEAVariable.hpp
 * \brief Main class for defining the variables of the adjoint FEA solver.
 * \author T. Albring, R. Sanchez.
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
 * \class CDiscAdjFEAVariable
 * \brief Main class for defining the variables of the adjoint solver.
 * \ingroup Discrete_Adjoint
 * \author T. Albring, R. Sanchez.
 * \version 7.0.0 "Blackbird"
 */
class CDiscAdjFEAVariable : public CVariable {
protected:
  MatrixType Sensitivity; /* Vector holding the derivative of target functional with respect to the coordinates at this node*/
  MatrixType Solution_Direct;

  MatrixType Dynamic_Derivative;
  MatrixType Dynamic_Derivative_n;
  MatrixType Dynamic_Derivative_Vel;
  MatrixType Dynamic_Derivative_Vel_n;
  MatrixType Dynamic_Derivative_Accel;
  MatrixType Dynamic_Derivative_Accel_n;

  MatrixType Solution_Vel;
  MatrixType Solution_Accel;

  MatrixType Solution_Vel_time_n;
  MatrixType Solution_Accel_time_n;

  MatrixType Solution_Old_Vel;
  MatrixType Solution_Old_Accel;

  MatrixType Solution_Direct_Vel;
  MatrixType Solution_Direct_Accel;

  MatrixType Cross_Term_Derivative;
  MatrixType Geometry_CrossTerm_Derivative;

  MatrixType Solution_BGS;

  /*!
   * \brief Constructor of the class.
   * \param[in] disp - Pointer to the adjoint value (initialization value).
   * \param[in] vel - Pointer to the adjoint value (initialization value).
   * \param[in] accel - Pointer to the adjoint value (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] unsteady - Allocate velocity and acceleration.
   * \param[in] config - Definition of the particular problem.
   */
  CDiscAdjFEAVariable(const su2double *disp, const su2double *vel, const su2double *accel,
                      unsigned long npoint, unsigned long ndim, unsigned long nvar, bool unsteady, CConfig *config);

public:
  /*!
   * \brief Destructor of the class.
   */
  ~CDiscAdjFEAVariable() = default;

  /*!
   * \brief Set the sensitivity at the node
   * \param[in] iDim - spacial component
   * \param[in] val - value of the Sensitivity
   */
  inline void SetSensitivity(unsigned long iPoint, unsigned long iDim, su2double val) final { Sensitivity(iPoint,iDim) = val; }

  /*!
   * \brief Get the Sensitivity at the node
   * \param[in] iDim - spacial component
   * \return value of the Sensitivity
   */
  inline su2double GetSensitivity(unsigned long iPoint, unsigned long iDim) const final { return Sensitivity(iPoint,iDim);}

  inline void SetDynamic_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative(iPoint,iVar) = der;
  }

  inline void SetDynamic_Derivative_n(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_n(iPoint,iVar) = der;
  }

  inline su2double GetDynamic_Derivative(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative(iPoint,iVar);
  }

  inline su2double GetDynamic_Derivative_n(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_n(iPoint,iVar);
  }

  inline void SetDynamic_Derivative_Vel(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_Vel(iPoint,iVar) = der;
  }

  inline void SetDynamic_Derivative_Vel_n(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_Vel_n(iPoint,iVar) = der;
  }

  inline su2double GetDynamic_Derivative_Vel(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_Vel(iPoint,iVar);
  }

  inline su2double GetDynamic_Derivative_Vel_n(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_Vel_n(iPoint,iVar);
  }

  inline void SetDynamic_Derivative_Accel(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_Accel(iPoint,iVar) = der;
  }

  inline void SetDynamic_Derivative_Accel_n(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Dynamic_Derivative_Accel_n(iPoint,iVar) = der;
  }

  inline su2double GetDynamic_Derivative_Accel(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_Accel(iPoint,iVar);
  }

  inline su2double GetDynamic_Derivative_Accel_n(unsigned long iPoint, unsigned long iVar) const final {
    return Dynamic_Derivative_Accel_n(iPoint,iVar);
  }

  inline void SetSolution_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline void SetSolution_Vel_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct_Vel(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline void SetSolution_Accel_Direct(unsigned long iPoint, const su2double *val_solution_direct) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Direct_Accel(iPoint,iVar) = val_solution_direct[iVar];
  }

  inline su2double* GetSolution_Direct(unsigned long iPoint) final { return Solution_Direct[iPoint]; }

  inline su2double* GetSolution_Vel_Direct(unsigned long iPoint) final { return Solution_Direct_Vel[iPoint]; }

  inline su2double* GetSolution_Accel_Direct(unsigned long iPoint) final { return Solution_Direct_Accel[iPoint]; }

  inline su2double GetSolution_Old_Vel(unsigned long iPoint, unsigned long iVar) const final { return Solution_Old_Vel(iPoint,iVar); }

  inline su2double GetSolution_Old_Accel(unsigned long iPoint, unsigned long iVar) const final { return Solution_Old_Accel(iPoint,iVar); }

  /*!
   * \brief Get the acceleration (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_Accel(unsigned long iPoint, unsigned long iVar) const final { return Solution_Accel(iPoint,iVar); }

  /*!
   * \brief Get the acceleration of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Accel_time_n(unsigned long iPoint, unsigned long iVar) const final { return Solution_Accel_time_n(iPoint,iVar); }

  /*!
   * \brief Get the velocity (Structural Analysis).
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetSolution_Vel(unsigned long iPoint, unsigned long iVar) const final { return Solution_Vel(iPoint,iVar); }

  /*!
   * \brief Get the velocity of the nodes (Structural Analysis) at time n.
   * \param[in] iVar - Index of the variable.
   * \return Pointer to the old solution vector.
   */
  inline su2double GetSolution_Vel_time_n(unsigned long iPoint, unsigned long iVar) const final { return Solution_Vel_time_n(iPoint,iVar); }

  /*!
   * \brief Set the value of the acceleration (Structural Analysis - adjoint).
   * \param[in] val_solution - Solution of the problem (acceleration).
   */
  inline void SetSolution_Accel(unsigned long iPoint, const su2double *val_solution_accel) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution_Accel(iPoint,iVar) = val_solution_accel[iVar];
  }

  /*!
   * \brief Set the value of the velocity (Structural Analysis - adjoint).
   * \param[in] val_solution - Solution of the problem (velocity).
   */
  inline void SetSolution_Vel(unsigned long iPoint, const su2double *val_solution_vel) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Vel(iPoint,iVar) = val_solution_vel[iVar];
  }

  /*!
   * \brief Set the value of the adjoint acceleration (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Accel_time_n(unsigned long iPoint, const su2double *val_solution_accel_time_n) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Accel_time_n(iPoint,iVar) = val_solution_accel_time_n[iVar];
  }

  /*!
   * \brief Set the value of the adjoint velocity (Structural Analysis) at time n.
   * \param[in] val_solution_old - Pointer to the residual vector.
   */
  inline void SetSolution_Vel_time_n(unsigned long iPoint, const su2double *val_solution_vel_time_n) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) Solution_Vel_time_n(iPoint,iVar) = val_solution_vel_time_n[iVar];
  }

  /*!
   * \brief Set the value of the old acceleration (Structural Analysis - adjoint).
   */
  void Set_OldSolution_Accel() final;

  /*!
   * \brief Set the value of the old velocity (Structural Analysis - adjoint).
   */
  void Set_OldSolution_Vel() final;

  /*!
   * \brief Set the contribution of crossed terms into the derivative.
   */
  inline void SetCross_Term_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Cross_Term_Derivative(iPoint,iVar) = der;
  }

  /*!
   * \brief Get the contribution of crossed terms into the derivative.
   */
  inline su2double GetCross_Term_Derivative(unsigned long iPoint, unsigned long iVar) const final { return Cross_Term_Derivative(iPoint,iVar); }

  /*!
   * \brief A virtual member. Get the geometry solution.
   * \param[in] iVar - Index of the variable.
   * \return Value of the solution for the index <i>iVar</i>.
   */
  inline su2double GetGeometry_CrossTerm_Derivative(unsigned long iPoint, unsigned long iVar) const final {
    return Geometry_CrossTerm_Derivative(iPoint,iVar);
  }

  /*!
   * \brief A virtual member. Set the value of the mesh solution (adjoint).
   * \param[in] der - cross term derivative.
   */
  inline void SetGeometry_CrossTerm_Derivative(unsigned long iPoint, unsigned long iVar, su2double der) final {
    Geometry_CrossTerm_Derivative(iPoint,iVar) = der;
  }

  /*!
   * \brief Set the value of the adjoint solution in the current BGS subiteration.
   */
  inline void Set_BGSSolution(unsigned long iPoint, unsigned long iDim, su2double val_solution) final { Solution_BGS(iPoint,iDim) = val_solution; }

  /*!
   * \brief Get the value of the adjoint solution in the previous BGS subiteration.
   * \param[out] val_solution - adjoint solution in the previous BGS subiteration.
   */
  inline su2double Get_BGSSolution(unsigned long iPoint, unsigned long iDim) const final { return Solution_BGS(iPoint,iDim); }

};
