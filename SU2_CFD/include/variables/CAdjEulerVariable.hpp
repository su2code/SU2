/*!
 * \file CAdjEulerVariable.hpp
 * \brief Main class for defining the variables of the adjoint Euler solver.
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

#include "CVariable.hpp"

/*!
 * \class CAdjEulerVariable
 * \brief Main class for defining the variables of the adjoint Euler solver.
 * \ingroup Euler_Equations
 * \author F. Palacios, T. Economon
 */
class CAdjEulerVariable : public CVariable {
protected:
  MatrixType Psi;                /*!< \brief Vector of the adjoint variables. */
  MatrixType ForceProj_Vector;   /*!< \brief Vector d. */
  MatrixType ObjFuncSource;      /*!< \brief Vector containing objective function sensitivity for discrete adjoint. */
  MatrixType HB_Source;          /*!< \brief Harmonic balance source term. */

  CVectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  CVectorOfMatrix Gradient_Aux;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] psirho - Value of the adjoint density (initialization value).
   * \param[in] phi - Value of the adjoint velocity (initialization value).
   * \param[in] psie - Value of the adjoint energy (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjEulerVariable(su2double psirho, const su2double *phi, su2double psie,
                    unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjEulerVariable() override = default;

  /*!
   * \brief Set all the primitive variables for compressible flows.
   */
  bool SetPrimVar(unsigned long iPoint, su2double SharpEdge_Distance, bool check, CConfig *config) final;

  /*!
   * \brief Set the value of the adjoint velocity.
   * \param[in] val_phi - Value of the adjoint velocity.
   */
  inline void SetPhi_Old(unsigned long iPoint, const su2double *val_phi) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution_Old(iPoint,iDim+1)=val_phi[iDim];
  }

  /*!
   * \brief Set the value of the force projection vector.
   * \param[in] val_ForceProj_Vector - Pointer to the force projection vector.
   */
  inline void SetForceProj_Vector(unsigned long iPoint, const su2double *val_ForceProj_Vector) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) ForceProj_Vector(iPoint,iDim) = val_ForceProj_Vector[iDim];
  }

  /*!
   * \brief Set the value of the objective function source.
   * \param[in] val_ObjFuncSource - Pointer to the objective function source.
   */
  inline void SetObjFuncSource(unsigned long iPoint, const su2double *val_ObjFuncSource) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) ObjFuncSource(iPoint,iVar) = val_ObjFuncSource[iVar];
  }

  /*!
   * \brief Set the velocity vector from the old solution.
   * \param[in] val_velocity - Pointer to the velocity.
   */
  inline void SetVelocity_Old(unsigned long iPoint, const su2double *val_velocity) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Solution_Old(iPoint,iDim+1) = val_velocity[iDim]*Solution(iPoint,0);
  }

  /*!
   * \brief Set the momentum part of the truncation error to zero.
   * \param[in] iPoint - Point index.
   */
  inline void SetVel_ResTruncError_Zero(unsigned long iPoint) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Res_TruncError(iPoint,iDim+1) = 0.0;
  }

  /*!
   * \brief Specify a vector to set the velocity components of the solution. Multiplied by density for compressible cases.
   * \param[in] iPoint - Point index.
   * \param[in] val_vector - Pointer to the vector.
   */
  inline void SetVelSolutionVector(unsigned long iPoint, const su2double *val_vector) final {
    for (unsigned long iDim = 0; iDim < nDim; iDim++) Solution(iPoint, iDim+1) = GetDensity(iPoint) * val_vector[iDim];
  }

  /*!
   * \brief Get the value of the force projection vector.
   * \return Pointer to the force projection vector.
   */
  inline su2double *GetForceProj_Vector(unsigned long iPoint) final { return ForceProj_Vector[iPoint]; }

  /*!
   * \brief Get the value of the objective function source.
   * \param[in] val_SetObjFuncSource - Pointer to the objective function source.
   */
  inline su2double *GetObjFuncSource(unsigned long iPoint) final { return ObjFuncSource[iPoint]; }

  /*!
   * \brief Set the harmonic balance source term.
   * \param[in] iVar - Index of the variable.
   * \param[in] val_solution - Value of the harmonic balance source term. for the index <i>iVar</i>.
   */
  inline void SetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar, su2double val_source) final {
    HB_Source(iPoint,iVar) = val_source;
  }

  /*!
   * \brief Get the harmonic balance source term.
   * \param[in] iVar - Index of the variable.
   * \return Value of the harmonic balance source term for the index <i>iVar</i>.
   */
  inline su2double GetHarmonicBalance_Source(unsigned long iPoint, unsigned long iVar) const final {
    return HB_Source(iPoint,iVar);
  }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline CMatrixView<su2double> GetGradient_Reconstruction(unsigned long iPoint) final {
    return Gradient_Reconstruction[iPoint];
  }

  /*!
   * \brief Get the reconstruction gradient for variables at all points.
   * \return Reference to reconstruction gradient.
   */
  inline CVectorOfMatrix& GetGradient_Reconstruction() final { return Gradient_Reconstruction; }
  inline const CVectorOfMatrix& GetGradient_Reconstruction() const final { return Gradient_Reconstruction; }

};
