/*!
 * \file CAdjEulerVariable.hpp
 * \brief Main class for defining the variables of the adjoint Euler solver.
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
  MatrixType IntBoundary_Jump;   /*!< \brief Interior boundary jump vector. */
  MatrixType HB_Source;          /*!< \brief Harmonic balance source term. */

  VectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  VectorOfMatrix Gradient_Aux;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

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
  virtual ~CAdjEulerVariable() = default;

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
   * \brief Set the value of the interior boundary jump vector vector.
   * \param[in] val_IntBoundary_Jump - Pointer to the interior boundary jump vector.
   */
  inline void SetIntBoundary_Jump(unsigned long iPoint, const su2double *val_IntBoundary_Jump) final {
    for (unsigned long iVar = 0; iVar < nVar; iVar++) IntBoundary_Jump(iPoint,iVar) = val_IntBoundary_Jump[iVar];
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
   * \brief Get the value of the force projection vector.
   * \return Pointer to the force projection vector.
   */
  inline su2double *GetIntBoundary_Jump(unsigned long iPoint) final { return IntBoundary_Jump[iPoint]; }

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
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \return Value of the reconstruction variables gradient at a node.
   */
  inline su2double GetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim) const final {
    return Gradient_Reconstruction(iPoint,iVar,iDim);
  }

  /*!
   * \brief Get the value of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \param[in] iVar   - Index of the variable.
   * \param[in] iDim   - Index of the dimension.
   * \param[in] value  - Value of the reconstruction gradient component.
   */
  inline void SetGradient_Reconstruction(unsigned long iPoint, unsigned long iVar, unsigned long iDim, su2double value) final {
    Gradient_Reconstruction(iPoint,iVar,iDim) = value;
  }

  /*!
   * \brief Get the array of the reconstruction variables gradient at a node.
   * \param[in] iPoint - Index of the current node.
   * \return Array of the reconstruction variables gradient at a node.
   */
  inline su2double **GetGradient_Reconstruction(unsigned long iPoint) final { return Gradient_Reconstruction[iPoint]; }

  /*!
   * \brief Get the reconstruction gradient for variables at all points.
   * \return Reference to reconstruction gradient.
   */
  inline VectorOfMatrix& GetGradient_Reconstruction(void) final { return Gradient_Reconstruction; }

};
