/*!
 * \file CTurbVariable.hpp
 * \brief Base class for defining the variables of the turbulence model.
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
 * \class CTurbVariable
 * \brief Base class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbVariable : public CVariable {
protected:
  VectorType muT;         /*!< \brief Eddy viscosity. */
  MatrixType HB_Source;   /*!< \brief Harmonic Balance source term. */

  VectorOfMatrix& Gradient_Reconstruction;  /*!< \brief Reference to the gradient of the primitive variables for MUSCL reconstruction for the convective term */
  VectorOfMatrix Gradient_Aux;              /*!< \brief Auxiliary structure to store a second gradient for reconstruction, if required. */

  VectorOfMatrix &Smatrix_Reconstruction;   /*!< \brief Geometry-based matrix for weighted least squares gradient calculations. */
  VectorOfMatrix Smatrix_Aux;               /*!< \brief Geometry-based matrix for weighted least squares gradient calculations. */

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbVariable() override = default;

  /*!
   * \brief Get the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \return the value of the eddy viscosity.
   */
  inline su2double GetmuT(unsigned long iPoint) const final { return muT(iPoint); }

  /*!
   * \brief Set the value of the eddy viscosity.
   * \param[in] iPoint - Point index.
   * \param[in] val_muT - Value of the eddy viscosity.
   */
  inline void SetmuT(unsigned long iPoint, su2double val_muT) final { muT(iPoint) = val_muT; }

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
   * \brief Get the reconstruction gradient for primitive variable at all points.
   * \return Reference to variable reconstruction gradient.
   */
  inline VectorOfMatrix& GetGradient_Reconstruction(void) final { return Gradient_Reconstruction; }

  /*!
   * \overload
   * \param[in] iPoint - Point index.
   * \param[in] iVar - Index of the variable.
   * \param[in] iDim - Index of the dimension.
   * \param[in] value - Value of the Smatrix.
   */
  inline void SetSmatrix_Reconstruction(unsigned long iPoint, unsigned long iDim, unsigned long jDim, su2double value) { Smatrix_Reconstruction(iPoint,iDim,jDim) = value; }

  /*!
   * \brief Get the value of the Smatrix entry for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \param[in] iDim - Index of the dimension.
   * \param[in] jDim - Index of the dimension.
   * \return Value of the Smatrix entry.
   */
  inline su2double GetSmatrix_Reconstruction(unsigned long iPoint, unsigned long iDim, unsigned long jDim) const { return Smatrix_Reconstruction(iPoint,iDim,jDim); }

  /*!
   * \brief Get the value of the Smatrix entry for least squares gradient calculations.
   * \param[in] iPoint - Point index.
   * \return Value of the Smatrix entry.
   */
  inline su2double **GetSmatrix_Reconstruction(unsigned long iPoint) { return Smatrix_Reconstruction[iPoint]; }

  /*!
   * \brief Get the value Smatrix for the entire domain.
   * \return Reference to the Smatrix.
   */
  inline VectorOfMatrix& GetSmatrix_Reconstruction(void) { return Smatrix_Reconstruction; }

  /*!
   * \brief A virtual member.
   */
  inline virtual const MatrixType& GetPrimitive(void) const { static const MatrixType dummy_mat; return dummy_mat; }

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double GetPrimitive(unsigned long iPoint, unsigned long iVar) const { return 0.0; }

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(unsigned long iPoint, unsigned long iVar, su2double val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual void SetPrimitive(unsigned long iPoint, const su2double *val_prim) {}

  /*!
   * \brief A virtual member.
   */
  inline virtual su2double *GetPrimitive(unsigned long iPoint) { return nullptr; }


  

};

