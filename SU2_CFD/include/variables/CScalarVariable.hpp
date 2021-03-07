/*!
 * \file CScalarVariable.hpp
 * \brief Main subroutines defining the variables for the  transported scalar model.
 * \author D. Mayer, T. Economon
 * \version 7.1.0 "Blackbird"
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
 * \class CScalarVariable
 * \brief Main class for defining the variables of scalar transport eqns.
 * \ingroup Scalar_Equations
 * \author T. Economon
 */
class CScalarVariable : public CVariable {
protected:
  MatrixType Diffusivity;  /*!< \brief Vector of mass diffusivities for scalar transport. */
  
public:
  /*!
   * \brief Constructor of the class.
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim   - Number of dimensions of the problem.
   * \param[in] nvar   - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CScalarVariable(unsigned long npoint,
                  unsigned long ndim,
                  unsigned long nvar,
                  CConfig       *config);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CScalarVariable() = default;
  
  /*!
   * \brief Set the value of the mass diffusivity
   * \param[in] val_diffusivity - the mass diffusivity.
   * \param[in] val_ivar        - eqn. index to the mass diffusivity.
   */
  inline void SetDiffusivity(unsigned long  iPoint,
                             su2double      val_diffusivity,
                             unsigned short val_ivar) {
    Diffusivity(iPoint,val_ivar) = val_diffusivity;
  }
  
  /*!
   * \brief Get the value of the mass diffusivity
   * \param[in] val_ivar - eqn. index to the mass diffusivity.
   * \return Value of the mass diffusivity
   */
  inline su2double GetDiffusivity(unsigned long  iPoint,
                                  unsigned short val_ivar) {
    return Diffusivity(iPoint,val_ivar);
  }
  
  /*!
   * \brief Get the value of the mass diffusivities
   * \return Pointer to the mass diffusivities
   */
  inline su2double* GetDiffusivity(unsigned long iPoint) {
    return Diffusivity[iPoint];
  }

};
