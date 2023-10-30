/*!
 * \file CSpeciesVariable.hpp
 * \brief Base class for defining the variables of the species transport model.
 * \author T. Kattmann
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

#include "CScalarVariable.hpp"

/*!
 * \class CSpeciesVariable
 * \brief Base class for defining the variables of the species transport.
 */
class CSpeciesVariable : public CScalarVariable {
 protected:
  MatrixType Diffusivity; /*!< \brief Matrix (nPoint,nVar) of mass diffusivities for scalar transport. */

 public:
  static constexpr size_t MAXNVAR = 20; /*!< \brief Max number of variables for static arrays. Increase, if necessary. */

  /*!
   * \brief Constructor of the class.
   * \param[in] species_inf - species variable values (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSpeciesVariable(const su2double* species_inf, unsigned long npoint, unsigned long ndim, unsigned long nvar,
                   const CConfig* config);

  /*!
   * \brief Set the value of the mass diffusivity
   * \param[in] val_diffusivity - the mass diffusivity.
   * \param[in] val_ivar        - eqn. index to the mass diffusivity.
   */
  inline void SetDiffusivity(unsigned long iPoint, su2double val_diffusivity, unsigned short val_ivar) {
    Diffusivity(iPoint, val_ivar) = val_diffusivity;
  }

  /*!
   * \brief Get the value of the mass diffusivity
   * \param[in] val_ivar - eqn. index to the mass diffusivity.
   * \return Value of the mass diffusivity
   */
  inline su2double GetDiffusivity(unsigned long iPoint, unsigned short val_ivar) const {
    return Diffusivity(iPoint, val_ivar);
  }

  /*!
   * \brief Get the value of the mass diffusivities
   * \return Pointer to the mass diffusivities
   */
  inline const su2double* GetDiffusivity(unsigned long iPoint) const { return Diffusivity[iPoint]; }
};
