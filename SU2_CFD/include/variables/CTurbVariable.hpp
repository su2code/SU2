/*!
 * \file CTurbVariable.hpp
 * \brief Base class for defining the variables of the turbulence model.
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

#include "CScalarVariable.hpp"

/*!
 * \class CTurbVariable
 * \brief Base class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CTurbVariable : public CScalarVariable {
protected:
  VectorType muT; /*!< \brief Eddy viscosity. */

public:
  static constexpr size_t MAXNVAR = 2;
  VectorType turb_index;
  VectorType intermittency;         /*!< \brief Value of the intermittency for the trans. model. */

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
    * \brief Set the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \param[in] val_turb_index - Value of the turbulence index.
   */
  inline void SetTurbIndex(unsigned long iPoint, su2double val_turb_index) final { turb_index(iPoint) = val_turb_index; }

  /*!
   * \brief Get the value of the turbulence index.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the turbulence index.
   */
  inline su2double GetTurbIndex(unsigned long iPoint) const final { return turb_index(iPoint); }

  /*!
   * \brief Get the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \return Value of the intermittency of the transition model.
   */
  inline su2double GetIntermittency(unsigned long iPoint) const final { return intermittency(iPoint); }

  /*!
   * \brief Set the intermittency of the transition model.
   * \param[in] iPoint - Point index.
   * \param[in] val_intermittency - New value of the intermittency.
   */
  inline void SetIntermittency(unsigned long iPoint, su2double val_intermittency) final { intermittency(iPoint) = val_intermittency; }

};

