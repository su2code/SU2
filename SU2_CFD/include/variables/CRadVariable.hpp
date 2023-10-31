/*!
 * \file CRadVariable.hpp
 * \brief Class for defining the variables of the radiation solver.
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

#include "CVariable.hpp"

class CRadVariable : public CVariable {

private:

protected:

  MatrixType Radiative_SourceTerm;
  su2vector<bool> Vol_HeatSource;    /*!< \brief Volumetric heat source. */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] npoint  - Number of points in the problem.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CRadVariable(unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CRadVariable(void) override = default;

  /*!
   * \brief Get the Radiative source term at the node
   * \return Radiative source term for the energy equation
   */
  inline const su2double *GetRadiative_SourceTerm(unsigned long iPoint) const final { return Radiative_SourceTerm[iPoint];}

  /*!
   * \brief Set the Radiative source term at the node
   * \param[in] val_RadSourceTerm - value of the radiative source term
   */
  inline void SetRadiative_SourceTerm(unsigned long iPoint, unsigned long iVar, su2double val_RadSourceTerm) final {
    Radiative_SourceTerm(iPoint, iVar) = val_RadSourceTerm;
  }

  /*!
   * \brief Get whether a volumetric heat source is to be introduced in point iPoint
   * \return Bool, determines if this point introduces volumetric heat
   */
  inline bool GetVol_HeatSource(unsigned long iPoint) const final { return Vol_HeatSource(iPoint);}

  /*!
   * \brief Set as true a volumetric heat source for point iPoint
   */
  inline void SetVol_HeatSource(unsigned long iPoint) { Vol_HeatSource(iPoint) = true;}

  /*!
   * \brief Reset as false a volumetric heat source for all points
   */
  inline void ResetVol_HeatSource(void) { Vol_HeatSource = false;}

};
