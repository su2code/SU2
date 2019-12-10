/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
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

#include "CTurbVariable.hpp"

/*!
 * \class CTransLMVariable
 * \brief Transition model variables.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTransLMVariable final : public CTurbVariable {
protected:
  VectorType gamma_sep;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_intermittency
   * \param[in] val_REth
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] constants -
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double intermittency, su2double REth, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransLMVariable() = default;

  /*!
   * \brief ________________.
   */
  inline su2double GetIntermittency(unsigned long iPoint) const override { return Solution(iPoint,0); }

  /*!
   * \brief ________________.
   * \param[in] gamma_sep_in
   */
  inline void SetGammaSep(unsigned long iPoint, su2double gamma_sep_in) override { gamma_sep(iPoint) = gamma_sep_in; }

  /*!
   * \brief Correction for separation-induced transition.
   */
  inline void SetGammaEff(unsigned long iPoint) override { Solution(iPoint,0) = max(Solution(iPoint,0), gamma_sep(iPoint)); }
};
