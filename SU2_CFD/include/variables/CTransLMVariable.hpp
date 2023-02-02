/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author F. Palacios, T. Economon
 * \version 7.5.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
 * \author A. Bueno, S. Kang.
 */

class CTransLMVariable final : public CTurbVariable {
protected:
  VectorType Intermittency_Eff;
  VectorType Intermittency_Sep;
  VectorType Corr_Rec;
  VectorType Re_t;

public:
  /*!
   * \brief Constructor of the class.
   * \param[in] Intermittency - intermittency(gamma) (initialization value).
   * \param[in] ReThetaT - momentum thickness Reynolds number(ReThetaT)(initialization value).
   * \param[in] gammaSep - separation intermittency(gamma) (initialization value).
   * \param[in] gammaEff - effective intermittency(gamma) (initialization value).
   * \param[in] npoint - Number of points/nodes/vertices in the domain.
   * \param[in] ndim - Number of dimensions of the problem.
   * \param[in] nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double Intermittency, su2double ReThetaT, su2double gammaSep, su2double gammaEff, unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransLMVariable() override = default;

  /*!
   * \brief Set Separation intermittency.
   */
  void SetIntermittencySep(unsigned long iPoint, su2double val_Intermittency_sep) override;
  
  /*!
   * \brief Set Effective intermittency.
   */
  void SetIntermittencyEff(unsigned long iPoint, su2double val_Intermittency_sep) override;

  /*!
   * \brief Set Value of Transition Momentum Thickness Reynolds number from correlations.
   */
  void SetCorr_Rec(unsigned long iPoint, su2double val_Corr_Rec) override;

  /*!
   * \brief Set Value of Momentum Thickness Reynolds number from correlations (substitute to the second equation of original LM model).
   */
  void SetRe_t(unsigned long iPoint, su2double val_Re_t) override;

  /*!
   * \brief Calculate effective intermittency.
   */
  inline su2double GetIntermittencyEff(unsigned long iPoint) const override { return Intermittency_Eff(iPoint); }

  /*!
   * \brief Value of separation intermittency.
   */
  inline su2double GetIntermittencySep(unsigned long iPoint) const override { return Intermittency_Sep(iPoint); }

  /*!
   * \brief Get Value of Transition Momentum Thickness Reynolds number from correlations.
   */
  inline su2double GetCorr_Rec(unsigned long iPoint) const override { return Corr_Rec(iPoint); }

  /*!
   * \brief Get Value of Momentum Thickness Reynolds number from correlations (substitute to the second equation of original LM model).
   */
  inline su2double GetRe_t(unsigned long iPoint) const override { return Re_t(iPoint); }
};
