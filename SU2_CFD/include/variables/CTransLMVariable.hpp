/*!
 * \file CTransLMVariable.hpp
 * \brief Declaration of the variables of the transition model.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

class CTransLMVariable : public CTurbVariable {
protected:
  su2double gamma_sep;

public:

  /*!
   * \brief Constructor of the class.
   */
  CTransLMVariable(void);

  /*!
   * \overload
   * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
   * \param[in] val_intermittency
   * \param[in] val_REth
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTransLMVariable(su2double val_nu_tilde, su2double val_intermittency, su2double val_REth, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTransLMVariable(void);

  /*!
   * \brief ________________.
   */
  inline su2double GetIntermittency(void) { return Solution[0]; }

  /*!
   * \brief ________________.
   * \param[in] gamma_sep_in
   */
  inline void SetGammaSep(su2double gamma_sep_in) {gamma_sep = gamma_sep_in;}

  /*!
   * \brief Correction for separation-induced transition.
   */
  inline void SetGammaEff(void) {Solution[0] = max(Solution[0], gamma_sep);}
};
