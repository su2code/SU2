/*!
 * \file CAdjTurbVariable.hpp
 * \brief Main class for defining the variables of the adjoint turbulence model.
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

#include "CVariable.hpp"

/*!
 * \class CAdjTurbVariable
 * \brief Main class for defining the variables of the adjoint turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */
class CAdjTurbVariable : public CVariable {
protected:
  su2double *dmuT_dUTvar;       /*!< \brief Sensitivity of eddy viscosity to mean flow and turbulence vars. */
  su2double **dRTstar_dUTvar;   /*!< \brief Sensitivity of modified turbulence residual (no boundary flux)
                                 to mean flow and turbulence vars. */
  su2double **dFT_dUTvar;   /*!< \brief Sensitivity of boundary flux
                             to mean flow and turbulence vars. */
  su2double *EddyViscSens;    /*!< \brief Eddy Viscosity Sensitivity. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CAdjTurbVariable(void);

  /*!
   * \overload
   * \param[in] val_psinu_inf - Value of the adjoint turbulence variable at the infinity (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAdjTurbVariable(su2double val_psinu_inf, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAdjTurbVariable(void);

  /*!
   * \brief Set the Eddy Viscosity Sensitivity of the problem.
   * \param[in] val_EddyViscSens - Eddy Viscosity Sensitivity.
   */
  inline void SetEddyViscSens(su2double *val_EddyViscSens, unsigned short numTotalVar) {
    for (unsigned short iVar = 0; iVar < numTotalVar; iVar++) EddyViscSens[iVar] = val_EddyViscSens[iVar];
  }

  /*!
   * \brief Get the Eddy Viscosity Sensitivity of the problem.
   * \return Pointer to the Eddy Viscosity Sensitivity.
   */
  inline su2double *GetEddyViscSens(void) {return EddyViscSens; }
};
