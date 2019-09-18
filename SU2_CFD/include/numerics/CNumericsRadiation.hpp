/*!
 * \file CNumericsRadiation.hpp
 * \brief Declaration and inlines of the class to compute
 *        residual terms in radiation problems.
 * \author Ruben Sanchez
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

#include "../numerics_structure.hpp"

class CNumericsRadiation : public CNumerics {
 private:

 protected:

  bool implicit, incompressible;
  su2double *RadVar_i,        /*!< \brief Vector of radiation variables at point i. */
  *RadVar_j;                  /*!< \brief Vector of radiation variables at point j. */
  su2double **RadVar_Grad_i,  /*!< \brief Gradient of turbulent variables at point i. */
  **RadVar_Grad_j;            /*!< \brief Gradient of turbulent variables at point j. */
  su2double Absorption_Coeff; /*!< \brief Absorption coefficient. */
  su2double Scattering_Coeff; /*!< \brief Scattering coefficient. */

  su2double Temperature_Ref;  /*!< \brief Reference temperature for redimensionalization of P1 solver. */

 public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CNumericsRadiation(unsigned short val_nDim, unsigned short val_nVar,
                     CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CNumericsRadiation(void);

  /*!
   * \brief Set the value of the radiation variable.
   * \param[in] val_radvar_i - Value of the turbulent variable at point i.
   * \param[in] val_radvar_j - Value of the turbulent variable at point j.
   */
  inline void SetRadVar(su2double *val_radvar_i, su2double *val_radvar_j){
    RadVar_i = val_radvar_i;
    RadVar_j = val_radvar_j;
  }


  /*!
   * \brief Set the gradient of the radiation variables.
   * \param[in] val_radvar_grad_i - Gradient of the turbulent variable at point i.
   * \param[in] val_radvar_grad_j - Gradient of the turbulent variable at point j.
   */
  inline void SetRadVarGradient(su2double **val_radvar_grad_i, su2double **val_radvar_grad_j){
    RadVar_Grad_i = val_radvar_grad_i;
    RadVar_Grad_j = val_radvar_grad_j;
  }

};
