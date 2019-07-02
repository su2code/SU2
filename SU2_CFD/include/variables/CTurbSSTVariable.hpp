/*!
 * \file CTurbSSTVariable.hpp
 * \brief Declaration of the variables of the SST turbulence model.
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
 * \class CTurbSSTVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSSTVariable : public CTurbVariable {
protected:
  su2double sigma_om2,
  beta_star;
  su2double F1,    /*!< \brief Menter blending function for blending of k-w and k-eps. */
  F2,            /*!< \brief Menter blending function for stress limiter. */
  CDkw;           /*!< \brief Cross-diffusion. */

public:
  /*!
   * \brief Constructor of the class.
   */
  CTurbSSTVariable(void);

  /*!
   * \overload
   * \param[in] val_rho_kine - Turbulent variable value (initialization value).
   * \param[in] val_rho_omega - Turbulent variable value (initialization value).
   * \param[in] val_muT - Turbulent variable value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] constants -
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSSTVariable(su2double val_rho_kine, su2double val_rho_omega, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar,
                   su2double *constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSSTVariable(void);

  /*!
   * \brief Set the blending function for the blending of k-w and k-eps.
   * \param[in] val_viscosity - Value of the vicosity.
   * \param[in] val_dist - Value of the distance to the wall.
   * \param[in] val_density - Value of the density.
   */
  void SetBlendingFunc(su2double val_viscosity, su2double val_dist, su2double val_density);

  /*!
   * \brief Get the first blending function.
   */
  inline su2double GetF1blending(void) { return F1; }

  /*!
   * \brief Get the second blending function.
   */
  inline su2double GetF2blending(void) { return F2; }

  /*!
   * \brief Get the value of the cross diffusion of tke and omega.
   */
  inline su2double GetCrossDiff(void) { return CDkw; }
};
