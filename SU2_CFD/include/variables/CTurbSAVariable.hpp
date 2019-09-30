/*!
 * \file CTurbSAVariable.hpp
 * \brief Declaration of the variables of the SA turbulence model.
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
 * \class CTurbSAVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author A. Bueno.
 */

class CTurbSAVariable : public CTurbVariable {

private:
  su2double gamma_BC; /*!< \brief Value of the intermittency for the BC trans. model. */
  su2double DES_LengthScale;
  su2double Vortex_Tilting;

public:
  /*!
   * \brief Constructor of the class.
   */
  CTurbSAVariable(void);

  /*!
   * \overload
   * \param[in] val_nu_tilde - Turbulent variable value (initialization value).
   * \param[in] val_muT  - The eddy viscosity
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CTurbSAVariable(su2double val_nu_tilde, su2double val_muT, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbSAVariable(void);

  /*!
   * \brief Set the harmonic balance source term.
   * \param[in] val_var - Index of the variable.
   * \param[in] val_source - Value of the harmonic balance source term. for the index <i>val_var</i>.
   */
  inline void SetHarmonicBalance_Source(unsigned short val_var, su2double val_source) {HB_Source[val_var] = val_source; }

  /*!
   * \brief Get the harmonic balance source term.
   * \param[in] val_var - Index of the variable.
   * \return Value of the harmonic balance source term for the index <i>val_var</i>.
   */
  inline su2double GetHarmonicBalance_Source(unsigned short val_var) {return HB_Source[val_var]; }

  /*!
   * \brief Get the intermittency of the BC transition model.
   * \return Value of the intermittency of the BC transition model.
   */
  inline su2double GetGammaBC(void) {return gamma_BC; }

  /*!
   * \brief Set the intermittency of the BC transition model.
   * \param[in] val_gamma - New value of the intermittency.
   */
  inline void SetGammaBC(su2double val_gamma) {gamma_BC = val_gamma; }

  /*!
   * \brief Get the DES length scale
   * \return Value of the DES length Scale.
   */
  inline su2double GetDES_LengthScale(void) {return DES_LengthScale; }

  /*!
   * \brief Set the DES Length Scale.
   */
  inline void SetDES_LengthScale(su2double val_des_lengthscale) {DES_LengthScale = val_des_lengthscale; }

  /*!
   * \brief Set the vortex tilting measure for computation of the EDDES length scale
   */
  void SetVortex_Tilting(su2double **PrimGrad_Flow, su2double* Vorticity, su2double LaminarViscosity);

  /*!
   * \brief Get the vortex tilting measure for computation of the EDDES length scale
   * \return Value of the DES length Scale
   */
  inline su2double GetVortex_Tilting() {return Vortex_Tilting; }

};
