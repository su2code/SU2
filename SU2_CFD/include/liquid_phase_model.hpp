/*!
 * \file liquid_phase_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#ifndef LIQUID_PHASE_MODEL_HPP_
#define LIQUID_PHASE_MODEL_HPP_
#endif /* LIQUID_PHASE_MODEL_HPP_ */
#pragma once

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

#include "../../Common/include/datatype_structure.hpp"

using namespace std;


/*!
 * \class CLiquidModel
 * \version 1.0
 */
class CLiquidModel {

protected:
	su2double rho_l, h_l;
	su2double dGibbs, Tsat, Psat, sigma;

public:


    /*!
     * \brief Constructor of the class.
     */
    CLiquidModel(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CLiquidModel(void);

    /*!
     * \brief return liquid density value.
     */
    su2double GetLiquidDensity_PT(void);

    /*!
     * \brief return liquid enthalpy value.
     */
    su2double GetLiquidEnthalpy_Prho(void);

    /*!
     * \brief return surface tension value.
     */
    su2double GetSurfaceTension_T(void);

    /*!
     * \brief return saturation temperature value at correspondent pressure.
     */
    su2double GetTsat_P(void);

    /*!
     * \brief return saturation pressure value at correspondent temperature.
     */
    su2double GetPsat_T(void);

    /*!
     * \brief return free gibbs energy variation.
     */
    su2double GetdGibbs_PT(void);


    /*!
     * \brief return liquid density value.
     */
    virtual   void SetLiquidDensity_PT(su2double P, su2double T, su2double Tstar);

    /*!
     * \brief return liquid enthalpy value.
     */
    virtual   void SetLiquidEnthalpy_Prho(su2double P, su2double T, su2double h_v, su2double Tstar);

    /*!
     * \brief return surface tension value.
     */
    virtual   void SetSurfaceTension_T(su2double T, su2double Tstar);

    /*!
     * \brief return saturation temperature value at correspondent pressure.
     */
    virtual   void SetTsat_P(su2double P);

    /*!
     * \brief return free gibbs energy variation and Psat at same temperature.
     */
    virtual   void SetdGibbs_PT(su2double P, su2double T, su2double Gas_Constant);


};


/*!
 * \class CWater
 * \version 1.0
 */
class CWater : public CLiquidModel {
protected:

  su2double rho_l, h_l, dGibbs, Tsat, Psat, sigma;
  
  su2double *coeff_saturation, *coeff_latent_heat;
  su2double Asat, Bsat,
            Csat, Dsat,
		    Esat, Fsat,
		    Gsat, Hsat;

public:
  

  /*!
   * \brief Constructor of the class.
   */
  CWater(void);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CWater(void);
  
  /*!
   * \brief return liquid density value.
   */
  void SetLiquidDensity_PT(su2double P, su2double T, su2double Tstar);

  /*!
   * \brief return liquid enthalpy value.
   */
  void SetLiquidEnthalpy_PT(su2double P, su2double T, su2double h_v, su2double Tstar);

  /*!
   * \brief return surface tension value.
   */
  void SetSurfaceTension_T(su2double T, su2double Tstar);

  /*!
   * \brief return saturation temperature value at correspondent pressure.
   */
  void SetTsat_P(su2double P);

  /*!
   * \brief return free gibbs energy variation and Psat at same temperature.
   */
  void SetdGibbs_PT(su2double P, su2double T, su2double Gas_Constant);


  
};




#include "liquid_phase_model.inl"
