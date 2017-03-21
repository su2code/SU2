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
#include "../../Common/include/config_structure.hpp"

using namespace std;


/*!
 * \class CLiquidModel
 * \version 1.0
 */
class CLiquidModel {

protected:

	su2double rho_l, rho_m, h_l, dGibbs, Tsat, Psat, y, R, T_l, sigma, Rc;
	su2double Tstar, Gas_Constant;

CConfig *config;


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
    su2double GetLiquidDensity(void);

    /*!
     * \brief return mixture density value.
     */
    su2double GetMixtureDensity(void);

    /*!
     * \brief return liquid density value.
     */
    su2double GetLiquidTemperature(void);

    /*!
     * \brief return liquid enthalpy value.
     */
    su2double GetLiquidEnthalpy(void);

    /*!
     * \brief return surface tension value.
     */
    su2double GetSurfaceTension(void);

    /*!
     * \brief return saturation temperature value at correspondent pressure.
     */
    su2double GetTsat(void);

    /*!
     * \brief return saturation pressure value at correspondent temperature.
     */
    su2double GetPsat(void);

    /*!
     * \brief return critical radius.
     */
    su2double GetCriticalRadius(void);

    /*!
	* \brief return radius.
	*/
	su2double GetRadius(void);


    /*!
     * \brief return liquid density value.
     */
    void SetLiquidProp(su2double P, su2double T, su2double rho, su2double h_v, su2double *Two_Phase_Var);

    void SetRadius(su2double *Two_Phase_Var);

    void SetRCritical (su2double P, su2double T);

    void SetDensity_Mixture (su2double rho, su2double *Two_Phase_Var);

    virtual void SetTsat(su2double P);

    virtual void SetPsat (su2double T);

    virtual void SetLiquidDensity();

    virtual void SetTLiquid( su2double T);

    virtual void SetLiquidEnthalpy(su2double h_v);

    virtual void SetSurfaceTension(su2double T);

};


/*!
 * \class CWater
 * \version 1.0
 */
class CWater : public CLiquidModel {
protected:
  
  su2double *coeff_saturation, *coeff_latent_heat;
  su2double Asat, Bsat,
            Csat, Dsat,
		    Esat, Fsat,
		    Gsat, Hsat;

public:
  

  /*!
   * \brief Constructor of the class.
   */
  CWater();
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CWater(void);
  
  /*!
   * \brief return liquid density value.
   */

  void SetTsat(su2double P);

  void SetPsat (su2double T);

  void SetLiquidDensity();

  void SetTLiquid( su2double T);

  void SetLiquidEnthalpy(su2double h_v);

  void SetSurfaceTension(su2double T);

  
};




#include "liquid_phase_model.inl"
