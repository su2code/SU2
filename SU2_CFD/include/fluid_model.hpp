/*!
 * \file fluid_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#define LEN_COMPONENTS 32

#include "stdio.h"
#include "math.h"

using namespace std;

#include "../include/transport_model.hpp"
#include "../../Common/include/config_structure.hpp"


/*!
 * \class CFluidModel
 * \brief Main class for defining the Thermo-Physical Model
 * a child class for each particular Model (Ideal-Gas, Van der Waals, etc.)
 * \author: S.Vitale, G.Gori, M.Pini
 * \version 4.2.0 "Cardinal"
 */
class CFluidModel {
protected:
su2double   	 StaticEnergy,			/*!< \brief Internal Energy. */
			 Entropy,  				/*!< \brief Entropy. */
			 Density,  				/*!< \brief Density. */
			 Pressure, 				/*!< \brief Pressure. */
			 SoundSpeed2, 			/*!< \brief SpeedSound. */
			 Temperature,			/*!< \brief Temperature. */
			 dPdrho_e, 				/*!< \brief DpDd_e. */
			 dPde_rho, 				/*!< \brief DpDe_d. */
			 dTdrho_e, 				/*!< \brief DTDd_e. */
			 dTde_rho, 				/*!< \brief DTDe_d. */
             Cp,                    /*!< \brief Specific Heat Capacity at constant pressure. */
			 Mu,					/*!< \brief Specific Heat Capacity at constant pressure. */
		     dmudrho_T, 			/*!< \brief Specific Heat Capacity at constant pressure. */
		     dmudT_rho,				/*!< \brief Specific Heat Capacity at constant pressure. */
		     Kt,					/*!< \brief Specific Heat Capacity at constant pressure. */
		     dktdrho_T, 			/*!< \brief Specific Heat Capacity at constant pressure. */
		     dktdT_rho;				/*!< \brief Specific Heat Capacity at constant pressure. */

CViscosityModel *LaminarViscosity;	          /*!< \brief Laminar Viscosity Model */
CConductivityModel *ThermalConductivity;	  /*!< \brief Thermal Conductivity Model */

public:

	/*!
		 * \brief Constructor of the class.
		 */
		CFluidModel(void);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CFluidModel(void);

		/*!
		 * \brief Get fluid pressure.
		 */
		su2double GetPressure ();

		/*!
		 * \brief Get fluid temperature.
		 */
		su2double GetTemperature ();

		/*!
		 * \brief Get fluid entropy.
		 */
		su2double GetEntropy ();

		/*!
		 * \brief Get fluid internal energy.
		 */
		su2double GetStaticEnergy ();

		/*!
		 * \brief Get fluid density.
		 */
		su2double GetDensity ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		su2double GetSoundSpeed ();

		/*!
		 * \brief Get fluid speed of sound squared.
		 */
		su2double GetSoundSpeed2 ();

		/*!
		 * \brief Get fluid specific heat at constant pressure.
		 */
		su2double GetCp ();

		/*!
		 * \brief Get fluid dynamic viscosity
		 */

		su2double GetLaminarViscosity ();

		/*!
		 * \brief Get fluid thermal conductivity
		 */

		su2double GetThermalConductivity ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double GetdPdrho_e ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		su2double GetdPde_rho ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double GetdTdrho_e ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		su2double GetdTde_rho ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		su2double Getdmudrho_T ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		su2double GetdmudT_rho ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		su2double Getdktdrho_T ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		su2double GetdktdT_rho ();

		/*!
		 * \brief Set viscosity model.
		 */
		void SetLaminarViscosityModel (CConfig *config);

		/*!
		 * \brief Set thermal conductivity model.
		 */
		void SetThermalConductivityModel (CConfig *config);

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("e, rho").
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */

		virtual void SetTDState_rhoe (su2double rho, su2double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("PT").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (T).
		 */

		virtual void SetTDState_PT (su2double P, su2double T );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		virtual void SetTDState_Prho (su2double P, su2double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 *
		 */

		virtual void SetEnergy_Prho (su2double P, su2double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		virtual void SetTDState_hs (su2double h, su2double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		virtual void SetTDState_rhoT (su2double rho, su2double T );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */

		virtual void SetTDState_Ps (su2double P, su2double s );

};


/*!
 * \class CIdealGas
 * \brief Child class for defining ideal gas model.
 * \author: S.Vitale, M.Pini.
 * \version 4.2.0 "Cardinal"
 */
class CIdealGas : public CFluidModel {

protected:
	su2double Gamma, 						/*!< \brief Heat Capacity Ratio. */
	        Gamma_Minus_One, 			/*!< \brief Heat Capacity Ratio Minus One. */
	        Gas_Constant;				/*!< \brief Gas Constant. */


public:

	   /*!
		 * \brief Constructor of the class.
		 */
		CIdealGas(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CIdealGas(su2double gamma, su2double R);


		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CIdealGas(void);

		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */

		void SetTDState_rhoe (su2double rho, su2double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure  and Temperature
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */

		void SetTDState_PT (su2double P, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetTDState_Prho (su2double P, su2double rho );

		/*!
		 * \brief Set the Dimensionless Internal Energy using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetEnergy_Prho (su2double P, su2double rho );

		/*!
		 * \brief Set the Dimensionless State using Enthalpy and Entropy
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (su2double h, su2double s );


		/*!
		 * \brief Set the Dimensionless State using Density and Temperature
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (su2double rho, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Entropy
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */

		void SetTDState_Ps (su2double P, su2double s );
};


/*!
 * derived class CVanDerWaalsGas
 * \brief Child class for defining the Van der Waals model.
 * \author: S.Vitale, M.Pini
 * \version 4.2.0 "Cardinal"
 */
class CVanDerWaalsGas : public CIdealGas {

protected:
	su2double
			a, b, Zed;   					/*!< \brief Parameters for the Dimensionless Equation. */

public:

	   /*!
		 * \brief Constructor of the class.
		 */
		CVanDerWaalsGas(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CVanDerWaalsGas(su2double gamma, su2double R, su2double Pstar, su2double Tstar);


		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CVanDerWaalsGas(void);

		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */
		void SetTDState_rhoe (su2double rho, su2double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */
		void SetTDState_PT (su2double P, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetTDState_Prho (su2double P, su2double rho );

		/*!
		 * \brief Set the Dimensionless Internal Energy using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetEnergy_Prho (su2double P, su2double rho );

		/*!
		 * \brief Set the Dimensionless state using Enthalpy and Entropy
		 * \param[in] h - first thermodynamic variable (h).
		 * \param[in] s - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (su2double h, su2double s );


		/*!
		 * \brief Set the Dimensionless state using Density and Temperature
		 * \param[in] rho - first thermodynamic variable (rho).
		 * \param[in] T - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (su2double rho, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Entropy
		 * \param[in] P - first thermodynamic variable (P).
		 * \param[in] s - second thermodynamic variable (s).
		 */

		void SetTDState_Ps (su2double P, su2double s );

};


/*!
 * \derived class CPengRobinson
 * \brief Child class for defining the Peng-Robinson model.
 * \author: S.Vitale, G. Gori
 * \version 4.2.0 "Cardinal"
 */
class CPengRobinson : public CIdealGas {

protected:
	su2double  a, 						/*!< \brief model parameter. */
    		b, 						/*!< \brief model parameter. */
    		k, 						/*!< \brief model parameter (computed with acentric factor). */
    		Zed, 						/*!< \brief compressibility factor. */
    		TstarCrit;				/*!< \brief Critical temperature. */

private:

       /*!
	    * \brief Internal model parameter.
	    */
	    su2double  alpha2 (su2double T);


	   /*!
		* \brief Internal function for the implicit call hs.
		*/
		su2double  T_v_h (su2double v, su2double h);
		/*!
		* \brief Internal function for the implicit call Ps.
		*/
		su2double T_P_rho(su2double P, su2double rho);



public:

	    /*!
		 * \brief Constructor of the class.
		 */
		CPengRobinson(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CPengRobinson(su2double gamma, su2double R, su2double Pstar, su2double Tstar, su2double w);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CPengRobinson(void);

		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */
		void SetTDState_rhoe (su2double rho, su2double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */
		void SetTDState_PT (su2double P, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetTDState_Prho (su2double P, su2double rho );

		/*!
		 * \brief Set the Dimensionless Energy using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetEnergy_Prho (su2double P, su2double rho );
		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (su2double h, su2double s );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (su2double rho, su2double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Entropy
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */

		void SetTDState_Ps (su2double P, su2double s );

};



#include "fluid_model.inl"

