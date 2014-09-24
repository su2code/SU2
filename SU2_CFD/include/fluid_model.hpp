/*!
 * gas_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna
 * \version 1.0.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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
 * \version 3.2.1 "eagle"
 */
class CFluidModel {
protected:
double   	 StaticEnergy,			/*!< \brief Internal Energy. */
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
		double GetPressure ();

		/*!
		 * \brief Get fluid temperature.
		 */
		double GetTemperature ();

		/*!
		 * \brief Get fluid entropy.
		 */
		double GetEntropy ();

		/*!
		 * \brief Get fluid internal energy.
		 */
		double GetStaticEnergy ();

		/*!
		 * \brief Get fluid density.
		 */
		double GetDensity ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		double GetSoundSpeed ();

		/*!
		 * \brief Get fluid speed of sound squared.
		 */
		double GetSoundSpeed2 ();

		/*!
		 * \brief Get fluid specific heat at constant pressure.
		 */
		double GetCp ();

		/*!
		 * \brief Get fluid dynamic viscosity
		 */

		double GetLaminarViscosity ();

		/*!
		 * \brief Get fluid thermal conductivity
		 */

		double GetThermalConductivity ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		double GetdPdrho_e ();

		/*!
		 * \brief Get fluid pressure partial derivative.
		 */
		double GetdPde_rho ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		double GetdTdrho_e ();

		/*!
		 * \brief Get fluid temperature partial derivative.
		 */
		double GetdTde_rho ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		double Getdmudrho_T ();

		/*!
		 * \brief Get fluid dynamic viscosity partial derivative.
		 */
		double GetdmudT_rho ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		double Getdktdrho_T ();

		/*!
		 * \brief Get fluid thermal conductivity partial derivative.
		 */
		double GetdktdT_rho ();

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
		 * \param[in] InputSpec - Input pair for FLP calls ("e,rho").
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */

		virtual void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("PT").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (T).
		 */

		virtual void SetTDState_PT (double P, double T );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		virtual void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 *
		 */

		virtual void SetEnergy_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		virtual void SetTDState_hs (double h, double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		virtual void SetTDState_rhoT (double rho, double T );

};


/*!
 * \class CIdealGas
 * \brief Child class for defining ideal gas model.
 * \author: S.Vitale, M.Pini.
 * \version 3.2.1 "eagle"
 */
class CIdealGas : public CFluidModel {

protected:
	double Gamma, 						/*!< \brief Heat Capacity Ratio. */
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
		CIdealGas(double gamma, double R);


		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CIdealGas(void);

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("e,rho").
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */

		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("PT").
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */

		void SetTDState_PT (double P, double T );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Prho").
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Prho").
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetEnergy_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (double h, double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (double rho, double T );
};


/*!
 * \derived class CVanDerWaalsGas
 * \brief Child class for defining the Van der Waals model.
 * \author: S.Vitale, M.Pini
 * \version 3.2.1 "eagle"
 */
class CVanDerWaalsGas : public CIdealGas {

protected:
	double
			a, b, Zed;   					/*!< \brief Parameters for the Dimensionless Equation. */

public:

	   /*!
		 * \brief Constructor of the class.
		 */
		CVanDerWaalsGas(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CVanDerWaalsGas(double gamma, double R, double Pstar, double Tstar);


		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CVanDerWaalsGas(void);

		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */
		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */
		void SetTDState_PT (double P, double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief Set the Dimensionless Energy using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */

		void SetEnergy_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (double h, double s );


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (double rho, double T );

};


/*!
 * \derived class CPengRobinson
 * \brief Child class for defining the Peng-Robinson model.
 * \author: S.Vitale, G. Gori
 * \version 3.2.1 "eagle"
 */
class CPengRobinson : public CIdealGas {

protected:
	double  a, 						/*!< \brief model parameter. */
    		b, 						/*!< \brief model parameter. */
    		k, 						/*!< \brief model parameter (computed with acentric factor). */
    		Zed, 						/*!< \brief compressibility factor. */
    		TstarCrit;				/*!< \brief Critical temperature. */

private:

       /*!
	    * \brief Internal model parameter.
	    */
	    double  alpha2 (double T);

public:

	    /*!
		 * \brief Constructor of the class.
		 */
		CPengRobinson(void);

		/*!
		 * \brief Constructor of the class.
		 */
		CPengRobinson(double gamma, double R, double Pstar, double Tstar, double w);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CPengRobinson(void);

		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] rho - first thermodynamic variable.
		 * \param[in] e - second thermodynamic variable.
		 */
		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] T - second thermodynamic variable.
		 */
		void SetTDState_PT (double P, double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief Set the Dimensionless Energy using Pressure and Density
		 * \param[in] P - first thermodynamic variable.
		 * \param[in] rho - second thermodynamic variable.
		 */
		void SetEnergy_Prho (double P, double rho );
		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("hs").
		 * \param[in] th1 - first thermodynamic variable (h).
		 * \param[in] th2 - second thermodynamic variable (s).
		 *
		 */
		void SetTDState_hs (double h, double s );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("rhoT").
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (T).
		 *
		 */
		void SetTDState_rhoT (double rho, double T );

};



#include "fluid_model.inl"

