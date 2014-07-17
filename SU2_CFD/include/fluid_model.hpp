/*!
 * gas_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author TUDelft Polimi
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
//#include <fluidprop.h>
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
 * \class CGModel
 * \brief Main class for defining the Thermo-Physical Model
 * a child class for each particular Model (Ideal-Gas, Van der Waals, etc.)
 * \author
 * \version 1.0
 */
class CFluidModel {
protected:
double   	 StaticEnergy,			/*!< \brief Internal Energy. */
			 Entropy,  				/*!< \brief Entropy. */
			 Density,  				/*!< \brief Density. */
			 Pressure, 				/*!< \brief Pressure. */
			 SoundSpeed2, 			/*!< \brief SpeedSound. */
			 Temperature,			/*!< \brief Temperature. */
			 ThermalConductivity, 	/*!< \brief ThermalConductivity. */
			 DpDd_e, 				/*!< \brief DpDd_e. */
			 DpDe_d, 				/*!< \brief DpDe_d. */
			 DTDd_e, 				/*!< \brief DTDd_e. */
			 DTDe_d; 				/*!< \brief DTDe_d. */

CViscosityModel *DynamicViscosity;	        /*!< \brief Viscosity */

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
		 * \brief Get fluid all sat props (p,T,d,h,s,c,..).
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
		 * \brief Get fluid entropy.
		 */
		double GetStaticEnergy ();

		/*!
		 * \brief Get fluid entropy.
		 */
		double GetDensity ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		double GetSoundSpeed ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		double GetSoundSpeed2 ();

		/*!
		 * \brief Get fluid speed of sound.
		 */

		double GetLaminarViscosity (double T, double rho);

		/*!
		 * \brief Get fluid speed of sound.
		 */
		double GetDpDd_e ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		double GetDpDe_d ();

		/*!
		 * \brief Get fluid speed of sound.
		 */
		void SetViscosityModel (CConfig *config);


		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		virtual void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
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
		 *brief virtual member that would be different for each gas model implemented
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */
//		virtual void SetTDState_Ps (double P, double s );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 *
		 */


		virtual void SetEnergy_Prho (double P, double rho );



};


/*!
 * \class CEulerSolver
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author F. Palacios.
 * \version 3.1.0 "eagle"
 */
class CIdealGas : public CFluidModel {

protected:
	double Gamma, 						/*!< \brief Heat Capacity Ratio. */
	        Gamma_Minus_One, 			/*!< \brief Heat Capacity Ratio Minus One. */
	        Gas_Constant;				/*!< \brief Enthalphy. */


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
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		void SetTDState_PT (double P, double T );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief virtual member that would be different for each gas model implemented
		 * \param[in] InputSpec - Input pair for FLP calls ("Pv").
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (v).
		 */

		void SetEnergy_Prho (double P, double rho );

};


/*!
 * \derived class CVanDerWaalsGas
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author M. Pini
 * \version 3.1.0 "eagle"
 */
class CVanDerWaalsGas : public CIdealGas {

protected:
	double
			a, b;   					/*!< \brief Parameters for the Dimensionless Equation. */

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
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (e).
		 */
		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (T).
		 */
		void SetTDState_PT (double P, double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (rho).
		 */
		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Entropy
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */
//		void SetTDState_Ps (double P, double s );

		/*!
		 * \brief Set the Dimensionless Energy using Pressure and Density
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (rho).
		 */
		void SetEnergy_Prho (double P, double rho );

		/*!
		 * \brief Set the Dimensionless Entropy using Pressure and Density
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (rho).
//		 */
//		void SetEntropy_Pd (double P, double rho );
//
//		/*!
//		 * \brief Set the Dimensional State using Pressure and Temperature
//		 * \param[in] th1 - first thermodynamic variable (P).
//		 * \param[in] th2 - second thermodynamic variable (T).
//		 */
//		void SetDimensionTDState_PT (double P, double T, double *params);
//
//		/*!
//		 * \brief Set the Dimensional State using Pressure and Density
//		 * \param[in] th1 - first thermodynamic variable (P).
//		 * \param[in] th2 - second thermodynamic variable (rho).
//		 */
//		void SetDimensionTDState_Pd (double P, double rho, double *params);

};


/*!
 * \derived class CPengRobinson
 * \brief Main class for defining the Euler's flow solver.
 * \ingroup Euler_Equations
 * \author G. Gori
 * \version 3.1.0 "eagle"
 */
class CPengRobinson : public CIdealGas {

protected:
	double  a, 						/*!< \brief Heat Capacity Ratio. */
    		b, 						/*!< \brief Heat Capacity Ratio Minus One. */
	        k, 						/*!< \brief Critical Temperature. */
			TstarCrit;					/*!< \brief Critical Temperature. */



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
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (e).
		 */
		double alpha2(double T);
		/*!
		 * \brief Set the Dimensionless State using Density and Internal Energy
		 * \param[in] th1 - first thermodynamic variable (rho).
		 * \param[in] th2 - second thermodynamic variable (e).
		 */
		void SetTDState_rhoe (double rho, double e );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Temperature
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (T).
		 */
		void SetTDState_PT (double P, double T );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Density
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (rho).
		 */
		void SetTDState_Prho (double P, double rho );

		/*!
		 * \brief Set the Dimensionless State using Pressure and Entropy
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (s).
		 */
//		void SetTDState_Ps (double P, double s );

		/*!
		 * \brief Set the Dimensionless Energy using Pressure and Density
		 * \param[in] th1 - first thermodynamic variable (P).
		 * \param[in] th2 - second thermodynamic variable (rho).
		 */
		void SetEnergy_Prho (double P, double rho );
};



#include "fluid_model.inl"

