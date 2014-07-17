/*!
 * transport_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S.Vitale, M.Pini, G.Gori, A.Guardone, P.Colonna
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


#ifndef TRANSPORT_MODEL_HPP_
#define TRANSPORT_MODEL_HPP_





#endif /* TRANSPORT_MODEL_HPP_ */
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




/*!
 * \class CGModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CViscosityModel {
protected:
double   	 Mu,			/*!< \brief Internal Energy. */
			 DmuDT_d, 				/*!< \brief DpDd_e. */
			 DmuDd_T; 				/*!< \brief DpDe_d. */
public:

		/*!
		 * \brief Constructor of the class.
		 */
		CViscosityModel(void);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CViscosityModel(void);

		/*!
		 * \brief return viscosity value.
		 */
		double GetViscosity(void);

		/*!
		 * \brief return viscosity value.
		 */
		double GetDmuDT_d(void);

		/*!
		 * \brief return viscosity value.
		 */
		double GetDmuDd_T(void);

		/*!
		 * \brief Set Viscosity.
		 */
		virtual	void SetViscosity(double T, double rho);


};


/*!
 * \class CConstantViscosity
 * \brief this class defines a constant viscosity
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantViscosity : public CViscosityModel {

public:

		/*!
		 * \brief Constructor of the class.
		 */
	    CConstantViscosity(void);

	    /*!
		 * \brief Constructor of the class.
		 */
		CConstantViscosity(double mu_const);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CConstantViscosity(void);


};


/*!
 * \class CSutherland
 * \brief this class defines a constant viscosity
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CSutherland : public CViscosityModel {
protected:
	double   	 Mu_ref,		/*!< \brief Internal Energy. */
				 T_ref, 		/*!< \brief DpDd_e. */
				 S; 			/*!< \brief DpDe_d. */

public:

		/*!
		 * \brief Constructor of the class.
		 */
	    CSutherland(void);

	    /*!
		 * \brief Constructor of the class.
		 */
		CSutherland(double mu_ref, double t_ref, double s);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CSutherland(void);

		/*!
		 * \brief Set Viscosity.
		 */
		void SetViscosity(double T, double rho);

};

#include "transport_model.inl"
