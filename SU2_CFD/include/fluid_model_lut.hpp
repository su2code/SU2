/*!
 * \file fluid_model.hpp
 * \brief Headers of the main thermodynamic subroutines of the SU2 solvers.
 * \author S. Vitale, G. Gori, M. Pini, A. Guardone, P. Colonna
 * \version 4.1.2 "Cardinal"
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

#include "../include/fluid_model.hpp"
#include "../../Common/include/config_structure.hpp"


struct KD_node
{
	int Branch_Splitting_Direction, Branch_Dimension, *Flattened_Point_Index;
	su2double * x_values, * y_values;
	KD_node* upper;
	KD_node* lower;
};

/*!
 * \class CFluidModel
 * \brief Main class for defining the Thermo-Physical Model to be stored in a table
 * \author: M. Kosec, A.Rubino, S.Vitale
 * \version 4.1.2 "Cardinal"
 */
class CThermoList {
public:
	su2double   	 StaticEnergy,			/*!< \brief Internal Energy. */
	Entropy,  				/*!< \brief Entropy. */
	Density,  				/*!< \brief Density. */
	Enthalpy,					/*!< \brief Enthalpy. */
	Pressure, 				/*!< \brief Pressure. */
	SoundSpeed2, 		/*!< \brief SpeedSound. */
	Temperature,			/*!< \brief Temperature. */
	dPdrho_e, 				/*!< \brief Fluid derivative DpDd_e. */
	dPde_rho, 				/*!< \brief Fluid derivative DpDe_d. */
	dTdrho_e, 				/*!< \brief Fluid derivative DTDd_e. */
	dTde_rho, 				/*!< \brief Fluid derivative DTDe_d. */
	Cp,              /*!< \brief Specific Heat Capacity at constant pressure. */
	Mu,					    /*!< \brief Laminar Viscosity. */
	dmudrho_T, 			/*!< \brief Fluid derivative DmuDrho_T */
	dmudT_rho,				/*!< \brief Fluid derivative DmuDT_rho. */
	Kt,					    /*!< \brief Thermal Conductivity. */
	dktdrho_T, 			/*!< \brief Fluid derivative DktDrho_T.  */
	dktdT_rho;				/*!< \brief Fluid derivative DktDT_rho. */




	/*!
	 * \brief Constructor of the class.
	 */
	CThermoList(void);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CThermoList(void);


//	void CTLprint(void);


};


/*!
 * \class CLookUpTable
 * \brief Child class for defining ideal gas model.
 * \author: A. Rubino, S.Vitale.
 * \version 4.1.2 "Cardinal"
 */
class CLookUpTable : public CFluidModel {

protected:
	CThermoList **ThermoTables;
	su2double Interpolation_Coeff[4][4]; /*!< \brief Fluid derivative DktDT_rho. */
	su2double Interpolation_Matrix[4][4]; /*!< \brief Fluid derivative DktDT_rho. */
	long iIndex, jIndex;
	int Table_Pressure_Stations, Table_Density_Stations; /*!< \brief The pressure and density dimensions of the table */
	KD_node *HS_tree; //KD tree for HS thermoPair
//	CThermoList interpolated;

	su2double StaticEnergy_Table_Limits[2];
	su2double Entropy_Table_Limits[2];
	su2double Enthalpy_Table_Limits[2];
	su2double Density_Table_Limits[2];
	su2double Pressure_Table_Limits[2];
	su2double SoundSpeed2_Table_Limits[2];
	su2double Temperature_Table_Limits[2];
	su2double dPdrho_e_Table_Limits[2];
	su2double dPde_rho_Table_Limits[2];
	su2double dTdrho_e_Table_Limits[2];
	su2double dTde_rho_Table_Limits[2];
	su2double Cp_Table_Limits[2];
	su2double Mu_Table_Limits[2];
	su2double dmudrho_T_Table_Limits[2];
	su2double dmudT_rho_Table_Limits[2];
	su2double Kt_Table_Limits[2];
	su2double dktdrho_T_Table_Limits[2];
	su2double dktdT_rho_Table_Limits[2];
	//Nearest neighbour's i and j indexes
	int* Nearest_Neighbour_iIndex;
	int* Nearest_Neighbour_jIndex;

public:

	/*!
	 * \brief default Constructor of the class.
	 */
	CLookUpTable(void);

	/*!
	 * \brief Constructor of the class.
	 */
	CLookUpTable(CConfig *config);

	/*!
	 * \brief Destructor of the class.
	 */
	virtual ~CLookUpTable(void);

	/*!
	 * \brief Search Thermo Pair
	 * \param[in] thermo1 - first thermodynamic variable.
	 * \param[in] thermo2 - second thermodynamic variable
	 * \param[in] input thermodynamic pair.
	 */
	void SearchThermoPair (su2double thermo1, su2double thermo2,  unsigned short thermoPair );


	struct KD_node* KD_Tree(su2double* x_values, su2double* y_values, int* Flattened_Point_Index, int Branch_Dimension, int Branch_Splitting_Direction);
	su2double Dist_KD_Tree (su2double x, su2double y, KD_node *branch);
	void free_KD_tree(KD_node* root);
	void NN_N_KD_Tree(int N, su2double thermo1, su2double thermo2,KD_node *root, su2double *best_dist);
	void SearchZigZag (su2double thermo1, su2double thermo2,  unsigned long thermoPair );
	void SearchThermoPair (su2double thermo1, su2double thermo2,  unsigned long thermoPair );

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


	void Interp2D_SingleSkewCoeff(std::string grid_var);
	su2double Interp2D_Inv_Dist(int N, std::string interpolant_var, su2double* dist);
	void Gaussian_Inverse(int nDim);
	void Interp2D_ArbitrarySkewCoeff(su2double x, su2double y, std::string grid_var);
	su2double Interp2D_lin(std::string interpolant_var);
	void TableLoadCFX(string filename);
	void LUTprint(void);
	void TableDump(char* filename);
	void RecordState(char* file);
};
