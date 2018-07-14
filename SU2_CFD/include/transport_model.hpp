/*!
 * \file transport_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 5.0.0 "Raven"
 *
 * SU2 Original Developers: Dr. Francisco D. Palacios.
 *                          Dr. Thomas D. Economon.
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

#ifndef TRANSPORT_MODEL_HPP_
#define TRANSPORT_MODEL_HPP_
#endif /* TRANSPORT_MODEL_HPP_ */
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
 * \class CViscosityModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Power law, Sutherland, Chung, etc.)
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CViscosityModel {
protected:
su2double      Mu,      /*!< \brief Dynamic viscosity. */
       dmudrho_T,   /*!< \brief DmuDrho_T. */
       dmudT_rho;   /*!< \brief DmuDT_rho. */
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
    su2double GetViscosity(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double Getdmudrho_T(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double GetdmudT_rho(void);

    /*!
     * \brief Set Viscosity.
     */
    virtual   void SetViscosity(su2double T, su2double rho);

    /*!
     * \brief Set Viscosity Derivatives.
     */
    virtual   void SetDerViscosity(su2double T, su2double rho);

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
  CConstantViscosity(su2double mu_const);
  
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
  su2double      Mu_ref,    /*!< \brief Internal Energy. */
  T_ref,     /*!< \brief DpDd_e. */
  S;       /*!< \brief DpDe_d. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CSutherland(void);
  
  /*!
   * \brief Constructor of the class.
   */
  CSutherland(su2double mu_ref, su2double t_ref, su2double s);
  
  /*!
   * \brief Destructor of the class.
   */
  virtual ~CSutherland(void);
  
  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho);
  
  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho);
  
};

class CLookUpTable_Viscosity: public CViscosityModel {

protected:
	int rank;
	bool skewed_linear_table;/*!< \brief Boolean to check for the type P-rho sample domain*/
	bool LUT_Debug_Mode;/*!< \brief If true, master node prints errors of points outside LUT*/
	su2double Density_Reference_Value;
	su2double Temperature_Reference_Value;
	su2double Viscosity_Reference_Value;
	su2double Pressure_Reference_Value;
	su2double Pressure, Density, Temperature;

	su2double
	**ThermoTables_Density, /*!< \brief Density look up table values. */
	**ThermoTables_Pressure, /*!< \brief Pressure look up table values. */
	**ThermoTables_Temperature, /*!< \brief Temperature look up table values. */
	**ThermoTables_Mu, /*!< \brief Laminar Viscosity look up table values. */
	**ThermoTables_dmudrho_T, /*!< \brief Fluid derivative DmuDrho_T look up table values. */
	**ThermoTables_dmudT_rho; /*!< \brief Fluid derivative DmuDT_rho look up table values. */

	su2double Interpolation_Matrix[4][4]; /*!< \brief The (Vandermonde) matrix for the interpolation (bilinear) */
	su2double Interpolation_Coeff[4][4]; /*!< \brief Used to hold inverse of Interpolation_Matrix, and solution vector */
	int LowerI, UpperI, middleI, LowerJ, UpperJ, middleJ;/*!< \brief The i,j indexes (rho, P) of the position of the table search. Can be used as a restart for next search.*/
	int Table_Pressure_Stations;/*!< \brief The pressure dimensions of the table */
	int Table_Density_Stations; /*!< \brief The density dimensions of the table */

	su2double Density_Table_Limits[2];/*!< \brief The [min,max] values of the Density values in the LUT */
	su2double Pressure_Table_Limits[2];/*!< \brief The [min,max] values of the Pressure values in the LUT */
	su2double SoundSpeed2_Table_Limits[2]; /*!< \brief The [min,max] values of the SoundSpeed squared values in the LUT */
	su2double Temperature_Table_Limits[2];/*!< \brief The [min,max] values of the Temperature values in the LUT */
	su2double Mu_Table_Limits[2];/*!< \brief The [min,max] values of the dPde_rho  values in the LUT */
	su2double dmudrho_T_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dmudrho_T  values in the LUT */
	su2double dmudT_rho_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dmudT_rho  values in the LUT */


public:

	/*!
	 * \brief default Constructor of the class.
	 */
	CLookUpTable_Viscosity(void);

	/*!
	 * \brief Constructor the LUT by reading it in from a file.
	 * \param[in] Filename - The name of the (.rgp) file from which to load the table
	 */
	CLookUpTable_Viscosity(CConfig *config, bool dimensional);

	/*!
	 * \brief Destructor of the class, primarily handling the dealloc of the KD_trees and LUT itself.
	 */
	virtual ~CLookUpTable_Viscosity(void);

	/*!
	 * \brief Set Viscosity using T and rho pair
	 */

	void SetViscosity(su2double T, su2double rho);

	/*!
	 * \brief Set Viscosity Derivatives using T, rho pair
	 */
	void SetDerViscosity(su2double T, su2double rho);

	void Search_NonEquispaced_Rho_Index(su2double rho);
	void Search_j_for_Y_given_i(su2double x, su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y );

	void Gaussian_Inverse(int nDim);

	/*!
	 * \brief Calculate the bilinear interpolation coefficients for a quad with arbitrary skew.
	 *  The entails building the Vandermonde matrix, inverting it, transposing it, and dot product with the search values of x, and y.
	 *  This formulation with the transpose means that the coefficients depend only on the x,y cooridinate of the search and
	 *  not on the thermodynamic variable being interpolated. Thus, the same coefficients can be used across
	 *  the interpolation of all desired thermodynamic properties.
	 * \param[in] x - the x value used to set the thermodynamic state. (e.g. rho in rhoe)
	 * \param[in] x - the y value used to set the thermodynamic state. (e.g. rho in e)
	 * \param[in] grid_var - the pair of thermodynamic variables which define the grid i.e. the interpolation quad. (e.g. RHOE for rhoe)
	 */

	void Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x, su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y, std::string grid_var);


	/*!
	 * \brief Use the interpolation coefficients to interpolate a given thermodynamic variable property. (Must calculate the interpolation coefficients first)
	 * \param[in] interpolant_var - the name of the variable to be interpolated e.g Density
	 */

	su2double Interpolate_2D_Bilinear(su2double ** ThermoTables_Z);
	void Check_Interpolated_PRHO_Limits(std::string interpolation_case);

	/*!
	 * \brief Load the LUT table from a CFX file format. X axis must be Density, and Y axis pressure. Equal spacing not required.
	 * \param[in] filename - the name of the CFX file containing the table
	 */

	void LookUpTable_Malloc();
	void LookUpTable_Load_CFX(std::string filename);
	void CFX_Import_Table_By_Number(ifstream *tab, su2double **ThermoTables_X, bool skip_prho);
	void LookUpTable_Load_DAT(std::string filename);
	void Find_Table_Limits();
	void NonDimensionalise_Table_Values();
};


/*!
 * \class CThermalConductivityModel
 * \brief Main class for defining the Transport-Physical Model
 * a child class for each particular Model (Prandtl-based, etc.)
 * \author S. Vitale, M. Pini
 * \version 1.0
 */
class CConductivityModel {
protected:
su2double      Kt,      /*!< \brief Thermal conductivity. */
       dktdrho_T,   /*!< \brief DktDrho_T. */
       dktdT_rho;   /*!< \brief DktDT_rho. */
public:

    /*!
     * \brief Constructor of the class.
     */
    CConductivityModel(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConductivityModel(void);

    /*!
     * \brief return viscosity value.
     */
    su2double GetConductivity(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double Getdktdrho_T(void);

    /*!
     * \brief return viscosity partial derivative value.
     */
    su2double GetdktdT_rho(void);

    /*!
     * \brief Set Thermal conductivity.
     */
    virtual   void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    virtual   void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};


/*!
 * \class CConstantPrandtl
 * \brief this class defines a constant thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantConductivity : public CConductivityModel {

public:

    /*!
     * \brief Constructor of the class.
     */
      CConstantConductivity(void);

    /*!
     * \brief Constructor of the class.
     */
      CConstantConductivity(su2double kt_const);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConstantConductivity(void);

};


/*!
 * \class CConstantPrandtl
 * \brief this class defines a non-constant thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantPrandtl : public CConductivityModel {
protected:
  su2double      Pr_const;    /*!< \brief Prandtl's number. */

public:

    /*!
     * \brief Constructor of the class.
     */
      CConstantPrandtl(void);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CConstantPrandtl(void);

    /*!
     * \brief Constructor of the class.
     */
      CConstantPrandtl(su2double pr_const);

    /*!
     * \brief Set Thermal conductivity.
     * \brief par1 -> Cp.
     * \brief par2 -> Mu.
     */
    void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};

/*!
 * \class CLookUpTable_Conductivity
 * \brief Enables reading conductivity values from a look up table
 * Values are returned based on an input of Rho and T
 * \author M.Kosec, A. Rubino, S.Vitale
 */
class CLookUpTable_Conductivity: public CConductivityModel {

protected:
	int rank;
	bool skewed_linear_table;/*!< \brief Boolean to check for the type P-rho sample domain*/
	bool LUT_Debug_Mode;/*!< \brief If true, master node prints errors of points outside LUT*/
	su2double Density_Reference_Value;
	su2double Pressure_Reference_Value;
	su2double Temperature_Reference_Value;
	su2double	Conductivity_Reference_Value;
	su2double Density, Temperature, Pressure;

	su2double
	**ThermoTables_Density, /*!< \brief Density look up table values. */
	**ThermoTables_Pressure, /*!< \brief Pressure look up table values. */
	**ThermoTables_Temperature, /*!< \brief Temperature look up table values. */
	**ThermoTables_Kt, /*!< \brief Thermal Conductivity look up table values. */
	**ThermoTables_dktdrho_T, /*!< \brief Fluid derivative DktDrho_T look up table values. */
	**ThermoTables_dktdT_rho; /*!< \brief Fluid derivative DktDT_rho look up table values. */

	su2double Interpolation_Matrix[4][4]; /*!< \brief The (Vandermonde) matrix for the interpolation (bilinear) */
	su2double Interpolation_Coeff[4][4]; /*!< \brief Used to hold inverse of Interpolation_Matrix, and solution vector */
	int LowerI, UpperI, middleI, LowerJ, UpperJ, middleJ;/*!< \brief The i,j indexes (rho, P) of the position of the table search. Can be used as a restart for next search.*/
	int Table_Pressure_Stations;/*!< \brief The pressure dimensions of the table */
	int Table_Density_Stations; /*!< \brief The density dimensions of the table */

	su2double Density_Table_Limits[2];/*!< \brief The [min,max] values of the Density values in the LUT */
	su2double Pressure_Table_Limits[2];/*!< \brief The [min,max] values of the Pressure values in the LUT */
	su2double Temperature_Table_Limits[2];/*!< \brief The [min,max] values of the Pressure values in the LUT */
	su2double Kt_Table_Limits[2];/*!< \brief The [min,max] values of the Kt values in the LUT */
	su2double dktdrho_T_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dktdrho_T values in the LUT */
	su2double dktdT_rho_Table_Limits[2];/*!< \brief (UNUSED) The [min,max] values of the dktdT_rho values in the LUT */


public:

	/*!
	 * \brief Constructor the LUT by reading it in from a file.
	 * \param[in] Filename - The name of the (.rgp) file from which to load the table
	 */
	CLookUpTable_Conductivity(CConfig *config);

	/*!
	 * \brief Destructor of the class, primarily handling the dealloc of the KD_trees and LUT itself.
	 */
	virtual ~CLookUpTable_Conductivity(void);

	/*!
	 * \brief Set Thermal conductivity.
	 * \brief par1 -> Cp.
	 * \brief par2 -> Mu.
	 */
	//WARNING only T and RHO are relevant inputs: mu and cp are UNUSED
	void SetConductivity(su2double T, su2double rho, su2double mu,
			su2double cp);

	//WARNING only T and RHO are relevant inputs: mu and cp are UNUSED
	void SetDerConductivity(su2double T, su2double rho, su2double mu,
			su2double cp);
	/*!
	 * \brief Calculate the inverse of a square matrix (e.g. the Vandermonde matrix) with pivoting Gaussian elimination
	 * \param[in] nDim - the dimension of the square block to invert
	 */

	void Search_NonEquispaced_Rho_Index(su2double rho);
	void Search_j_for_Y_given_i(su2double x, su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y );

	void Gaussian_Inverse(int nDim);

	/*!
	 * \brief Calculate the bilinear interpolation coefficients for a quad with arbitrary skew.
	 *  The entails building the Vandermonde matrix, inverting it, transposing it, and dot product with the search values of x, and y.
	 *  This formulation with the transpose means that the coefficients depend only on the x,y cooridinate of the search and
	 *  not on the thermodynamic variable being interpolated. Thus, the same coefficients can be used across
	 *  the interpolation of all desired thermodynamic properties.
	 * \param[in] x - the x value used to set the thermodynamic state. (e.g. rho in rhoe)
	 * \param[in] x - the y value used to set the thermodynamic state. (e.g. rho in e)
	 * \param[in] grid_var - the pair of thermodynamic variables which define the grid i.e. the interpolation quad. (e.g. RHOE for rhoe)
	 */

	void Interpolate_2D_Bilinear_Arbitrary_Skew_Coeff(su2double x, su2double y, su2double **ThermoTables_X, su2double **ThermoTables_Y, std::string grid_var);


	/*!
	 * \brief Use the interpolation coefficients to interpolate a given thermodynamic variable property. (Must calculate the interpolation coefficients first)
	 * \param[in] interpolant_var - the name of the variable to be interpolated e.g Density
	 */

	su2double Interpolate_2D_Bilinear(su2double ** ThermoTables_Z);
	void Check_Interpolated_PRHO_Limits(std::string interpolation_case);

	/*!
	 * \brief Load the LUT table from a CFX file format. X axis must be Density, and Y axis pressure. Equal spacing not required.
	 * \param[in] filename - the name of the CFX file containing the table
	 */

	void LookUpTable_Malloc();
	void LookUpTable_Load_CFX(std::string filename);
	void CFX_Import_Table_By_Number(ifstream *tab, su2double **ThermoTables_X, bool skip_prho);
	void LookUpTable_Load_DAT(std::string filename);
	void Find_Table_Limits();
	void NonDimensionalise_Table_Values();
};

#include "transport_model.inl"
