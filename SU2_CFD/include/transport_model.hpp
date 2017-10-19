/*!
 * \file transport_model.hpp
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


/*!
 * \class CFluidPropViscosity
 * \brief this class defines a viscosity according the Chung method implemented
 * in FluidProp (Power law, Sutherland, Chung, etc.)
 * \author T.P. van der Stelt
 * \version 1.0
 */
class CFluidPropViscosity : public CViscosityModel {
protected:

public:
		/*!
		 * \brief Constructor of the class.
		 */
	         CFluidPropViscosity(void);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CFluidPropViscosity(void);

		/*!
		 * \brief Set Viscosity.
		 */
		void SetViscosity(double T, double rho);

		/*!
		 * \brief Set Viscosity Derivatives.
		 */
		void SetDerViscosity(double T, double rho);
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
 * \class CFluidPropConductivity
 * \brief this class defines a thermal conductivity according the Chung method implemented in FluidProp
 * \author T.P. van der Stelt
 */
class CFluidPropConductivity : public CConductivityModel {
protected:
		double Pr_const;		/*!< \brief Prandtl's number. */

public:

		/*!
		 * \brief Constructor of the class.
		 */
	    	CFluidPropConductivity(void);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CFluidPropConductivity(void);

		/*!
		 * \brief Constructor of the class.
		 */
	    	CFluidPropConductivity(double pr_const);

		/*!
		 * \brief Set Thermal conductivity.
		 * \brief par1 -> Cp.
		 * \brief par2 -> Mu.
		 */
		void SetConductivity(double T, double rho, double mu, double cp);

		/*!
		 * \brief Set Thermal conductivity derivatives.
		 */
		void SetDerConductivity(double T, double rho, double dmudrho_T, double dmudT_rho, double cp);
};


class CCO2_Viscosity : public CViscosityModel {
protected:
  su2double      Mu_ref,    /*!< \brief Internal Energy. */
  T_ref, T_crit,    /*!< \brief DpDd_e. */
  S;       /*!< \brief DpDe_d. */
  su2double *visc_coeff;

public:

  /*!
   * \bief Constructor of the class.
   */
  CCO2_Viscosity(void);

  /*!
   * \brief Constructor of the class.
   */
  CCO2_Viscosity(su2double mu_ref, su2double t_ref, su2double t_crit);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CCO2_Viscosity(void);

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho);

  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho);

};

class CR22_Viscosity : public CViscosityModel {
protected:
  su2double      Mu_ref,    /*!< \brief Internal Energy. */
  T_ref,     /*!< \brief DpDd_e. */
  S;       /*!< \brief DpDe_d. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CR22_Viscosity(void);

  /*!
   * \brief Constructor of the class.
   */
  CR22_Viscosity(su2double mu_ref, su2double t_ref);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CR22_Viscosity(void);

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho);

  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho);

};

class CR12_Viscosity : public CViscosityModel {
protected:
  su2double      Mu_ref,    /*!< \brief Internal Energy. */
  T_ref,     /*!< \brief DpDd_e. */
  S;       /*!< \brief DpDe_d. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CR12_Viscosity(void);

  /*!
   * \brief Constructor of the class.
   */
  CR12_Viscosity(su2double mu_ref, su2double t_ref);

  /*!
   * \brief Destructor of the class.
   */
  virtual ~CR12_Viscosity(void);

  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho);

  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho);

};

class CCO2_Conductivity: public CConductivityModel {
protected:

su2double K_ref, T_ref, T_crit;
su2double *cond_coeff;
public:

    /*!
     * \brief Constructor of the class.
     */
    CCO2_Conductivity(void);

    CCO2_Conductivity(su2double k_ref, su2double t_ref, su2double t_crit);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CCO2_Conductivity(void);

    /*!
     * \brief Set Thermal conductivity.
     */
    void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};

class CR22_Conductivity: public CConductivityModel {
protected:

su2double K_ref, T_ref;
public:

    /*!
     * \brief Constructor of the class.
     */
    CR22_Conductivity(void);

    CR22_Conductivity(su2double k_ref, su2double t_ref);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CR22_Conductivity(void);

    /*!
     * \brief Set Thermal conductivity.
     */
    void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};

class CR12_Conductivity: public CConductivityModel {
protected:

su2double K_ref, T_ref;
public:

    /*!
     * \brief Constructor of the class.
     */
    CR12_Conductivity(void);

    CR12_Conductivity(su2double k_ref, su2double t_ref);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~CR12_Conductivity(void);

    /*!
     * \brief Set Thermal conductivity.
     */
    void SetConductivity(su2double T, su2double rho, su2double mu, su2double cp);

    /*!
     * \brief Set Thermal conductivity derivatives.
     */
    void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};


#include "transport_model.inl"
