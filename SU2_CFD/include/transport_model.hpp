/*!
 * \file transport_model.hpp
 * \brief Headers of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna
 * \version 7.0.4 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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
 * \version 6.2.0
 */
class CViscosityModel {
protected:
  su2double Mu,  /*!< \brief Dynamic viscosity. */
  dmudrho_T,     /*!< \brief DmuDrho_T. */
  dmudT_rho;     /*!< \brief DmuDT_rho. */
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
  ~CConstantViscosity(void) override;
   
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
  su2double Mu_ref,  /*!< \brief Internal Energy. */
  T_ref,             /*!< \brief DpDd_e. */
  S;                 /*!< \brief DpDe_d. */
  
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
  ~CSutherland(void) override;
  
  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho) override;
  
  /*!
   * \brief Set Viscosity Derivatives.
   */
  void SetDerViscosity(su2double T, su2double rho) override;
  
};

/*!
 * \class CPolynomialViscosity
 * \brief Defines viscosity as a polynomial function of temperature.
 * \author T. Economon
 */
class CPolynomialViscosity : public CViscosityModel {
protected:
  unsigned short nPolyCoeffs; /*!< \brief Number of coefficients in the temperature polynomial. */
  su2double *b;               /*!< \brief Polynomial coefficients for viscosity as a function of temperature. */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialViscosity(void);
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialViscosity(unsigned short val_nCoeffs, su2double* val_b);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPolynomialViscosity(void) override;
  
  /*!
   * \brief Set Viscosity.
   */
  void SetViscosity(su2double T, su2double rho) override;
  
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
  su2double Kt,  /*!< \brief Thermal conductivity. */
  dktdrho_T,     /*!< \brief DktDrho_T. */
  dktdT_rho;     /*!< \brief DktDT_rho. */
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
  virtual void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp);

  /*!
   * \brief Set Thermal conductivity derivatives.
   */
  virtual void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp);

};

/*!
 * \class CConstantConductivity
 * \brief this class defines a constant thermal conductivity.
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
  ~CConstantConductivity(void) override;

};

/*!
 * \class CConstantConductivityRANS
 * \brief this class defines a constant thermal conductivity with a turbulent Prandtl number for including effects of turbulent heat transfer.
 * \author T. Economon
 * \version 1.0
 */
class CConstantConductivityRANS : public CConductivityModel {
  
protected:
  su2double Kt_Lam;         /*!< \brief Constant laminar conductivity. */
  su2double Prandtl_Turb;   /*!< \brief Turbulent Prandtl number. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CConstantConductivityRANS(void);
  
  /*!
   * \brief Constructor of the class.
   */
  CConstantConductivityRANS(su2double kt_const, su2double pr_turb);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CConstantConductivityRANS(void) override;
  
  /*!
   * \brief Set effective thermal conductivity.
   */
  void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override;
  
};

/*!
 * \class CConstantPrandtl
 * \brief this class defines a non-constant thermal conductivity using a constant Prandtl's number
 * \author S.Vitale, M.Pini
 * \version 1.0
 */
class CConstantPrandtl : public CConductivityModel {
protected:
  su2double Pr_const;    /*!< \brief Prandtl's number. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CConstantPrandtl(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CConstantPrandtl(void) override;

  /*!
   * \brief Constructor of the class.
   */
  CConstantPrandtl(su2double pr_const);

  /*!
   * \brief Set Thermal conductivity.
   * \brief par1 -> Cp.
   * \brief par2 -> Mu.
   */
  void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override;

  /*!
   * \brief Set Thermal conductivity derivatives.
   */
  void SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) override;

};

/*!
 * \class CConstantPrandtlRANS
 * \brief Defines a non-constant effective thermal conductivity for RANS problems using Prandtl numbers.
 * \author T. Economon
 * \version 1.0
 */
class CConstantPrandtlRANS : public CConductivityModel {

protected:
  su2double Prandtl_Lam;    /*!< \brief Laminar Prandtl number. */
  su2double Prandtl_Turb;   /*!< \brief Turbulent Prandtl number. */

public:

  /*!
   * \brief Constructor of the class.
   */
  CConstantPrandtlRANS(void);

  /*!
   * \brief Destructor of the class.
   */
  ~CConstantPrandtlRANS(void) override;

  /*!
   * \brief Constructor of the class.
   */
  CConstantPrandtlRANS(su2double pr_lam, su2double pr_turb);

  /*!
   * \brief Set effective thermal conductivity.
   */
  void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override;

};

/*!
 * \class CPolynomialConductivity
 * \brief Defines a non-constant thermal conductivity using polynomial function of temperature.
 * \author T. Economon
 */
class CPolynomialConductivity : public CConductivityModel {
protected:
  unsigned short nPolyCoeffs; /*!< \brief Number of coefficients in the temperature polynomial. */
  su2double *b;               /*!< \brief Polynomial coefficients for thermal conductivity as a function of temperature. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialConductivity(void);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPolynomialConductivity(void) override;
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialConductivity(unsigned short val_nCoeffs, su2double* val_b);
  
  /*!
   * \brief Set Thermal conductivity.
   */
  void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override;
  
};

/*!
 * \class CPolynomialConductivityRANS
 * \brief Defines a non-constant thermal conductivity using polynomial function of temperature for RANS problems wioth the addition of a turbulent component based on a turbulent Prandtl number.
 * \author T. Economon
 */
class CPolynomialConductivityRANS : public CConductivityModel {
protected:
  unsigned short nPolyCoeffs; /*!< \brief Number of coefficients in the temperature polynomial. */
  su2double *b;               /*!< \brief Polynomial coefficients for thermal conductivity as a function of temperature. */
  su2double Prandtl_Turb;     /*!< \brief Turbulent Prandtl number. */

public:
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialConductivityRANS(void);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CPolynomialConductivityRANS(void) override;
  
  /*!
   * \brief Constructor of the class.
   */
  CPolynomialConductivityRANS(unsigned short val_nCoeffs, su2double* val_b, su2double pr_turb);
  
  /*!
   * \brief Set Thermal conductivity.
   */
  void SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) override;
  
};

#include "transport_model.inl"
