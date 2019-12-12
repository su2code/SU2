/*!
 * transport_model.cpp
 * \brief Source of the main transport properties subroutines of the SU2 solvers.
 * \author S. Vitale, M. Pini, G. Gori, A. Guardone, P. Colonna, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation 
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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


#include "../include/transport_model.hpp"

CViscosityModel::CViscosityModel(void) {

  /*--- Attributes initialization ---*/

  Mu = 0.0;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CViscosityModel::~CViscosityModel(void) { }


CConstantViscosity::CConstantViscosity(void) : CViscosityModel() { }

CConstantViscosity::CConstantViscosity(su2double mu_const) : CViscosityModel() {

  /*--- Attributes initialization ---*/

  Mu = mu_const;
  dmudrho_T = 0.0;
  dmudT_rho = 0.0;

}

CConstantViscosity::~CConstantViscosity(void) { }




CSutherland::CSutherland(void) : CViscosityModel() {
  Mu_ref = 0.0;
  T_ref = 0.0;
  S = 0.0;

}

CSutherland::CSutherland(su2double mu_ref, su2double t_ref, su2double s) : CViscosityModel() {

  Mu_ref = mu_ref;
  T_ref = t_ref;
  S = s;
}

CSutherland::~CSutherland(void) { }


void CSutherland::SetViscosity(su2double T, su2double rho) {

  const su2double TnonDim = T/T_ref;
  Mu = Mu_ref*TnonDim*sqrt(TnonDim)*((T_ref + S)/(T + S));

}

void CSutherland::SetDerViscosity(su2double T, su2double rho) {

  dmudrho_T = 0.0;

  const su2double T_refInv = 1.0/T_ref;
  const su2double TnonDim  = T_refInv*T;
  const su2double TSInv    = 1.0/(T + S);
 
  dmudT_rho = Mu_ref*(T_ref + S)*TSInv*sqrt(TnonDim)
            * (1.5*T_refInv - TnonDim*TSInv);

}

CPolynomialViscosity::CPolynomialViscosity(void) : CViscosityModel() {
  nPolyCoeffs = 0;
  b           = NULL;
}

CPolynomialViscosity::CPolynomialViscosity(unsigned short val_nCoeffs, su2double* val_b) : CViscosityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
}

CPolynomialViscosity::~CPolynomialViscosity(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialViscosity::SetViscosity(su2double T, su2double rho) {
  
  /*--- Evaluate the new Mu from the coefficients and temperature. ---*/
  
  Mu = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Mu += b[iVar]*pow(T,iVar);
  
}

/*-------------------------------------------------*/
/*---------- Thermal Conductivity Models ----------*/
/*-------------------------------------------------*/

CConductivityModel::CConductivityModel(void) {

  /*--- Attributes initialization ---*/

  Kt = 0.0;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConductivityModel::~CConductivityModel(void) { }


CConstantConductivity::CConstantConductivity(void) : CConductivityModel() { }

CConstantConductivity::CConstantConductivity(su2double kt_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Kt = kt_const;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;

}

CConstantConductivity::~CConstantConductivity(void) { }

CConstantConductivityRANS::CConstantConductivityRANS(void) : CConductivityModel() { }

CConstantConductivityRANS::CConstantConductivityRANS(su2double kt_const, su2double pr_turb) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  Kt_Lam    = kt_const;
  dktdrho_T = 0.0;
  dktdT_rho = 0.0;
  
  Prandtl_Turb = pr_turb;

}

void CConstantConductivityRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  Kt = Kt_Lam + cp*mu_turb/Prandtl_Turb;
  
}

CConstantConductivityRANS::~CConstantConductivityRANS(void) { }

CConstantPrandtl::CConstantPrandtl(void) : CConductivityModel() { }

CConstantPrandtl::CConstantPrandtl(su2double pr_const) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Pr_const = pr_const;

}

void CConstantPrandtl::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {

  Kt = mu_lam*cp/Pr_const;

}

void CConstantPrandtl::SetDerConductivity(su2double T, su2double rho, su2double dmudrho_T, su2double dmudT_rho, su2double cp) {

  dktdrho_T = dmudrho_T*cp/Pr_const;
  dktdT_rho = dmudT_rho*cp/Pr_const;

}

CConstantPrandtl::~CConstantPrandtl(void) { }

CConstantPrandtlRANS::CConstantPrandtlRANS(void) : CConductivityModel() { }

CConstantPrandtlRANS::CConstantPrandtlRANS(su2double pr_lam, su2double pr_turb) : CConductivityModel() {

  /*--- Attributes initialization ---*/

  Prandtl_Lam  = pr_lam;
  Prandtl_Turb = pr_turb;
}

void CConstantPrandtlRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {

  Kt = cp * ((mu_lam/Prandtl_Lam) + (mu_turb/Prandtl_Turb));
  
}

CConstantPrandtlRANS::~CConstantPrandtlRANS(void) { }

CPolynomialConductivity::CPolynomialConductivity(void) : CConductivityModel() {
  nPolyCoeffs = 0;
  b           = NULL;
}

CPolynomialConductivity::CPolynomialConductivity(unsigned short val_nCoeffs, su2double* val_b) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
}

CPolynomialConductivity::~CPolynomialConductivity(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialConductivity::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  /*--- Evaluate the new Kt from the coefficients and temperature. ---*/
  
  Kt = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Kt += b[iVar]*pow(T,iVar);
  
}

CPolynomialConductivityRANS::CPolynomialConductivityRANS(unsigned short val_nCoeffs, su2double* val_b, su2double pr_turb) : CConductivityModel() {
  
  /*--- Attributes initialization ---*/
  
  nPolyCoeffs = val_nCoeffs;
  b = new su2double[nPolyCoeffs];
  
  for (unsigned short iVar = 0; iVar < nPolyCoeffs; iVar++)
    b[iVar] = val_b[iVar];
  
  Prandtl_Turb = pr_turb;
  
}

CPolynomialConductivityRANS::~CPolynomialConductivityRANS(void) {
  if (b != NULL) delete [] b;
}

void CPolynomialConductivityRANS::SetConductivity(su2double T, su2double rho, su2double mu_lam, su2double mu_turb, su2double cp) {
  
  /*--- Evaluate the new Kt from the coefficients and temperature. ---*/
  
  Kt = b[0];
  for (unsigned short iVar = 1; iVar < nPolyCoeffs; iVar++)
    Kt += b[iVar]*pow(T,iVar);
  
  /*--- Add a component due to turbulence to compute effective conductivity. ---*/
  
  Kt += cp*mu_turb/Prandtl_Turb;
  
}

