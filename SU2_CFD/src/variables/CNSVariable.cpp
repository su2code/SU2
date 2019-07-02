/*!
 * \file CNSVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../include/variables/CNSVariable.hpp"


CNSVariable::CNSVariable(void) : CEulerVariable() { }

CNSVariable::CNSVariable(su2double val_density, su2double *val_velocity, su2double val_energy,
                         unsigned short val_nDim, unsigned short val_nvar, CConfig *config) :
                         CEulerVariable(val_density, val_velocity, val_energy, val_nDim, val_nvar, config) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();

  inv_TimeScale   = config->GetModVel_FreeStream() / config->GetRefLength();
  Roe_Dissipation = 0.0;
  Vortex_Tilting  = 0.0;
  Tau_Wall        = -1.0;

}

CNSVariable::CNSVariable(su2double *val_solution, unsigned short val_nDim,
                         unsigned short val_nvar, CConfig *config) :
                         CEulerVariable(val_solution, val_nDim, val_nvar, config) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();
  Prandtl_Turb    = config->GetPrandtl_Turb();

  inv_TimeScale   = config->GetModVel_FreeStream() / config->GetRefLength();
  Roe_Dissipation = 0.0;
  Vortex_Tilting  = 0.0;
  Tau_Wall        = -1.0;

}

CNSVariable::~CNSVariable(void) { }

bool CNSVariable::SetVorticity(void) {

  Vorticity[0] = 0.0; Vorticity[1] = 0.0;

  Vorticity[2] = Gradient_Primitive[2][0]-Gradient_Primitive[1][1];

  if (nDim == 3) {
    Vorticity[0] = Gradient_Primitive[3][1]-Gradient_Primitive[2][2];
    Vorticity[1] = -(Gradient_Primitive[3][0]-Gradient_Primitive[1][2]);
  }

  return false;

}

bool CNSVariable::SetStrainMag(void) {

  su2double Div;
  unsigned short iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive, nDim+1, nDim);

  Div = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Div += Gradient_Primitive[iDim+1][iDim];
  }

  StrainMag = 0.0;

  /*--- Add diagonal part ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    StrainMag += pow(Gradient_Primitive[iDim+1][iDim] - 1.0/3.0*Div, 2.0);
  }

  /*--- Add off diagonals ---*/

  StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][1] + Gradient_Primitive[2][0]), 2.0);

  if (nDim == 3) {
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[1][2] + Gradient_Primitive[3][0]), 2.0);
    StrainMag += 2.0*pow(0.5*(Gradient_Primitive[2][2] + Gradient_Primitive[3][1]), 2.0);
  }

  StrainMag = sqrt(2.0*StrainMag);

  AD::SetPreaccOut(StrainMag);
  AD::EndPreacc();

  return false;

}

void CNSVariable::SetRoe_Dissipation_NTS(su2double val_delta,
                                         su2double val_const_DES){

  static const su2double cnu = pow(0.09, 1.5),
                         ch1 = 3.0,
                         ch2 = 1.0,
                         ch3 = 2.0,
                         sigma_max = 1.0;

  unsigned short iDim;
  su2double Omega, Omega_2 = 0, Baux, Gaux, Lturb, Kaux, Aaux;

  AD::StartPreacc();
  AD::SetPreaccIn(Vorticity, 3);
  AD::SetPreaccIn(StrainMag);
  AD::SetPreaccIn(val_delta);
  AD::SetPreaccIn(val_const_DES);
  /*--- Density ---*/
  AD::SetPreaccIn(Solution[0]);
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive[nDim+5]);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive[nDim+6]);

  /*--- Central/upwind blending based on:
   * Zhixiang Xiao, Jian Liu, Jingbo Huang, and Song Fu.  "Numerical
   * Dissipation Effects on Massive Separation Around Tandem Cylinders",
   * AIAA Journal, Vol. 50, No. 5 (2012), pp. 1119-1136.
   * https://doi.org/10.2514/1.J051299
   * ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega_2 += Vorticity[iDim]*Vorticity[iDim];
  }
  Omega = sqrt(Omega_2);

  Baux = (ch3 * Omega * max(StrainMag, Omega)) /
      max((pow(StrainMag,2)+Omega_2)*0.5, 1E-20);
  Gaux = tanh(pow(Baux,4.0));

  Kaux = max(sqrt((Omega_2 + pow(StrainMag, 2))*0.5), 0.1 * inv_TimeScale);

  const su2double nu = GetLaminarViscosity()/GetDensity();
  const su2double nu_t = GetEddyViscosity()/GetDensity();
  Lturb = sqrt((nu + nu_t)/(cnu*Kaux));

  Aaux = ch2*max((val_const_DES*val_delta/Lturb)/Gaux -  0.5, 0.0);

  Roe_Dissipation = sigma_max * tanh(pow(Aaux, ch1));

  AD::SetPreaccOut(Roe_Dissipation);
  AD::EndPreacc();

}

void CNSVariable::SetRoe_Dissipation_FD(su2double val_wall_dist){

  /*--- Constants for Roe Dissipation ---*/

  static const su2double k2 = pow(0.41,2.0);

  su2double uijuij = 0;
  unsigned short iDim, jDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive, nVar, nDim);
  AD::SetPreaccIn(val_wall_dist);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive[nDim+5]);
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive[nDim+6]);

  for(iDim=0;iDim<nDim;++iDim){
    for(jDim=0;jDim<nDim;++jDim){
      uijuij+= Gradient_Primitive[1+iDim][jDim]*Gradient_Primitive[1+iDim][jDim];
    }
  }

  uijuij=sqrt(fabs(uijuij));
  uijuij=max(uijuij,1e-10);

  const su2double nu = GetLaminarViscosity()/GetDensity();
  const su2double nu_t = GetEddyViscosity()/GetDensity();
  const su2double r_d = (nu + nu_t)/(uijuij*k2*pow(val_wall_dist, 2.0));

  Roe_Dissipation = 1.0-tanh(pow(8.0*r_d,3.0));

  AD::SetPreaccOut(Roe_Dissipation);
  AD::EndPreacc();

}

bool CNSVariable::SetPrimVar(su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {

    unsigned short iVar;
  su2double density, staticEnergy;
  bool check_dens = false, check_press = false, check_sos = false,
  check_temp = false, RightVol = true;


  SetVelocity(); // Computes velocity and velocity^2
  density = GetDensity();
  staticEnergy = GetEnergy()-0.5*Velocity2 - turb_ke;

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  check_dens  = SetDensity();
  check_press = SetPressure(FluidModel->GetPressure());
  check_sos   = SetSoundSpeed(FluidModel->GetSoundSpeed2());
  check_temp  = SetTemperature(FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos  || check_temp) {

    /*--- Copy the old solution ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      Solution[iVar] = Solution_Old[iVar];

    /*--- Recompute the primitive variables ---*/

    SetVelocity(); // Computes velocity and velocity^2
    density = GetDensity();
    staticEnergy = GetEnergy()-0.5*Velocity2 - turb_ke;

    /*--- Check will be moved inside fluid model plus error description strings ---*/

    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity();
    SetPressure(FluidModel->GetPressure());
    SetSoundSpeed(FluidModel->GetSoundSpeed2());
    SetTemperature(FluidModel->GetTemperature());

    RightVol = false;

  }

  /*--- Set enthalpy ---*/

  SetEnthalpy();                                  // Requires pressure computation.

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity ---*/

  SetEddyViscosity(eddy_visc);

  /*--- Set thermal conductivity ---*/

  SetThermalConductivity(FluidModel->GetThermalConductivity());

  /*--- Set specific heat ---*/

  SetSpecificHeatCp(FluidModel->GetCp());

  return RightVol;

}

void CNSVariable::SetSecondaryVar(CFluidModel *FluidModel) {

    /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/

    SetdPdrho_e( FluidModel->GetdPdrho_e() );
    SetdPde_rho( FluidModel->GetdPde_rho() );

    SetdTdrho_e( FluidModel->GetdTdrho_e() );
    SetdTde_rho( FluidModel->GetdTde_rho() );

    /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

    Setdmudrho_T( FluidModel->Getdmudrho_T() );
    SetdmudT_rho( FluidModel->GetdmudT_rho() );

    Setdktdrho_T( FluidModel->Getdktdrho_T() );
    SetdktdT_rho( FluidModel->GetdktdT_rho() );

}

