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


CNSVariable::CNSVariable(su2double density, const su2double *velocity, su2double energy,
                         unsigned long npoint, unsigned long ndim, unsigned long nvar, CConfig *config) :
                         CEulerVariable(density,velocity,energy,npoint,ndim,nvar,config) {

  inv_TimeScale = config->GetModVel_FreeStream() / config->GetRefLength();

  Vorticity.resize(nPoint,3) = su2double(0.0);
  StrainMag.resize(nPoint) = su2double(0.0);
  Tau_Wall.resize(nPoint) = su2double(-1.0);
  DES_LengthScale.resize(nPoint) = su2double(0.0);
  Roe_Dissipation.resize(nPoint) = su2double(0.0);
  Vortex_Tilting.resize(nPoint) = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint) = su2double(0.0);
}

bool CNSVariable::SetVorticity_StrainMag() {

  for (unsigned long iPoint = 0; iPoint < nPoint; ++iPoint) {

    /*--- Vorticity ---*/

    Vorticity(iPoint,0) = 0.0; Vorticity(iPoint,1) = 0.0;

    Vorticity(iPoint,2) = Gradient_Primitive(iPoint,2,0)-Gradient_Primitive(iPoint,1,1);

    if (nDim == 3) {
      Vorticity(iPoint,0) = Gradient_Primitive(iPoint,3,1)-Gradient_Primitive(iPoint,2,2);
      Vorticity(iPoint,1) = -(Gradient_Primitive(iPoint,3,0)-Gradient_Primitive(iPoint,1,2));
    }

    /*--- Strain Magnitude ---*/

    AD::StartPreacc();
    AD::SetPreaccIn(Gradient_Primitive[iPoint], nDim+1, nDim);

    su2double Div = 0.0;
    for (unsigned long iDim = 0; iDim < nDim; iDim++)
      Div += Gradient_Primitive(iPoint,iDim+1,iDim);

    StrainMag(iPoint) = 0.0;

    /*--- Add diagonal part ---*/

    for (unsigned long iDim = 0; iDim < nDim; iDim++) {
      StrainMag(iPoint) += pow(Gradient_Primitive(iPoint,iDim+1,iDim) - 1.0/3.0*Div, 2.0);
    }
    if (nDim == 2) {
      StrainMag(iPoint) += pow(1.0/3.0*Div, 2.0);
    }

    /*--- Add off diagonals ---*/

    StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,1,1) + Gradient_Primitive(iPoint,2,0)), 2);

    if (nDim == 3) {
      StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,1,2) + Gradient_Primitive(iPoint,3,0)), 2);
      StrainMag(iPoint) += 2.0*pow(0.5*(Gradient_Primitive(iPoint,2,2) + Gradient_Primitive(iPoint,3,1)), 2);
    }

    StrainMag(iPoint) = sqrt(2.0*StrainMag(iPoint));

    AD::SetPreaccOut(StrainMag(iPoint));
    AD::EndPreacc();
  }
  return false;
}

void CNSVariable::SetRoe_Dissipation_NTS(unsigned long iPoint,
                                         su2double val_delta,
                                         su2double val_const_DES){

  static const su2double cnu = pow(0.09, 1.5),
                         ch1 = 3.0,
                         ch2 = 1.0,
                         ch3 = 2.0,
                         sigma_max = 1.0;

  unsigned long iDim;
  su2double Omega, Omega_2 = 0.0, Baux, Gaux, Lturb, Kaux, Aaux;

  AD::StartPreacc();
  AD::SetPreaccIn(Vorticity[iPoint], 3);
  AD::SetPreaccIn(StrainMag(iPoint));
  AD::SetPreaccIn(val_delta);
  AD::SetPreaccIn(val_const_DES);
  /*--- Density ---*/
  AD::SetPreaccIn(Solution(iPoint,0));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive(iPoint,nDim+5));
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive(iPoint,nDim+6));

  /*--- Central/upwind blending based on:
   * Zhixiang Xiao, Jian Liu, Jingbo Huang, and Song Fu.  "Numerical
   * Dissipation Effects on Massive Separation Around Tandem Cylinders",
   * AIAA Journal, Vol. 50, No. 5 (2012), pp. 1119-1136.
   * https://doi.org/10.2514/1.J051299
   * ---*/

  for (iDim = 0; iDim < 3; iDim++){
    Omega_2 += pow(Vorticity(iPoint,iDim),2);
  }
  Omega = sqrt(Omega_2);

  Baux = (ch3 * Omega * max(StrainMag(iPoint), Omega)) /
      max((pow(StrainMag(iPoint),2)+Omega_2)*0.5, 1E-20);
  Gaux = tanh(pow(Baux,4.0));

  Kaux = max(sqrt((Omega_2 + pow(StrainMag(iPoint), 2))*0.5), 0.1 * inv_TimeScale);

  const su2double nu = GetLaminarViscosity(iPoint)/GetDensity(iPoint);
  const su2double nu_t = GetEddyViscosity(iPoint)/GetDensity(iPoint);
  Lturb = sqrt((nu + nu_t)/(cnu*Kaux));

  Aaux = ch2*max((val_const_DES*val_delta/Lturb)/Gaux -  0.5, 0.0);

  Roe_Dissipation(iPoint) = sigma_max * tanh(pow(Aaux, ch1));

  AD::SetPreaccOut(Roe_Dissipation(iPoint));
  AD::EndPreacc();

}

void CNSVariable::SetRoe_Dissipation_FD(unsigned long iPoint, su2double val_wall_dist){

  /*--- Constants for Roe Dissipation ---*/

  const passivedouble k2 = pow(0.41,2.0);

  AD::StartPreacc();
  AD::SetPreaccIn(Gradient_Primitive[iPoint], nVar, nDim);
  AD::SetPreaccIn(val_wall_dist);
  /*--- Eddy viscosity ---*/
  AD::SetPreaccIn(Primitive(iPoint,nDim+5));
  /*--- Laminar viscosity --- */
  AD::SetPreaccIn(Primitive(iPoint,nDim+6));

  su2double uijuij = 0.0;

  for(unsigned long iDim = 0; iDim < nDim; ++iDim)
    for(unsigned long jDim = 0; jDim < nDim; ++jDim)
      uijuij += pow(Gradient_Primitive(iPoint,1+iDim,jDim),2);

  uijuij = max(sqrt(uijuij),1e-10);

  const su2double nu = GetLaminarViscosity(iPoint)/GetDensity(iPoint);
  const su2double nu_t = GetEddyViscosity(iPoint)/GetDensity(iPoint);
  const su2double r_d = (nu + nu_t)/(uijuij*k2*pow(val_wall_dist,2));

  Roe_Dissipation(iPoint) = 1.0-tanh(pow(8.0*r_d,3.0));

  AD::SetPreaccOut(Roe_Dissipation(iPoint));
  AD::EndPreacc();
}

bool CNSVariable::SetPrimVar(unsigned long iPoint, su2double eddy_visc, su2double turb_ke, CFluidModel *FluidModel) {

  bool RightVol = true;

  SetVelocity(iPoint); // Computes velocity and velocity^2
  su2double density      = GetDensity(iPoint);
  su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint) - turb_ke;

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos  || check_temp) {

    /*--- Copy the old solution ---*/

    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    /*--- Recompute the primitive variables ---*/

    SetVelocity(iPoint); // Computes velocity and velocity^2
    density      = GetDensity(iPoint);
    staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint) - turb_ke;

    /*--- Check will be moved inside fluid model plus error description strings ---*/

    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity(iPoint);
    SetPressure(iPoint, FluidModel->GetPressure());
    SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint, FluidModel->GetTemperature());

    RightVol = false;

  }

  /*--- Set enthalpy ---*/

  SetEnthalpy(iPoint); // Requires pressure computation.

  /*--- Set laminar viscosity ---*/

  SetLaminarViscosity(iPoint, FluidModel->GetLaminarViscosity());

  /*--- Set eddy viscosity ---*/

  SetEddyViscosity(iPoint, eddy_visc);

  /*--- Set thermal conductivity ---*/

  SetThermalConductivity(iPoint, FluidModel->GetThermalConductivity());

  /*--- Set specific heat ---*/

  SetSpecificHeatCp(iPoint, FluidModel->GetCp());

  return RightVol;
}

void CNSVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {

    /*--- Compute secondary thermodynamic properties (partial derivatives...) ---*/

    SetdPdrho_e( iPoint, FluidModel->GetdPdrho_e() );
    SetdPde_rho( iPoint, FluidModel->GetdPde_rho() );

    SetdTdrho_e( iPoint, FluidModel->GetdTdrho_e() );
    SetdTde_rho( iPoint, FluidModel->GetdTde_rho() );

    /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

    Setdmudrho_T( iPoint, FluidModel->Getdmudrho_T() );
    SetdmudT_rho( iPoint, FluidModel->GetdmudT_rho() );

    Setdktdrho_T( iPoint, FluidModel->Getdktdrho_T() );
    SetdktdT_rho( iPoint, FluidModel->GetdktdT_rho() );

}

