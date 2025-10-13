/*!
 * \file CEulerVariable.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, T. Economon
 * \version 8.3.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2025, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/variables/CEulerVariable.hpp"
#include "../../include/fluid/CFluidModel.hpp"
#include "../../../Common/include/toolboxes/random_toolbox.hpp"

unsigned long EulerNPrimVarGrad(const CConfig *config, unsigned long ndim) {
  if (config->GetContinuous_Adjoint()) return ndim + 4;
  if (config->GetKind_ConvNumScheme_Flow() == SPACE_CENTERED) return ndim + 1;

  const bool ideal_gas = config->GetKind_FluidModel() == STANDARD_AIR ||
                         config->GetKind_FluidModel() == IDEAL_GAS;
  if (ideal_gas && config->GetKind_Upwind_Flow() == UPWIND::ROE && !config->Low_Mach_Correction()) {
    // Based on CRoeBase (numerics_simd).
    return ndim + 2;
  }
  return ndim + 4;
}

CEulerVariable::CEulerVariable(su2double density, const su2double *velocity, su2double energy, unsigned long npoint,
                               unsigned long ndim, unsigned long nvar, const CConfig *config)
  : CFlowVariable(npoint, ndim, nvar, ndim + 9, EulerNPrimVarGrad(config, ndim), config),
    indices(ndim, 0) {

  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);
  const bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);

  nSecondaryVar = config->GetViscous() ? 8 : 2,
  nSecondaryVarGrad = 2;

  /*--- Solution initialization ---*/

  su2double val_solution[5] = {su2double(1.0), velocity[0], velocity[1], energy, energy};
  if(nDim==3) val_solution[3] = velocity[2];

  for(unsigned long iPoint = 0; iPoint < nPoint; ++iPoint)
    for (unsigned long iVar = 0; iVar < nVar; ++iVar)
      Solution(iPoint,iVar) = density*val_solution[iVar];

  Solution_Old = Solution;

  if (classical_rk4) Solution_New = Solution;

  /*--- Allocate and initializate solution for dual time strategy ---*/

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

  Secondary.resize(nPoint,nSecondaryVar) = su2double(0.0);

  if (config->GetAxisymmetric()){
    nAuxVar = 3;
    Grad_AuxVar.resize(nPoint,nAuxVar,nDim,0.0);
    AuxVar.resize(nPoint,nAuxVar) = su2double(0.0);
  }

  if (config->GetWind_Gust()) {
    WindGust.resize(nPoint,nDim);
  }

  if (config->GetVorticityConfinement()) {
    nAuxVar = 1;
    Grad_AuxVar.resize(nPoint, nAuxVar, nDim, 0.0);
    AuxVar.resize(nPoint, nAuxVar) = su2double(0.0);
  }

  if (config->GetKind_FluidModel() == ENUM_FLUIDMODEL::DATADRIVEN_FLUID){
    DataDrivenFluid = true;
    DatasetExtrapolation.resize(nPoint) = 0;
    NIterNewtonsolver.resize(nPoint) = 0;
    FluidEntropy.resize(nPoint) = su2double(0.0);
  }
}

bool CEulerVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  bool RightVol = true;

  SetVelocity(iPoint);   // Computes velocity and velocity^2
  su2double density      = GetDensity(iPoint);
  su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint);

  /*--- Check will be moved inside fluid model plus error description strings ---*/

  FluidModel->SetTDState_rhoe(density, staticEnergy);

  bool check_dens  = SetDensity(iPoint);
  bool check_press = SetPressure(iPoint, FluidModel->GetPressure());
  bool check_sos   = SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
  bool check_temp  = SetTemperature(iPoint, FluidModel->GetTemperature());

  /*--- Check that the solution has a physical meaning ---*/

  if (check_dens || check_press || check_sos || check_temp) {

    /*--- Copy the old solution ---*/

    for (unsigned long iVar = 0; iVar < nVar; iVar++)
      Solution(iPoint, iVar) = Solution_Old(iPoint, iVar);

    /*--- Recompute the primitive variables ---*/

    SetVelocity(iPoint);   // Computes velocity and velocity^2
    su2double density = GetDensity(iPoint);
    su2double staticEnergy = GetEnergy(iPoint)-0.5*Velocity2(iPoint);
    /* check will be moved inside fluid model plus error description strings*/
    FluidModel->SetTDState_rhoe(density, staticEnergy);

    SetDensity(iPoint);
    SetPressure(iPoint, FluidModel->GetPressure());
    SetSoundSpeed(iPoint, FluidModel->GetSoundSpeed2());
    SetTemperature(iPoint, FluidModel->GetTemperature());

    RightVol = false;

  }

  SetEnthalpy(iPoint); // Requires pressure computation.

  /*--- Set look-up variables in case of data-driven fluid model ---*/
  if (DataDrivenFluid) {
    SetDataExtrapolation(iPoint, FluidModel->GetExtrapolation());
    SetEntropy(iPoint, FluidModel->GetEntropy());
  }

  return RightVol;
}

void CEulerVariable::SetSecondaryVar(unsigned long iPoint, CFluidModel *FluidModel) {

   /*--- Compute secondary thermo-physical properties (partial derivatives...) ---*/

   SetdPdrho_e(iPoint, FluidModel->GetdPdrho_e());
   SetdPde_rho(iPoint, FluidModel->GetdPde_rho());

}

void CEulerVariable::SetVelDIHT(unsigned long npoint, unsigned long ndim, unsigned long nvar, const CConfig *config, 
                                const CGeometry *geometry) {

  const unsigned short MAXnModes = 5000;
  const unsigned short nModes = config->GetDIHT_nModes();
  if (nModes > MAXnModes) SU2_MPI::Error("The maximum number of Fourier modes for DIHT is 5000.", CURRENT_FUNCTION);
  if (nModes <= 0) SU2_MPI::Error("Assign a valid value for the number of Fourier modes. (DIHT)", CURRENT_FUNCTION);

  const su2double Lx = config->GetDIHT_DomainLength(0);
  const su2double Ly = config->GetDIHT_DomainLength(1);
  const su2double Lz = config->GetDIHT_DomainLength(2);
  if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0) SU2_MPI::Error("Assign a valid value for the computational domain size. (DIHT)", CURRENT_FUNCTION);

  const unsigned long nx = static_cast<unsigned long>(config->GetDIHT_nPoint(0));
  const unsigned long ny = static_cast<unsigned long>(config->GetDIHT_nPoint(1));
  const unsigned long nz = static_cast<unsigned long>(config->GetDIHT_nPoint(2));
  if (nx <= 0 || ny <= 0 || nz <= 0) SU2_MPI::Error("Assign a valid value for the number of nodes. (DIHT)", CURRENT_FUNCTION);

  const su2double pi = 4.0 * atan(1.0);
  const su2double vref = config->GetVelocity_Ref();
  su2double Density_Inf  = config->GetDensity_FreeStreamND();
  su2double ModVel_Freestream = config->GetModVel_FreeStream();

  su2double k_cbc[39] = {0.11, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.70, 1.00,            
                         1.50, 2.00, 2.50, 3.00, 4.00, 6.00, 8.00, 10.0, 12.5,
                         15.0, 17.5, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                         90.0, 100., 110., 120., 130., 140., 150., 160., 170.,
                         180., 190., 200.};

  su2double espec_cbc[39] = {30.0, 60.0, 129., 230., 322., 435., 457., 380.,
                             270., 168., 120., 89.0, 70.3, 47.0, 24.7, 12.6,
                             7.42, 3.96, 2.33, 1.34, 0.80,  0.0,  0.0,  0.0,
                              0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
                              0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};

  for (unsigned short ind = 0; ind < 39; ind++) {
    k_cbc[ind] *= 100.0;
    espec_cbc[ind] *= 1.0e-6;
  }

  su2double dx = Lx/nx;
  su2double dy = Ly/ny;
  su2double dz = Lz/nz;

  if (ndim < 3) {
    SU2_MPI::Error("DIHT: Decay of turbulence must be three-dimensional.", CURRENT_FUNCTION);
  }

  su2double wmin = min(2.0*pi/Lx, min(2.0*pi/Ly, 2.0*pi/Lz));
  su2double wmax = max(pi/dx, max(pi/dy, pi/dz));
  su2double dk = (wmax-wmin)/nModes;

  auto espec_interpol = [&] (su2double k_) {

    if (k_<k_cbc[0]) return espec_cbc[0];
    if (k_>k_cbc[38]) return espec_cbc[38];

    unsigned short ind = 1;
    while (ind < 39 && k_cbc[ind] < k_) ind++;
    
    su2double km = k_cbc[ind-1];
    su2double kp = k_cbc[ind];
    su2double em = espec_cbc[ind-1];
    su2double ep = espec_cbc[ind];

    su2double de_dk = (ep-em)/(kp-km);
    su2double e_interp = em + de_dk * (k_-km);

    return e_interp;

  };

  unsigned long seed = RandomToolbox::GetSeed(npoint, ndim + nvar);

  su2double phi[MAXnModes], theta[MAXnModes], psi[MAXnModes], phi1[MAXnModes], theta1[MAXnModes];

  for (unsigned long iMode = 0; iMode < nModes; iMode++) {
    std::mt19937 gen(seed + iMode*nModes);
    phi[iMode] = RandomToolbox::GetRandomUniform(gen, 0.0, 2.0*pi);
    theta[iMode] = acos(RandomToolbox::GetRandomUniform(gen, -1.0, 1.0));
    psi[iMode] = RandomToolbox::GetRandomUniform(gen, -0.5*pi, 0.5*pi);
    phi1[iMode] = RandomToolbox::GetRandomUniform(gen, 0.0, 2.0*pi);
    theta1[iMode] = acos(RandomToolbox::GetRandomUniform(gen, -1.0, 1.0));
  }

  for (unsigned long iPoint = 0; iPoint < npoint; iPoint++) {
    su2double u = 0.0, v = 0.0, w = 0.0;
    su2double xp = geometry->nodes->GetCoord(iPoint,0);
    su2double yp = geometry->nodes->GetCoord(iPoint,1);
    su2double zp = geometry->nodes->GetCoord(iPoint,2);
    for (unsigned long iMode = 0; iMode < nModes; iMode++) {
      su2double wn = wmin + 0.5*dk + iMode*dk;
      su2double kx = sin(theta[iMode]) * cos(phi[iMode]) * wn;
      su2double ky = sin(theta[iMode]) * sin(phi[iMode]) * wn;
      su2double kz = cos(theta[iMode]) * wn;
      su2double ktx = sin(kx * dx * 0.5) / dx;
      su2double kty = sin(ky * dy * 0.5) / dy;
      su2double ktz = sin(kz * dz * 0.5) / dz;
      su2double zetax = sin(theta1[iMode]) * cos(phi1[iMode]);
      su2double zetay = sin(theta1[iMode]) * sin(phi1[iMode]);
      su2double zetaz = cos(theta1[iMode]);
      su2double sxm = zetay*ktz - zetaz*kty;
      su2double sym = zetaz*ktx - zetax*ktz;
      su2double szm = zetax*kty - zetay*ktx;
      su2double smag = sqrt(sxm*sxm + sym*sym + szm*szm + 1.0e-16);
      sxm /= smag; sym /= smag; szm /= smag;
      su2double espec = espec_interpol(wn);
      su2double qm = sqrt(espec*dk);
      su2double arg = kx * xp + ky * yp + kz * zp - psi[iMode];
      su2double facx = 2.0 * qm * cos(arg - kx*dx*0.5);
      su2double facy = 2.0 * qm * cos(arg - ky*dy*0.5);
      su2double facz = 2.0 * qm * cos(arg - kz*dz*0.5);
      u += facx * sxm;
      v += facy * sym;
      w += facz * szm;
    }
    Solution(iPoint, 1) = Density_Inf * u/vref;
    Solution(iPoint, 2) = Density_Inf * v/vref;
    Solution(iPoint, 3) = Density_Inf * w/vref;
    su2double q2 = Solution(iPoint, 1)*Solution(iPoint, 1) +
                   Solution(iPoint, 2)*Solution(iPoint, 2) +
                   Solution(iPoint, 3)*Solution(iPoint, 3);
    Solution(iPoint, 4) += 0.5 * (q2/Density_Inf - Density_Inf*ModVel_Freestream*ModVel_Freestream);
  }
  
  const bool dual_time = (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_1ST) ||
                         (config->GetTime_Marching() == TIME_MARCHING::DT_STEPPING_2ND);

  const bool classical_rk4 = (config->GetKind_TimeIntScheme_Flow() == CLASSICAL_RK4_EXPLICIT);

  Solution_Old = Solution;

  if (classical_rk4) Solution_New = Solution;

  if (dual_time) {
    Solution_time_n = Solution;
    Solution_time_n1 = Solution;
  }

}