/*!
 * \file CNEMONSVariable.cpp
 * \brief Definition of the solution fields.
 * \author C. Garbacz, W. Maier, S.R. Copeland
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../include/variables/CNEMONSVariable.hpp"
#include <cmath>

CNEMONSVariable::CNEMONSVariable(su2double val_pressure,
                                 const su2double *val_massfrac,
                                 const su2double *val_mach,
                                 su2double val_temperature,
                                 su2double val_temperature_ve,
                                 unsigned long npoint,
                                 unsigned long val_ndim,
                                 unsigned long val_nvar,
                                 unsigned long val_nvarprim,
                                 unsigned long val_nvarprimgrad,
                                 const CConfig *config,
                                 CNEMOGas *fluidmodel) : CNEMOEulerVariable(val_pressure,
                                                                       val_massfrac,
                                                                       val_mach,
                                                                       val_temperature,
                                                                       val_temperature_ve,
                                                                       npoint,
                                                                       val_ndim,
                                                                       val_nvar,
                                                                       val_nvarprim,
                                                                       val_nvarprimgrad,
                                                                       config,
                                                                       fluidmodel) {

  Temperature_Ref = config->GetTemperature_Ref();
  Viscosity_Ref   = config->GetViscosity_Ref();
  Viscosity_Inf   = config->GetViscosity_FreeStreamND();
  Prandtl_Lam     = config->GetPrandtl_Lam();

  DiffusionCoeff.resize(nPoint, nSpecies)  = su2double(0.0);
  LaminarViscosity.resize(nPoint)          = su2double(0.0);
  ThermalCond.resize(nPoint)               = su2double(0.0);
  ThermalCond_ve.resize(nPoint)            = su2double(0.0);
  Enthalpys.resize(nPoint, nSpecies)       = su2double(0.0);

  Max_Lambda_Visc.resize(nPoint) = su2double(0.0);
  inv_TimeScale = config->GetModVel_FreeStream() / config->GetRefLength();

  Vorticity.resize(nPoint,3)     = su2double(0.0);
  StrainMag.resize(nPoint)       = su2double(0.0);
  Tau_Wall.resize(nPoint)        = su2double(-1.0);
  DES_LengthScale.resize(nPoint) = su2double(0.0);
  Roe_Dissipation.resize(nPoint) = su2double(0.0);
  Vortex_Tilting.resize(nPoint)  = su2double(0.0);
  Max_Lambda_Visc.resize(nPoint) = su2double(0.0);

}

bool CNEMONSVariable::SetPrimVar(unsigned long iPoint, CFluidModel *FluidModel) {

  fluidmodel = static_cast<CNEMOGas*>(FluidModel);

  /*--- Convert conserved to primitive variables ---*/
  bool nonPhys = Cons2PrimVar(Solution[iPoint], Primitive[iPoint], dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);

  /*--- Reset solution to previous one, if nonphys ---*/
  if (nonPhys) {
    for (auto iVar = 0u; iVar < nVar; iVar++)
      Solution(iPoint,iVar) = Solution_Old(iPoint,iVar);

    /*--- Recompute Primitive from previous solution ---*/
    Cons2PrimVar(Solution[iPoint], Primitive[iPoint], dPdU[iPoint], dTdU[iPoint], dTvedU[iPoint], eves[iPoint], Cvves[iPoint]);
  }

  /*--- Set additional point quantities ---*/
  Gamma(iPoint) = fluidmodel->ComputeGamma();

  SetVelocity2(iPoint);

  const auto& Ds = fluidmodel->GetDiffusionCoeff();
  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    DiffusionCoeff(iPoint, iSpecies) = Ds[iSpecies];

  su2double T   =  Primitive(iPoint,nSpecies);
  su2double Tve =  Primitive(iPoint,nSpecies+1);

  su2double* val_eves = GetEve(iPoint);
  const auto& hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, val_eves);
  for (auto iSpecies = 0u; iSpecies < nSpecies; iSpecies++)
    Enthalpys(iPoint, iSpecies) = hs[iSpecies];

  LaminarViscosity(iPoint) = fluidmodel->GetViscosity();

  const auto& thermalconductivities = fluidmodel->GetThermalConductivities();
  ThermalCond(iPoint)      = thermalconductivities[0];
  ThermalCond_ve(iPoint)   = thermalconductivities[1];

  Primitive(iPoint, LAM_VISC_INDEX) = LaminarViscosity(iPoint);

  return nonPhys;
}
