/*!
 * \file CNEMOGas.cpp
 * \brief Source of the nonequilibrium gas model.
 * \author C. Garbacz, W. Maier, S. R. Copeland
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

#include "../../include/fluid/CNEMOGas.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

CNEMOGas::CNEMOGas(const CConfig* config, unsigned short val_nDim): CFluidModel(){

  nSpecies = config->GetnSpecies();
  nDim     = val_nDim;

  MassFrac.resize(nSpecies,0.0);
  MolarMass.resize(nSpecies,0.0);
  MolarFractions.resize(nSpecies,0.0);
  rhos.resize(nSpecies,0.0);
  Cvtrs.resize(nSpecies,0.0);
  Cvves.resize(nSpecies,0.0);
  eves.resize(nSpecies,0.0);
  hs.resize(2*nSpecies,0.0);
  ws.resize(nSpecies,0.0);
  DiffusionCoeff.resize(nSpecies,0.0);
  Enthalpy_Formation.resize(nSpecies,0.0);
  Ref_Temperature.resize(nSpecies,0.0);
  temperatures.resize(nEnergyEq,0.0);
  energies.resize(nEnergyEq,0.0);
  ThermalConductivities.resize(nEnergyEq,0.0);

  gas_model            = config->GetGasModel();
  Kind_TransCoeffModel = config->GetKind_TransCoeffModel();

  frozen               = config->GetFrozen();
  ionization           = config->GetIonization();
}

void CNEMOGas::SetTDStatePTTv(su2double val_pressure, const su2double *val_massfrac,
                              su2double val_temperature, su2double val_temperature_ve){

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    MassFrac[iSpecies] = val_massfrac[iSpecies];
  Pressure = val_pressure;
  T        = val_temperature;
  Tve      = val_temperature_ve;

  su2double denom = 0.0;

  /*--- Calculate mixture density from supplied primitive quantities ---*/
  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    denom += MassFrac[iSpecies] * (Ru/MolarMass[iSpecies]) * Tve;
  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++)
    denom += MassFrac[iSpecies] * (Ru/MolarMass[iSpecies]) * T;
  Density = Pressure / denom;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++){
    rhos[iSpecies]     = MassFrac[iSpecies]*Density;
    MassFrac[iSpecies] = rhos[iSpecies]/Density;
  }
}

su2double CNEMOGas::ComputeSoundSpeed(){

  su2double conc    = 0.0;
  su2double rhoCvtr = 0.0;
  Density = 0.0;

  auto& Cvtrs = GetSpeciesCvTraRot();

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Density+=rhos[iSpecies];

  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++){
    conc += rhos[iSpecies]/MolarMass[iSpecies];
    rhoCvtr += rhos[iSpecies] * Cvtrs[iSpecies];
  }
  SoundSpeed2 = (1.0 + Ru/rhoCvtr*conc) * Pressure/Density;

 return(sqrt(SoundSpeed2));
}

su2double CNEMOGas::ComputePressure(){

  su2double P = 0.0;

  for (iSpecies = 0; iSpecies < nEl; iSpecies++)
    P += rhos[iSpecies] * Ru/MolarMass[iSpecies] * Tve;
  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++)
    P += rhos[iSpecies] * Ru/MolarMass[iSpecies] * T;

  Pressure = P;

  return P;

}

su2double CNEMOGas::ComputeGasConstant(){

  su2double Mass = 0.0;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mass += MassFrac[iSpecies] * MolarMass[iSpecies];
  GasConstant = Ru / Mass;

  return GasConstant;
}

su2double CNEMOGas::ComputeGamma(){

  /*--- Extract Values ---*/
  rhoCvtr = ComputerhoCvtr();

  /*--- Gamma Computation ---*/
  su2double rhoR = 0.0;
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhoR += rhos[iSpecies]*Ru/MolarMass[iSpecies];

  gamma = rhoR/rhoCvtr+1;

  return gamma;

}

su2double CNEMOGas::ComputerhoCvve() {

    Cvves = ComputeSpeciesCvVibEle(Tve);

    rhoCvve = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      rhoCvve += rhos[iSpecies]*Cvves[iSpecies];

    return rhoCvve;
}

void CNEMOGas::ComputedPdU(const su2double *V, const vector<su2double>& val_eves, su2double *val_dPdU){

  // Note: Electron energy not included properly.

  if (val_dPdU == nullptr) {
    SU2_MPI::Error("Array dPdU not allocated!", CURRENT_FUNCTION);
  }

  /*--- Determine the electron density (if ionized) ---*/
  su2double rho_el = 0.0;
  if (ionization) { rho_el = rhos[0]; }

  /*--- Necessary indexes to assess primitive variables ---*/
  unsigned long RHOS_INDEX    = 0;
  unsigned long T_INDEX       = nSpecies;
  unsigned long VEL_INDEX     = nSpecies+2;
  unsigned long RHOCVTR_INDEX = nSpecies+nDim+6;
  unsigned long RHOCVVE_INDEX = nSpecies+nDim+7;

  /*--- Extract variables ---*/
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies] = V[RHOS_INDEX+iSpecies];

  Cvtrs              = GetSpeciesCvTraRot();
  Enthalpy_Formation = GetSpeciesFormationEnthalpy();
  Ref_Temperature    = GetRefTemperature();

  /*--- Rename for convenience ---*/
  su2double rhoCvtr = V[RHOCVTR_INDEX];
  su2double rhoCvve = V[RHOCVVE_INDEX];
  T = V[T_INDEX];

  /*--- Pre-compute useful quantities ---*/
  su2double CvtrBAR = 0.0;
  su2double sqvel   = 0.0;
  su2double conc    = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    CvtrBAR += rhos[iSpecies]*Cvtrs[iSpecies];
    conc    += rhos[iSpecies]/MolarMass[iSpecies];
  }

  /*--- Species density derivatives ---*/
  su2double ef = 0.0;
  for (iSpecies = nEl; iSpecies < nHeavy; iSpecies++) {
    ef = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies];
    val_dPdU[iSpecies] = T*Ru/MolarMass[iSpecies] + Ru*conc/rhoCvtr *
                         (-Cvtrs[iSpecies]*(T-Ref_Temperature[iSpecies]) -
                         ef + 0.5*sqvel);
  }

  if (ionization) {
    for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
      //      evibs = Ru/MolarMass[iSpecies] * thetav[iSpecies]/(exp(thetav[iSpecies]/Tve)-1.0);
      //      num = 0.0;
      //      denom = g[iSpecies][0] * exp(-thetae[iSpecies][0]/Tve);
      //      for (iEl = 1; iEl < nElStates[iSpecies]; iEl++) {
      //        num   += g[iSpecies][iEl] * thetae[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //        denom += g[iSpecies][iEl] * exp(-thetae[iSpecies][iEl]/Tve);
      //      }
      //      eels = Ru/MolarMass[iSpecies] * (num/denom);

      val_dPdU[iSpecies] -= rho_el * Ru/MolarMass[0] * (val_eves[iSpecies])/rhoCvve;
    }
    ef = Enthalpy_Formation[0] - Ru/MolarMass[0]*Ref_Temperature[0];
    val_dPdU[0] = Ru*conc/rhoCvtr * (-ef + 0.5*sqvel) + Ru/MolarMass[0]*Tve
                - rho_el*Ru/MolarMass[0] * (-3.0/2.0*Ru/MolarMass[0]*Tve)/rhoCvve;
  }

  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dPdU[nSpecies+iDim] = -conc*Ru*V[VEL_INDEX+iDim]/rhoCvtr;

  /*--- Total energy derivative ---*/
  val_dPdU[nSpecies+nDim]   = conc*Ru / rhoCvtr;

  /*--- Vib.-el energy derivative ---*/
  val_dPdU[nSpecies+nDim+1] = -val_dPdU[nSpecies+nDim] +
                               rho_el*Ru/MolarMass[0]*1.0/rhoCvve;

}

void CNEMOGas::ComputedTdU(const su2double *V, su2double *val_dTdU){

  /*--- Necessary indexes to assess primitive variables ---*/
  const unsigned long T_INDEX       = nSpecies;
  const unsigned long VEL_INDEX     = nSpecies+2;
  const unsigned long RHOCVTR_INDEX = nSpecies+nDim+6;

  /*--- Rename for convenience ---*/
  T                 = V[T_INDEX];
  su2double rhoCvtr = V[RHOCVTR_INDEX];

  Cvtrs              = GetSpeciesCvTraRot();
  Enthalpy_Formation = GetSpeciesFormationEnthalpy();
  Ref_Temperature    = GetRefTemperature();

  /*--- Calculate supporting quantities ---*/
  const su2double v2 = GeometryToolbox::SquaredNorm(nDim, &V[VEL_INDEX]);

  /*--- Species density derivatives ---*/
  for (iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
    su2double ef    = Enthalpy_Formation[iSpecies] - Ru/MolarMass[iSpecies]*Ref_Temperature[iSpecies];
    val_dTdU[iSpecies]   = (-ef + 0.5*v2 + Cvtrs[iSpecies]*(Ref_Temperature[iSpecies]-T)) / rhoCvtr;
  }

  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTdU[nSpecies+iDim] = -V[VEL_INDEX+iDim] / V[RHOCVTR_INDEX];

  /*--- Energy derivatives ---*/
  val_dTdU[nSpecies+nDim]   =  1.0 / V[RHOCVTR_INDEX];
  val_dTdU[nSpecies+nDim+1] = -1.0 / V[RHOCVTR_INDEX];

}

void CNEMOGas::ComputedTvedU(const su2double *V, const vector<su2double>& val_eves, su2double *val_dTvedU){

  /*--- Necessary indexes to assess primitive variables ---*/
  unsigned long RHOCVVE_INDEX = nSpecies+nDim+7;

  /*--- Rename for convenience ---*/
  su2double rhoCvve = V[RHOCVVE_INDEX];

  /*--- Species density derivatives ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    val_dTvedU[iSpecies] = -val_eves[iSpecies]/rhoCvve;
  }
  /*--- Momentum derivatives ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    val_dTvedU[nSpecies+iDim] = 0.0;

  /*--- Energy derivatives ---*/
  val_dTvedU[nSpecies+nDim]   = 0.0;
  val_dTvedU[nSpecies+nDim+1] = 1.0 / rhoCvve;

}
