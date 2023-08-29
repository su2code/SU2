/*!
 * \file NEMO_sources.cpp
 * \brief Implementation of numerics classes for integration
 *        of source terms in fluid flow NEMO problems.
 * \author C. Garbacz, W. Maier, S. Copeland.
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

#include "../../../include/numerics/NEMO/NEMO_sources.hpp"

CSource_NEMO::CSource_NEMO(unsigned short val_nDim,
                           unsigned short val_nVar,
                           unsigned short val_nPrimVar,
                           unsigned short val_nPrimVarGrad,
                           CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  /*--- Allocate arrays ---*/
  Y = new su2double[nSpecies];

  dYdr = new su2double*[nSpecies];
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    dYdr[iSpecies] = new su2double[nSpecies];
  }

  residual = new su2double[nVar]();
  jacobian = new su2double* [nVar];
  for(auto iVar = 0ul; iVar < nVar; ++iVar)
    jacobian[iVar] = new su2double [nVar]();
}

CSource_NEMO::~CSource_NEMO() {

  /*--- Deallocate arrays ---*/

  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    delete [] dYdr[iSpecies];
  delete [] dYdr;
  delete [] Y;

  delete [] residual;
  if (jacobian) {
    for(auto iVar = 0ul; iVar < nVar; ++iVar)
      delete [] jacobian[iVar];
    delete [] jacobian;
  }
}

CNumerics::ResidualType<> CSource_NEMO::ComputeChemistry(const CConfig *config) {

  /*--- Nonequilibrium chemistry ---*/
  vector<su2double> rhos;
  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    residual[iVar] = 0.0;

  if (implicit)
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;

  /*--- Rename for convenience ---*/
  su2double T = V_i[T_INDEX];
  su2double Tve = V_i[TVE_INDEX];
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  /*--- Set mixture state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  /*---Compute Prodcution/destruction terms ---*/
  const auto& ws = fluidmodel->ComputeNetProductionRates(implicit, V_i, eve_i, Cvve_i,
                                                         dTdU_i, dTvedU_i, jacobian);

  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++){
    residual[iSpecies] = ws[iSpecies] * Volume;}

  if (implicit) {
    for (auto iVar = 0ul; iVar<nVar; iVar++) {
      for (auto jVar = 0ul; jVar<nVar; jVar++) {
        jacobian[iVar][jVar] *= Volume;
      }
    }
  }

  return ResidualType<>(residual, jacobian, nullptr);

}

CNumerics::ResidualType<> CSource_NEMO::ComputeVibRelaxation(const CConfig *config) {

  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section
  const su2double res_min = -1E6;
  const su2double res_max = 1E6;

  vector<su2double> rhos;
  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }
  if (implicit) {
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;
  }

  /*--- Rename for convenience ---*/
  const su2double T = V_i[T_INDEX];
  const su2double Tve = V_i[TVE_INDEX];
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  /*--- Set fluid state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  /*--- Compute residual and jacobians ---*/
  const su2double VTterm = fluidmodel -> ComputeEveSourceTerm();
  if (implicit) {
    fluidmodel->GetEveSourceTermJacobian(V_i, eve_i, Cvve_i, dTdU_i,
                                         dTvedU_i, jacobian);
  }

  residual[nSpecies+nDim+1] = VTterm * Volume;

  if (implicit) {
    for (auto iVar = 0ul; iVar<nVar; iVar++) {
      for (auto jVar = 0ul; jVar<nVar; jVar++) {
        jacobian[iVar][jVar] *= Volume;
      }
    }
  }

  /*--- Relax/limit vt transfer ---*/
  if (config->GetVTTransferResidualLimiting()) {
    if (residual[nSpecies+nDim+1]>res_max) residual[nSpecies+nDim+1]=res_max;
    if (residual[nSpecies+nDim+1]<res_min) residual[nSpecies+nDim+1]=res_min;
  }

  return ResidualType<>(residual, jacobian, nullptr);
}

CNumerics::ResidualType<> CSource_NEMO::ComputeAxisymmetric(const CConfig *config) {

  /*--- Rename for convenience ---*/
  const auto Ds = Diffusion_Coeff_i;
  const auto GV = PrimVar_Grad_i;
  const su2double ktr = Thermal_Conductivity_i;
  const su2double kve = Thermal_Conductivity_ve_i;
  const su2double rho = V_i[RHO_INDEX];
  const su2double RuSI = UNIVERSAL_GAS_CONSTANT;
  const su2double Ru = 1000.0*RuSI;
  const su2double rhou = U_i[nSpecies];
  const su2double rhov = U_i[nSpecies+1];
  const su2double H = V_i[H_INDEX];
  const su2double rhoEve = U_i[nVar-1];
  const auto& Ms = fluidmodel->GetSpeciesMolarMass();
  const auto& hs = fluidmodel->ComputeSpeciesEnthalpy(V_i[T_INDEX], V_i[TVE_INDEX], eve_i );

  const bool viscous = config->GetViscous();
  const bool rans = (config->GetKind_Turb_Model() != TURB_MODEL::NONE);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) residual[iVar] = 0.0;

  /*--- Calculate inverse of y coordinate ---*/
  su2double yinv = 0.0;
  if (Coord_i[1]!= 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;

  /*--- Rename mass flux for convenience ---*/
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    Y[iSpecies] = V_i[RHOS_INDEX+iSpecies] / rho;

  /*--- Compute residual for inviscid axisym flow---*/
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
    residual[iSpecies] = yinv*rhov*Y[iSpecies]*Volume;
  residual[nSpecies]   = yinv*rhov*rhou/rho*Volume;
  residual[nSpecies+1] = yinv*rhov*rhov/rho*Volume;
  residual[nSpecies+2] = yinv*rhov*H*Volume;
  residual[nSpecies+3] = yinv*rhov*rhoEve/rho*Volume;

  /*---Compute Jacobian for inviscid axisym flow ---*/
  if (implicit) {

    /*--- Initialize matrices ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;

    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++)
        dYdr[iSpecies][jSpecies] = 0.0;

    /*--- Calculate additional quantities ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
        dYdr[iSpecies][jSpecies] += -1/rho*Y[iSpecies];
      }
      dYdr[iSpecies][iSpecies] += 1/rho;
    }

    /*--- Populate Jacobian ---*/

    // Species density
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
        jacobian[iSpecies][jSpecies] = dYdr[iSpecies][jSpecies]*rhov;
      }
      jacobian[iSpecies][nSpecies+1] = Y[iSpecies];
    }

    // X-momentum
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      jacobian[nSpecies][iSpecies] = -rhou*rhov/(rho*rho);
    jacobian[nSpecies][nSpecies]   = rhov/rho;
    jacobian[nSpecies][nSpecies+1] = rhou/rho;

    // Y-momentum
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
     jacobian[nSpecies+1][iSpecies] = -rhov*rhov/(rho*rho);
    jacobian[nSpecies+1][nSpecies+1] = 2*rhov/rho;

    // Energy
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      jacobian[nSpecies+nDim][iSpecies]      = -H*rhov/rho + dPdU_i[iSpecies]*rhov/rho;
    jacobian[nSpecies+nDim][nSpecies]        = dPdU_i[nSpecies]*rhov/rho;
    jacobian[nSpecies+nDim][nSpecies+1]      = H + dPdU_i[nSpecies+1]*rhov/rho;
    jacobian[nSpecies+nDim][nSpecies+nDim]   = (1+dPdU_i[nSpecies+nDim])*rhov/rho;
    jacobian[nSpecies+nDim][nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1]*rhov/rho;

    // Vib-el energy
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      jacobian[nSpecies+nDim+1][iSpecies] = -rhoEve*rhov/(rho*rho);
    jacobian[nSpecies+nDim+1][nSpecies+1] = rhoEve/rho;
    jacobian[nSpecies+nDim+1][nSpecies+nDim+1] = rhov/rho;

    for (auto iVar = 0ul; iVar < nVar; iVar++)
      for (auto jVar = 0ul; jVar < nVar; jVar++)
        jacobian[iVar][jVar] *= yinv*Volume;
  }

  /*--- Compute residual for viscous portion of axisym flow ---*/
  if (viscous) {
    if (!rans){ turb_ke_i = 0.0; }

    su2double Vector = 0.0;
    for (auto iSpecies = 0ul; iSpecies < nHeavy; iSpecies++)
      Vector += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1];

    su2double Mass = 0.0;
    for (auto iSpecies=0ul; iSpecies<nSpecies; iSpecies++)
      Mass += V_i[iSpecies]/rho*Ms[iSpecies];

    const su2double heat_capacity_cp_i   = V_i[RHOCVTR_INDEX]/rho + Ru/Mass;
    const su2double total_viscosity_i    = Laminar_Viscosity_i + Eddy_Viscosity_i;
    const su2double total_conductivity_i = ktr + kve + heat_capacity_cp_i*Eddy_Viscosity_i/Prandtl_Turb;
    const su2double u                    = V_i[VEL_INDEX];
    const su2double v                    = V_i[VEL_INDEX+1];
    const su2double qy_t                 = -total_conductivity_i*GV[T_INDEX][1];
    const su2double qy_ve                = -kve*GV[TVE_INDEX][1];

    /*--- Enthalpy and vib-el energy transport due to y-direction diffusion---*/
    su2double sumJhs_y, sumJeve_y;
    sumJhs_y = sumJeve_y = 0.0;
    for (auto iSpecies = 0ul; iSpecies < nHeavy; iSpecies++) {
      sumJhs_y  += -(rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1] - V_i[RHOS_INDEX+iSpecies]*Vector) * hs[iSpecies];
      sumJeve_y += -(rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][1] - V_i[RHOS_INDEX+iSpecies]*Vector) * eve_i[iSpecies];
    }

    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      residual[iSpecies] -= 0.0;
    residual[nSpecies] -= Volume*(yinv*total_viscosity_i*(GV[nSpecies+2][1]+GV[nSpecies+3][0])
                                                         -TWO3*AuxVar_Grad_i[0][0]);
    residual[nSpecies+1] -= Volume*(yinv*total_viscosity_i*2*(GV[nSpecies+3][1]-v*yinv)
                                                             -TWO3*AuxVar_Grad_i[0][1]);
    residual[nSpecies+2] -= Volume*(yinv*(-sumJhs_y + total_viscosity_i*(u*(GV[nSpecies+3][0]+GV[nSpecies+2][1])
                                                                        +v*TWO3*(2*GV[nSpecies+3][1]-GV[nSpecies+2][0]
                                                                        -v*yinv+rho*turb_ke_i))-qy_t)
                                                                        -TWO3*(AuxVar_Grad_i[1][1]+AuxVar_Grad_i[2][0]));
    residual[nSpecies+3] -= Volume*(yinv*(-sumJeve_y -qy_ve));
  }

  return ResidualType<>(residual, jacobian, nullptr);
}
