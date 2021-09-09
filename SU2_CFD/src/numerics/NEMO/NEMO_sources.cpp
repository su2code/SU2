/*!
 * \file NEMO_sources.cpp
 * \brief Implementation of numerics classes for integration
 *        of source terms in fluid flow NEMO problems.
 * \author C. Garbacz, W. Maier, S. Copeland.
 * \version 7.2.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2021, SU2 Contributors (cf. AUTHORS.md)
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

  unsigned short iSpecies;

  ws.resize(nSpecies,0.0);

  /*--- Allocate arrays ---*/
  Y      = new su2double[nSpecies];
  
  dYdr = new su2double*[nSpecies];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dYdr[iSpecies] = new su2double[nSpecies];
  }

  residual = new su2double[nVar]();
  jacobian = new su2double* [nVar];
  for(unsigned short iVar = 0; iVar < nVar; ++iVar)
    jacobian[iVar] = new su2double [nVar]();

}

CSource_NEMO::~CSource_NEMO(void) {
  unsigned short iSpecies;

  /*--- Deallocate arrays ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    delete [] dYdr[iSpecies];
  delete [] dYdr;
  delete [] Y;

  delete [] residual;
  if (jacobian){
    for(unsigned short iVar = 0; iVar < nVar; ++iVar)
      delete [] jacobian[iVar];
    delete [] jacobian;
  }

}

CNumerics::ResidualType<> CSource_NEMO::ComputeChemistry(const CConfig *config) {

  /*--- Nonequilibrium chemistry ---*/
  unsigned short iSpecies, iVar;
  unsigned short jSpecies, jVar;
  su2double T, Tve;
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }
  if (implicit)
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;

  /*--- Rename for convenience ---*/
  T   = V_i[T_INDEX];
  Tve = V_i[TVE_INDEX];
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  /*--- Set mixture state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  /*---Compute Prodcution/destruction terms ---*/
  ws = fluidmodel->ComputeNetProductionRates(implicit, V_i, eve_i, Cvve_i,
                                             dTdU_i, dTvedU_i, jacobian);

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) 
    residual[iSpecies] = ws[iSpecies] *Volume;

  if (implicit)
    for (iVar = 0; iVar<nVar; iVar++)
      for (jVar = 0; jVar<nVar; jVar++)
        jacobian[iVar][jVar] = jacobian[iVar][jVar] * Volume;

  return ResidualType<>(residual, jacobian, nullptr);

}

CNumerics::ResidualType<> CSource_NEMO::ComputeVibRelaxation(const CConfig *config) {

  /*--- Trans.-rot. & vibrational energy exchange via inelastic collisions ---*/
  // Note: Electronic energy not implemented
  // Note: Landau-Teller formulation
  // Note: Millikan & White relaxation time (requires P in Atm.)
  // Note: Park limiting cross section
  unsigned short iSpecies, iVar;
  unsigned short jSpecies, jVar;
  su2double T, Tve;
  su2double VTterm;
  su2double res_min = -1E6;
  su2double res_max = 1E6;
  vector<su2double> rhos;

  rhos.resize(nSpecies,0.0);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }
  if (implicit) {
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        jacobian[iVar][jVar] = 0.0;
  }

  /*--- Rename for convenience ---*/
  T   = V_i[T_INDEX];
  Tve = V_i[TVE_INDEX];
  for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    rhos[iSpecies]=V_i[RHOS_INDEX+iSpecies];

  /*--- Set fluid state ---*/
  fluidmodel->SetTDStateRhosTTv(rhos, T, Tve);

  /*--- Compute residual and jacobians ---*/
  VTterm = fluidmodel -> ComputeEveSourceTerm();
  if (implicit) 
    fluidmodel->GetEveSourceTermJacobian(V_i, eve_i, Cvve_i, dTdU_i,
                                         dTvedU_i, jacobian);

  residual[nSpecies+nDim+1] = VTterm * Volume;
  
  if (implicit)
    for (iVar = 0; iVar<nVar; iVar++)
      for (jVar = 0; jVar<nVar; jVar++)
        jacobian[iVar][jVar] = jacobian[iVar][jVar] * Volume; 

  /*--- Relax/limit vt transfer ---*/
  if(config->GetVTTransferResidualLimiting()){
    if(residual[nSpecies+nDim+1]>res_max) residual[nSpecies+nDim+1]=res_max;
    if(residual[nSpecies+nDim+1]<res_min) residual[nSpecies+nDim+1]=res_min;
  }

  return ResidualType<>(residual, jacobian, nullptr);
}

CNumerics::ResidualType<> CSource_NEMO::ComputeAxisymmetric(const CConfig *config) {

  unsigned short iDim, iSpecies, iVar, jVar,jSpecies;
  su2double rho, rhou, rhov, rhoEve, vel2, H, yinv, T, Tve, Ru, RuSI;
  su2double *Ds, **GV, ktr, kve;

  /*--- Rename for convenience ---*/
  Ds = Diffusion_Coeff_i;
  ktr = Thermal_Conductivity_i;
  kve = Thermal_Conductivity_ve_i;
  rho = V_i[RHO_INDEX];
  T   = V_i[T_INDEX];
  Tve  = V_i[TVE_INDEX];
  GV  = PrimVar_Grad_i;
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;

  auto& Ms = fluidmodel->GetSpeciesMolarMass();

  bool viscous = config->GetViscous();
  bool rans = (config->GetKind_Turb_Model() != NONE);
  hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, eve_i);

  /*--- Initialize residual and Jacobian arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    residual[iVar] = 0.0;
  }

  /*--- Calculate inverse of y coordinate ---*/
  if (Coord_i[1]!= 0.0) yinv = 1.0/Coord_i[1];
  else yinv = 0.0;

  /*--- Rename for convenience ---*/
  rho    = V_i[RHO_INDEX];
  rhou   = U_i[nSpecies];
  rhov   = U_i[nSpecies+1];
  H      = V_i[H_INDEX];
  rhoEve = U_i[nVar-1];
  vel2   = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    vel2 += V_i[VEL_INDEX+iDim]*V_i[VEL_INDEX+iDim];
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Y[iSpecies] = V_i[RHOS_INDEX+iSpecies] / rho;

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    residual[iSpecies] = yinv*rhov*Y[iSpecies]*Volume;
  residual[nSpecies]   = yinv*rhov*U_i[nSpecies]/rho*Volume;
  residual[nSpecies+1] = yinv*rhov*U_i[nSpecies+1]/rho*Volume;
  residual[nSpecies+2] = yinv*rhov*H*Volume;
  residual[nSpecies+3] = yinv*rhov*U_i[nSpecies+nDim+1]/rho*Volume;

 if (implicit) {

   /*--- Initialize ---*/
   for (iVar = 0; iVar < nVar; iVar++)
     for (jVar = 0; jVar < nVar; jVar++)
       jacobian[iVar][jVar] = 0.0;
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
     for (jSpecies = 0; jSpecies < nSpecies; jSpecies++)
       dYdr[iSpecies][jSpecies] = 0.0;

   /*--- Calculate additional quantities ---*/
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
     for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
       dYdr[iSpecies][jSpecies] += -1/rho*Y[iSpecies];
     }
     dYdr[iSpecies][iSpecies] += 1/rho;
   }

   /*--- Populate Jacobian ---*/

   // Species density
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
     for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
       jacobian[iSpecies][jSpecies] = dYdr[iSpecies][jSpecies]*rhov;
     }
     jacobian[iSpecies][nSpecies+1] = Y[iSpecies];
   }

   // X-momentum
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
     jacobian[nSpecies][iSpecies] = -rhou*rhov/(rho*rho);
   jacobian[nSpecies][nSpecies]   = rhov/rho;
   jacobian[nSpecies][nSpecies+1] = rhou/rho;

   // Y-momentum
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
     jacobian[nSpecies+1][iSpecies] = -rhov*rhov/(rho*rho);
   jacobian[nSpecies+1][nSpecies+1] = 2*rhov/rho;

   // Energy
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
     jacobian[nSpecies+nDim][iSpecies]      = -H*rhov/rho + dPdU_i[iSpecies]*rhov/rho;
   jacobian[nSpecies+nDim][nSpecies]        = dPdU_i[nSpecies]*rhov/rho;
   jacobian[nSpecies+nDim][nSpecies+1]      = H + dPdU_i[nSpecies+1]*rhov/rho;
   jacobian[nSpecies+nDim][nSpecies+nDim]   = (1+dPdU_i[nSpecies+nDim])*rhov/rho;
   jacobian[nSpecies+nDim][nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1]*rhov/rho;

   // Vib-el energy
   for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
     jacobian[nSpecies+nDim+1][iSpecies] = -rhoEve*rhov/(rho*rho);
   jacobian[nSpecies+nDim+1][nSpecies+1] = rhoEve/rho;
   jacobian[nSpecies+nDim+1][nSpecies+nDim+1] = rhov/rho;

   for (iVar = 0; iVar < nVar; iVar++)
     for (jVar = 0; jVar < nVar; jVar++)
       jacobian[iVar][jVar] *= yinv*Volume;
 }

  return ResidualType<>(residual, jacobian, nullptr);

}
