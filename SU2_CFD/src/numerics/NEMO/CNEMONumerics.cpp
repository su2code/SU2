/*!
 * \file CNEMONumerics.cpp
 * \brief Implementation of the base for NEMO numerics classes.
 *        Contains methods for common tasks, e.g. compute flux
 *        Jacobians.
 * \author S.R. Copeland, W. Maier, C. Garbacz
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

#include "../../../include/numerics/NEMO/CNEMONumerics.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"
CNEMONumerics::CNEMONumerics(unsigned short val_nDim, unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             const CConfig* config) :
                             CNumerics(val_nDim, val_nVar, config) {

    nSpecies     = nVar - nDim - 2;
    nPrimVar     = val_nPrimVar;
    nPrimVarGrad = val_nPrimVarGrad;

    RHOS_INDEX      = 0;
    T_INDEX         = nSpecies;
    TVE_INDEX       = nSpecies+1;
    VEL_INDEX       = nSpecies+2;
    P_INDEX         = nSpecies+nDim+2;
    RHO_INDEX       = nSpecies+nDim+3;
    H_INDEX         = nSpecies+nDim+4;
    A_INDEX         = nSpecies+nDim+5;
    RHOCVTR_INDEX   = nSpecies+nDim+6;
    RHOCVVE_INDEX   = nSpecies+nDim+7;
    LAM_VISC_INDEX  = nSpecies+nDim+8;
    EDDY_VISC_INDEX = nSpecies+nDim+9;

    /*--- Read from CConfig ---*/
    implicit   = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

    ionization = config->GetIonization();
    if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
    else            { nHeavy = nSpecies;   nEl = 0; }

    /*--- Instatiate the correct fluid model ---*/
    switch (config->GetKind_FluidModel()) {
      case MUTATIONPP:
        #if defined(HAVE_MPP) && !defined(CODI_REVERSE_TYPE) && !defined(CODI_FORWARD_TYPE)
          fluidmodel = new CMutationTCLib(config, nDim);
        #else
          SU2_MPI::Error(string("Mutation++ has not been configured/compiled. Add 1) '-Denable-mpp=true' to your meson string or 2) '-DHAVE_MPP' to the CXX FLAGS of your configure string, and recompile."),
          CURRENT_FUNCTION);
        #endif
      break;
      case SU2_NONEQ:
        fluidmodel = new CSU2TCLib(config, nDim, false);
      break;
    }
}

CNEMONumerics::~CNEMONumerics() {

  delete fluidmodel;
}

void CNEMONumerics::GetInviscidProjFlux(const su2double *val_U,
                                        const su2double *val_V,
                                        const su2double *val_normal,
                                        su2double *val_Proj_Flux) {

  const su2double *rhos;

  /*--- Initialize vectors ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    val_Proj_Flux[iVar] = 0.0;

  /*--- Rename for convienience ---*/
  const su2double rho = val_V[RHO_INDEX];
  const su2double u = val_V[VEL_INDEX];
  const su2double v = val_V[VEL_INDEX+1];
  const su2double w = val_V[VEL_INDEX+2];
  const su2double P = val_V[P_INDEX];
  const su2double H = val_V[H_INDEX];
  const su2double rhoEve = val_U[nSpecies+nDim+1];
  rhos   = &val_V[RHOS_INDEX];

  if (nDim == 2) {

    /*--- iDim = 0 (x-direction) ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies]  = (rhos[iSpecies]*u) * val_normal[0];
    val_Proj_Flux[nSpecies]    = (rho*u*u + P)      * val_normal[0];
    val_Proj_Flux[nSpecies+1]  = (rho*u*v)          * val_normal[0];
    val_Proj_Flux[nSpecies+2]  = (rho*u*H)          * val_normal[0];
    val_Proj_Flux[nSpecies+3]  = (rhoEve*u)         * val_normal[0];

    /*---- iDim = 1 (y-direction) ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies] += (rhos[iSpecies]*v) * val_normal[1];
    val_Proj_Flux[nSpecies]   += (rho*v*u)          * val_normal[1];
    val_Proj_Flux[nSpecies+1] += (rho*v*v + P)      * val_normal[1];
    val_Proj_Flux[nSpecies+2] += (rho*v*H)          * val_normal[1];
    val_Proj_Flux[nSpecies+3] += (rhoEve*v)         * val_normal[1];
  }
  else {

    /*--- iDim = 0 (x-direction) ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies]  = (rhos[iSpecies]*u) * val_normal[0];
    val_Proj_Flux[nSpecies]    = (rho*u*u + P)      * val_normal[0];
    val_Proj_Flux[nSpecies+1]  = (rho*u*v)          * val_normal[0];
    val_Proj_Flux[nSpecies+2]  = (rho*u*w)          * val_normal[0];
    val_Proj_Flux[nSpecies+3]  = (rho*u*H)          * val_normal[0];
    val_Proj_Flux[nSpecies+4]  = (rhoEve*u)         * val_normal[0];

    /*--- iDim = 0 (y-direction) ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies] += (rhos[iSpecies]*v) * val_normal[1];
    val_Proj_Flux[nSpecies]   += (rho*v*u)          * val_normal[1];
    val_Proj_Flux[nSpecies+1] += (rho*v*v + P)      * val_normal[1];
    val_Proj_Flux[nSpecies+2] += (rho*v*w)          * val_normal[1];
    val_Proj_Flux[nSpecies+3] += (rho*v*H)          * val_normal[1];
    val_Proj_Flux[nSpecies+4] += (rhoEve*v)         * val_normal[1];

    /*--- iDim = 0 (z-direction) ---*/
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies] += (rhos[iSpecies]*w) * val_normal[2];
    val_Proj_Flux[nSpecies]   += (rho*w*u)          * val_normal[2];
    val_Proj_Flux[nSpecies+1] += (rho*w*v)          * val_normal[2];
    val_Proj_Flux[nSpecies+2] += (rho*w*w + P)      * val_normal[2];
    val_Proj_Flux[nSpecies+3] += (rho*w*H)          * val_normal[2];
    val_Proj_Flux[nSpecies+4] += (rhoEve*w)         * val_normal[2];
  }
}

void CNEMONumerics::GetInviscidProjJac(const su2double *val_U,    const su2double *val_V,
                                       const su2double *val_dPdU, const su2double *val_normal,
                                       const su2double val_scale,
                                       su2double **val_Proj_Jac_Tensor) {

  const su2double *rhos;
  rhos = &val_V[RHOS_INDEX];

  /*--- Initialize the Jacobian tensor ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++)
      val_Proj_Jac_Tensor[iVar][jVar] = 0.0;

  /*--- Rename for convenience ---*/
  const su2double rho = val_V[RHO_INDEX];
  const su2double H = val_V[H_INDEX];
  const su2double rhoEve = val_U[nSpecies+nDim+1];

  su2double u[MAXNDIM];
  for (auto iDim = 0ul; iDim < nDim; iDim++)
    u[iDim] = val_V[VEL_INDEX+iDim];

  /*--- Calculate projected velocity ---*/
  su2double proj_vel = GeometryToolbox::DotProduct(nDim, u, val_normal);

  /*--- Species density rows ---*/
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
      val_Proj_Jac_Tensor[iSpecies][jSpecies] += -(rhos[iSpecies]/rho) * proj_vel;
    }
    val_Proj_Jac_Tensor[iSpecies][iSpecies]   += proj_vel;
    for (auto iDim  = 0ul; iDim < nDim; iDim++) {
      val_Proj_Jac_Tensor[iSpecies][nSpecies+iDim] += (rhos[iSpecies]/rho) * val_normal[iDim];
      val_Proj_Jac_Tensor[nSpecies+iDim][iSpecies] += val_dPdU[iSpecies]*val_normal[iDim] - proj_vel*u[iDim];
    }
    val_Proj_Jac_Tensor[nSpecies+nDim][iSpecies]   += (val_dPdU[iSpecies]-H) * proj_vel;
    val_Proj_Jac_Tensor[nSpecies+nDim+1][iSpecies] += -proj_vel * rhoEve/rho;
  }

  /*--- Momentum rows ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    for (auto jDim = 0ul; jDim < nDim; jDim++) {
      val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+jDim] += val_dPdU[nSpecies+jDim]*val_normal[iDim] + u[iDim]*val_normal[jDim];
    }
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+iDim]   += proj_vel;
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+nDim]   += val_dPdU[nSpecies+nDim]*val_normal[iDim];
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+nDim+1] += val_dPdU[nSpecies+nDim+1]*val_normal[iDim];
  }

  /*--- Total energy row ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+iDim] += val_dPdU[nSpecies+iDim]*proj_vel + H*val_normal[iDim];
  val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+nDim]   += (1+val_dPdU[nSpecies+nDim])*proj_vel;
  val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+nDim+1] +=  val_dPdU[nSpecies+nDim+1] *proj_vel;

  /*--- Vib.-el. energy row ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nSpecies+nDim+1][nSpecies+iDim] = rhoEve/rho*val_normal[iDim];
  val_Proj_Jac_Tensor[nSpecies+nDim+1][nSpecies+nDim+1] = proj_vel;

  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++)
      val_Proj_Jac_Tensor[iVar][jVar] = val_scale * val_Proj_Jac_Tensor[iVar][jVar];
}

void CNEMONumerics::GetViscousProjFlux(const su2double *val_primvar,
                                       const su2double* const* val_gradprimvar,
                                       su2double *val_eve,
                                       const su2double *val_normal,
                                       const su2double *val_diffusioncoeff,
                                       su2double val_lam_viscosity,
                                       su2double val_eddy_viscosity,
                                       su2double val_therm_conductivity,
                                       su2double val_therm_conductivity_ve,
                                       const CConfig *config) {

  // Requires a slightly non-standard primitive vector:
  // Assumes -     V = [Y1, ... , Yn, T, Tve, ... ]
  // and gradient GV = [GY1, ... , GYn, GT, GTve, ... ]
  // rather than the standard V = [r1, ... , rn, T, Tve, ... ]

  su2activematrix Flux_Tensor(nVar,nDim);

  /*--- Initialize ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (auto iDim = 0ul; iDim < nDim; iDim++)
      Flux_Tensor[iVar][iDim] = 0.0;
  }

  /*--- Rename variables for convenience ---*/
  const auto& Ds  = val_diffusioncoeff;
  const su2double mu  = val_lam_viscosity+val_eddy_viscosity;
  const su2double ktr = val_therm_conductivity;
  const su2double kve = val_therm_conductivity_ve;
  const su2double rho = val_primvar[RHO_INDEX];
  const su2double T = val_primvar[T_INDEX];
  const su2double Tve = val_primvar[TVE_INDEX];
  const auto& V   = val_primvar;
  const auto& GV  = val_gradprimvar;
  const auto& hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, val_eve);

  /*--- Pre-compute mixture quantities ---*/  //TODO
  su2double Vector[MAXNDIM] = {0.0};
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    for (auto iSpecies = 0ul; iSpecies < nHeavy; iSpecies++) {
      Vector[iDim] += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim];
    }
  }

  /*--- Compute the viscous stress tensor ---*/
  ComputeStressTensor(nDim,tau,val_gradprimvar+VEL_INDEX, mu);

  /*--- Populate entries in the viscous flux vector ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++) {

    /*--- Species diffusion velocity ---*/
    for (auto iSpecies = 0ul; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[iSpecies][iDim] = rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim]
          - V[RHOS_INDEX+iSpecies]*Vector[iDim];
    }

    /*--- Shear-stress/momentum related terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim] = 0.0;
    for (auto jDim = 0ul; jDim < nDim; jDim++) {
      Flux_Tensor[nSpecies+jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[nSpecies+nDim][iDim] += tau[iDim][jDim]*val_primvar[VEL_INDEX+jDim];
    }

    /*--- Diffusion terms ---*/
    for (auto iSpecies = 0ul; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[nSpecies+nDim][iDim]   += Flux_Tensor[iSpecies][iDim] * hs[iSpecies];
      Flux_Tensor[nSpecies+nDim+1][iDim] += Flux_Tensor[iSpecies][iDim] * val_eve[iSpecies];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim] += ktr*GV[T_INDEX][iDim] + kve*GV[TVE_INDEX][iDim];
    Flux_Tensor[nSpecies+nDim+1][iDim] += kve*GV[TVE_INDEX][iDim];
  }

  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
    }
  }
}

#define COMPUTE_VISCOUS_JACS(NVAR, NSPECIES) \
ComputeViscousJacs_impl<NVAR, NSPECIES>(val_Mean_PrimVar, \
                                        val_Mean_Eve, val_Mean_Cvve, val_diffusion_coeff, \
                                        val_laminar_viscosity, val_eddy_viscosity, \
                                        val_thermal_conductivity, val_thermal_conductivity_ve, \
                                        val_dist_ij, val_normal, val_dS, val_Fv, \
                                        val_Jac_i, val_Jac_j, config);


void CNEMONumerics::GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                                       su2double *val_Mean_Eve, const su2double *val_Mean_Cvve,
                                       const su2double *val_diffusion_coeff, su2double val_laminar_viscosity,
                                       su2double val_eddy_viscosity, su2double val_thermal_conductivity,
                                       su2double val_thermal_conductivity_ve,
                                       su2double val_dist_ij,
                                       const su2double *val_normal,
                                       su2double val_dS, const su2double *val_Fv,
                                       su2double **val_Jac_i, su2double **val_Jac_j,
                                       const CConfig *config){

  switch (nVar) {
    case 5:
      return COMPUTE_VISCOUS_JACS(5, 1);

    case 6:
      switch (nSpecies) {
        case 1: return COMPUTE_VISCOUS_JACS(6, 1);

        case 2: return COMPUTE_VISCOUS_JACS(6, 2);

        default: SU2_MPI::Error("nVar and nSpecies mismatch.", CURRENT_FUNCTION);
      }
    break;
    case 7:
      return COMPUTE_VISCOUS_JACS(7, 2);

     case 9:
      return COMPUTE_VISCOUS_JACS(9, 5);

    case 10:
      return COMPUTE_VISCOUS_JACS(10, 5);

    case 11:
      return COMPUTE_VISCOUS_JACS(11, 7);

    case 12:
      return COMPUTE_VISCOUS_JACS(12, 7);

    default:
      return COMPUTE_VISCOUS_JACS(DynamicSize,DynamicSize);

  }
}
#undef COMPUTE_VISCOUS_JACS

void CNEMONumerics::GetPMatrix(const su2double *U, const su2double *V, const su2double *val_dPdU,
                               const su2double *val_normal, const su2double *l, const su2double *m,
                               su2double **val_p_tensor) const {

  /*--- Initialize the P matrix to zero ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++)
      val_p_tensor[iVar][jVar] = 0.0;

  /*--- Pre-compute useful quantities ---*/
  su2double sqvel = 0.0;
  const su2double rho = V[RHO_INDEX];
  const su2double eve = U[nSpecies+nDim+1]/rho;
  su2double vU, vV, vW;
  vU = vV = vW = 0.0;
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    vU    += V[VEL_INDEX+iDim] * val_normal[iDim];
    vV    += V[VEL_INDEX+iDim] * l[iDim];
    vW    += V[VEL_INDEX+iDim] * m[iDim];
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  }
  const su2double a  = V[A_INDEX];
  const su2double a2 = V[A_INDEX]*V[A_INDEX];

  if(nDim == 2) {
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      val_p_tensor[iSpecies][iSpecies]   += 1.0/a2;
      val_p_tensor[iSpecies][nSpecies]   += 0.0;
      val_p_tensor[iSpecies][nSpecies+1] += V[RHOS_INDEX+iSpecies] / (2.0*rho*a2);
      val_p_tensor[iSpecies][nSpecies+2] += V[RHOS_INDEX+iSpecies] / (2.0*rho*a2);
      val_p_tensor[iSpecies][nSpecies+3] += 0.0;

      val_p_tensor[nSpecies][iSpecies]   += V[VEL_INDEX]   / a2;
      val_p_tensor[nSpecies+1][iSpecies] += V[VEL_INDEX+1] / a2;
      val_p_tensor[nSpecies+2][iSpecies] += (val_dPdU[nSpecies+nDim]*sqvel-val_dPdU[iSpecies]) /
                                            (val_dPdU[nSpecies+nDim]*a2);
      val_p_tensor[nSpecies+3][iSpecies] += 0.0;
    }

    for (auto iDim = 0ul; iDim < nDim; iDim++){
      val_p_tensor[nSpecies+iDim][nSpecies]     += l[iDim];
      val_p_tensor[nSpecies+iDim][nSpecies+1]   += (V[VEL_INDEX+iDim]+a*val_normal[iDim]) / (2.0*a2);
      val_p_tensor[nSpecies+iDim][nSpecies+2]   += (V[VEL_INDEX+iDim]-a*val_normal[iDim]) / (2.0*a2);
      val_p_tensor[nSpecies+iDim][nSpecies+3]   += 0.0;
    }

    val_p_tensor[nSpecies+2][nSpecies]   += vV;
    val_p_tensor[nSpecies+2][nSpecies+1] += ((V[H_INDEX])+a*vU) / (2.0*a2);
    val_p_tensor[nSpecies+2][nSpecies+2] += ((V[H_INDEX])-a*vU) / (2.0*a2);
    val_p_tensor[nSpecies+2][nSpecies+3] += -val_dPdU[nSpecies+nDim+1] / (val_dPdU[nSpecies+nDim]*a2);

    val_p_tensor[nSpecies+3][nSpecies]   += 0.0;
    val_p_tensor[nSpecies+3][nSpecies+1] += eve / (2.0*a2);
    val_p_tensor[nSpecies+3][nSpecies+2] += eve / (2.0*a2);
    val_p_tensor[nSpecies+3][nSpecies+3] += 1.0 / a2;

  } else {

    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      val_p_tensor[iSpecies][iSpecies]   = 1.0/a2;
      val_p_tensor[iSpecies][nSpecies]   = 0.0;
      val_p_tensor[iSpecies][nSpecies+1] = 0.0;
      val_p_tensor[iSpecies][nSpecies+2] = V[RHOS_INDEX+iSpecies] / (2.0*rho*a2);
      val_p_tensor[iSpecies][nSpecies+3] = V[RHOS_INDEX+iSpecies] / (2.0*rho*a2);
      val_p_tensor[iSpecies][nSpecies+4] = 0.0;

      val_p_tensor[nSpecies][iSpecies]   = V[VEL_INDEX]   / a2;
      val_p_tensor[nSpecies+1][iSpecies] = V[VEL_INDEX+1] / a2;
      val_p_tensor[nSpecies+2][iSpecies] = V[VEL_INDEX+2] / a2;
      val_p_tensor[nSpecies+3][iSpecies] = (val_dPdU[nSpecies+3]*sqvel-val_dPdU[iSpecies])
          / (val_dPdU[nSpecies+3]*a2);
      val_p_tensor[nSpecies+4][iSpecies] = 0.0;
    }

    for (auto iDim = 0ul; iDim < nDim; iDim++){
      val_p_tensor[nSpecies+iDim][nSpecies]     = l[iDim];
      val_p_tensor[nSpecies+iDim][nSpecies+1]   = m[iDim];
      val_p_tensor[nSpecies+iDim][nSpecies+2]   = (V[VEL_INDEX+iDim]+a*val_normal[iDim]) / (2.0*a2);
      val_p_tensor[nSpecies+iDim][nSpecies+3]   = (V[VEL_INDEX+iDim]-a*val_normal[iDim]) / (2.0*a2);
      val_p_tensor[nSpecies+iDim][nSpecies+4]   = 0.0;
    }

    val_p_tensor[nSpecies+3][nSpecies]   = vV;
    val_p_tensor[nSpecies+3][nSpecies+1] = vW;
    val_p_tensor[nSpecies+3][nSpecies+2] = ((V[H_INDEX])+a*vU) / (2.0*a2);
    val_p_tensor[nSpecies+3][nSpecies+3] = ((V[H_INDEX])-a*vU) / (2.0*a2);
    val_p_tensor[nSpecies+3][nSpecies+4] = -val_dPdU[nSpecies+nDim+1] / (val_dPdU[nSpecies+nDim]*a2);

    val_p_tensor[nSpecies+4][nSpecies]   = 0.0;
    val_p_tensor[nSpecies+4][nSpecies+1] = 0.0;
    val_p_tensor[nSpecies+4][nSpecies+2] = eve / (2.0*a2);
    val_p_tensor[nSpecies+4][nSpecies+3] = eve / (2.0*a2);
    val_p_tensor[nSpecies+4][nSpecies+4] = 1.0 / a2;
  }
}

void CNEMONumerics::GetPMatrix_inv(const su2double *U, const su2double *V, const su2double *val_dPdU,
                                   const su2double *val_normal, const su2double *l, const su2double *m,
                                   su2double **val_invp_tensor) const {

  for (auto iVar = 0ul; iVar < nVar; iVar++)
    for (auto jVar = 0ul; jVar < nVar; jVar++)
      val_invp_tensor[iVar][jVar] = 0.0;

  /*--- Pre-compute useful quantities ---*/
  const su2double rho = V[RHO_INDEX];
  const su2double eve = U[nSpecies+nDim+1]/rho;
  su2double vU, vV, vW;
  vU = vV = vW = 0.0;
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    vU += V[VEL_INDEX+iDim] * val_normal[iDim];
    vV += V[VEL_INDEX+iDim] * l[iDim];
    vW += V[VEL_INDEX+iDim] * m[iDim];
  }
  const su2double a  = V[A_INDEX];
  const su2double a2 = V[A_INDEX]*V[A_INDEX];

  if (nDim == 2) {

    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
        val_invp_tensor[iSpecies][jSpecies] += -(V[RHOS_INDEX+iSpecies]/rho) * val_dPdU[jSpecies];
      }
      val_invp_tensor[iSpecies][iSpecies]   += a2;
      val_invp_tensor[iSpecies][nSpecies]   += val_dPdU[nSpecies+nDim] * V[VEL_INDEX] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+1] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+1] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+2] += -val_dPdU[nSpecies+nDim] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+3] += -val_dPdU[nSpecies+nDim+1] * (V[RHOS_INDEX+iSpecies]/rho);

      val_invp_tensor[nSpecies][iSpecies]   += -vV;
      val_invp_tensor[nSpecies+1][iSpecies] += val_dPdU[iSpecies] - vU*a;
      val_invp_tensor[nSpecies+2][iSpecies] += val_dPdU[iSpecies] + vU*a;
      val_invp_tensor[nSpecies+3][iSpecies] += -eve * val_dPdU[iSpecies];
    }

    val_invp_tensor[nSpecies][nSpecies]     += l[0];
    val_invp_tensor[nSpecies][nSpecies+1]   += l[1];
    val_invp_tensor[nSpecies][nSpecies+2]   += 0.0;
    val_invp_tensor[nSpecies][nSpecies+3]   += 0.0;

    val_invp_tensor[nSpecies+1][nSpecies]   += a*val_normal[0] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX];
    val_invp_tensor[nSpecies+1][nSpecies+1] += a*val_normal[1] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+1];
    val_invp_tensor[nSpecies+1][nSpecies+2] += val_dPdU[nSpecies+nDim];
    val_invp_tensor[nSpecies+1][nSpecies+3] += val_dPdU[nSpecies+nDim+1];

    val_invp_tensor[nSpecies+2][nSpecies]   += -a*val_normal[0] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX];
    val_invp_tensor[nSpecies+2][nSpecies+1] += -a*val_normal[1] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+1];
    val_invp_tensor[nSpecies+2][nSpecies+2] += val_dPdU[nSpecies+nDim];
    val_invp_tensor[nSpecies+2][nSpecies+3] += val_dPdU[nSpecies+nDim+1];

    val_invp_tensor[nSpecies+3][nSpecies]   += val_dPdU[nSpecies+nDim] * V[VEL_INDEX] * eve;
    val_invp_tensor[nSpecies+3][nSpecies+1] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+1] * eve;
    val_invp_tensor[nSpecies+3][nSpecies+2] += -val_dPdU[nSpecies+nDim] * eve;
    val_invp_tensor[nSpecies+3][nSpecies+3] += a2 - val_dPdU[nSpecies+nDim+1]*eve;

  } else {

    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      for (auto jSpecies = 0ul; jSpecies < nSpecies; jSpecies++) {
        val_invp_tensor[iSpecies][jSpecies] += -(V[RHOS_INDEX+iSpecies]/rho) * val_dPdU[jSpecies];
      }
      val_invp_tensor[iSpecies][iSpecies]   += a2;
      val_invp_tensor[iSpecies][nSpecies]   += val_dPdU[nSpecies+nDim] * V[VEL_INDEX] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+1] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+1] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+2] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+2] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+3] += -val_dPdU[nSpecies+nDim] * (V[RHOS_INDEX+iSpecies]/rho);
      val_invp_tensor[iSpecies][nSpecies+4] += -val_dPdU[nSpecies+nDim+1] * (V[RHOS_INDEX+iSpecies]/rho);

      val_invp_tensor[nSpecies][iSpecies]   += -vV;
      val_invp_tensor[nSpecies+1][iSpecies] += -vW;
      val_invp_tensor[nSpecies+2][iSpecies] += val_dPdU[iSpecies] - vU*a;
      val_invp_tensor[nSpecies+3][iSpecies] += val_dPdU[iSpecies] + vU*a;
      val_invp_tensor[nSpecies+4][iSpecies] += -eve * val_dPdU[iSpecies];
    }

    val_invp_tensor[nSpecies][nSpecies]     += l[0];
    val_invp_tensor[nSpecies][nSpecies+1]   += l[1];
    val_invp_tensor[nSpecies][nSpecies+2]   += l[2];
    val_invp_tensor[nSpecies][nSpecies+3]   += 0.0;
    val_invp_tensor[nSpecies][nSpecies+4]   += 0.0;

    val_invp_tensor[nSpecies+1][nSpecies]   += m[0];
    val_invp_tensor[nSpecies+1][nSpecies+1] += m[1];
    val_invp_tensor[nSpecies+1][nSpecies+2] += m[2];
    val_invp_tensor[nSpecies+1][nSpecies+3] += 0.0;
    val_invp_tensor[nSpecies+1][nSpecies+4] += 0.0;

    val_invp_tensor[nSpecies+2][nSpecies]   += a*val_normal[0] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX];
    val_invp_tensor[nSpecies+2][nSpecies+1] += a*val_normal[1] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+1];
    val_invp_tensor[nSpecies+2][nSpecies+2] += a*val_normal[2] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+2];
    val_invp_tensor[nSpecies+2][nSpecies+3] += val_dPdU[nSpecies+nDim];
    val_invp_tensor[nSpecies+2][nSpecies+4] += val_dPdU[nSpecies+nDim+1];

    val_invp_tensor[nSpecies+3][nSpecies]   += -a*val_normal[0] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX];
    val_invp_tensor[nSpecies+3][nSpecies+1] += -a*val_normal[1] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+1];
    val_invp_tensor[nSpecies+3][nSpecies+2] += -a*val_normal[2] - val_dPdU[nSpecies+nDim]*V[VEL_INDEX+2];
    val_invp_tensor[nSpecies+3][nSpecies+3] += val_dPdU[nSpecies+nDim];
    val_invp_tensor[nSpecies+3][nSpecies+4] += val_dPdU[nSpecies+nDim+1];

    val_invp_tensor[nSpecies+4][nSpecies]   += val_dPdU[nSpecies+nDim] * V[VEL_INDEX] * eve;
    val_invp_tensor[nSpecies+4][nSpecies+1] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+1] * eve;
    val_invp_tensor[nSpecies+4][nSpecies+2] += val_dPdU[nSpecies+nDim] * V[VEL_INDEX+2] * eve;
    val_invp_tensor[nSpecies+4][nSpecies+3] += -val_dPdU[nSpecies+nDim] * eve;
    val_invp_tensor[nSpecies+4][nSpecies+4] += a2 - val_dPdU[nSpecies+nDim+1]*eve;
  }
}
