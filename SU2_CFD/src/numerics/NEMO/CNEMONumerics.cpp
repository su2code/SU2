/*!
 * \file CNEMONumerics.cpp
 * \brief Implementation of the base for NEMO numerics classes.
 *        Contains methods for common tasks, e.g. compute flux
 *        Jacobians.
 * \author S.R. Copeland, W. Maier, C. Garbacz
 * \version 7.2.1 "Blackbird"
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

    hs.resize(nSpecies,0.0);

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

CNEMONumerics::~CNEMONumerics(void) {

  delete fluidmodel;
}

void CNEMONumerics::GetInviscidProjFlux(const su2double *val_U,
                                        const su2double *val_V,
                                        const su2double *val_normal,
                                        su2double *val_Proj_Flux) {
  unsigned short iSpecies, iVar;
  su2double rho, u, v, w, rhoEve, P, H;
  const su2double *rhos;

  /*--- Initialize vectors ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_Proj_Flux[iVar] = 0.0;

  /*--- Rename for convienience ---*/
  rho    = val_V[RHO_INDEX];
  u      = val_V[VEL_INDEX];
  v      = val_V[VEL_INDEX+1];
  w      = val_V[VEL_INDEX+2];
  P      = val_V[P_INDEX];
  H      = val_V[H_INDEX];
  rhoEve = val_U[nSpecies+nDim+1];
  rhos   = &val_V[RHOS_INDEX];

  if (nDim == 2) {

    /*--- iDim = 0 (x-direction) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies]  = (rhos[iSpecies]*u) * val_normal[0];
    val_Proj_Flux[nSpecies]    = (rho*u*u + P)      * val_normal[0];
    val_Proj_Flux[nSpecies+1]  = (rho*u*v)          * val_normal[0];
    val_Proj_Flux[nSpecies+2]  = (rho*u*H)          * val_normal[0];
    val_Proj_Flux[nSpecies+3]  = (rhoEve*u)         * val_normal[0];

    /*---- iDim = 1 (y-direction) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies] += (rhos[iSpecies]*v) * val_normal[1];
    val_Proj_Flux[nSpecies]   += (rho*v*u)          * val_normal[1];
    val_Proj_Flux[nSpecies+1] += (rho*v*v + P)      * val_normal[1];
    val_Proj_Flux[nSpecies+2] += (rho*v*H)          * val_normal[1];
    val_Proj_Flux[nSpecies+3] += (rhoEve*v)         * val_normal[1];
  }
  else {

    /*--- iDim = 0 (x-direction) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies]  = (rhos[iSpecies]*u) * val_normal[0];
    val_Proj_Flux[nSpecies]    = (rho*u*u + P)      * val_normal[0];
    val_Proj_Flux[nSpecies+1]  = (rho*u*v)          * val_normal[0];
    val_Proj_Flux[nSpecies+2]  = (rho*u*w)          * val_normal[0];
    val_Proj_Flux[nSpecies+3]  = (rho*u*H)          * val_normal[0];
    val_Proj_Flux[nSpecies+4]  = (rhoEve*u)         * val_normal[0];

    /*--- iDim = 0 (y-direction) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
      val_Proj_Flux[iSpecies] += (rhos[iSpecies]*v) * val_normal[1];
    val_Proj_Flux[nSpecies]   += (rho*v*u)          * val_normal[1];
    val_Proj_Flux[nSpecies+1] += (rho*v*v + P)      * val_normal[1];
    val_Proj_Flux[nSpecies+2] += (rho*v*w)          * val_normal[1];
    val_Proj_Flux[nSpecies+3] += (rho*v*H)          * val_normal[1];
    val_Proj_Flux[nSpecies+4] += (rhoEve*v)         * val_normal[1];

    /*--- iDim = 0 (z-direction) ---*/
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
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
  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    for (unsigned short jVar = 0; jVar < nVar; jVar++)
      val_Proj_Jac_Tensor[iVar][jVar] = 0.0;

  /*--- Rename for convenience ---*/
  su2double rho    = val_V[RHO_INDEX];
  su2double H      = val_V[H_INDEX];
  su2double rhoEve = val_U[nSpecies+nDim+1];

  su2double u[MAXNDIM];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    u[iDim] = val_V[VEL_INDEX+iDim];

  /*--- Calculate projected velocity ---*/
  su2double proj_vel = GeometryToolbox::DotProduct(nDim, u, val_normal);

  /*--- Species density rows ---*/
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (unsigned short jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
      val_Proj_Jac_Tensor[iSpecies][jSpecies] += -(rhos[iSpecies]/rho) * proj_vel;
    }
    val_Proj_Jac_Tensor[iSpecies][iSpecies]   += proj_vel;
    for (unsigned short iDim  = 0; iDim < nDim; iDim++) {
      val_Proj_Jac_Tensor[iSpecies][nSpecies+iDim] += (rhos[iSpecies]/rho) * val_normal[iDim];
      val_Proj_Jac_Tensor[nSpecies+iDim][iSpecies] += val_dPdU[iSpecies]*val_normal[iDim] - proj_vel*u[iDim];
    }
    val_Proj_Jac_Tensor[nSpecies+nDim][iSpecies]   += (val_dPdU[iSpecies]-H) * proj_vel;
    val_Proj_Jac_Tensor[nSpecies+nDim+1][iSpecies] += -proj_vel * rhoEve/rho;
  }

  /*--- Momentum rows ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+jDim] += val_dPdU[nSpecies+jDim]*val_normal[iDim] + u[iDim]*val_normal[jDim];
    }
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+iDim]   += proj_vel;
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+nDim]   += val_dPdU[nSpecies+nDim]*val_normal[iDim];
    val_Proj_Jac_Tensor[nSpecies+iDim][nSpecies+nDim+1] += val_dPdU[nSpecies+nDim+1]*val_normal[iDim];
  }

  /*--- Total energy row ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+iDim] += val_dPdU[nSpecies+iDim]*proj_vel + H*val_normal[iDim];
  val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+nDim]   += (1+val_dPdU[nSpecies+nDim])*proj_vel;
  val_Proj_Jac_Tensor[nSpecies+nDim][nSpecies+nDim+1] +=  val_dPdU[nSpecies+nDim+1] *proj_vel;

  /*--- Vib.-el. energy row ---*/
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    val_Proj_Jac_Tensor[nSpecies+nDim+1][nSpecies+iDim] = rhoEve/rho*val_normal[iDim];
  val_Proj_Jac_Tensor[nSpecies+nDim+1][nSpecies+nDim+1] = proj_vel;

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    for (unsigned short jVar = 0; jVar < nVar; jVar++)
      val_Proj_Jac_Tensor[iVar][jVar] = val_scale * val_Proj_Jac_Tensor[iVar][jVar];
}

void CNEMONumerics::GetViscousProjFlux(su2double *val_primvar,
                                       su2double **val_gradprimvar,
                                       su2double *val_eve,
                                       const su2double *val_normal,
                                       su2double *val_diffusioncoeff,
                                       su2double val_lam_viscosity,
                                       su2double val_eddy_viscosity,
                                       su2double val_therm_conductivity,
                                       su2double val_therm_conductivity_ve,
                                       const CConfig *config) {

  // Requires a slightly non-standard primitive vector:
  // Assumes -     V = [Y1, ... , Yn, T, Tve, ... ]
  // and gradient GV = [GY1, ... , GYn, GT, GTve, ... ]
  // rather than the standard V = [r1, ... , rn, T, Tve, ... ]

  unsigned short iSpecies, iVar, iDim, jDim;
  su2double *Ds, *V, **GV, mu, ktr, kve;
  su2double rho, T, Tve, RuSI, Ru;
  auto& Ms = fluidmodel->GetSpeciesMolarMass();

  su2activematrix Flux_Tensor(nVar,nDim);

  /*--- Initialize ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Flux_Tensor[iVar][iDim] = 0.0;
  }

  /*--- Rename for convenience ---*/
  Ds  = val_diffusioncoeff;
  mu  = val_lam_viscosity+val_eddy_viscosity;
  ktr = val_therm_conductivity;
  kve = val_therm_conductivity_ve;
  rho = val_primvar[RHO_INDEX];
  T   = val_primvar[T_INDEX];
  Tve = val_primvar[TVE_INDEX];
  V   = val_primvar;
  GV  = val_gradprimvar;
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;

  hs = fluidmodel->ComputeSpeciesEnthalpy(T, Tve, val_eve);

  /*--- Scale thermal conductivity with turb visc ---*/
  // TODO: Need to determine proper way to incorporate eddy viscosity
  // This is only scaling Kve by same factor as ktr
  // NOTE: V[iSpecies] is == Ys.
  su2double Mass = 0.0;
  su2double tmp1, scl, Cptr;
  for (iSpecies=0;iSpecies<nSpecies;iSpecies++)
    Mass += V[iSpecies]*Ms[iSpecies];
  Cptr = V[RHOCVTR_INDEX]/V[RHO_INDEX]+Ru/Mass;
  tmp1 = Cptr*(val_eddy_viscosity/Prandtl_Turb);
  scl  = tmp1/ktr;
  ktr += Cptr*(val_eddy_viscosity/Prandtl_Turb);
  kve  = kve*(1.0+scl);
  //Cpve = V[RHOCVVE_INDEX]+Ru/Mass;
  //kve += Cpve*(val_eddy_viscosity/Prandtl_Turb);

  /*--- Pre-compute mixture quantities ---*/

  su2double Vector[MAXNDIM] = {0.0};

  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Vector[iDim] += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim];
    }
  }

  /*--- Compute the viscous stress tensor ---*/
  ComputeStressTensor(nDim,tau,val_gradprimvar+VEL_INDEX, mu);

  /*--- Populate entries in the viscous flux vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    /*--- Species diffusion velocity ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[iSpecies][iDim] = rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim]
          - V[RHOS_INDEX+iSpecies]*Vector[iDim];
    }
    if (ionization) {
      SU2_MPI::Error("NEED TO IMPLEMENT IONIZED FUNCTIONALITY!!!",CURRENT_FUNCTION);
    }

    /*--- Shear stress related terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim] = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      Flux_Tensor[nSpecies+jDim][iDim]  = tau[iDim][jDim];
      Flux_Tensor[nSpecies+nDim][iDim] += tau[iDim][jDim]*val_primvar[VEL_INDEX+jDim];
    }

    /*--- Diffusion terms ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[nSpecies+nDim][iDim]   += Flux_Tensor[iSpecies][iDim] * hs[iSpecies];
      Flux_Tensor[nSpecies+nDim+1][iDim] += Flux_Tensor[iSpecies][iDim] * val_eve[iSpecies];
    }

    /*--- Heat transfer terms ---*/
    Flux_Tensor[nSpecies+nDim][iDim]   += ktr*GV[T_INDEX][iDim] +
        kve*GV[TVE_INDEX][iDim];
    Flux_Tensor[nSpecies+nDim+1][iDim] += kve*GV[TVE_INDEX][iDim];
  }

  for (iVar = 0; iVar < nVar; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Proj_Flux_Tensor[iVar] += Flux_Tensor[iVar][iDim]*val_normal[iDim];
    }
  }
}

void CNEMONumerics::GetViscousProjJacs(su2double *val_Mean_PrimVar,
                                       su2double **val_Mean_GradPrimVar,
                                       su2double *val_Mean_Eve,
                                       su2double *val_Mean_Cvve,
                                       su2double *val_diffusion_coeff,
                                       su2double val_laminar_viscosity,
                                       su2double val_eddy_viscosity,
                                       su2double val_thermal_conductivity,
                                       su2double val_thermal_conductivity_ve,
                                       su2double val_dist_ij, su2double *val_normal,
                                       su2double val_dS, su2double *val_Fv,
                                       su2double **val_Jac_i, su2double **val_Jac_j,
                                       const CConfig *config) {


  //TODO UPDATE WITH EDDY VISC
//  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar;
//  su2double rho, rho_i, rho_j, vel[3], T, Tve;
//  su2double mu, ktr, kve, *Ds, dij, Ru, RuSI;
//  su2double theta, thetax, thetay, thetaz;
//  su2double etax, etay, etaz;
//  su2double pix, piy, piz;
//  su2double sumY, sumY_i, sumY_j;
//
//  /*--- Initialize arrays ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      dFdVi[iVar][jVar] = 0.0;
//      dFdVj[iVar][jVar] = 0.0;
//      dVdUi[iVar][jVar] = 0.0;
//      dVdUj[iVar][jVar] = 0.0;
//    }
//  }
//
//  /*--- Initialize the Jacobian matrices ---*/
//  for (iVar = 0; iVar < nVar; iVar++) {
//    for (jVar = 0; jVar < nVar; jVar++) {
//      val_Jac_i[iVar][jVar] = 0.0;
//      val_Jac_j[iVar][jVar] = 0.0;
//    }
//  }
//
//  /*--- Initialize storage vectors & matrices ---*/
//  for (iVar = 0; iVar < nSpecies; iVar++) {
//    sumdFdYjh[iVar]   = 0.0;
//    sumdFdYjeve[iVar] = 0.0;
//    for (jVar = 0; jVar < nSpecies; jVar++) {
//      dFdYi[iVar][jVar] = 0.0;
//      dFdYj[iVar][jVar] = 0.0;
//      dJdr_i[iVar][jVar] = 0.0;
//      dJdr_j[iVar][jVar] = 0.0;
//    }
//  }
//
//
//  /*--- Calculate preliminary geometrical quantities ---*/
//  dij = val_dist_ij;
//  theta = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    theta += val_normal[iDim]*val_normal[iDim];
//  }
//
//
//  /*--- Rename for convenience ---*/
//  rho = val_Mean_PrimVar[RHO_INDEX];
//  rho_i = V_i[RHO_INDEX];
//  rho_j = V_j[RHO_INDEX];
//  T   = val_Mean_PrimVar[T_INDEX];
//  Tve = val_Mean_PrimVar[TVE_INDEX];
//  Ds  = val_diffusion_coeff;
//  mu  = val_laminar_viscosity;
//  ktr = val_thermal_conductivity;
//  kve = val_thermal_conductivity_ve;
//  RuSI= UNIVERSAL_GAS_CONSTANT;
//  Ru  = 1000.0*RuSI;
//
//  hs = fluidmodel->ComputeSpeciesEnthalpy(T, val_Mean_Eve);
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Ys[iSpecies]   = val_Mean_PrimVar[RHOS_INDEX+iSpecies];
//    Ys_i[iSpecies] = V_i[RHOS_INDEX+iSpecies]/V_i[RHO_INDEX];
//    Ys_j[iSpecies] = V_j[RHOS_INDEX+iSpecies]/V_j[RHO_INDEX];
//    Cvtr[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
//  }
//  for (iDim = 0; iDim < nDim; iDim++)
//    vel[iDim] = val_Mean_PrimVar[VEL_INDEX+iDim];
//
//  /*--- Calculate useful diffusion parameters ---*/
//  // Summation term of the diffusion fluxes
//  sumY = 0.0;
//  sumY_i = 0.0;
//  sumY_j = 0.0;
//  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
//    sumY_i += Ds[iSpecies]*theta/dij*Ys_i[iSpecies];
//    sumY_j += Ds[iSpecies]*theta/dij*Ys_j[iSpecies];
//    sumY   += Ds[iSpecies]*theta/dij*(Ys_j[iSpecies]-Ys_i[iSpecies]);
//  }
//
//
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    for (jSpecies  = 0; jSpecies < nSpecies; jSpecies++) {
//
//      // first term
//      dJdr_j[iSpecies][jSpecies] +=  0.5*(Ds[iSpecies]*theta/dij *
//                                          (Ys_j[iSpecies]*rho_i/rho_j +
//                                           Ys_i[iSpecies]));
//      dJdr_i[iSpecies][jSpecies] += -0.5*(Ds[iSpecies]*theta/dij *
//                                          (Ys_j[iSpecies] +
//                                           Ys_i[iSpecies]*rho_j/rho_i));
//
//      // second term
//      dJdr_j[iSpecies][jSpecies] +=
//          0.25*(Ys_i[iSpecies] - rho_i/rho_j*Ys_j[iSpecies])*sumY
//          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_j)
//          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_j/rho_j;
//
//      dJdr_i[iSpecies][jSpecies] +=
//          0.25*(-rho_j/rho_i*Ys_i[iSpecies]+Ys_j[iSpecies])*sumY
//          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_i)
//          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_i/rho_i;
//    }
//
//    // first term
//    dJdr_j[iSpecies][iSpecies] += -0.5*Ds[iSpecies]*theta/dij*(1+rho_i/rho_j);
//    dJdr_i[iSpecies][iSpecies] +=  0.5*Ds[iSpecies]*theta/dij*(1+rho_j/rho_i);
//
//    // second term
//    dJdr_j[iSpecies][iSpecies] += 0.25*(1.0+rho_i/rho_j)*sumY;
//    dJdr_i[iSpecies][iSpecies] += 0.25*(1.0+rho_j/rho_i)*sumY;
//  }
//
//  /*--- Calculate transformation matrix ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    dVdUi[iSpecies][iSpecies] = 1.0;
//    dVdUj[iSpecies][iSpecies] = 1.0;
//  }
//  for (iDim = 0; iDim < nDim; iDim++) {
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      dVdUi[nSpecies+iDim][iSpecies] = -V_i[VEL_INDEX+iDim]/V_i[RHO_INDEX];
//      dVdUj[nSpecies+iDim][iSpecies] = -V_j[VEL_INDEX+iDim]/V_j[RHO_INDEX];
//    }
//    dVdUi[nSpecies+iDim][nSpecies+iDim] = 1.0/V_i[RHO_INDEX];
//    dVdUj[nSpecies+iDim][nSpecies+iDim] = 1.0/V_j[RHO_INDEX];
//  }
//  for (iVar = 0; iVar < nVar; iVar++) {
//    dVdUi[nSpecies+nDim][iVar]   = dTdU_i[iVar];
//    dVdUj[nSpecies+nDim][iVar]   = dTdU_j[iVar];
//    dVdUi[nSpecies+nDim+1][iVar] = dTvedU_i[iVar];
//    dVdUj[nSpecies+nDim+1][iVar] = dTvedU_j[iVar];
//  }
//
//
//  if (nDim == 2) {
//
//    /*--- Geometry parameters ---*/
//    thetax = theta + val_normal[0]*val_normal[0]/3.0;
//    thetay = theta + val_normal[1]*val_normal[1]/3.0;
//    etaz   = val_normal[0]*val_normal[1]/3.0;
//    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]  );
//    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1]);
//
//    /*--- Populate primitive Jacobian ---*/
//
//    // X-momentum
//    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
//    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;
//
//    // Y-momentum
//    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
//    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;
//
//    // Energy
//    dFdVj[nSpecies+2][nSpecies]   = pix*val_dS;
//    dFdVj[nSpecies+2][nSpecies+1] = piy*val_dS;
//    dFdVj[nSpecies+2][nSpecies+2] = ktr*theta/dij*val_dS;
//    dFdVj[nSpecies+2][nSpecies+3] = kve*theta/dij*val_dS;
//
//    // Vib-el Energy
//    dFdVj[nSpecies+3][nSpecies+3] = kve*theta/dij*val_dS;
//
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++)
//        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];
//
//    // Common terms
//    dFdVi[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
//    dFdVj[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
//    dFdVi[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
//    dFdVj[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      dFdVi[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
//                                                             Cvtr[iSpecies]   );
//      dFdVj[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
//                                                             Cvtr[iSpecies]   );
//      dFdVi[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVj[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//    }
//
//    // Unique terms
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
//        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
//        dFdVj[nSpecies+2][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
//        dFdVi[nSpecies+2][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
//        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
//        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
//      }
//    }
//
//  } //nDim == 2
//  else {
//
//    /*--- Geometry parameters ---*/
//    thetax = theta + val_normal[0]*val_normal[0]/3.0;
//    thetay = theta + val_normal[1]*val_normal[1]/3.0;
//    thetaz = theta + val_normal[2]*val_normal[2]/3.0;
//    etax   = val_normal[1]*val_normal[2]/3.0;
//    etay   = val_normal[0]*val_normal[2]/3.0;
//    etaz   = val_normal[0]*val_normal[1]/3.0;
//    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]   + etay*vel[2]  );
//    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1] + etax*vel[2]  );
//    piz    = mu/dij * (etay*vel[0]   + etax*vel[1]   + thetaz*vel[2]);
//
//    /*--- Populate primitive Jacobian ---*/
//
//    // X-momentum
//    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
//    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;
//    dFdVj[nSpecies][nSpecies+2]   = mu*etay/dij*val_dS;
//
//    // Y-momentum
//    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
//    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;
//    dFdVj[nSpecies+1][nSpecies+2] = mu*etax/dij*val_dS;
//
//    // Z-momentum
//    dFdVj[nSpecies+2][nSpecies]   = mu*etay/dij*val_dS;
//    dFdVj[nSpecies+2][nSpecies+1] = mu*etax/dij*val_dS;
//    dFdVj[nSpecies+2][nSpecies+2] = mu*thetaz/dij*val_dS;
//
//    // Energy
//    dFdVj[nSpecies+3][nSpecies]   = pix*val_dS;
//    dFdVj[nSpecies+3][nSpecies+1] = piy*val_dS;
//    dFdVj[nSpecies+3][nSpecies+2] = piz*val_dS;
//    dFdVj[nSpecies+3][nSpecies+3] = ktr*theta/dij*val_dS;
//    dFdVj[nSpecies+3][nSpecies+4] = kve*theta/dij*val_dS;
//
//    // Vib.-el energy
//    dFdVj[nSpecies+4][nSpecies+4] = kve*theta/dij*val_dS;
//
//    for (iVar = 0; iVar < nVar; iVar++)
//      for (jVar = 0; jVar < nVar; jVar++)
//        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];
//
//    // Common terms
//    for (iDim = 0; iDim < nDim; iDim++) {
//      dFdVi[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
//      dFdVj[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
//    }
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
//                                                             Cvtr[iSpecies]   );
//      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
//                                                             Cvtr[iSpecies]   );
//      dFdVi[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVj[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVi[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//      dFdVj[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
//    }
//
//    // Unique terms
//    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
//        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
//        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
//        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
//        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
//        dFdVj[nSpecies+4][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
//        dFdVi[nSpecies+4][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
//      }
//    }
//
//  } // nDim == 3
//
//  /*--- dFv/dUij = dFv/dVij * dVij/dUij ---*/
//  for (iVar = 0; iVar < nVar; iVar++)
//    for (jVar = 0; jVar < nVar; jVar++)
//      for (kVar = 0; kVar < nVar; kVar++) {
//        val_Jac_i[iVar][jVar] += dFdVi[iVar][kVar]*dVdUi[kVar][jVar];
//        val_Jac_j[iVar][jVar] += dFdVj[iVar][kVar]*dVdUj[kVar][jVar];
//      }
}

void CNEMONumerics::GetPMatrix(const su2double *U, const su2double *V, const su2double *val_dPdU,
                               const su2double *val_normal, const su2double *l, const su2double *m,
                               su2double **val_p_tensor) const {

  // P matrix is equivalent to the L matrix in Gnoffo
  unsigned short iSpecies, iDim, iVar, jVar;
  su2double sqvel, rho, a, a2, eve;
  su2double vU, vV, vW;

  /*--- Initialize the P matrix to zero ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      val_p_tensor[iVar][jVar] = 0.0;

  /*--- Pre-compute useful quantities ---*/
  sqvel = 0.0;
  rho = V[RHO_INDEX];
  eve = U[nSpecies+nDim+1]/rho;
  vU = 0.0;  vV = 0.0;  vW = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vU    += V[VEL_INDEX+iDim] * val_normal[iDim];
    vV    += V[VEL_INDEX+iDim] * l[iDim];
    vW    += V[VEL_INDEX+iDim] * m[iDim];
    sqvel += V[VEL_INDEX+iDim] * V[VEL_INDEX+iDim];
  }
  a  = V[A_INDEX];
  a2 = V[A_INDEX]*V[A_INDEX];

  if(nDim == 2) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
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

    for (iDim = 0; iDim < nDim; iDim++){
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

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
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

    for (iDim = 0; iDim < nDim; iDim++){
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

  unsigned short iSpecies, jSpecies, iDim, iVar, jVar;
  su2double rho, a, a2, eve;
  su2double vU, vV, vW;

  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      val_invp_tensor[iVar][jVar] = 0.0;

  /*--- Pre-compute useful quantities ---*/
  rho = V[RHO_INDEX];
  eve = U[nSpecies+nDim+1]/rho;
  vU = 0.0;  vV = 0.0;  vW = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    vU += V[VEL_INDEX+iDim] * val_normal[iDim];
    vV += V[VEL_INDEX+iDim] * l[iDim];
    vW += V[VEL_INDEX+iDim] * m[iDim];
  }
  a  = V[A_INDEX];
  a2 = V[A_INDEX]*V[A_INDEX];

  if (nDim == 2) {

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
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

    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
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
