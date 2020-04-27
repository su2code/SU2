/*!
 * \file NEMO_slau.hpp
 * \brief Declaration of numerics classes for the NEMO family of schemes,
 *        including SLAU. The implementation is in NEMO.cpp.
 * \author F. Palacios, T. Economon
 * \version 7.0.3 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/NEMO/NEMO_diffusion.hpp"

CAvgGrad_NEMO::CAvgGrad_NEMO(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             CConfig *config) : CNumerics(val_nDim,
                                                          val_nVar,
                                                          config) {

  implicit = (config->GetKind_TimeIntScheme_NEMO() == EULER_IMPLICIT);

  /*--- Rename for convenience ---*/
  nDim         = val_nDim;
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;

  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
  PrimVar_i    = new su2double [nPrimVar];
  PrimVar_j    = new su2double [nPrimVar];
  Mean_PrimVar = new su2double [nPrimVar];

  Mean_U      = new su2double[nVar];
  Mean_dPdU   = new su2double[nVar];
  Mean_dTdU   = new su2double[nVar];
  Mean_dTvedU = new su2double[nVar];
  Mean_Eve    = new su2double[nSpecies];
  Mean_Cvve   = new su2double[nSpecies];
  Mean_GU     = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GU[iVar] = new su2double[nDim];

  Mean_Diffusion_Coeff = new su2double[nSpecies];

  /*--- Compressible flow, primitive gradient variables nDim+3, (T,vx,vy,vz) ---*/
  Mean_GradPrimVar = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];
}

CAvgGrad_NEMO::~CAvgGrad_NEMO(void) {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;
  delete [] Mean_Diffusion_Coeff;

  delete [] Mean_U;
  delete [] Mean_dPdU;
  delete [] Mean_dTdU;
  delete [] Mean_dTvedU;
  delete [] Mean_Eve;
  delete [] Mean_Cvve;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GU[iVar];
  delete [] Mean_GU;

  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;
}

void CAvgGrad_NEMO::ComputeResidual(su2double *val_residual,
                                    su2double **val_Jacobian_i,
                                    su2double **val_Jacobian_j,
                                    CConfig *config) {

  unsigned short iSpecies, iVar, iDim;

  /*--- Normalized normal vector ---*/
  Area = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean transport coefficients ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coeff[iSpecies] = 0.5*(Diffusion_Coeff_i[iSpecies] +
                                          Diffusion_Coeff_j[iSpecies]);
  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i +
                                Laminar_Viscosity_j);
  Mean_Thermal_Conductivity = 0.5*(Thermal_Conductivity_i +
                                   Thermal_Conductivity_j);
  Mean_Thermal_Conductivity_ve = 0.5*(Thermal_Conductivity_ve_i +
                                      Thermal_Conductivity_ve_j);

  /*--- Mean gradient approximation ---*/
  // Mass fraction
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    PrimVar_i[iSpecies] = V_i[iSpecies]/V_i[RHO_INDEX];
    PrimVar_j[iSpecies] = V_j[iSpecies]/V_j[RHO_INDEX];
    Mean_PrimVar[iSpecies] = 0.5*(PrimVar_i[iSpecies] + PrimVar_j[iSpecies]);
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iSpecies][iDim] = 0.5*(1.0/V_i[RHO_INDEX] *
                                              (PrimVar_Grad_i[iSpecies][iDim] -
                                               PrimVar_i[iSpecies] *
                                               PrimVar_Grad_i[RHO_INDEX][iDim]) +
                                              1.0/V_j[RHO_INDEX] *
                                              (PrimVar_Grad_j[iSpecies][iDim] -
                                               PrimVar_j[iSpecies] *
                                               PrimVar_Grad_j[RHO_INDEX][iDim]));

    }
  }

  for (iVar = nSpecies; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  for (iVar = nSpecies; iVar < nPrimVarGrad; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] +
                                          PrimVar_Grad_j[iVar][iDim]);
    }
  }
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Mean_Eve[iSpecies]  = 0.5*(eve_i[iSpecies]  + eve_j[iSpecies]);
    Mean_Cvve[iSpecies] = 0.5*(Cvve_i[iSpecies] + Cvve_j[iSpecies]);
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve, Normal,
                     Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,
                     Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                     config);

  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);

    GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve, Mean_Cvve,
                       Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,
                       Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                       dist_ij, UnitNormal, Area, Proj_Flux_Tensor,
                       val_Jacobian_i, val_Jacobian_j, config);
  }
}


void CAvgGrad_NEMO::GetViscousProjFlux(su2double *val_primvar,
                                       su2double **val_gradprimvar,
                                       su2double *val_eve,
                                       su2double *val_normal,
                                       su2double *val_diffusioncoeff,
                                       su2double val_viscosity,
                                       su2double val_therm_conductivity,
                                       su2double val_therm_conductivity_ve,
                                       CConfig *config) {

  // Requires a slightly non-standard primitive vector:
  // Assumes -     V = [Y1, ... , Yn, T, Tve, ... ]
  // and gradient GV = [GY1, ... , GYn, GT, GTve, ... ]
  // rather than the standard V = [r1, ... , rn, T, Tve, ... ]

  bool ionization;
  unsigned short iSpecies, iVar, iDim, jDim, nHeavy, nEl;
  su2double *Ds, *V, **GV, mu, ktr, kve, div_vel;
  su2double Ru, RuSI;
  su2double rho, T, Tve;

  /*--- Initialize ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Flux_Tensor[iVar][iDim] = 0.0;
  }

  /*--- Read from CConfig ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Rename for convenience ---*/
  Ds  = val_diffusioncoeff;
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  kve = val_therm_conductivity_ve;
  rho = val_primvar[RHO_INDEX];
  T   = val_primvar[T_INDEX];
  Tve = val_primvar[TVE_INDEX];
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;
  V   = val_primvar;
  GV  = val_gradprimvar;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    hs[iSpecies]  = variable->CalcHs(config, T, val_eve[iSpecies], iSpecies);

  /*--- Calculate the velocity divergence ---*/
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += GV[VEL_INDEX+iDim][iDim];

  /*--- Pre-compute mixture quantities ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Vector[iDim] = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Vector[iDim] += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim];
    }
  }

  /*--- Compute the viscous stress tensor ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      tau[iDim][jDim] = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] += mu * (val_gradprimvar[VEL_INDEX+jDim][iDim] +
          val_gradprimvar[VEL_INDEX+iDim][jDim]);
    }
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    /*--- Species diffusion velocity ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[iSpecies][iDim] = rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim]
          - V[RHOS_INDEX+iSpecies]*Vector[iDim];
    }
    if (ionization) {
      cout << "GetViscProjFlux -- NEED TO IMPLEMENT IONIZED FUNCTIONALITY!!!" << endl;
      exit(1);
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

void CAvgGrad_NEMO::GetViscousProjJacs(su2double *val_Mean_PrimVar,
                                       su2double **val_Mean_GradPrimVar,
                                       su2double *val_Mean_Eve,
                                       su2double *val_Mean_Cvve,
                                       su2double *val_diffusion_coeff,
                                       su2double val_laminar_viscosity,
                                       su2double val_thermal_conductivity,
                                       su2double val_thermal_conductivity_ve,
                                       su2double val_dist_ij, su2double *val_normal,
                                       su2double val_dS, su2double *val_Fv,
                                       su2double **val_Jac_i, su2double **val_Jac_j,
                                       CConfig *config) {

  bool ionization;
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar, nHeavy, nEl;
  su2double rho, rho_i, rho_j, vel[3], T, Tve, *xi, *Ms;
  su2double mu, ktr, kve, *Ds, dij, Ru, RuSI;
  su2double theta, thetax, thetay, thetaz;
  su2double etax, etay, etaz;
  su2double pix, piy, piz;
  su2double sumY, sumY_i, sumY_j;

  /*--- Initialize arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      dFdVi[iVar][jVar] = 0.0;
      dFdVj[iVar][jVar] = 0.0;
      dVdUi[iVar][jVar] = 0.0;
      dVdUj[iVar][jVar] = 0.0;
    }
  }

  /*--- Initialize the Jacobian matrices ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jac_i[iVar][jVar] = 0.0;
      val_Jac_j[iVar][jVar] = 0.0;
    }
  }

  /*--- Initialize storage vectors & matrices ---*/
  for (iVar = 0; iVar < nSpecies; iVar++) {
    sumdFdYjh[iVar]   = 0.0;
    sumdFdYjeve[iVar] = 0.0;
    for (jVar = 0; jVar < nSpecies; jVar++) {
      dFdYi[iVar][jVar] = 0.0;
      dFdYj[iVar][jVar] = 0.0;
      dJdr_i[iVar][jVar] = 0.0;
      dJdr_j[iVar][jVar] = 0.0;
    }
  }

  /*--- Assign booleans from CConfig ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Calculate preliminary geometrical quantities ---*/
  dij = val_dist_ij;
  theta = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    theta += val_normal[iDim]*val_normal[iDim];
  }


  /*--- Rename for convenience ---*/
  rho = val_Mean_PrimVar[RHO_INDEX];
  rho_i = V_i[RHO_INDEX];
  rho_j = V_j[RHO_INDEX];
  T   = val_Mean_PrimVar[T_INDEX];
  Tve = val_Mean_PrimVar[TVE_INDEX];
  Ds  = val_diffusion_coeff;
  mu  = val_laminar_viscosity;
  ktr = val_thermal_conductivity;
  kve = val_thermal_conductivity_ve;
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;
  Ms  = config->GetMolar_Mass();
  xi  = config->GetRotationModes();
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys[iSpecies]   = val_Mean_PrimVar[RHOS_INDEX+iSpecies];
    Ys_i[iSpecies] = V_i[RHOS_INDEX+iSpecies]/V_i[RHO_INDEX];
    Ys_j[iSpecies] = V_j[RHOS_INDEX+iSpecies]/V_j[RHO_INDEX];
    hs[iSpecies]   = variable->CalcHs(config, T, val_Mean_Eve[iSpecies], iSpecies);
    Cvtr[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++)
    vel[iDim] = val_Mean_PrimVar[VEL_INDEX+iDim];

  /*--- Calculate useful diffusion parameters ---*/
  // Summation term of the diffusion fluxes
  sumY = 0.0;
  sumY_i = 0.0;
  sumY_j = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    sumY_i += Ds[iSpecies]*theta/dij*Ys_i[iSpecies];
    sumY_j += Ds[iSpecies]*theta/dij*Ys_j[iSpecies];
    sumY   += Ds[iSpecies]*theta/dij*(Ys_j[iSpecies]-Ys_i[iSpecies]);
  }


  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies  = 0; jSpecies < nSpecies; jSpecies++) {

      // first term
      dJdr_j[iSpecies][jSpecies] +=  0.5*(Ds[iSpecies]*theta/dij *
                                          (Ys_j[iSpecies]*rho_i/rho_j +
                                           Ys_i[iSpecies]));
      dJdr_i[iSpecies][jSpecies] += -0.5*(Ds[iSpecies]*theta/dij *
                                          (Ys_j[iSpecies] +
                                           Ys_i[iSpecies]*rho_j/rho_i));

      // second term
      dJdr_j[iSpecies][jSpecies] +=
          0.25*(Ys_i[iSpecies] - rho_i/rho_j*Ys_j[iSpecies])*sumY
          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_j)
          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_j/rho_j;

      dJdr_i[iSpecies][jSpecies] +=
          0.25*(-rho_j/rho_i*Ys_i[iSpecies]+Ys_j[iSpecies])*sumY
          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_i)
          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_i/rho_i;
    }

    // first term
    dJdr_j[iSpecies][iSpecies] += -0.5*Ds[iSpecies]*theta/dij*(1+rho_i/rho_j);
    dJdr_i[iSpecies][iSpecies] +=  0.5*Ds[iSpecies]*theta/dij*(1+rho_j/rho_i);

    // second term
    dJdr_j[iSpecies][iSpecies] += 0.25*(1.0+rho_i/rho_j)*sumY;
    dJdr_i[iSpecies][iSpecies] += 0.25*(1.0+rho_j/rho_i)*sumY;
  }

  /*--- Calculate transformation matrix ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dVdUi[iSpecies][iSpecies] = 1.0;
    dVdUj[iSpecies][iSpecies] = 1.0;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dVdUi[nSpecies+iDim][iSpecies] = -V_i[VEL_INDEX+iDim]/V_i[RHO_INDEX];
      dVdUj[nSpecies+iDim][iSpecies] = -V_j[VEL_INDEX+iDim]/V_j[RHO_INDEX];
    }
    dVdUi[nSpecies+iDim][nSpecies+iDim] = 1.0/V_i[RHO_INDEX];
    dVdUj[nSpecies+iDim][nSpecies+iDim] = 1.0/V_j[RHO_INDEX];
  }
  for (iVar = 0; iVar < nVar; iVar++) {
    dVdUi[nSpecies+nDim][iVar]   = dTdU_i[iVar];
    dVdUj[nSpecies+nDim][iVar]   = dTdU_j[iVar];
    dVdUi[nSpecies+nDim+1][iVar] = dTvedU_i[iVar];
    dVdUj[nSpecies+nDim+1][iVar] = dTvedU_j[iVar];
  }


  if (nDim == 2) {

    /*--- Geometry parameters ---*/
    thetax = theta + val_normal[0]*val_normal[0]/3.0;
    thetay = theta + val_normal[1]*val_normal[1]/3.0;
    etaz   = val_normal[0]*val_normal[1]/3.0;
    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]  );
    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1]);

    /*--- Populate primitive Jacobian ---*/

    // X-momentum
    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;

    // Y-momentum
    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;

    // Energy
    dFdVj[nSpecies+2][nSpecies]   = pix*val_dS;
    dFdVj[nSpecies+2][nSpecies+1] = piy*val_dS;
    dFdVj[nSpecies+2][nSpecies+2] = ktr*theta/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+3] = kve*theta/dij*val_dS;

    // Vib-el Energy
    dFdVj[nSpecies+3][nSpecies+3] = kve*theta/dij*val_dS;

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

    // Common terms
    dFdVi[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
    dFdVj[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
    dFdVi[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
    dFdVj[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dFdVi[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVj[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVi[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
    }

    // Unique terms
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
        dFdVj[nSpecies+2][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVi[nSpecies+2][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
      }
    }

  } //nDim == 2
  else {

    /*--- Geometry parameters ---*/
    thetax = theta + val_normal[0]*val_normal[0]/3.0;
    thetay = theta + val_normal[1]*val_normal[1]/3.0;
    thetaz = theta + val_normal[2]*val_normal[2]/3.0;
    etax   = val_normal[1]*val_normal[2]/3.0;
    etay   = val_normal[0]*val_normal[2]/3.0;
    etaz   = val_normal[0]*val_normal[1]/3.0;
    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]   + etay*vel[2]  );
    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1] + etax*vel[2]  );
    piz    = mu/dij * (etay*vel[0]   + etax*vel[1]   + thetaz*vel[2]);

    /*--- Populate primitive Jacobian ---*/

    // X-momentum
    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies][nSpecies+2]   = mu*etay/dij*val_dS;

    // Y-momentum
    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+2] = mu*etax/dij*val_dS;

    // Z-momentum
    dFdVj[nSpecies+2][nSpecies]   = mu*etay/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+1] = mu*etax/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+2] = mu*thetaz/dij*val_dS;

    // Energy
    dFdVj[nSpecies+3][nSpecies]   = pix*val_dS;
    dFdVj[nSpecies+3][nSpecies+1] = piy*val_dS;
    dFdVj[nSpecies+3][nSpecies+2] = piz*val_dS;
    dFdVj[nSpecies+3][nSpecies+3] = ktr*theta/dij*val_dS;
    dFdVj[nSpecies+3][nSpecies+4] = kve*theta/dij*val_dS;

    // Vib.-el energy
    dFdVj[nSpecies+4][nSpecies+4] = kve*theta/dij*val_dS;

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

    // Common terms
    for (iDim = 0; iDim < nDim; iDim++) {
      dFdVi[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
      dFdVj[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVi[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVi[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
    }

    // Unique terms
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVj[nSpecies+4][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        dFdVi[nSpecies+4][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
      }
    }

  } // nDim == 3

  /*--- dFv/dUij = dFv/dVij * dVij/dUij ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      for (kVar = 0; kVar < nVar; kVar++) {
        val_Jac_i[iVar][jVar] += dFdVi[iVar][kVar]*dVdUi[kVar][jVar];
        val_Jac_j[iVar][jVar] += dFdVj[iVar][kVar]*dVdUj[kVar][jVar];
      }
}

CAvgGradCorrected_NEMO::CAvgGradCorrected_NEMO(unsigned short val_nDim,
                                               unsigned short val_nVar,
                                               unsigned short val_nPrimVar,
                                               unsigned short val_nPrimVarGrad,
                                               CConfig *config) : CNumerics(val_nDim,
                                                                            val_nVar,
                                                                            config) {

  implicit = (config->GetKind_TimeIntScheme_NEMO() == EULER_IMPLICIT);

  /*--- Rename for convenience ---*/
  nDim         = val_nDim;
  nVar         = val_nVar;
  nPrimVar     = val_nPrimVar;
  nPrimVarGrad = val_nPrimVarGrad;

  /*--- Compressible flow, primitive variables nDim+3, (T,vx,vy,vz,P,rho) ---*/
  PrimVar_i    = new su2double [nPrimVar];
  PrimVar_j    = new su2double [nPrimVar];
  Mean_PrimVar = new su2double [nPrimVar];

  Mean_Eve  = new su2double[nSpecies];
  Mean_Cvve = new su2double[nSpecies];

  Mean_Diffusion_Coeff = new su2double[nSpecies];

  /*--- Compressible flow, primitive gradient variables nDim+3, (T,vx,vy,vz) ---*/
  Mean_GradPrimVar = new su2double* [nPrimVarGrad];
  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    Mean_GradPrimVar[iVar] = new su2double [nDim];

  Proj_Mean_GradPrimVar_Edge = new su2double[nPrimVarGrad];
  Edge_Vector = new su2double[3];
}

CAvgGradCorrected_NEMO::~CAvgGradCorrected_NEMO(void) {

  delete [] PrimVar_i;
  delete [] PrimVar_j;
  delete [] Mean_PrimVar;

  delete [] Mean_Eve;
  delete [] Mean_Cvve;

  delete [] Mean_Diffusion_Coeff;

  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
    delete [] Mean_GradPrimVar[iVar];
  delete [] Mean_GradPrimVar;

  delete [] Proj_Mean_GradPrimVar_Edge;
  delete [] Edge_Vector;
}

void CAvgGradCorrected_NEMO::GetViscousProjFlux(su2double *val_primvar,
                                                su2double **val_gradprimvar,
                                                su2double *val_eve,
                                                su2double *val_normal,
                                                su2double *val_diffusioncoeff,
                                                su2double val_viscosity,
                                                su2double val_therm_conductivity,
                                                su2double val_therm_conductivity_ve,
                                                CConfig *config) {

  // Requires a slightly non-standard primitive vector:
  // Assumes -     V = [Y1, ... , Yn, T, Tve, ... ]
  // and gradient GV = [GY1, ... , GYn, GT, GTve, ... ]
  // rather than the standard V = [r1, ... , rn, T, Tve, ... ]

  bool ionization;
  unsigned short iSpecies, iVar, iDim, jDim, nHeavy, nEl;
  su2double *Ds, *V, **GV, mu, ktr, kve, div_vel;
  su2double Ru, RuSI;
  su2double rho, T, Tve;

  /*--- Initialize ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Flux_Tensor[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Flux_Tensor[iVar][iDim] = 0.0;
  }

  /*--- Read from CConfig ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Rename for convenience ---*/
  Ds  = val_diffusioncoeff;
  mu  = val_viscosity;
  ktr = val_therm_conductivity;
  kve = val_therm_conductivity_ve;
  rho = val_primvar[RHO_INDEX];
  T   = val_primvar[T_INDEX];
  Tve = val_primvar[TVE_INDEX];
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;
  V   = val_primvar;
  GV  = val_gradprimvar;
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    hs[iSpecies]  = variable->CalcHs(config, T, val_eve[iSpecies], iSpecies);

  /*--- Calculate the velocity divergence ---*/
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += GV[VEL_INDEX+iDim][iDim];

  /*--- Pre-compute mixture quantities ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    Vector[iDim] = 0.0;
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Vector[iDim] += rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim];
    }
  }

  /*--- Compute the viscous stress tensor ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      tau[iDim][jDim] = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] += mu * (val_gradprimvar[VEL_INDEX+jDim][iDim] +
          val_gradprimvar[VEL_INDEX+iDim][jDim]);
    }
    tau[iDim][iDim] -= TWO3*mu*div_vel;
  }

  /*--- Populate entries in the viscous flux vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    /*--- Species diffusion velocity ---*/
    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
      Flux_Tensor[iSpecies][iDim] = rho*Ds[iSpecies]*GV[RHOS_INDEX+iSpecies][iDim]
          - V[RHOS_INDEX+iSpecies]*Vector[iDim];
    }
    if (ionization) {
      cout << "GetViscProjFlux -- NEED TO IMPLEMENT IONIZED FUNCTIONALITY!!!" << endl;
      exit(1);
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

void CAvgGradCorrected_NEMO::GetViscousProjJacs(su2double *val_Mean_PrimVar,
                                                su2double **val_Mean_GradPrimVar,
                                                su2double *val_Mean_Eve,
                                                su2double *val_Mean_Cvve,
                                                su2double *val_diffusion_coeff,
                                                su2double val_laminar_viscosity,
                                                su2double val_thermal_conductivity,
                                                su2double val_thermal_conductivity_ve,
                                                su2double val_dist_ij, su2double *val_normal,
                                                su2double val_dS, su2double *val_Fv,
                                                su2double **val_Jac_i, su2double **val_Jac_j,
                                                CConfig *config) {

  bool ionization;
  unsigned short iDim, iSpecies, jSpecies, iVar, jVar, kVar, nHeavy, nEl;
  su2double rho, rho_i, rho_j, vel[3], T, Tve, *xi, *Ms;
  su2double mu, ktr, kve, *Ds, dij, Ru, RuSI;
  su2double theta, thetax, thetay, thetaz;
  su2double etax, etay, etaz;
  su2double pix, piy, piz;
  su2double sumY, sumY_i, sumY_j;

  /*--- Initialize arrays ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      dFdVi[iVar][jVar] = 0.0;
      dFdVj[iVar][jVar] = 0.0;
      dVdUi[iVar][jVar] = 0.0;
      dVdUj[iVar][jVar] = 0.0;
    }
  }

  /*--- Initialize the Jacobian matrices ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    for (jVar = 0; jVar < nVar; jVar++) {
      val_Jac_i[iVar][jVar] = 0.0;
      val_Jac_j[iVar][jVar] = 0.0;
    }
  }

  /*--- Initialize storage vectors & matrices ---*/
  for (iVar = 0; iVar < nSpecies; iVar++) {
    sumdFdYjh[iVar]   = 0.0;
    sumdFdYjeve[iVar] = 0.0;
    for (jVar = 0; jVar < nSpecies; jVar++) {
      dFdYi[iVar][jVar] = 0.0;
      dFdYj[iVar][jVar] = 0.0;
      dJdr_i[iVar][jVar] = 0.0;
      dJdr_j[iVar][jVar] = 0.0;
    }
  }

  /*--- Assign booleans from CConfig ---*/
  ionization = config->GetIonization();
  if (ionization) { nHeavy = nSpecies-1; nEl = 1; }
  else            { nHeavy = nSpecies;   nEl = 0; }

  /*--- Calculate preliminary geometrical quantities ---*/
  dij = val_dist_ij;
  theta = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    theta += val_normal[iDim]*val_normal[iDim];
  }


  /*--- Rename for convenience ---*/
  rho = val_Mean_PrimVar[RHO_INDEX];
  rho_i = V_i[RHO_INDEX];
  rho_j = V_j[RHO_INDEX];
  T   = val_Mean_PrimVar[T_INDEX];
  Tve = val_Mean_PrimVar[TVE_INDEX];
  Ds  = val_diffusion_coeff;
  mu  = val_laminar_viscosity;
  ktr = val_thermal_conductivity;
  kve = val_thermal_conductivity_ve;
  RuSI= UNIVERSAL_GAS_CONSTANT;
  Ru  = 1000.0*RuSI;
  Ms  = config->GetMolar_Mass();
  xi  = config->GetRotationModes();
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Ys[iSpecies]   = val_Mean_PrimVar[RHOS_INDEX+iSpecies];
    Ys_i[iSpecies] = V_i[RHOS_INDEX+iSpecies]/V_i[RHO_INDEX];
    Ys_j[iSpecies] = V_j[RHOS_INDEX+iSpecies]/V_j[RHO_INDEX];
    hs[iSpecies]   = variable->CalcHs(config, T, val_Mean_Eve[iSpecies], iSpecies);
    Cvtr[iSpecies] = (3.0/2.0 + xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
  }
  for (iDim = 0; iDim < nDim; iDim++)
    vel[iDim] = val_Mean_PrimVar[VEL_INDEX+iDim];

  /*--- Calculate useful diffusion parameters ---*/
  // Summation term of the diffusion fluxes
  sumY = 0.0;
  sumY_i = 0.0;
  sumY_j = 0.0;
  for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
    sumY_i += Ds[iSpecies]*theta/dij*Ys_i[iSpecies];
    sumY_j += Ds[iSpecies]*theta/dij*Ys_j[iSpecies];
    sumY   += Ds[iSpecies]*theta/dij*(Ys_j[iSpecies]-Ys_i[iSpecies]);
  }


  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    for (jSpecies  = 0; jSpecies < nSpecies; jSpecies++) {

      // first term
      dJdr_j[iSpecies][jSpecies] +=  0.5*(Ds[iSpecies]*theta/dij *
                                          (Ys_j[iSpecies]*rho_i/rho_j +
                                           Ys_i[iSpecies]));
      dJdr_i[iSpecies][jSpecies] += -0.5*(Ds[iSpecies]*theta/dij *
                                          (Ys_j[iSpecies] +
                                           Ys_i[iSpecies]*rho_j/rho_i));

      // second term
      dJdr_j[iSpecies][jSpecies] +=
          0.25*(Ys_i[iSpecies] - rho_i/rho_j*Ys_j[iSpecies])*sumY
          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_j)
          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_j/rho_j;

      dJdr_i[iSpecies][jSpecies] +=
          0.25*(-rho_j/rho_i*Ys_i[iSpecies]+Ys_j[iSpecies])*sumY
          - 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*Ds[jSpecies]*theta/(dij*rho_i)
          + 0.25*(Ys_i[iSpecies]+Ys_j[iSpecies])*(rho_i+rho_j)*sumY_i/rho_i;
    }

    // first term
    dJdr_j[iSpecies][iSpecies] += -0.5*Ds[iSpecies]*theta/dij*(1+rho_i/rho_j);
    dJdr_i[iSpecies][iSpecies] +=  0.5*Ds[iSpecies]*theta/dij*(1+rho_j/rho_i);

    // second term
    dJdr_j[iSpecies][iSpecies] += 0.25*(1.0+rho_i/rho_j)*sumY;
    dJdr_i[iSpecies][iSpecies] += 0.25*(1.0+rho_j/rho_i)*sumY;
  }

  /*--- Calculate transformation matrix ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    dVdUi[iSpecies][iSpecies] = 1.0;
    dVdUj[iSpecies][iSpecies] = 1.0;
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dVdUi[nSpecies+iDim][iSpecies] = -V_i[VEL_INDEX+iDim]/V_i[RHO_INDEX];
      dVdUj[nSpecies+iDim][iSpecies] = -V_j[VEL_INDEX+iDim]/V_j[RHO_INDEX];
    }
    dVdUi[nSpecies+iDim][nSpecies+iDim] = 1.0/V_i[RHO_INDEX];
    dVdUj[nSpecies+iDim][nSpecies+iDim] = 1.0/V_j[RHO_INDEX];
  }
  for (iVar = 0; iVar < nVar; iVar++) {
    dVdUi[nSpecies+nDim][iVar]   = dTdU_i[iVar];
    dVdUj[nSpecies+nDim][iVar]   = dTdU_j[iVar];
    dVdUi[nSpecies+nDim+1][iVar] = dTvedU_i[iVar];
    dVdUj[nSpecies+nDim+1][iVar] = dTvedU_j[iVar];
  }


  if (nDim == 2) {

    /*--- Geometry parameters ---*/
    thetax = theta + val_normal[0]*val_normal[0]/3.0;
    thetay = theta + val_normal[1]*val_normal[1]/3.0;
    etaz   = val_normal[0]*val_normal[1]/3.0;
    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]  );
    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1]);

    /*--- Populate primitive Jacobian ---*/

    // X-momentum
    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;

    // Y-momentum
    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;

    // Energy
    dFdVj[nSpecies+2][nSpecies]   = pix*val_dS;
    dFdVj[nSpecies+2][nSpecies+1] = piy*val_dS;
    dFdVj[nSpecies+2][nSpecies+2] = ktr*theta/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+3] = kve*theta/dij*val_dS;

    // Vib-el Energy
    dFdVj[nSpecies+3][nSpecies+3] = kve*theta/dij*val_dS;

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

    // Common terms
    dFdVi[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
    dFdVj[nSpecies+2][nSpecies]   += 0.5*val_Fv[nSpecies];
    dFdVi[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
    dFdVj[nSpecies+2][nSpecies+1] += 0.5*val_Fv[nSpecies+1];
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dFdVi[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVj[nSpecies+2][nSpecies+2] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVi[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+2][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
    }

    // Unique terms
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
        dFdVj[nSpecies+2][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVi[nSpecies+2][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
      }
    }

  } //nDim == 2
  else {

    /*--- Geometry parameters ---*/
    thetax = theta + val_normal[0]*val_normal[0]/3.0;
    thetay = theta + val_normal[1]*val_normal[1]/3.0;
    thetaz = theta + val_normal[2]*val_normal[2]/3.0;
    etax   = val_normal[1]*val_normal[2]/3.0;
    etay   = val_normal[0]*val_normal[2]/3.0;
    etaz   = val_normal[0]*val_normal[1]/3.0;
    pix    = mu/dij * (thetax*vel[0] + etaz*vel[1]   + etay*vel[2]  );
    piy    = mu/dij * (etaz*vel[0]   + thetay*vel[1] + etax*vel[2]  );
    piz    = mu/dij * (etay*vel[0]   + etax*vel[1]   + thetaz*vel[2]);

    /*--- Populate primitive Jacobian ---*/

    // X-momentum
    dFdVj[nSpecies][nSpecies]     = mu*thetax/dij*val_dS;
    dFdVj[nSpecies][nSpecies+1]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies][nSpecies+2]   = mu*etay/dij*val_dS;

    // Y-momentum
    dFdVj[nSpecies+1][nSpecies]   = mu*etaz/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+1] = mu*thetay/dij*val_dS;
    dFdVj[nSpecies+1][nSpecies+2] = mu*etax/dij*val_dS;

    // Z-momentum
    dFdVj[nSpecies+2][nSpecies]   = mu*etay/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+1] = mu*etax/dij*val_dS;
    dFdVj[nSpecies+2][nSpecies+2] = mu*thetaz/dij*val_dS;

    // Energy
    dFdVj[nSpecies+3][nSpecies]   = pix*val_dS;
    dFdVj[nSpecies+3][nSpecies+1] = piy*val_dS;
    dFdVj[nSpecies+3][nSpecies+2] = piz*val_dS;
    dFdVj[nSpecies+3][nSpecies+3] = ktr*theta/dij*val_dS;
    dFdVj[nSpecies+3][nSpecies+4] = kve*theta/dij*val_dS;

    // Vib.-el energy
    dFdVj[nSpecies+4][nSpecies+4] = kve*theta/dij*val_dS;

    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        dFdVi[iVar][jVar] = -dFdVj[iVar][jVar];

    // Common terms
    for (iDim = 0; iDim < nDim; iDim++) {
      dFdVi[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
      dFdVj[nSpecies+3][nSpecies+iDim]   += 0.5*val_Fv[nSpecies+iDim];
    }
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      dFdVi[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVj[nSpecies+3][nSpecies+3] += 0.5*val_Fv[iSpecies]*(Ru/Ms[iSpecies] +
                                                             Cvtr[iSpecies]   );
      dFdVi[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+3][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVi[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
      dFdVj[nSpecies+4][nSpecies+4] += 0.5*val_Fv[iSpecies]*val_Mean_Cvve[iSpecies];
    }

    // Unique terms
    for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      for (jSpecies = 0; jSpecies < nSpecies; jSpecies++) {
        dFdVj[iSpecies][jSpecies]   += -dJdr_j[iSpecies][jSpecies]*val_dS;
        dFdVi[iSpecies][jSpecies]   += -dJdr_i[iSpecies][jSpecies]*val_dS;
        dFdVj[nSpecies+3][iSpecies] += -dJdr_j[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVi[nSpecies+3][iSpecies] += -dJdr_i[jSpecies][iSpecies]*hs[jSpecies]*val_dS;
        dFdVj[nSpecies+4][iSpecies] += -dJdr_j[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
        dFdVi[nSpecies+4][iSpecies] += -dJdr_i[jSpecies][iSpecies]*val_Mean_Eve[jSpecies]*val_dS;
      }
    }

  } // nDim == 3

  /*--- dFv/dUij = dFv/dVij * dVij/dUij ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    for (jVar = 0; jVar < nVar; jVar++)
      for (kVar = 0; kVar < nVar; kVar++) {
        val_Jac_i[iVar][jVar] += dFdVi[iVar][kVar]*dVdUi[kVar][jVar];
        val_Jac_j[iVar][jVar] += dFdVj[iVar][kVar]*dVdUj[kVar][jVar];
      }
}

void CAvgGradCorrected_NEMO::ComputeResidual(su2double *val_residual,
                                             su2double **val_Jacobian_i,
                                             su2double **val_Jacobian_j,
                                             CConfig *config) {

  unsigned short iSpecies;
  su2double dist_ij_2;

  /*--- Normalized normal vector ---*/
  Area = 0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }

  /*--- Make a local copy of the primitive variables ---*/
  // NOTE: We are transforming the species density terms to species mass fractions
  // Mass fraction
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    PrimVar_i[iSpecies] = V_i[iSpecies]/V_i[RHO_INDEX];
    PrimVar_j[iSpecies] = V_j[iSpecies]/V_j[RHO_INDEX];
    Mean_PrimVar[iSpecies] = 0.5*(PrimVar_i[iSpecies] + PrimVar_j[iSpecies]);
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iSpecies][iDim] = 0.5*(1.0/V_i[RHO_INDEX] *
                                              (PrimVar_Grad_i[iSpecies][iDim] -
                                               PrimVar_i[iSpecies] *
                                               PrimVar_Grad_i[RHO_INDEX][iDim]) +
                                              1.0/V_j[RHO_INDEX] *
                                              (PrimVar_Grad_j[iSpecies][iDim] -
                                               PrimVar_j[iSpecies] *
                                               PrimVar_Grad_j[RHO_INDEX][iDim]));

    }
  }
  for (iVar = nSpecies; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
  }
  for (iVar = nSpecies; iVar < nPrimVarGrad; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] +
                                          PrimVar_Grad_j[iVar][iDim]);
    }
  }

  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Mean_Eve[iSpecies] = 0.5*(eve_i[iSpecies] + eve_j[iSpecies]);
    Mean_Cvve[iSpecies] = 0.5*(Cvve_i[iSpecies] + Cvve_j[iSpecies]);
  }

  /*--- Mean transport coefficients ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coeff[iSpecies] = 0.5*(Diffusion_Coeff_i[iSpecies] +
                                          Diffusion_Coeff_j[iSpecies]);
  Mean_Laminar_Viscosity           = 0.5*(Laminar_Viscosity_i +
                                          Laminar_Viscosity_j);
  Mean_Thermal_Conductivity        = 0.5*(Thermal_Conductivity_i +
                                          Thermal_Conductivity_j);
  Mean_Thermal_Conductivity_ve     = 0.5*(Thermal_Conductivity_ve_i +
                                          Thermal_Conductivity_ve_j);


  /*--- Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
    Proj_Mean_GradPrimVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Proj_Mean_GradPrimVar_Edge[iVar] += Mean_GradPrimVar[iVar][iDim]*Edge_Vector[iDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] -= (Proj_Mean_GradPrimVar_Edge[iVar] -
                                       (PrimVar_j[iVar]-PrimVar_i[iVar]))*Edge_Vector[iDim] / dist_ij_2;
    }
  }

  /*--- Get projected flux tensor ---*/
  GetViscousProjFlux(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve,
                     Normal, Mean_Diffusion_Coeff,
                     Mean_Laminar_Viscosity,
                     Mean_Thermal_Conductivity,
                     Mean_Thermal_Conductivity_ve,
                     config);

  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);

    GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve, Mean_Cvve,
                       Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,
                       Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                       dist_ij, UnitNormal, Area, Proj_Flux_Tensor,
                       val_Jacobian_i, val_Jacobian_j, config);

  }
}

