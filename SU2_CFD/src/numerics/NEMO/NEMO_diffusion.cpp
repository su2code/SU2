/*!
 * \file NEMO_diffusion.cpp
 * \brief Implementation of numerics classes for discretization
 *        of viscous fluxes in fluid flow NEMO problems.
 * \author S.R. Copeland, W. Maier, C. Garbacz
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/numerics/NEMO/NEMO_diffusion.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CAvgGrad_NEMO::CAvgGrad_NEMO(unsigned short val_nDim,
                             unsigned short val_nVar,
                             unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad,
                             CConfig *config) : CNEMONumerics(val_nDim, val_nVar,
                                                              val_nPrimVar, val_nPrimVarGrad,
                                                              config) {

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

  Flux   = new su2double[nVar];

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
  delete [] Flux;
}

CNumerics::ResidualType<> CAvgGrad_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iSpecies, iVar, iDim;

  /*--- Normalized normal vector ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  /*--- Mean transport coefficients ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coeff[iSpecies] = 0.5*(Diffusion_Coeff_i[iSpecies] +
                                          Diffusion_Coeff_j[iSpecies]);
  Mean_Laminar_Viscosity           = 0.5*(Laminar_Viscosity_i +
                                          Laminar_Viscosity_j);
  Mean_Eddy_Viscosity              = 0.5*(Eddy_Viscosity_i +
                                          Eddy_Viscosity_j);
  Mean_Thermal_Conductivity        = 0.5*(Thermal_Conductivity_i +
                                          Thermal_Conductivity_j);
  Mean_Thermal_Conductivity_ve     = 0.5*(Thermal_Conductivity_ve_i +
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
                     Mean_Diffusion_Coeff, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                     Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                     config);


  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);

    GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve, Mean_Cvve,
                       Mean_Diffusion_Coeff, Mean_Laminar_Viscosity,Mean_Eddy_Viscosity,
                       Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                       dist_ij, UnitNormal, Area, Proj_Flux_Tensor,
                       Jacobian_i, Jacobian_j, config);
  }

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

CAvgGradCorrected_NEMO::CAvgGradCorrected_NEMO(unsigned short val_nDim,
                                               unsigned short val_nVar,
                                               unsigned short val_nPrimVar,
                                               unsigned short val_nPrimVarGrad,
                                               CConfig *config) : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad,
                                                          config) {

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

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

  Flux   = new su2double[nVar];

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

  delete [] Flux;

}

CNumerics::ResidualType<> CAvgGradCorrected_NEMO::ComputeResidual(const CConfig *config) {

  unsigned short iSpecies;
  su2double dist_ij_2;

  /*--- Normalized normal vector ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

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
    Mean_Eve[iSpecies]  = 0.5*(eve_i[iSpecies]  + eve_j[iSpecies]);
    Mean_Cvve[iSpecies] = 0.5*(Cvve_i[iSpecies] + Cvve_j[iSpecies]);
  }

  /*--- Mean transport coefficients ---*/
  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
    Mean_Diffusion_Coeff[iSpecies] = 0.5*(Diffusion_Coeff_i[iSpecies] +
                                          Diffusion_Coeff_j[iSpecies]);
  Mean_Laminar_Viscosity           = 0.5*(Laminar_Viscosity_i +
                                          Laminar_Viscosity_j);
  Mean_Eddy_Viscosity              = 0.5*(Eddy_Viscosity_i +
                                          Eddy_Viscosity_j);
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
                     Mean_Eddy_Viscosity,
                     Mean_Thermal_Conductivity,
                     Mean_Thermal_Conductivity_ve,
                     config);

  /*--- Update viscous residual ---*/
  for (iVar = 0; iVar < nVar; iVar++)
    Flux[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/
  if (implicit) {
    dist_ij = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
    dist_ij = sqrt(dist_ij);

    GetViscousProjJacs(Mean_PrimVar, Mean_GradPrimVar, Mean_Eve, Mean_Cvve,
                       Mean_Diffusion_Coeff, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                       Mean_Thermal_Conductivity, Mean_Thermal_Conductivity_ve,
                       dist_ij, UnitNormal, Area, Proj_Flux_Tensor,
                       Jacobian_i, Jacobian_j, config);

  }

  return ResidualType<>(Flux, Jacobian_j, Jacobian_j);
}

