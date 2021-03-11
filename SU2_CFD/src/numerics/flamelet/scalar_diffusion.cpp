/*!
 * \file scalar_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in turbulence problems.
 * \author F. Palacios, T. Economon
 * \version 7.1.1 "Blackbird"
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

#include "../../../include/numerics/flamelet/scalar_diffusion.hpp"

CAvgGrad_transportedScalar::CAvgGrad_transportedScalar(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 bool correct_grad,
                                 const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  correct_gradient(correct_grad),
  implicit(config->GetKind_TimeIntScheme_Scalar() == EULER_IMPLICIT),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE)
{
  Proj_Mean_GradScalarVar_Normal = new su2double [nVar] ();
  Proj_Mean_GradScalarVar_Edge = new su2double [nVar] ();
  Proj_Mean_GradScalarVar = new su2double [nVar] ();

  Flux = new su2double [nVar] ();
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar] ();
    Jacobian_j[iVar] = new su2double [nVar] ();
  }
}

CAvgGrad_transportedScalar::~CAvgGrad_transportedScalar(void) {

  delete [] Proj_Mean_GradScalarVar_Normal;
  delete [] Proj_Mean_GradScalarVar_Edge;
  delete [] Proj_Mean_GradScalarVar;

  delete [] Flux;
  if (Jacobian_i != nullptr) {
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete [] Jacobian_i[iVar];
      delete [] Jacobian_j[iVar];
    }
    delete [] Jacobian_i;
    delete [] Jacobian_j;
  }
}

CNumerics::ResidualType<> CAvgGrad_transportedScalar::ComputeResidual(const CConfig* config) {

  unsigned short iVar, iDim;

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(ScalarVar_Grad_j, nVar, nDim);
  if (correct_gradient) {
    AD::SetPreaccIn(ScalarVar_i, nVar); AD::SetPreaccIn(ScalarVar_j ,nVar);
  }
  ExtraADPreaccIn();

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6); AD::SetPreaccIn(V_j, nDim+6);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+4];  Laminar_Viscosity_j = V_j[nDim+4];
    Eddy_Viscosity_i = V_i[nDim+5];     Eddy_Viscosity_j = V_j[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradScalarVar_Normal[iVar] = 0.0;
    Proj_Mean_GradScalarVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      su2double Mean_GradScalarVar = 0.5*(ScalarVar_Grad_i[iVar][iDim] +
                                        ScalarVar_Grad_j[iVar][iDim]);

      Proj_Mean_GradScalarVar_Normal[iVar] += Mean_GradScalarVar * Normal[iDim];

      if (correct_gradient)
        Proj_Mean_GradScalarVar_Edge[iVar] += Mean_GradScalarVar * Edge_Vector[iDim];
    }
    Proj_Mean_GradScalarVar[iVar] = Proj_Mean_GradScalarVar_Normal[iVar];
    if (correct_gradient) {
      Proj_Mean_GradScalarVar[iVar] -= Proj_Mean_GradScalarVar_Edge[iVar]*proj_vector_ij -
      (ScalarVar_j[iVar]-ScalarVar_i[iVar])*proj_vector_ij;
    }
  }

  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);

}

CAvgGrad_transportedScalar_general::CAvgGrad_transportedScalar_general(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   bool correct_grad,
                                   const CConfig* config) :
  CAvgGrad_transportedScalar(val_nDim, val_nVar, correct_grad, config){ }


void CAvgGrad_transportedScalar_general::ExtraADPreaccIn() {
  //AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);
}

void CAvgGrad_transportedScalar_general::FinishResidualCalc(const CConfig* config) {

  for (auto iVar = 0u; iVar < nVar; iVar++) {
    
    /*--- Get the diffusion coefficient(s). ---*/
    
    //Mean_Diffusivity[iVar] = 0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar]);

    /*--- Compute the viscous residual. ---*/

    //val_residual[iVar] = Mean_Diffusivity[iVar]*Proj_Mean_GradScalarVar[iVar];
    Flux[iVar] = 0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar])*Proj_Mean_GradScalarVar[iVar];
    
    /*--- Use TSL approx. to compute derivatives of the gradients. ---*/

    if (implicit) {
      for (auto jVar = 0u; jVar < nVar; jVar++) {
        if (iVar == jVar) {
          //Jacobian_i[iVar][jVar] = -Mean_Diffusivity[iVar]*proj_vector_ij;
          //Jacobian_j[iVar][jVar] =  Mean_Diffusivity[iVar]*proj_vector_ij;
          Jacobian_i[iVar][jVar] = -0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar])*proj_vector_ij;
          Jacobian_j[iVar][jVar] =  0.5*(Diffusion_Coeff_i[iVar] + Diffusion_Coeff_j[iVar])*proj_vector_ij;
        } else {
          Jacobian_i[iVar][jVar] = 0.0;
          Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    }
  }


}
