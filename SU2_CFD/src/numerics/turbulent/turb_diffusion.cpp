/*!
 * \file turb_diffusion.cpp
 * \brief Implementation of numerics classes to compute viscous
 *        fluxes in turbulence problems.
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

#include "../../../include/numerics/turbulent/turb_diffusion.hpp"
#include "../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CAvgGrad_Scalar::CAvgGrad_Scalar(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 bool correct_grad,
                                 const CConfig* config) :
  CNumerics(val_nDim, val_nVar, config),
  correct_gradient(correct_grad),
  incompressible(config->GetKind_Regime() == INCOMPRESSIBLE),
  exact_jacobian(config->GetUse_Accurate_Visc_Jacobians())
{
  Proj_Mean_GradTurbVar = new su2double [nVar];

  Flux = new su2double [nVar];
  Jacobian_i = new su2double* [nVar];
  Jacobian_j = new su2double* [nVar];
  jacobianWeights_i = new su2double* [nVar];
  jacobianWeights_j = new su2double* [nVar];
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double [nVar] ();
    Jacobian_j[iVar] = new su2double [nVar] ();
    jacobianWeights_i[iVar] = new su2double [nDim] ();
    jacobianWeights_j[iVar] = new su2double [nDim] ();
  }
}

CAvgGrad_Scalar::~CAvgGrad_Scalar(void) {

  delete [] Proj_Mean_GradTurbVar;

  delete [] Flux;
  
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    delete [] Jacobian_i[iVar];
    delete [] Jacobian_j[iVar];
    delete [] jacobianWeights_i[iVar];
    delete [] jacobianWeights_j[iVar];
  }
  delete [] Jacobian_i;
  delete [] Jacobian_j;
  delete [] jacobianWeights_i;
  delete [] jacobianWeights_j;
}

CNumerics::ResidualType<> CAvgGrad_Scalar::ComputeResidual(const CConfig* config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Coord_i, nDim); 
  AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(TurbVar_j, nVar); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  ExtraADPreaccIn();

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+6); 
    AD::SetPreaccIn(V_j, nDim+6);

    Density_i = V_i[nDim+2]; Laminar_Viscosity_i = V_i[nDim+4]; Eddy_Viscosity_i = V_i[nDim+5];
    Density_j = V_j[nDim+2]; Laminar_Viscosity_j = V_j[nDim+4]; Eddy_Viscosity_j = V_j[nDim+5];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); 
    AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2]; Laminar_Viscosity_i = V_i[nDim+5]; Eddy_Viscosity_i = V_i[nDim+6];
    Density_j = V_j[nDim+2]; Laminar_Viscosity_j = V_j[nDim+5]; Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Flux[iVar] = 0.0;
    for (auto jVar = 0u; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] = 0.0;
      Jacobian_j[iVar][jVar] = 0.0;
    }
    for (auto iDim = 0u; iDim < nDim; iDim++) {
      jacobianWeights_i[iVar][iDim] = 0.0;
      jacobianWeights_j[iVar][iDim] = 0.0;
    }
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  GeometryToolbox::Distance(nDim,Coord_j,Coord_i,Edge_Vector);
  dist_ij_2 = GeometryToolbox::SquaredNorm(nDim,Edge_Vector);
  proj_vector_ij = correct_gradient? GeometryToolbox::DotProduct(nDim,Normal,Edge_Vector)/dist_ij_2 : su2double(1.0);

  /*--- Mean gradient approximation ---*/
  for (auto iVar = 0u; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar[iVar] = 0.5*(GeometryToolbox::DotProduct(nDim,TurbVar_Grad_i[iVar],Normal)
                                     + GeometryToolbox::DotProduct(nDim,TurbVar_Grad_j[iVar],Normal));      
    if (correct_gradient) {
      Proj_Mean_GradTurbVar_Edge = 0.5*(GeometryToolbox::DotProduct(nDim,TurbVar_Grad_i[iVar],Edge_Vector)
                                      + GeometryToolbox::DotProduct(nDim,TurbVar_Grad_j[iVar],Edge_Vector));
      Proj_Mean_GradTurbVar[iVar] -= Proj_Mean_GradTurbVar_Edge*proj_vector_ij -
                                    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
    }
  }

  FinishResidualCalc(config);

  AD::SetPreaccOut(Flux, nVar);
  AD::EndPreacc();

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j, jacobianWeights_i, jacobianWeights_j);

}

void CAvgGrad_Scalar::CorrectJacobian(const CConfig *config) {
  
  /*--- Add contributions of nodal gradients ---*/
 
  for (auto iDim = 0u; iDim < nDim; iDim++) {
    const su2double weight = 0.5*(Normal[iDim]/proj_vector_ij - Edge_Vector[iDim]);
    for (auto iVar= 0; iVar < nVar; iVar++) {
      jacobianWeights_i[iVar][iDim] = -weight*Jacobian_i[iVar][iVar];      
      jacobianWeights_j[iVar][iDim] =  weight*Jacobian_j[iVar][iVar];
    }
  }
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                 bool correct_grad, const CConfig* config) :
                 CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) { }

void CAvgGrad_TurbSA::ExtraADPreaccIn() { }

void CAvgGrad_TurbSA::FinishResidualCalc(const CConfig* config) {

  /*--- Compute mean effective viscosity ---*/

  su2double nu_i = Laminar_Viscosity_i/Density_i;
  su2double nu_j = Laminar_Viscosity_j/Density_j;
  su2double nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);

  Flux[0] = nu_e*Proj_Mean_GradTurbVar[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  Jacobian_i[0][0] = -nu_e*proj_vector_ij/sigma;
  Jacobian_j[0][0] =  nu_e*proj_vector_ij/sigma;
  
  if (exact_jacobian) CorrectJacobian(config);

  Jacobian_i[0][0] += 0.5*Proj_Mean_GradTurbVar[0]/sigma;
  Jacobian_j[0][0] += 0.5*Proj_Mean_GradTurbVar[0]/sigma;

}

CAvgGrad_TurbSA_Neg::CAvgGrad_TurbSA_Neg(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         bool correct_grad,
                                         const CConfig* config) :
                     CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) { }

void CAvgGrad_TurbSA_Neg::ExtraADPreaccIn() { }

void CAvgGrad_TurbSA_Neg::FinishResidualCalc(const CConfig* config) {

  /*--- Compute mean effective viscosity ---*/

  su2double nu_i = Laminar_Viscosity_i/Density_i;
  su2double nu_j = Laminar_Viscosity_j/Density_j;

  su2double nu_ij = 0.5*(nu_i+nu_j);
  su2double nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);

  su2double nu_e;

  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    su2double Xi = nu_tilde_ij/nu_ij;
    su2double fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
  }

  Flux[0] = nu_e*Proj_Mean_GradTurbVar[0]/sigma;

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar[0]-nu_e*proj_vector_ij)/sigma;
  Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar[0]+nu_e*proj_vector_ij)/sigma;

}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim,
                                   unsigned short val_nVar,
                                   const su2double *constants,
                                   bool correct_grad,
                                   const CConfig* config) :
  CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config),
  sigma_k1(constants[0]),
  sigma_k2(constants[1]),
  sigma_om1(constants[2]),
  sigma_om2(constants[3]),
  a1(constants[7]){

}

void CAvgGrad_TurbSST::ExtraADPreaccIn() {
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);
}

void CAvgGrad_TurbSST::FinishResidualCalc(const CConfig* config) {

  /*--- Compute the blended constant for the viscous terms ---*/
  const su2double sigma_kine_i  = F1_i*sigma_k1  + (1.0 - F1_i)*sigma_k2;
  const su2double sigma_kine_j  = F1_j*sigma_k1  + (1.0 - F1_j)*sigma_k2;
  const su2double sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  const su2double sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;

  /*--- Compute mean effective viscosity ---*/
  const su2double diff_i_kine  = Laminar_Viscosity_i + Eddy_Viscosity_i*sigma_kine_i;
  const su2double diff_j_kine  = Laminar_Viscosity_j + Eddy_Viscosity_j*sigma_kine_j;
  const su2double diff_i_omega = Laminar_Viscosity_i + Eddy_Viscosity_i*sigma_omega_i;
  const su2double diff_j_omega = Laminar_Viscosity_j + Eddy_Viscosity_j*sigma_omega_j;

  const su2double diff_kine  = 0.5*(diff_i_kine  + diff_j_kine);
  const su2double diff_omega = 0.5*(diff_i_omega + diff_j_omega);

  Flux[0] = Proj_Mean_GradTurbVar[0]*diff_kine;
  Flux[1] = Proj_Mean_GradTurbVar[1]*diff_omega;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  const bool wasActive = AD::BeginPassive();

  const su2double proj_on_rho_i = proj_vector_ij/Density_i;
  const su2double proj_on_rho_j = proj_vector_ij/Density_j;
      
  Jacobian_i[0][0] = -proj_on_rho_i*diff_kine; Jacobian_i[1][1] = -proj_on_rho_i*diff_omega;
  Jacobian_j[0][0] =  proj_on_rho_j*diff_kine; Jacobian_j[1][1] =  proj_on_rho_j*diff_omega;

  if (exact_jacobian) CorrectJacobian(config);
  if (!correct_gradient) {
    Jacobian_i[0][0] = Jacobian_i[1][1] = 0.0;
    Jacobian_j[0][0] = Jacobian_j[1][1] = 0.0;
  }

  if (exact_jacobian) {
    
    /*--- Jacobian wrt eddy viscosity ---*/
    
    const su2double zeta_i = max(TurbVar_i[1], VorticityMag_i*F2_i/a1);
    const su2double zeta_j = max(TurbVar_j[1], VorticityMag_j*F2_j/a1);

    Jacobian_i[0][0] += 0.5/zeta_i*Proj_Mean_GradTurbVar[0]*sigma_kine_i;
    Jacobian_i[1][0] += 0.5/zeta_i*Proj_Mean_GradTurbVar[1]*sigma_omega_i;
    if (TurbVar_i[1] > VorticityMag_i*F2_i/a1) {
      Jacobian_i[0][1] -= 0.5*TurbVar_i[0]/pow(TurbVar_i[1],2.0)*Proj_Mean_GradTurbVar[0]*sigma_kine_i;
      Jacobian_i[1][1] -= 0.5*TurbVar_i[0]/pow(TurbVar_i[1],2.0)*Proj_Mean_GradTurbVar[1]*sigma_omega_i;
    }

    Jacobian_j[0][0] += 0.5/zeta_j*Proj_Mean_GradTurbVar[0]*sigma_kine_j;
    Jacobian_j[1][0] += 0.5/zeta_j*Proj_Mean_GradTurbVar[1]*sigma_omega_j;
    if (TurbVar_j[1] > VorticityMag_j*F2_j/a1) {
      Jacobian_j[0][1] -= 0.5*TurbVar_j[0]/pow(TurbVar_j[1],2.0)*Proj_Mean_GradTurbVar[0]*sigma_kine_j;
      Jacobian_j[1][1] -= 0.5*TurbVar_j[0]/pow(TurbVar_j[1],2.0)*Proj_Mean_GradTurbVar[1]*sigma_omega_j;
    }
    
    /*--- Jacobian wrt laminar viscosity ---*/

    const su2double Cv    = Gas_Constant/Gamma_Minus_One;
    const su2double muref = config->GetMu_RefND();
    const su2double Tref  = config->GetMu_Temperature_RefND();
    const su2double Sref  = config->GetMu_SND();
    
    const su2double T_i      = V_i[0];
    const su2double T_j      = V_j[0];
    const su2double dmudT_i  = muref*(Tref+Sref)/pow(Tref,1.5) 
                             * (3.*Sref*sqrt(T_i) + pow(T_i,1.5))
                             / (2.*pow((T_i+Sref),2.));
    const su2double dmudT_j  = muref*(Tref+Sref)/pow(Tref,1.5) 
                             * (3.*Sref*sqrt(T_j) + pow(T_j,1.5))
                             / (2.*pow((T_j+Sref),2.));
    const su2double factor_i = dmudT_i/(Density_i*Cv);
    const su2double factor_j = dmudT_j/(Density_j*Cv);
    
    for (auto iVar = 0u; iVar < nVar; iVar++) {
      Jacobian_i[iVar][0] -= 0.5*factor_i*Proj_Mean_GradTurbVar[iVar];
      Jacobian_j[iVar][0] -= 0.5*factor_j*Proj_Mean_GradTurbVar[iVar];
    }
  }

  AD::EndPassive(wasActive);

}
