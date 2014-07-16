/*!
 * \file numerics_machine_learning_direct_turbulent.cpp
 * \brief This file contains all the convective term discretization.
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 3.2.0 "eagle"
 *
 * SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).
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

#include "../include/numerics_machine_learning_turbulent.hpp"


SpalartAllmarasOtherOutputs::SpalartAllmarasOtherOutputs(){}

SpalartAllmarasOtherOutputs::~SpalartAllmarasOtherOutputs(){}

SpalartAllmarasConstants::SpalartAllmarasConstants(){
  /*--- Spalart-Allmaras closure constants ---*/
  cv1_3 = pow(7.1,3.0);
  k2 = pow(0.41,2.0);
  cb1 = 0.1355;
  cw2 = 0.3;
  cw3_6 = pow(2.0,6.0);
  sigma = 2./3.;
  cb2 = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1+cb2)/sigma;
}

SpalartAllmarasConstants::~SpalartAllmarasConstants(){}

SpalartAllmarasInputs::SpalartAllmarasInputs(int nDim){
  init(nDim, 1e-10);
}

void SpalartAllmarasInputs::init(int nDim, double limiter){
  this->nDim = nDim;
  this->limiter = limiter;
  DUiDXj = new double*[nDim];
  for (int i = 0; i < nDim; i++){
    DUiDXj[i] = new double[nDim];
  }
  DTurb_Kin_Visc_DXj = new double[nDim];
  return;
}

SpalartAllmarasInputs::SpalartAllmarasInputs(int nDim, double limiter){
  init(nDim, limiter);
}

SpalartAllmarasInputs::~SpalartAllmarasInputs(){
  for (int i = 0; i < nDim; i++){
    delete [] DUiDXj[i];
  }
  delete [] DUiDXj;
  delete [] DTurb_Kin_Visc_DXj;
}

int SpalartAllmarasInputs::GetNumDim(){
  return nDim;
}

double SpalartAllmarasInputs::GetLimiter(){
  return limiter;
}

double** SpalartAllmarasInputs::GetMeanFlowGradient(){
  return DUiDXj;
}
double* SpalartAllmarasInputs::GetTurbKinViscGradient(){
  return DTurb_Kin_Visc_DXj;
}

void SpalartAllmarasInputs::Set(double** DUiDXj, double* DTurb_Kin_Visc_DXj, bool rotating_frame, bool transition, double dist, double Laminar_Viscosity, double Density, double Turbulent_Kinematic_Viscosity, double intermittency){
  for (int i=0; i <this->nDim; i++){
    for (int j = 0; j < this->nDim; j++){
      this->DUiDXj[i][j] = DUiDXj[i][j];
    }
    this->DTurb_Kin_Visc_DXj[i] = DTurb_Kin_Visc_DXj[i];
  }
  this->rotating_frame = rotating_frame;
  this->transition = transition;
  this->dist = dist;
  this->Laminar_Viscosity = Laminar_Viscosity;
  this -> Density = Density;
  this->Turbulent_Kinematic_Viscosity = Turbulent_Kinematic_Viscosity;
  this->intermittency = intermittency;
  return;
}

/* Computes the spalart-allmaras source term.
 the outputs are
 (Production, Destruction, CrossProduction, Total)
 jacobian is dSourceDNuHat (Turbulent_Kinematic_Viscosity)
    This ignores the contribution of the cross-production term
 
 Does not include the volume term
 */
void SpalartAllmarasSourceTerm(SpalartAllmarasInputs* inputs, SpalartAllmarasConstants* constants, double* output_residual, double* output_jacobian, SpalartAllmarasOtherOutputs* otherOutput){
  double dist = inputs->dist; // Wall distance
  int nDim = inputs->GetNumDim();
  // Limit if too close to the wall
  double limiter = inputs->GetLimiter();
  int nResidOutputs = 4;
  int nJacOutputs = 1;
  if (dist < limiter){
    for (int i = 0; i < nResidOutputs; i++){
      output_residual[i] = 0;
    }
    for (int i = 0; i < nJacOutputs; i++){
      output_jacobian[i] = 0;
    }
    return;
  }
  
  double **DUiDXj = inputs->GetMeanFlowGradient();
  double Vorticity = ComputeVorticity(nDim,DUiDXj);
  double Omega = sqrt(Vorticity);
  
  double StrainMag;
  double Laminar_Viscosity = inputs->Laminar_Viscosity;
  double Density = inputs->Density;
  double Turbulent_Kinematic_Viscosity = inputs->Turbulent_Kinematic_Viscosity;
  bool transition = inputs->transition;
  double intermittency = inputs->intermittency;
  
  double div,dist_2, Laminar_Kinematic_Viscosity, J, J_2, J_3,
  fv1, fv2, S, inv_k2_d2, Shat,inv_Shat, r, g, g_6, glim, fw, norm2_Grad,
  dfv1,dfv2, dr, dg, dfw;
  
  double Production, Destruction, CrossProduction;
  double *DTurb_Kin_Visc_DXj = inputs->GetTurbKinViscGradient();
  double dShat;
  
  // Correction for rotating frame
  if (inputs->rotating_frame) {
    div = DUiDXj[0][0] + DUiDXj[1][1];
    if (nDim == 3) div += DUiDXj[2][2];
    StrainMag = 0.0;
    // add diagonals
    StrainMag += pow(DUiDXj[0][0] - 1.0/3.0*div,2.0);
    StrainMag += pow(DUiDXj[1][1] - 1.0/3.0*div,2.0);
    if (nDim == 3) StrainMag += pow(DUiDXj[2][2] - 1.0/3.0*div,2.0);
    // add off diagonals
    StrainMag += 2.0*pow(0.5*(DUiDXj[0][1]+DUiDXj[1][0]),2.0);
    if (nDim == 3) {
      StrainMag += 2.0*pow(0.5*(DUiDXj[0][2]+DUiDXj[2][0]),2.0);
      StrainMag += 2.0*pow(0.5*(DUiDXj[1][2]+DUiDXj[2][1]),2.0);
    }
    StrainMag = sqrt(2.0*StrainMag);
    Omega += 2.0*min(0.0,StrainMag-Omega);
  }
  /*--- Production term ---*/
  
  dist_2 = dist*dist;
  Laminar_Kinematic_Viscosity = Laminar_Viscosity/Density;
  J = Turbulent_Kinematic_Viscosity/Laminar_Kinematic_Viscosity;
  J_2 = J*J;
  J_3 = J_2*J;
  fv1 = J_3/(J_3+constants->cv1_3);
  fv2 = 1.0 - J/(1.0+J*fv1);
  S = Omega;
  inv_k2_d2 = 1.0/(constants->k2*dist_2);
  Shat = S + Turbulent_Kinematic_Viscosity*fv2*inv_k2_d2;
  inv_Shat = 1.0/max(Shat, 1.0e-10);
  Production = constants->cb1*Shat*Turbulent_Kinematic_Viscosity;
  if (transition){
    Production *= intermittency;
  }
  
  double limitOmega = max(Omega, 1.0e-10);
  otherOutput->mul_production = constants->cb1 * (1 + Turbulent_Kinematic_Viscosity * inv_k2_d2 * fv2 / limitOmega);
//  cout << "mul prod = " << otherOutput->mul_production << endl;
//  cout << "omega = " << Omega << endl;
//  cout << "production = " << Production << endl;
//  cout << "mul = " << otherOutput->mul_production * Omega * Turbulent_Kinematic_Viscosity << endl;
//  cout << "turb kin visc = " << Turbulent_Kinematic_Viscosity << endl;
//  cout << "wall dist = " << dist << endl;
//  cout << "fv2 = " << fv2 << endl;
//  cout << "inv k2 d2 = " << inv_k2_d2 << endl;
//  if (otherOutput->mul_production < -1e+5){
//    cout << endl;
//    for (int i = 0; i < 1; i++){
//  throw "rah";
//      cout << "rah" << endl;
//    }
//  }
  
  /*--- Destruction term ---*/
  
  r = min(Turbulent_Kinematic_Viscosity*inv_Shat*inv_k2_d2,10.0);
  g = r + constants->cw2*(pow(r,6.0)-r);
  g_6 =	pow(g,6.0);
  
  double cw3_6 = constants->cw3_6;
  glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
  fw = g*glim;
  
  otherOutput->fw = fw;
  
  otherOutput->mul_destruction = constants->cw1 * fw;
  
  Destruction = constants->cw1*fw*Turbulent_Kinematic_Viscosity*Turbulent_Kinematic_Viscosity/dist_2;
  if (transition){
    Destruction *= min(max(intermittency,0.1),1.0);
  }
  
  /*--- Diffusion term ---*/
  norm2_Grad = 0.0;
  for (int iDim = 0; iDim < nDim; iDim++){
    norm2_Grad += DTurb_Kin_Visc_DXj[iDim]*DTurb_Kin_Visc_DXj[iDim];
  }
  CrossProduction = constants->cb2_sigma*norm2_Grad;
  
  otherOutput->mul_crossproduction = constants->cb2_sigma;
  
  output_residual[0] = Production;
  output_residual[1] = Destruction;
  output_residual[2] = CrossProduction;
  output_residual[3] = Production - Destruction + CrossProduction;
  
  /*--- Implicit part ---*/
  
  /*--- Production term ---*/
  
  dfv1 = 3.0*J_2*constants->cv1_3/(Laminar_Kinematic_Viscosity*pow(J_3+constants->cv1_3,2.));
  dfv2 = -(1/Laminar_Kinematic_Viscosity-J_2*dfv1)/pow(1.+J*fv1,2.);
  if ( Shat <= 1.0e-10 ){
    dShat = 0.0;
  }
  else{
    dShat = (fv2+Turbulent_Kinematic_Viscosity*dfv2)*inv_k2_d2;
  }
  output_jacobian[0] = constants->cb1*(Turbulent_Kinematic_Viscosity*dShat+Shat);
  
  /*--- Destruction term ---*/
  
  dr = (Shat-Turbulent_Kinematic_Viscosity*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
  if (r == 10.0){
    dr = 0.0;
  }
  dg = dr*(1.+constants->cw2*(6.*pow(r,5.)-1.));
  dfw = dg*glim*(1.-g_6/(g_6+constants->cw3_6));
  output_jacobian[0] -= constants->cw1*(dfw*Turbulent_Kinematic_Viscosity + 2.*fw)*Turbulent_Kinematic_Viscosity/dist_2;
  
  // NOTE: Do not have derivative with respect to the cross production term
  
  inputs->Omega = Omega;
  
  return;
  
};

double ComputeVorticity(int nDim, double** DUiDXj){
  double Vorticity = (DUiDXj[1][0]-DUiDXj[0][1])*(DUiDXj[1][0]-DUiDXj[0][1]);
  if (nDim == 3){
    Vorticity += ( (DUiDXj[2][1]-DUiDXj[1][2])*(DUiDXj[2][1]-DUiDXj[1][2]) + (DUiDXj[0][2]-DUiDXj[2][0])*(DUiDXj[0][2]-DUiDXj[2][0]));
  }
  return Vorticity;
};
