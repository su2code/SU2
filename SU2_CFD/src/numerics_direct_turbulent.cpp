/*!
 * \file numerics_direct_turbulent.cpp
 * \brief This file contains all the convective term discretization.
 * \author F. Palacios, A. Bueno
 * \version 4.3.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"
#include <limits>

CUpwSca_TurbSA::CUpwSca_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbSA::~CUpwSca_TurbSA(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSA::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  q_ij = 0.0;
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(TurbVar_i[0]); AD::SetPreaccIn(TurbVar_j[0]);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement){
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();
}

CAvgGrad_TurbSA::CAvgGrad_TurbSA(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSA::~CAvgGrad_TurbSA(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSA::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {    
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
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
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGrad_TurbSA_Neg::CAvgGrad_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  cn1   = 16.0;
  fn    = 0.0;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSA_Neg::~CAvgGrad_TurbSA_Neg(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);   AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);   AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  
  nu_ij = 0.5*(nu_i+nu_j);
  nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);

  Xi = nu_tilde_ij/nu_ij;
  
  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
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
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CAvgGradCorrected_TurbSA::CAvgGradCorrected_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSA::~CAvgGradCorrected_TurbSA(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSA::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);   AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);   AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

CAvgGradCorrected_TurbSA_Neg::CAvgGradCorrected_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  cn1   = 16.0;
  fn    = 0.0;

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSA_Neg::~CAvgGradCorrected_TurbSA_Neg(void) {
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim);  AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  
  nu_ij = 0.5*(nu_i+nu_j);
  nu_tilde_ij = 0.5*(TurbVar_i[0]+TurbVar_j[0]);
  
  Xi = nu_tilde_ij/nu_ij;
  
  if (nu_tilde_ij > 0.0) {
    nu_e = nu_ij + nu_tilde_ij;
  }
  else {
    fn = (cn1 + Xi*Xi*Xi)/(cn1 - Xi*Xi*Xi);
    nu_e = nu_ij + fn*nu_tilde_ij;
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
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbSA::CSourcePieceWise_TurbSA(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  
  /*--- Spalart-Allmaras closure constants ---*/
  
  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  
}

CSourcePieceWise_TurbSA::~CSourcePieceWise_TurbSA(void) { }

void CSourcePieceWise_TurbSA::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  val_residual[0] = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  /*--- Evaluate Omega ---*/
  
  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);
  
  /*--- Rotational correction term ---*/
  
  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }
  
  if (dist_i > 1e-10) {
    
    /*--- Production term ---*/
    
    dist_i_2 = dist_i*dist_i;
    nu = Laminar_Viscosity_i/Density_i;
    Ji = TurbVar_i[0]/nu;
    Ji_2 = Ji*Ji;
    Ji_3 = Ji_2*Ji;
    fv1 = Ji_3/(Ji_3+cv1_3);
    fv2 = 1.0 - Ji/(1.0+Ji*fv1);
    ft2 = ct3*exp(-ct4*Ji_2);
    S = Omega;
    inv_k2_d2 = 1.0/(k2*dist_i_2);
    
    Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
    Shat = max(Shat, 1.0e-10);
    inv_Shat = 1.0/Shat;
    
    /*--- Production term ---*/;

//    Original SA model
//    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;
    
    Production = cb1*Shat*TurbVar_i[0]*Volume;

    /*--- Destruction term ---*/
    
    r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
    g = r + cw2*(pow(r,6.0)-r);
    g_6 = pow(g,6.0);
    glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
    fw = g*glim;
    
//    Original SA model
//    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
    
    Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;

    /*--- Diffusion term ---*/
    
    norm2_Grad = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
    
    CrossProduction = cb2_sigma*norm2_Grad*Volume;
    
    val_residual[0] = Production - Destruction + CrossProduction;
    
    /*--- Implicit part, production term ---*/
    
    dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
    dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
    if ( Shat <= 1.0e-10 ) dShat = 0.0;
    else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
    val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
    
    /*--- Implicit part, destruction term ---*/
    
    dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
    if (r == 10.0) dr = 0.0;
    dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
    dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
    val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;
    
  }

//  AD::SetPreaccOut(val_residual[0]);
//  AD::EndPreacc();
  
}

CSourcePieceWise_TurbSA_Neg::CSourcePieceWise_TurbSA_Neg(unsigned short val_nDim, unsigned short val_nVar,
                                                         CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  rotating_frame = config->GetRotating_Frame();
  
  /*--- Negative Spalart-Allmaras closure constants ---*/
  
  cv1_3 = pow(7.1, 3.0);
  k2    = pow(0.41, 2.0);
  cb1   = 0.1355;
  cw2   = 0.3;
  ct3   = 1.2;
  ct4   = 0.5;
  cw3_6 = pow(2.0, 6.0);
  sigma = 2./3.;
  cb2   = 0.622;
  cb2_sigma = cb2/sigma;
  cw1 = cb1/k2+(1.0+cb2)/sigma;
  
}

CSourcePieceWise_TurbSA_Neg::~CSourcePieceWise_TurbSA_Neg(void) {
  
}

void CSourcePieceWise_TurbSA_Neg::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
//  AD::StartPreacc();
//  AD::SetPreaccIn(V_i, nDim+6);
//  AD::SetPreaccIn(Vorticity_i, nDim);
//  AD::SetPreaccIn(StrainMag_i);
//  AD::SetPreaccIn(TurbVar_i[0]);
//  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
//  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
  }
  
  val_residual[0] = 0.0;
  Production      = 0.0;
  Destruction     = 0.0;
  CrossProduction = 0.0;
  val_Jacobian_i[0][0] = 0.0;
  
  /*--- Evaluate Omega ---*/
  
  Omega = sqrt(Vorticity_i[0]*Vorticity_i[0] + Vorticity_i[1]*Vorticity_i[1] + Vorticity_i[2]*Vorticity_i[2]);

  /*--- Rotational correction term ---*/
  
  if (rotating_frame) { Omega += 2.0*min(0.0, StrainMag_i-Omega); }
  
  if (dist_i > 1e-10) {
    
    if (TurbVar_i[0] > 0.0) {
      
      /*--- Production term ---*/
      
      dist_i_2 = dist_i*dist_i;
      nu = Laminar_Viscosity_i/Density_i;
      Ji = TurbVar_i[0]/nu;
      Ji_2 = Ji*Ji;
      Ji_3 = Ji_2*Ji;
      fv1 = Ji_3/(Ji_3+cv1_3);
      fv2 = 1.0 - Ji/(1.0+Ji*fv1);
      ft2 = ct3*exp(-ct4*Ji_2);
      S = Omega;
      inv_k2_d2 = 1.0/(k2*dist_i_2);
      
      Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
      Shat = max(Shat, 1.0e-10);
      inv_Shat = 1.0/Shat;
      
      /*--- Production term ---*/;
      
      //    Original SA model
      //    Production = cb1*(1.0-ft2)*Shat*TurbVar_i[0]*Volume;
      
      Production = cb1*Shat*TurbVar_i[0]*Volume;
      
      /*--- Destruction term ---*/
      
      r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
      g = r + cw2*(pow(r,6.0)-r);
      g_6 =	pow(g,6.0);
      glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
      fw = g*glim;
      
      //    Original SA model
      //    Destruction = (cw1*fw-cb1*ft2/k2)*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      Destruction = cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      /*--- Diffusion term ---*/
      
      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
      
      CrossProduction = cb2_sigma*norm2_Grad*Volume;
      
      val_residual[0] = Production - Destruction + CrossProduction;
      
      /*--- Implicit part, production term ---*/
      
      dfv1 = 3.0*Ji_2*cv1_3/(nu*pow(Ji_3+cv1_3,2.));
      dfv2 = -(1/nu-Ji_2*dfv1)/pow(1.+Ji*fv1,2.);
      if ( Shat <= 1.0e-10 ) dShat = 0.0;
      else dShat = (fv2+TurbVar_i[0]*dfv2)*inv_k2_d2;
      val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;
      
      /*--- Implicit part, destruction term ---*/
      
      dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
      if (r == 10.0) dr = 0.0;
      dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
      dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
      val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] +	2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;
      
    }
    
    else {
      
      /*--- Production term ---*/
      
      dist_i_2 = dist_i*dist_i;
      
      /*--- Production term ---*/;
      
      Production = cb1*(1.0-ct3)*Omega*TurbVar_i[0]*Volume;
      
      /*--- Destruction term ---*/
      
      Destruction = cw1*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
      
      /*--- Diffusion term ---*/
      
      norm2_Grad = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
      
      CrossProduction = cb2_sigma*norm2_Grad*Volume;
      
      val_residual[0] = Production + Destruction + CrossProduction;
      
      /*--- Implicit part, production term ---*/
      
      val_Jacobian_i[0][0] += cb1*(1.0-ct3)*Omega*Volume;
      
      /*--- Implicit part, destruction term ---*/
      
      val_Jacobian_i[0][0] += 2.0*cw1*TurbVar_i[0]/dist_i_2*Volume;
      
    }
    
  }

//  AD::SetPreaccOut(val_residual, nVar);
//  AD::EndPreacc();
}

CUpwSca_TurbSST::CUpwSca_TurbSST(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbSST::~CUpwSca_TurbSST(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbSST::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(TurbVar_i,2);
  AD::SetPreaccIn(TurbVar_j,2);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement){
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+2);
    AD::SetPreaccIn(V_j, nDim+2);

    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+3);
    AD::SetPreaccIn(V_j, nDim+3);

    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }
  
  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  
  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;		val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;		val_Jacobian_i[1][1] = a0;
    
    val_Jacobian_j[0][0] = a1;		val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[1][0] = 0.0;		val_Jacobian_j[1][1] = a1;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CAvgGrad_TurbSST::CAvgGrad_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbSST::~CAvgGrad_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbSST::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}


CAvgGradCorrected_TurbSST::CAvgGradCorrected_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k1  = constants[0];
  sigma_om1 = constants[2];
  sigma_k2  = constants[1];
  sigma_om2 = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbSST::~CAvgGradCorrected_TurbSST(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbSST::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_omega_i, sigma_omega_j;
  su2double diff_i_kine, diff_i_omega, diff_j_kine, diff_j_omega;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j ,nVar);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F1_j);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i  = F1_i*sigma_k1 + (1.0 - F1_i)*sigma_k2;
  sigma_kine_j  = F1_j*sigma_k1 + (1.0 - F1_j)*sigma_k2;
  sigma_omega_i = F1_i*sigma_om1 + (1.0 - F1_i)*sigma_om2;
  sigma_omega_j = F1_j*sigma_om1 + (1.0 - F1_j)*sigma_om2;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine  = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine  = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_omega = Laminar_Viscosity_i + sigma_omega_i*Eddy_Viscosity_i;
  diff_j_omega = Laminar_Viscosity_j + sigma_omega_j*Eddy_Viscosity_j;
  
  diff_kine  = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_omega = 0.5*(diff_i_omega + diff_j_omega);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0; proj_vector_ij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_omega*Proj_Mean_GradTurbVar_Corrected[1];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;		Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;									    Jacobian_i[1][1] = -diff_omega*proj_vector_ij/Density_i;
    
    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 		Jacobian_j[0][1] = 0.0;
    Jacobian_j[1][0] = 0.0;									    Jacobian_j[1][1] = diff_omega*proj_vector_ij/Density_j;
  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbSST::CSourcePieceWise_TurbSST(unsigned short val_nDim, unsigned short val_nVar, su2double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Closure constants ---*/
  beta_star     = constants[6];
  sigma_omega_1 = constants[2];
  sigma_omega_2 = constants[3];
  beta_1        = constants[4];
  beta_2        = constants[5];
  alfa_1        = constants[8];
  alfa_2        = constants[9];
  a1            = constants[7];
}

CSourcePieceWise_TurbSST::~CSourcePieceWise_TurbSST(void) { }

void CSourcePieceWise_TurbSST::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(F1_i); AD::SetPreaccIn(F2_i); AD::SetPreaccIn(CDkw_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pw, zeta;
  
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5);

    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7);

    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }
  
  val_residual[0] = 0.0;        val_residual[1] = 0.0;
  val_Jacobian_i[0][0] = 0.0;	val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[1][0] = 0.0;	val_Jacobian_i[1][1] = 0.0;
  
  //  cout<<" F1_i: "<<F1_i<<"\n";

  /*--- Computation of blended constants for the source terms---*/
  
  alfa_blended = F1_i*alfa_1 + (1.0 - F1_i)*alfa_2;
  beta_blended = F1_i*beta_1 + (1.0 - F1_i)*beta_2;
  
  if (dist_i > 1e-10) {
    
    /*--- Production ---*/
    diverg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      diverg += PrimVar_Grad_i[iDim+1][iDim];
    
    pk = Eddy_Viscosity_i*StrainMag_i*StrainMag_i - 2.0/3.0*Density_i*TurbVar_i[0]*diverg;
    pk = min(pk,20.0*beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]);
    pk = max(pk,0.0);

    zeta = max(TurbVar_i[1], StrainMag_i*F2_i/a1);
    pw = StrainMag_i*StrainMag_i - 2.0/3.0*zeta*diverg;
    pw = max(pw,0.0);
    
    val_residual[0] += pk*Volume;
    val_residual[1] += alfa_blended*Density_i*pw*Volume;
    
    /*--- Dissipation ---*/
    val_residual[0] -= beta_star*Density_i*TurbVar_i[1]*TurbVar_i[0]*Volume;
    val_residual[1] -= beta_blended*Density_i*TurbVar_i[1]*TurbVar_i[1]*Volume;
    
    /*--- Cross diffusion ---*/
    val_residual[1] += (1.0 - F1_i)*CDkw_i*Volume;
    
    /*--- Implicit part ---*/

    // original:
    /*
    val_Jacobian_i[0][0] = -beta_star*TurbVar_i[1]*Volume;    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][0] = 0.0;                               val_Jacobian_i[1][1] = -2.0*beta_blended*TurbVar_i[1]*Volume;
    */


    // swh:
    val_Jacobian_i[0][0] = -beta_star*Density_i*TurbVar_i[1]*Volume;
    val_Jacobian_i[0][1] = -beta_star*Density_i*TurbVar_i[0]*Volume;
    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[1][1] = -2.0*beta_blended*Density_i*TurbVar_i[1]*Volume;

  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}



//swh
CUpwSca_TurbKE::CUpwSca_TurbKE(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbKE::~CUpwSca_TurbKE(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbKE::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+3);
  AD::SetPreaccIn(V_j, nDim+3);
  AD::SetPreaccIn(TurbVar_i,nVar);
  AD::SetPreaccIn(TurbVar_j,nVar);
  AD::SetPreaccIn(Normal, nDim);

  if (incompressible) {
    Density_i = V_i[nDim+1];
    Density_j = V_j[nDim+1];
  }
  else {
    Density_i = V_i[nDim+2];
    Density_j = V_j[nDim+2];
  }
  
  q_ij = 0.0;
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  
  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  val_residual[2] = a0*Density_i*TurbVar_i[2]+a1*Density_j*TurbVar_j[2];
  val_residual[3] = 0.0; // no convection in f scalar


  if (implicit) {
    val_Jacobian_i[0][0] = a0;   val_Jacobian_i[0][1] = 0.0;  val_Jacobian_i[0][2] = 0.0;  val_Jacobian_i[0][3] = 0.0;
    val_Jacobian_i[1][0] = 0.0;	 val_Jacobian_i[1][1] = a0;   val_Jacobian_i[1][2] = 0.0;  val_Jacobian_i[1][3] = 0.0;
    val_Jacobian_i[2][0] = 0.0;	 val_Jacobian_i[2][1] = 0.0;  val_Jacobian_i[2][2] = a0;   val_Jacobian_i[2][3] = 0.0;
    val_Jacobian_i[3][0] = 0.0;	 val_Jacobian_i[3][1] = 0.0;  val_Jacobian_i[3][2] = 0.0;  val_Jacobian_i[3][3] = 0.0;

    val_Jacobian_j[0][0] = a1;   val_Jacobian_j[0][1] = 0.0;  val_Jacobian_j[0][2] = 0.0;  val_Jacobian_j[0][3] = 0.0;
    val_Jacobian_j[1][0] = 0.0;	 val_Jacobian_j[1][1] = a1;   val_Jacobian_j[1][2] = 0.0;  val_Jacobian_j[1][3] = 0.0;
    val_Jacobian_j[2][0] = 0.0;	 val_Jacobian_j[2][1] = 0.0;  val_Jacobian_j[2][2] = a1;   val_Jacobian_j[2][3] = 0.0;
    val_Jacobian_j[3][0] = 0.0;	 val_Jacobian_j[3][1] = 0.0;  val_Jacobian_j[3][2] = 0.0;  val_Jacobian_j[3][3] = 0.0;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}

CAvgGrad_TurbKE::CAvgGrad_TurbKE(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbKE::~CAvgGrad_TurbKE(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
  delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbKE::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_epsi_i, sigma_epsi_j;
  su2double sigma_zeta_i, sigma_zeta_j;
  su2double diff_i_kine, diff_j_kine;
  su2double diff_i_epsi, diff_j_epsi;
  su2double diff_i_zeta, diff_j_zeta;
  su2double diff_i_f, diff_j_f;
  //  su2double Lm_i, Lm_j;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(Lm_i); AD::SetPreaccIn(Lm_j);

  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  // there are already stored as inverses
  sigma_kine_i = sigma_k;
  sigma_kine_j = sigma_k;
  sigma_epsi_i = sigma_e;
  sigma_epsi_j = sigma_e;
  sigma_zeta_i = sigma_z;
  sigma_zeta_j = sigma_z;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_epsi = Laminar_Viscosity_i + sigma_epsi_i*Eddy_Viscosity_i;
  diff_j_epsi = Laminar_Viscosity_j + sigma_epsi_j*Eddy_Viscosity_j;
  diff_i_zeta = Laminar_Viscosity_i + sigma_zeta_i*Eddy_Viscosity_i;
  diff_j_zeta = Laminar_Viscosity_j + sigma_zeta_j*Eddy_Viscosity_j;
  
  diff_kine = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_epsi = 0.5*(diff_i_epsi + diff_j_epsi);
  diff_zeta = 0.5*(diff_i_zeta + diff_j_zeta);
  //  diff_f = Lm_i*Lm_i; //here
  diff_f = Density_i;
  //  cout << "lm_i:" << Lm_i << "\n";
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_epsi*Proj_Mean_GradTurbVar_Corrected[1];
  val_residual[2] = diff_zeta*Proj_Mean_GradTurbVar_Corrected[2];
  val_residual[3] = diff_f*Proj_Mean_GradTurbVar_Corrected[3];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/ //here
  if (implicit) {
    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i; 
    Jacobian_i[0][1] = 0.0; 
    Jacobian_i[0][2] = 0.0; 
    Jacobian_i[0][3] = 0.0;

    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -diff_epsi*proj_vector_ij/Density_i;
    Jacobian_i[1][2] = 0.0;
    Jacobian_i[1][3] = 0.0;

    Jacobian_i[2][0] = 0.0;
    Jacobian_i[2][1] = 0.0;
    Jacobian_i[2][2] = -diff_zeta*proj_vector_ij/Density_i;
    Jacobian_i[2][3] = 0.0;
    
    Jacobian_i[3][0] = 0.0;
    Jacobian_i[3][1] = 0.0;
    Jacobian_i[3][2] = 0.0;
    Jacobian_i[3][3] = -diff_f*proj_vector_ij/Density_i;


    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 
    Jacobian_j[0][1] = 0.0; 
    Jacobian_j[0][2] = 0.0; 
    Jacobian_j[0][3] = 0.0;

    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_epsi*proj_vector_ij/Density_j;
    Jacobian_j[1][2] = 0.0;
    Jacobian_j[1][3] = 0.0;

    Jacobian_j[2][0] = 0.0;
    Jacobian_j[2][1] = 0.0;
    Jacobian_j[2][2] = diff_zeta*proj_vector_ij/Density_j;
    Jacobian_j[2][3] = 0.0;
    
    Jacobian_j[3][0] = 0.0;
    Jacobian_j[3][1] = 0.0;
    Jacobian_j[3][2] = 0.0;
    Jacobian_j[3][3] = diff_f*proj_vector_ij/Density_i;

  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
  
}


CAvgGradCorrected_TurbKE::CAvgGradCorrected_TurbKE(unsigned short val_nDim, unsigned short val_nVar, su2double *constants, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Normal = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbKE::~CAvgGradCorrected_TurbKE(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Normal;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbKE::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  su2double sigma_kine_i, sigma_kine_j, sigma_epsi_i, sigma_epsi_j, sigma_zeta_i, sigma_zeta_j;
  su2double diff_i_kine, diff_i_epsi, diff_j_kine, diff_j_epsi;
  su2double diff_i_zeta, diff_j_zeta, diff_i_f, diff_j_f;
  //  su2double Lm_i, Lm_j;
  
  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim); AD::SetPreaccIn(TurbVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(TurbVar_i, nVar); AD::SetPreaccIn(TurbVar_j ,nVar);
  AD::SetPreaccIn(Lm_i); AD::SetPreaccIn(Lm_j);

  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i = sigma_k;
  sigma_kine_j = sigma_k;
  sigma_epsi_i = sigma_e;
  sigma_epsi_j = sigma_e;
  sigma_zeta_i = sigma_z;
  sigma_zeta_j = sigma_z;
  
  /*--- Compute mean effective viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_epsi = Laminar_Viscosity_i + sigma_epsi_i*Eddy_Viscosity_i;
  diff_j_epsi = Laminar_Viscosity_j + sigma_epsi_j*Eddy_Viscosity_j;
  diff_i_zeta = Laminar_Viscosity_i + sigma_zeta_i*Eddy_Viscosity_i;
  diff_j_zeta = Laminar_Viscosity_j + sigma_zeta_j*Eddy_Viscosity_j;
  
  diff_kine = 0.5*(diff_i_kine + diff_j_kine);    // Could instead use weighted average!
  diff_epsi = 0.5*(diff_i_epsi + diff_j_epsi);
  diff_zeta = 0.5*(diff_i_zeta + diff_j_zeta);
  //  diff_f = Lm_i*Lm_i; //here
  diff_f = Density_i;
  //  cout << "lm_i:" << Lm_i << "\n";
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0.0; proj_vector_ij = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Normal[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Normal[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Normal[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar_Corrected[0];
  val_residual[1] = diff_epsi*Proj_Mean_GradTurbVar_Corrected[1];
  val_residual[2] = diff_zeta*Proj_Mean_GradTurbVar_Corrected[2];
  val_residual[3] = diff_f*Proj_Mean_GradTurbVar_Corrected[3];
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {

    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i; 
    Jacobian_i[0][1] = 0.0; 
    Jacobian_i[0][2] = 0.0; 
    Jacobian_i[0][3] = 0.0;

    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -diff_epsi*proj_vector_ij/Density_i;
    Jacobian_i[1][2] = 0.0;
    Jacobian_i[1][3] = 0.0;

    Jacobian_i[2][0] = 0.0;
    Jacobian_i[2][1] = 0.0;
    Jacobian_i[2][2] = -diff_zeta*proj_vector_ij/Density_i;
    Jacobian_i[2][3] = 0.0;
    
    Jacobian_i[3][0] = 0.0;
    Jacobian_i[3][1] = 0.0;
    Jacobian_i[3][2] = 0.0;
    Jacobian_i[3][3] = -diff_f*proj_vector_ij/Density_i;


    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j; 
    Jacobian_j[0][1] = 0.0; 
    Jacobian_j[0][2] = 0.0; 
    Jacobian_j[0][3] = 0.0;

    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_epsi*proj_vector_ij/Density_j;
    Jacobian_j[1][2] = 0.0;
    Jacobian_j[1][3] = 0.0;

    Jacobian_j[2][0] = 0.0;
    Jacobian_j[2][1] = 0.0;
    Jacobian_j[2][2] = diff_zeta*proj_vector_ij/Density_j;
    Jacobian_j[2][3] = 0.0;
    
    Jacobian_j[3][0] = 0.0;
    Jacobian_j[3][1] = 0.0;
    Jacobian_j[3][2] = 0.0;
    Jacobian_j[3][3] = diff_f*proj_vector_ij/Density_i;

  }
  
  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}

CSourcePieceWise_TurbKE::CSourcePieceWise_TurbKE(unsigned short val_nDim, unsigned short val_nVar, su2double *constants,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  /*--- Closure constants ---*/
  /*sigma_k = constants[0];
  sigma_e = constants[1];
  C_mu    = constants[2];
  C_e1    = constants[3];
  C_e2    = constants[4];*/

  /* RNG
  C_mu    = constants[0];
  sigma_k = constants[1];
  sigma_e = constants[2];
  C_e1    = constants[4];
  C_e2    = constants[5];
  eta_o   = constants[6];
  beta    = constants[7];
  */

  /* zeta f */
  C_mu    = constants[0];
  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
  C_e1o   = constants[4];
  C_e2    = constants[5];
  C_1     = constants[6];
  C_2p    = constants[7];
  C_T     = constants[8];
  C_L     = constants[9];
  C_eta   = constants[10];

}

CSourcePieceWise_TurbKE::~CSourcePieceWise_TurbKE(void) { }

void CSourcePieceWise_TurbKE::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+7);
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);
  //  AD::SetPreaccIn(Lm_i); AD::SetPreaccIn(Tm_i);

  unsigned short iDim;
  su2double alfa_blended, beta_blended;
  su2double diverg, pk, pe, pz, pf, pv2;
  su2double dk, de, dz, df, dv2, S, Vol;
  su2double T1, T2, T3, L1, L2, L3, R1, R2, R3, R;
  su2double yplus, C_e2star, eta, C_2f, C_e1;
  su2double tke, tdr, zeta, v2, f, L, T, mu, nu, nuT, muT, rho;
  su2double tke_d, tdr_d, zeta_d, v2_d, f_d;
  su2double dTdrk, dTdre, dTdrz, dLdrk, dLdre, dLdrz;
  su2double dTdk, dTde, dTdz, dLdk, dLde, dLdz, dCe1dz, dPedT, dPedC, dDzdz, dDzdT,
    dDede, dDedT, dPzdf, dDzdk, dPfdT, dPfdL, dPfdz, dPfde, dDfdL, dDfdf, dDkde;
  su2double dTdrv2, dLdrv2, dTdv2, dLdv2, dCe1dv2, dDv2dv2, dDv2dT, dPv2df, dDv2dk, 
    dPfdv2, dCe1dk, dPv2dk, dPfdk, dDv2de, dPede, dPkde;


  //yplus = CNSSolver.YPlus;
  
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  // for readability...
  tke  = TurbVar_i[0]; ///Density_i;
  tdr  = TurbVar_i[1]; ///Density_i;
  //  zeta = TurbVar_i[2];
  v2   = TurbVar_i[2];
  f    = TurbVar_i[3];
  mu   = Laminar_Viscosity_i;
  muT  = Eddy_Viscosity_i;
  rho  = Density_i;
  S    = StrainMag_i; //*sqrt(2.0) already included
  Vol  = Volume;
  //v2 = zeta*tke;
  //  zeta = v2/tke;
  nu = mu/rho;
  nuT = muT/rho;

  //  if (tke>=20.0) {cout << "TKE: " << tke << "\n";}

  // limits?
  //  tke  = max(tke,0.0);
  //  tdr  = max(tdr,0.0);
  //  zeta = min(zeta,2.0/3.0);
  //  zeta = max(zeta,0.0);
  //  v2 = min(v2,2.0/3.0*tke);
  //  v2 = max(v2,0.0);
  //  f  = max(f,0.0);

  // from previous timestep
  //  tke0 = solver_container[TURB_SOL]->node[iPoint]->GetSolution_Old(6);
  //  tdr0 = solver_container[TURB_SOL]->node[iPoint]->GetSolution_Old(7);
  //  zeta0 = solver_container[TURB_SOL]->node[iPoint]->GetSolution_Old(8);
  //  f0 = solver_container[TURB_SOL]->node[iPoint]->GetSolution_Old(9);

  // denominator floors <jump>
  //  su2double scalar_min = 1.0E-14;
  su2double VelMag, *VelInf, L_Inf, scalar_min;
  su2double tke_raw, tdr_raw;
  tke_raw = tke;
  tdr_raw = tdr;
  VelInf = config->GetVelocity_FreeStreamND();
  L_Inf = config->GetLength_Reynolds();
  for (iDim = 0; iDim < nDim; iDim++)
  VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);
  //  su2double solve_tol = config->GetLinear_Solver_Error()
  scalar_min = 1.0E-8/(VelMag*VelMag); // setting based on tke min being 1e-8
  //  scalar_min = 1.0E-12;
  tke = max(tke, scalar_min*VelMag*VelMag);
  tdr = max(tdr, scalar_min*VelMag*VelMag*VelMag/L_Inf);
  v2 = max(v2, 2.0/3.0*scalar_min*VelMag*VelMag);
  f = max(f, scalar_min*VelMag/L_Inf);
  zeta = v2/tke;
  //  S = max(S,scalar_min*VelMag/L_Inf); // no checked...
  S = max(S,1.0E-14);
  /*
  tke_d = max(tke,1.0E-8);
  tdr_d = max(tdr,1.0E-8);
  zeta_d = max(zeta,1.0E-8);
  v2_d = max(v2,1.0E-8);
  f_d = max(f,1.0E-8);
  */

  // must find T&L here due to Jacobian branching...
  //L = Lm_i;
  //T = Tm_i;

  //--- Model time scale ---//
  T1 = tke/tdr;
  T2 = 0.6/(sqrt(6.0)*C_mu*S*zeta);
  T3 = C_T*sqrt(nu/tdr);
  T = max(min(T1,T2),T3); 
  T = max(T,1.0E-8);
  //  L = min(L,10.0);

  //--- Model rate ---//
  //  R1 = max(tdr,0.0)/tke_d;
  //  R2 = (sqrt(6.0)*C_mu*S*zeta)/0.6;
  //  R3 = 1.0/C_T*sqrt(max(tdr,0.0)/nu);
  //  R = min(max(R1,R2),R3); 

  //--- Model length scale ---//
  L1 = pow(tke,1.5)/tdr;
  L2 = sqrt(tke)/(sqrt(6.0)*C_mu*S*zeta);
  L3 = C_eta*pow(pow(nu,3.0)/tdr,0.25);
  L = C_L * max(min(L1,L2),L3);
  L = max(L,1.0E-6);
  //  L = min(L,1.0E6);
  //  cout << "L:" << L << "\n";

  //--- Initial Jacobian ---//
  val_Jacobian_i[0][0] = 0.0;
  val_Jacobian_i[0][1] = 0.0;
  val_Jacobian_i[0][2] = 0.0;
  val_Jacobian_i[0][3] = 0.0;

  val_Jacobian_i[1][0] = 0.0;
  val_Jacobian_i[1][1] = 0.0;
  val_Jacobian_i[1][2] = 0.0;
  val_Jacobian_i[1][3] = 0.0;

  val_Jacobian_i[2][0] = 0.0;
  val_Jacobian_i[2][1] = 0.0;
  val_Jacobian_i[2][2] = 0.0;
  val_Jacobian_i[2][3] = 0.0;

  val_Jacobian_i[3][0] = 0.0;
  val_Jacobian_i[3][1] = 0.0;
  val_Jacobian_i[3][2] = 0.0;
  val_Jacobian_i[3][3] = 0.0;


  /*
    //--- zeta-f ---//
    C_e1 = C_e1o*(1.0+0.012/zeta_d); // error in paper?
    //  C_e1 = C_e1o*(1.0+0.012/sqrt(zeta_d));
    //  C_e1 = C_e1o*(1.0+0.012/sqrt(2.0/3.0));

    //--- divergence of velocity ---//
    diverg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      diverg += PrimVar_Grad_i[iDim+1][iDim];
 
    //--- Production ---// //<warp>//
    //pk = muT*(S*S - 2.0/3.0*diverg*diverg) - 2.0/3.0*rho*tke*diverg;
    pk = muT*S*S - 2.0/3.0*rho*tke*diverg;
    pe = C_e1*pk/T;
    pz = rho*f;
    pf = (C_1-1.0 + C_2p*pk/(rho*tdr_d)) * (2.0/3.0-zeta) * rho/T * 1.0/(L*L);

    pk = max(pk,0.0);
    pe = max(pe,0.0);
    pz = max(pz,0.0);
    pf = max(pf,0.0);

    //--- Dissipation ---//
    dk = rho*tdr;
    de = C_e2*rho*tdr/T;
    dz = min(zeta,2.0/3.0)/tke_d*pk;
    df = rho*f/(L*L);

    //-- Store in residual --//
    val_residual[0] = (pk-dk) * Vol;
    val_residual[1] = (pe-de) * Vol;
    val_residual[2] = (pz-dz) * Vol; // I HAZ STOOPID?
    val_residual[3] = (pf-df) * Vol;

        
    //--- Implicit part ---//

    // precompute T&L portions of the jacobian
    if (T==T3) {

      dTdk = 0.0;
      dTde = -0.5*C_T*sqrt(nu)*pow(tdr_d,-1.5);
      dTdz = 0.0;
      dLdk = 0.0;
      dLde = -C_L*0.25*C_eta*pow(nu,0.75)*pow(tdr_d,-5.0/4.0);
      dLdz = 0.0;

      dTdrk = 0.0;
      dTdre = -0.5*C_T*sqrt(nu)*pow(rho*tdr_d,-1.5) * sqrt(rho) ;
      dTdrz = 0.0;
      dLdrk = 0.0;
      dLdre = -C_L*0.25*C_eta*pow(nu,0.75)*pow(rho*tdr_d,-5.0/4.0) * pow(rho,0.25);
      dLdrz = 0.0;

    }
    else if (T==T2) {

      dTdk = 0.0;
      dTde = 0.0;
      dTdz = -0.6/(sqrt(6.0)*C_mu*S*zeta_d*zeta_d);
      dLdk = C_L*0.5/(sqrt(6.0)*C_mu*S*zeta_d*sqrt(tke_d));
      dLde = -C_L*sqrt(tke_d)/(sqrt(6.0)*C_mu*S*zeta_d*zeta_d);
      dLdz = 0.0;

      dTdrk = 0.0;
      dTdre = 0.0;
      dTdrz = -0.6/(sqrt(6.0)*C_mu*S*zeta_d*zeta_d*rho);
      dLdrk = C_L*0.5/(sqrt(6.0)*C_mu*S*zeta_d*sqrt(rho*tke_d))*1.0/sqrt(rho) ;
      dLdre = -C_L*sqrt(tke_d)/(sqrt(6.0)*C_mu*S*zeta_d*zeta_d*rho);
      dLdrz = 0.0;

    }
    else {

      dTdk = 1.0/tdr_d;
      dTde = -tke/(tdr_d*tdr_d);
      dTdz = 0.0;
      dLdk = C_L*1.5*sqrt(tke_d)/tdr_d;
      dLde = -C_L*pow(tke_d,1.5)/(tdr_d*tdr_d);
      dLdz = 0.0;

      dTdrk = 1.0/(rho*tdr_d);
      dTdre = -tke/(tdr_d*tdr_d*rho);
      dTdrz = 0.0;
      dLdrk = C_L*1.5*sqrt(rho*tke_d)/tdr_d*1.0/pow(rho,1.3);
      dLdre = -C_L*pow(tke_d,1.5)/(tdr_d*tdr_d*rho);
      dLdrz = 0.0;

    }

    // other jacobian portions
    dDkde = 1.0; //rho;
    dCe1dz = -C_e1o*0.012/(zeta_d*zeta_d*rho); //-C_e1o*0.012/(zeta_d*zeta_d);
    //    dCe1dz = -C_e1o*0.012/(pow(zeta_d,1.5)*rho);
    dPedT = -C_e1*pk/(T*T);
    dPedC = pk/T;
    dDede = C_e2/T; //C_e2*rho/T;
    dDedT = -C_e2*rho*tdr/(T*T);
    dPzdf = 1.0; //rho;

    dDzdk = -zeta*pk/(tke_d*tke_d*rho); //-zeta*pk/(tke_d*tke_d);
    dDzdz = pk/(rho*tke_d); //pk/tke_d;

    dPfdT = -(C_1-1.0+C_2p*pk/(rho*tdr_d)) * (2.0/3.0-zeta) * rho/(T*T) * 1.0/(L*L);
    dPfdL = -(C_1-1.0+C_2p*pk/(rho*tdr_d)) * (2.0/3.0-zeta) * rho/T * 2.0/(L*L*L);
    dPfde = -C_2p*pk/(rho*tdr_d*rho*tdr_d) * (2.0/3.0-zeta) * rho/T * 1.0/(L*L); //-C_2p*pk/(rho*tdr_d*tdr_d) * (2.0/3.0-zeta) * 1.0/T * 1.0/(L*L);
    dPfdz = -(C_1-1.0+C_2p*pk/(rho*tdr_d)) * 1.0/T * 1.0/(L*L);
    dDfdL = -2.0*rho*f/(L*L*L);
    dDfdf = 1.0/(L*L);

  */
    // production...
    /*
    val_Jacobian_i[0][0] += 0.0;
    val_Jacobian_i[0][1] += 0.0;
    val_Jacobian_i[0][2] += 0.0;
    val_Jacobian_i[0][3] += 0.0;

    val_Jacobian_i[1][0] += (dPedT*dTdrk) * Vol;
    val_Jacobian_i[1][1] += (dPedT*dTdre) * Vol;
    val_Jacobian_i[1][2] += (dPedC*dCe1dz + dPedT*dTdrz) * Vol;
    val_Jacobian_i[1][3] += 0.0;

    val_Jacobian_i[2][0] += 0.0;
    val_Jacobian_i[2][1] += 0.0;
    val_Jacobian_i[2][2] += 0.0;
    val_Jacobian_i[2][3] += dPzdf*Vol;

    val_Jacobian_i[3][0] += (dPfdT*dTdrk + dPfdL*dLdrk) * Vol;
    val_Jacobian_i[3][1] += (dPfdT*dTdre + dPfdL*dLdre + dPfde) * Vol;
    val_Jacobian_i[3][2] += (dPfdT*dTdrz + dPfdL*dLdrz + dPfdz) * Vol;
    val_Jacobian_i[3][3] += 0.0;
    */

    /*
    // destruction...
    val_Jacobian_i[0][0] -= 0.0;
    val_Jacobian_i[0][1] -= dDkde * Vol;
    val_Jacobian_i[0][2] -= 0.0;
    val_Jacobian_i[0][3] -= 0.0;

    val_Jacobian_i[1][0] -= dDedT*dTdrk * Vol;
    val_Jacobian_i[1][1] -= (dDedT*dTdre + dDede) * Vol;
    val_Jacobian_i[1][2] -= dDedT*dTdrz * Vol;
    val_Jacobian_i[1][3] -= 0.0;

    val_Jacobian_i[2][0] -= dDzdk * Vol;
    val_Jacobian_i[2][1] -= 0.0;
    val_Jacobian_i[2][2] -= dDzdz * Vol;
    val_Jacobian_i[2][3] -= 0.0;

    val_Jacobian_i[3][0] -= dDfdL*dLdrk * Vol;
    val_Jacobian_i[3][1] -= dDfdL*dLdre * Vol;
    val_Jacobian_i[3][2] -= dDfdL*dLdrz * Vol;
    val_Jacobian_i[3][3] -= dDfdf * Vol;
  */

  /*
  cout << "C_2p: " << C_2p << "\n";
  cout << "C_1: " << C_1 << "\n";
  cout << "C_e1o: " << C_e1o << "\n";
  cout << "C_e2: " << C_e2 << "\n";
  cout << "C_T: " << C_T << "\n";
  cout << "C_L: " << C_L << "\n";
  cout << "C_eta: " << C_eta << "\n";
  cout << "C_mu: " << C_mu << "\n";
  */

    //--- v2-f ---//
    //    C_2f = C_2p;
    //    f = (C_2f*pk - ((C_1-6.0)*v2 - 2.0/3.0*(C_1-1.0)*tke)*rho/T)/tke_d;
    //    su2double kf = (C_2f*pk - ((C_1-6.0)*v2 - 2.0/3.0*(C_1-1.0)*tke)*rho/T);
    //    kf = max(kf,0.0);

    C_e1 = C_e1o*(1.0+0.045*sqrt(tke/v2));

    //--- divergence of velocity ---//
    diverg = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      diverg += PrimVar_Grad_i[iDim+1][iDim];
 
    //--- Production ---// //<warp>//
    //pk = muT*(S*S - 2.0/3.0*diverg*diverg) - 2.0/3.0*rho*tke*diverg;
    pk = muT*S*S - 2.0/3.0*rho*max(tke_raw,0.0)*diverg;
    //    pk = 0.0;
    pk = max(pk,0.0);
    pe = C_e1*pk/T;
    pv2 = rho*tke*f;
    pv2 = max(pv2,0.0);
    pv2 = min(pv2,2.0/3.0*pk+5.0*rho*v2/tke*tdr);
    //    pv2 = min(pv2,2.0/3.0*(pk+5.0*rho*tdr));
    //    pv2 = rho*f; //f=kf
    C_2f = C_2p + 0.5*(2.0/3.0-C_2p)*(1.0+tanh(50.0*(v2/tke-0.55)));
    //    C_2f = C_2p;
    pf = (C_2f*pk/tke - ((C_1-6.0)*v2/tke - 2.0/3.0*(C_1-1.0))*rho/T) * 1.0/(L*L); // jee C1=1.4
    //    pf = (C_2f*pk - (v2*(C_1-6.0) - 2.0/3.0*tke*(C_1-1.0))*rho/T) * 1.0/(L*L); //f=kf
    //    pf = (C_2f*pk + C_1*R*(2.0/3.0*tke-v2)*rho + 5.0*v2*R*rho)/tke_d * 1.0/(L*L); C1=0.4

    //    pe = max(pe,0.0);
    //    pf = max(pf,0.0);

    //--- Dissipation ---//
    dk = rho*max(tdr_raw,0.0);
    de = C_e2*rho*tdr_raw/T;
    dv2 = 6.0*(v2/tke)*rho*max(tdr_raw,0.0);
    df = rho*max(f,0.0)/(L*L);

    //-- Store in residual --//
    val_residual[0] = (pk-dk) * Vol;
    val_residual[1] = (pe-de) * Vol;
    val_residual[2] = (pv2-dv2) * Vol;
    val_residual[3] = (pf-df) * Vol;

        
    //--- Implicit part ---//

    // precompute T&L portions of the jacobian
    if (T==T3) {

      dTdk = 0.0;
      dTde = -0.5*C_T*sqrt(nu)*pow(tdr,-1.5);
      dTdv2 = 0.0;
      dLdk = 0.0;
      dLde = -C_L*0.25*C_eta*pow(nu,0.75)*pow(tdr,-5.0/4.0);
      dLdv2 = 0.0;

      dTdrk = 0.0;
      dTdre = -0.5*C_T*sqrt(nu)*pow(rho*tdr,-1.5) * sqrt(rho) ;
      dTdrv2 = 0.0;
      dLdrk = 0.0;
      dLdre = -C_L*0.25*C_eta*pow(nu,0.75)*pow(rho*tdr,-5.0/4.0) * pow(rho,0.25);
      dLdrv2 = 0.0;

    }
    else if (T==T2) {

      dTdk = 0.6/(sqrt(6.0)*C_mu*S*v2);
      dTde = 0.0;
      dTdv2 = -0.6*tke/(sqrt(6.0)*C_mu*S*v2*v2);
      dLdk = C_L*1.5*sqrt(tke)/(sqrt(6.0)*C_mu*S*v2);
      dLde = -C_L*pow(tke,1.5)/(sqrt(6.0)*C_mu*S*v2*v2);
      dLdv2 = 0.0;

      dTdrk = 0.6/(sqrt(6.0)*C_mu*S*v2*rho);
      dTdre = 0.0;
      dTdrv2 = -0.6/(sqrt(6.0)*C_mu*S*v2*v2*rho);
      dLdrk = C_L*1.5*sqrt(rho*tke)/(sqrt(6.0)*C_mu*S*v2)*1.0/pow(rho,1.5) ;
      dLdre = -C_L*pow(tke,1.5)/(sqrt(6.0)*C_mu*S*v2*v2*rho);
      dLdrv2 = 0.0;

    }
    else {

      dTdk = 1.0/tdr;
      dTde = -tke/(tdr*tdr);
      dTdv2 = 0.0;
      dLdk = C_L*1.5*sqrt(tke)/tdr;
      dLde = -C_L*pow(tke,1.5)/(tdr*tdr);
      dLdv2 = 0.0;

      dTdrk = 1.0/(rho*tdr);
      dTdre = -tke/(tdr*tdr*rho);
      dTdrv2 = 0.0;
      dLdrk = C_L*1.5*sqrt(rho*tke)/tdr*1.0/pow(rho,1.3);
      dLdre = -C_L*pow(tke,1.5)/(tdr*tdr*rho);
      dLdrv2 = 0.0;

    }

    // other jacobian portions
    dPkde = 0.0;
    dDkde = 1.0;
    if(tdr_raw<0.0) {
    //       dPkde = -1.0;
       dDkde = 0.0;
    }

    dCe1dk = C_e1o*0.045*0.5/sqrt(rho*tke*rho*v2);
    dCe1dv2 = -C_e1o*0.045*0.5*sqrt(tke/(v2*v2*rho));
    dPedT = -C_e1*pk/(T*T);
    dPedC = pk/T;
    dPede = 0.0;
    dDede = C_e2/T;
    if(tdr_raw<0.0) {
       dPede = -C_e2/T;
       dDede = 0.0;
    }
    dDedT = -C_e2*rho*tdr/(T*T);

    dPv2dk = f;//    dPv2dk = 0.0;
    dPv2df = tke;//    dPv2df = 1.0;
    dDv2dk = -6.0*v2/(tke*tke*rho)*rho*tdr;
    dDv2de = 6.0*v2/tke;
    dDv2dv2 = 6.0/(tke*rho)*rho*tdr;

    dPfdT = ((C_1-6.0)*v2/tke - 2.0/3.0*(C_1-1.0))*rho/(T*T) * 1.0/(L*L);
    dPfdL = -(C_2f*pk/tke - ((C_1-6.0)*v2/tke - 2.0/3.0*(C_1-1.0))*rho/T) * 2.0/(L*L*L);
    dPfdk = (-C_2f*pk/(tke*tke*rho) - (-(C_1-6.0)*v2/(tke*tke*rho))*rho/T) * 1.0/(L*L);
    dPfde = 0.0;
    dPfdv2 = -(C_1-6.0)*1.0/(rho*tke) * rho/T * 1.0/(L*L);
    dDfdL = -2.0*rho*f/(L*L*L);
    dDfdf = 1.0/(L*L);

    /*
    dPfdT = ((C_1-6.0)*v2 - 2.0/3.0*tke*(C_1-1.0))*rho/(T*T) * 1.0/(L*L);
    dPfdL = -(C_2f*pk - ((C_1-6.0)*v2 - 2.0/3.0*tke*(C_1-1.0))*rho/T) * 2.0/(L*L*L);
    dPfdk = 2.0/3.0*(C_1-1.0)/T * 1.0/(L*L);
    dPfde = 0.0;
    dPfdv2 = -(C_1-6.0)/T * 1.0/(L*L);
    dDfdL = -2.0*rho*f/(L*L*L);
    dDfdf = 1.0/(L*L);
    */

    // production...
    /*
    val_Jacobian_i[0][0] += 0.0;
    val_Jacobian_i[0][1] += 0.0;
    val_Jacobian_i[0][2] += 0.0;
    val_Jacobian_i[0][3] += 0.0;

    val_Jacobian_i[1][0] += 0.0; //(dPedC*dCe1dk + dPedT*dTdrk) * Vol;
    val_Jacobian_i[1][1] += 0.0; //(dPedT*dTdre) * Vol;
    val_Jacobian_i[1][2] += 0.0; //(dPedC*dCe1dv2 + dPedT*dTdrv2) * Vol;
    val_Jacobian_i[1][3] += 0.0;

    val_Jacobian_i[2][0] += dPv2dk*Vol;
    val_Jacobian_i[2][1] += 0.0;
    val_Jacobian_i[2][2] += 0.0;
    val_Jacobian_i[2][3] += dPv2df*Vol;

    val_Jacobian_i[3][0] += dPfdk * Vol  + (dPfdT*dTdrk + dPfdL*dLdrk + dPfdk) * Vol;
    val_Jacobian_i[3][1] += dPfde * Vol + (dPfdT*dTdre + dPfdL*dLdre + dPfde) * Vol;
    val_Jacobian_i[3][2] += dPfdv2 * Vol + (dPfdT*dTdrv2 + dPfdL*dLdrv2 + dPfdv2) * Vol;
    val_Jacobian_i[3][3] += 0.0;
    */
    //    val_Jacobian_i[0][1] += dPkde * Vol;
    val_Jacobian_i[1][1] += dPede * Vol;


    // destruction...
    val_Jacobian_i[0][0] -= 0.0;
    val_Jacobian_i[0][1] -= 0.0;// dDkde * Vol;
    val_Jacobian_i[0][2] -= 0.0;
    val_Jacobian_i[0][3] -= 0.0;

    val_Jacobian_i[1][0] -= 0.0; //dDedT*dTdrk * Vol;
    val_Jacobian_i[1][1] -= dDede * Vol; // + dDedT*dTdre * Vol;
    val_Jacobian_i[1][2] -= 0.0; //dDedT*dTdrv2 * Vol;
    val_Jacobian_i[1][3] -= 0.0;

    val_Jacobian_i[2][0] -= 0.0; //dDv2dk * Vol;
    val_Jacobian_i[2][1] -= 0.0; //dDv2de * Vol;
    val_Jacobian_i[2][2] -= dDv2dv2 * Vol;
    val_Jacobian_i[2][3] -= 0.0;

    val_Jacobian_i[3][0] -= 0.0; //dDfdL*dLdrk * Vol;
    val_Jacobian_i[3][1] -= 0.0; //dDfdL*dLdre * Vol;
    val_Jacobian_i[3][2] -= 0.0; //dDfdL*dLdrv2 * Vol;
    val_Jacobian_i[3][3] -= dDfdf * Vol;


  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}



CUpwSca_TurbML::CUpwSca_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                               CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwSca_TurbML::~CUpwSca_TurbML(void) {
  
  delete [] Velocity_i;
  delete [] Velocity_j;
  
}

void CUpwSca_TurbML::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  
  q_ij = 0.0;
  
  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TurbVar_i[0]+a1*TurbVar_j[0];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }
  
  
}

CAvgGrad_TurbML::CAvgGrad_TurbML(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGrad_TurbML::~CAvgGrad_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGrad_TurbML::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
    }
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Kappa[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Kappa[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}

CAvgGradCorrected_TurbML::CAvgGradCorrected_TurbML(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
  
}

CAvgGradCorrected_TurbML::~CAvgGradCorrected_TurbML(void) {
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
  
}

void CAvgGradCorrected_TurbML::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  if (incompressible) {
    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }
  
  /*--- Compute mean effective viscosity ---*/
  
  nu_i = Laminar_Viscosity_i/Density_i;
  nu_j = Laminar_Viscosity_j/Density_j;
  nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  
  /*--- Compute vector going from iPoint to jPoint ---*/
  
  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  proj_vector_ij = proj_vector_ij/dist_ij_2;
  
  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
    Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
    Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
    (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  }
  
  val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  
  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  
  if (implicit) {
    Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
    Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  }
  
}
