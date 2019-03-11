/*!
 * \file numerics_direct_transition.cpp
 * \brief This file contains all the convective term discretization.
 * \author A. Aranake
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

CUpwLin_TransLM::CUpwLin_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  grid_movement  = config->GetGrid_Movement();
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  
}

CUpwLin_TransLM::~CUpwLin_TransLM(void) {
  delete [] Velocity_i;
  delete [] Velocity_j;
}

void CUpwLin_TransLM::ComputeResidual (su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  
  Density_i = U_i[0];
  Density_j = U_j[0];
  
  q_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = U_i[iDim+1]/Density_i;
    Velocity_j[iDim] = U_j[iDim+1]/Density_j;
    q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
  val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];
  //  cout << "Velicity x: " << Velocity_i[0] << ", " << Velocity_j[0] << endl;
  //  cout << "Velicity y: " << Velocity_i[1] << ", " << Velocity_j[1] << endl;
  //  cout << "val_resid: " << val_residual[0] << ", " << val_residual[1] << endl;
  
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_i[1][1] = a0;
  }
}

CUpwSca_TransLM::CUpwSca_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
}

CUpwSca_TransLM::~CUpwSca_TransLM(void) {
  delete [] Velocity_i;
  delete [] Velocity_j;
}

void CUpwSca_TransLM::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {
  
  q_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    q_ij += 0.5*(U_i[iDim]+U_j[iDim])*Normal[iDim];
  }
  
  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*TransVar_i[0]+a1*TransVar_j[0];
  val_residual[1] = a0*TransVar_i[1]+a1*TransVar_j[1];
  
  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
    val_Jacobian_i[1][1] = a0;
    val_Jacobian_j[1][1] = a1;
    
    /*--- Zero out off-diagonal terms just in case ---*/
    val_Jacobian_i[0][1] = 0;
    val_Jacobian_j[0][1] = 0;
    val_Jacobian_i[1][0] = 0;
    val_Jacobian_j[1][0] = 0;
  }
  
}

CAvgGrad_TransLM::CAvgGrad_TransLM(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTransVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTransVar_Edge = new su2double [nVar];
  Mean_GradTransVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTransVar[iVar] = new su2double [nDim];
}

CAvgGrad_TransLM::~CAvgGrad_TransLM(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTransVar_Kappa;
  delete [] Proj_Mean_GradTransVar_Edge;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTransVar[iVar];
  delete [] Mean_GradTransVar;
}

void CAvgGrad_TransLM::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
 /*--- This section is commented out on 04/11/2016 
       after review of the static scan ---*/
 // su2double *Density_Grad_i      = new su2double[nDim];
 // su2double *Density_Grad_j      = new su2double[nDim];
 // su2double *Conservative_Grad_i = new su2double[nDim];
 // su2double *Conservative_Grad_j = new su2double[nDim];
 // su2double *Primitive_Grad_i    = new su2double[nDim];
 // su2double *Primitive_Grad_j    = new su2double[nDim];
 // 
 // /*--- Intermediate values for combining viscosities ---*/
 // su2double Inter_Viscosity_i, Inter_Viscosity_j, REth_Viscosity_i, REth_Viscosity_j, Inter_Viscosity_Mean, REth_Viscosity_Mean;
 // 
 // /*--- Model constants---*/
 // su2double sigmaf       = 1.0;
 // su2double sigma_thetat = 2.0;
 // 
 // /*--- Get density ---*/
 // Density_i = U_i[0];
 // Density_j = U_j[0];
 // 
 // /*--- Construct combinations of viscosity ---*/
 // Inter_Viscosity_i    = (Laminar_Viscosity_i+Eddy_Viscosity_i/sigmaf);
 // Inter_Viscosity_j    = (Laminar_Viscosity_j+Eddy_Viscosity_j/sigmaf);
 // Inter_Viscosity_Mean = 0.5*(Inter_Viscosity_i+Inter_Viscosity_j);
 // REth_Viscosity_i     = sigma_thetat*(Laminar_Viscosity_i+Eddy_Viscosity_i);
 // REth_Viscosity_j     = sigma_thetat*(Laminar_Viscosity_j+Eddy_Viscosity_j);
 // REth_Viscosity_Mean  = 0.5*(REth_Viscosity_i+REth_Viscosity_j);
 // 
  ///*--- Compute vector going from iPoint to jPoint ---*/
  //dist_ij_2 = 0; proj_vector_ij = 0;
  //for (iDim = 0; iDim < nDim; iDim++) {
  //  Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
  //  dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  //  proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  //}
  //proj_vector_ij = proj_vector_ij/dist_ij_2; // to normalize vectors
 // 
  ///*--- Mean gradient approximation ---*/
  //for (iVar = 0; iVar < nVar; iVar++) {
  //  Proj_Mean_GradTransVar_Kappa[iVar] = 0.0;
  //  // Proj_Mean_GradTransVar_Edge[iVar] = 0.0;
  //  for (iDim = 0; iDim < nDim; iDim++) {
 //     
 //     /* -- Compute primitive grad using chain rule -- */
 //     Density_Grad_i[iDim]      = ConsVar_Grad_i[0][iDim];
 //     Density_Grad_j[iDim]      = ConsVar_Grad_j[0][iDim];
 //     Conservative_Grad_i[iDim] = TransVar_Grad_i[iVar][iDim];
 //     Conservative_Grad_j[iDim] = TransVar_Grad_j[iVar][iDim];
 //     Primitive_Grad_i[iDim]    = 1./Density_i*(Conservative_Grad_i[iDim]-TransVar_i[iVar]*Density_Grad_i[iDim]);
 //     Primitive_Grad_j[iDim]    = 1./Density_j*(Conservative_Grad_j[iDim]-TransVar_j[iVar]*Density_Grad_j[iDim]);
 //     
 //     /*--- Compute the average primitive gradient and project it in the normal direction ---*/
 //     Mean_GradTransVar[iVar][iDim] = 0.5*(Primitive_Grad_i[iDim] + Primitive_Grad_j[iDim]);
  //    Proj_Mean_GradTransVar_Kappa[iVar] += Mean_GradTransVar[iVar][iDim]*Normal[iDim];
  //  }
  //}
 // 
  //val_residual[0] = Inter_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[0];
  //val_residual[1] = REth_Viscosity_Mean*Proj_Mean_GradTransVar_Kappa[1];
 // 
  ///*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  //if (implicit) {
  //  Jacobian_i[0][0] = (0.5*Proj_Mean_GradTransVar_Kappa[0]-Inter_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_j[0][0] = (0.5*Proj_Mean_GradTransVar_Kappa[0]+Inter_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_i[1][1] = (0.5*Proj_Mean_GradTransVar_Kappa[1]-REth_Viscosity_Mean*proj_vector_ij);
  //  Jacobian_j[1][1] = (0.5*Proj_Mean_GradTransVar_Kappa[1]+REth_Viscosity_Mean*proj_vector_ij);
  //}
 // 
 // /*--- Free locally allocated memory. For efficiency, these arrays
 //  should really be allocated/deallocated in the constructor/destructor. ---*/
 // delete [] Density_Grad_i;
 // delete [] Density_Grad_j;
 // delete [] Conservative_Grad_i;
 // delete [] Conservative_Grad_j;
 // delete [] Primitive_Grad_i;
 // delete [] Primitive_Grad_j;
 // 
}

CAvgGradCorrected_TransLM::CAvgGradCorrected_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                     CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  unsigned short iVar;
  
  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  sigma = 2./3.;
  
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradTurbVar_Kappa = new su2double [nVar];
  Proj_Mean_GradTurbVar_Edge = new su2double [nVar];
  Proj_Mean_GradTurbVar_Corrected = new su2double [nVar];
  Mean_GradTurbVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradTurbVar[iVar] = new su2double [nDim];
}

CAvgGradCorrected_TransLM::~CAvgGradCorrected_TransLM(void) {
  
  unsigned short iVar;
  
  delete [] Edge_Vector;
  delete [] Proj_Mean_GradTurbVar_Kappa;
  delete [] Proj_Mean_GradTurbVar_Edge;
  delete [] Proj_Mean_GradTurbVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradTurbVar[iVar];
  delete [] Mean_GradTurbVar;
}

void CAvgGradCorrected_TransLM::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {
  
  //  switch (config->GetKind_Turb_Model()) {
  //  case SA :
  //    /*--- Compute mean effective viscosity ---*/
  //    nu_i = Laminar_Viscosity_i/U_i[0];
  //    nu_j = Laminar_Viscosity_j/U_j[0];
  //    nu_e = 0.5*(nu_i+nu_j+TurbVar_i[0]+TurbVar_j[0]);
  //
  //    /*--- Compute vector going from iPoint to jPoint ---*/
  //    dist_ij_2 = 0; proj_vector_ij = 0;
  //    for (iDim = 0; iDim < nDim; iDim++) {
  //      Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
  //      dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  //      proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  //    }
  //    proj_vector_ij = proj_vector_ij/dist_ij_2;
  //
  //    /*--- Mean gradient approximation. Projection of the mean gradient
  //       in the direction of the edge ---*/
  //    for (iVar = 0; iVar < nVar; iVar++) {
  //      Proj_Mean_GradTurbVar_Kappa[iVar] = 0.0;
  //      Proj_Mean_GradTurbVar_Edge[iVar] = 0.0;
  //      for (iDim = 0; iDim < nDim; iDim++) {
  //        Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
  //        Proj_Mean_GradTurbVar_Kappa[iVar] += Mean_GradTurbVar[iVar][iDim]*Normal[iDim];
  //        Proj_Mean_GradTurbVar_Edge[iVar] += Mean_GradTurbVar[iVar][iDim]*Edge_Vector[iDim];
  //      }
  //      Proj_Mean_GradTurbVar_Corrected[iVar] = Proj_Mean_GradTurbVar_Kappa[iVar];
  //      Proj_Mean_GradTurbVar_Corrected[iVar] -= Proj_Mean_GradTurbVar_Edge[iVar]*proj_vector_ij -
  //          (TurbVar_j[iVar]-TurbVar_i[iVar])*proj_vector_ij;
  //    }
  //
  //    val_residual[0] = nu_e*Proj_Mean_GradTurbVar_Corrected[0]/sigma;
  //
  //    /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  //    if (implicit) {
  //      Jacobian_i[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]-nu_e*proj_vector_ij)/sigma;
  //      Jacobian_j[0][0] = (0.5*Proj_Mean_GradTurbVar_Corrected[0]+nu_e*proj_vector_ij)/sigma;
  //    }
  //    break;
  //
  //  }
}

CSourcePieceWise_TransLM::CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {
  
  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;
  
  /*--- Spalart-Allmaras closure constants ---*/
  cv1_3 = pow(7.1,3.0);
  k2 = pow(0.41,2.0);
  cb1 = 0.1355;
  cw2 = 0.3;
  cw3_6 = pow(2.0,6.0);
  sigma = 2./3.;
  cb2 = 0.622;
  cw1 = cb1/k2+(1+cb2)/sigma;
  
  /*-- Gamma-theta closure constants --*/
  c_e1    = 1.0;
  c_a1    = 2.0;
  c_e2    = 50.0;
  c_a2    = 0.06;
  sigmaf  = 1.0;
  s1      = 2.0;
  c_theta = 0.03;
  sigmat  = 2.0;
  
  /*-- Correlation constants --*/
  flen_global  = 12.0;
  alpha_global = 0.85;
  
  /*-- For debugging -AA --*/
  debugme = 0;
}

CSourcePieceWise_TransLM::~CSourcePieceWise_TransLM(void) { }

void CSourcePieceWise_TransLM::ComputeResidual_TransLM(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config, su2double &gamma_sep) {
  //************************************************//
  // Please do not delete //SU2_CPP2C comment lines //
  //************************************************//
  
  //SU2_CPP2C START CSourcePieceWise_TransLM::ComputeResidual_TransLM
  //SU2_CPP2C CALL_LIST START
  //SU2_CPP2C INVARS *TransVar_i
  //SU2_CPP2C OUTVARS *val_residual
  //SU2_CPP2C VARS DOUBLE *U_i **PrimVar_Grad_i Laminar_Viscosity_i Eddy_Viscosity_i dist_i
  //SU2_CPP2C VARS DOUBLE SCALAR c_a1 c_e1 c_a2 c_e2 c_theta alpha_global flen_global
  //SU2_CPP2C CALL_LIST END
  
  //SU2_CPP2C DEFINE nDim
  
  //SU2_CPP2C DECL_LIST START
  //SU2_CPP2C VARS DOUBLE SCALAR Vorticity
  //SU2_CPP2C DECL_LIST END
  
  /*-- Local intermediate variables --*/
  su2double rey_tc, flen, re_v, strain, f_onset1, f_onset2, f_onset3, f_onset, f_turb, tu;
  
  su2double prod, des;
  su2double f_lambda, re_theta = 0.0, re_theta_lim, r_t;
  su2double Velocity_Mag = 0.0, du_ds, theta, lambda, time_scale, var1, f_theta;
  su2double f_reattach;
  su2double dU_dx, dU_dy, dU_dz = 0.0;
  
  //SU2_CPP2C COMMENT START
  su2double val_residuald[2], TransVar_id[2];
  
  //SU2_CPP2C COMMENT END
  
  val_residual[0] = 0.0;
  val_residual[1] = 0.0;
  
  //SU2_CPP2C COMMENT START
  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  if (implicit) {
    val_Jacobian_i[0][0] = 0.0;
    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[1][1] = 0.0;
  }
  //SU2_CPP2C COMMENT END
  
  /* -- These lines included just so Tapenade doesn't complain --*/
//  rey  = 0.0;
//  mach = 0.0;
//  tu   = 0.0;
  //SU2_CPP2C COMMENT START
  /* -- These lines must be manually reinserted into the differentiated routine! --*/
//  rey  = config->GetReynolds();
//  mach = config->GetMach();
  tu   = config->GetTurbulenceIntensity_FreeStream();
  //SU2_CPP2C COMMENT END
  
  /*--- Compute vorticity and strain (TODO: Update for 3D) ---*/
  Vorticity = fabs(PrimVar_Grad_i[1][1]-PrimVar_Grad_i[2][0]);
  
  /*-- Strain = sqrt(2*Sij*Sij) --*/
  strain = sqrt(2.*(    PrimVar_Grad_i[1][0]*PrimVar_Grad_i[1][0]
                    +  0.5*pow(PrimVar_Grad_i[1][1]+PrimVar_Grad_i[2][0],2)
                    +  PrimVar_Grad_i[2][1]*PrimVar_Grad_i[2][1]  ));
  
  /*-- Note: no incompressible for now! --*/
  
  if (dist_i > 0.0) {   // Only operate away from wall
    
    /*-- Intermittency eq.: --*/
    
    rey_tc = (4.45*pow(tu,3) - 5.7*pow(tu,2) + 1.37*tu + 0.585)*TransVar_i[1];
    flen   = 0.171*pow(tu,2) - 0.0083*tu + 0.0306;
    
    re_v   = U_i[0]*pow(dist_i,2.)/Laminar_Viscosity_i*strain;  // Vorticity Reynolds number
    
    /*-- f_onset controls transition onset location --*/
    r_t      = Eddy_Viscosity_i/Laminar_Viscosity_i;
    f_onset1 = re_v / (2.193*rey_tc);
    f_onset2 = min(max(f_onset1, pow(f_onset1,4.)), 2.);
    f_onset3 = max(1. - pow(0.4*r_t,3),0.);
    f_onset  = max(f_onset2 - f_onset3, 0.);
    
    f_turb = exp(-pow(0.25*r_t,4));  // Medida eq. 10
    
    prod = flen*c_a1*U_i[0]*strain*sqrt(f_onset*TransVar_i[0]);
    prod = prod*(1. - c_e1*TransVar_i[0]);
    
    des = c_a2*U_i[0]*Vorticity*TransVar_i[0]*f_turb;
    des = des*(c_e2*TransVar_i[0] - 1.);
    
    val_residual[0] = prod - des;
    
    /*-- REtheta eq: --*/
    if (nDim==2) {
      Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2])/U_i[0];
    } else if (nDim==3) {
      Velocity_Mag = sqrt(U_i[1]*U_i[1]+U_i[2]*U_i[2]+U_i[3]*U_i[3])/U_i[0];
    }
    
    /*-- Gradient of velocity magnitude ---*/
    dU_dx = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][0]
                              +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][0]);
    if (nDim==3)
      dU_dx += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][0]);
    
    dU_dy = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][1]
                              +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][1]);
    if (nDim==3)
      dU_dy += 0.5*Velocity_Mag*( 2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][1]);
    
    if (nDim==3)
      dU_dz = 0.5*Velocity_Mag*( 2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][2]
                                +2*U_i[2]/U_i[0]*PrimVar_Grad_i[2][2]
                                +2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][2]);
    
    du_ds = U_i[1]/(U_i[0]*Velocity_Mag) * dU_dx +  // Streamwise velocity derivative
    U_i[2]/(U_i[0]*Velocity_Mag) * dU_dy;
    if (nDim==3)
      du_ds += U_i[3]/(U_i[0]*Velocity_Mag) * dU_dz;
    
    re_theta_lim = 20.;
    
    /*-- Fixed-point iterations to solve REth correlation --*/
    f_lambda = 1.;
    for (int iter=0; iter<10; iter++) {
      if (tu <= 1.3) {
        re_theta = f_lambda * (1173.51-589.428*tu+0.2196/(tu*tu));
      } else {
        re_theta = 331.5 * f_lambda*pow(tu-0.5658,-0.671);
      }
      re_theta = max(re_theta, re_theta_lim);
      
      theta  = re_theta * Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
      
      lambda = U_i[0]*theta*theta*du_ds / Laminar_Viscosity_i;
      lambda = min(max(-0.1, lambda),0.1);
      
      if (lambda<=0.0) {
        f_lambda = 1. - (-12.986*lambda - 123.66*lambda*lambda -
                         405.689*lambda*lambda*lambda)*exp(-pow(2./3*tu,1.5));
      } else {
        f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu);
      }
    }
    
    /*-- Calculate blending function f_theta --*/
    time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*Velocity_Mag);
    
    // Deactivated the f_wake parameter...
    //theta_bl   = TransVar_i[1]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
    //delta_bl   = 7.5*theta_bl;
    //delta      = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;
    //
    //f_wake = 1.;
    
    var1 = (TransVar_i[0]-1./c_e2)/(1.0-1./c_e2);
    var1 = 1. - pow(var1,2);
    
    //f_theta = min(max(f_wake*exp(-pow(dist_i/delta,4)), var1),1.0);
    f_theta = min(var1,1.0);
    
    val_residual[1] = c_theta*U_i[0]/time_scale *  (1.-f_theta) * (re_theta-TransVar_i[1]);
    
    //SU2_CPP2C COMMENT START
    cout << "val_res0: "  << val_residual[0]      << endl;
    cout << "val_res1: "  << val_residual[1]      << endl;
    cout << "dist_i: "    << dist_i               << endl;
    cout << "re_v: "      << re_v                 << endl;
    cout << "c_a1: "      << c_a1                 << endl;
    cout << "strain: "    << strain               << endl;
    cout << "primgrad10: "<< PrimVar_Grad_i[1][0] << endl;
    cout << "primgrad11: "<< PrimVar_Grad_i[1][1] << endl;
    cout << "primgrad20: "<< PrimVar_Grad_i[2][0] << endl;
    cout << "primgrad21: "<< PrimVar_Grad_i[2][1] << endl;
    cout << "f_onset: "   << f_onset              << endl;
    cout << "TransVar0: " << TransVar_i[0]        << endl;
    cout << "prod: "      << prod                 << endl;
    cout << "c_a2: "      << c_a2                 << endl;
    cout << "Vorticity: " << Vorticity            << endl;
    cout << "f_turb: "    << f_turb               << endl;
    cout << "des: "       << des                  << endl;
    cout << "du_ds: "     << du_ds                << endl;
    cout << "r_t:    "    << r_t                  << endl;
    cout << "rey_tc: "    << rey_tc               << endl;
    cout << "re_theta: "  << re_theta             << endl;
    
    /*-- Calculate term for separation correction --*/
    f_reattach = exp(-pow(0.05*r_t,4));
    gamma_sep = s1*max(0., re_v/(3.235*rey_tc)-1.)*f_reattach;
    gamma_sep = min(gamma_sep,2.0)*f_theta;
    
    /*--- Implicit part ---*/
    TransVar_id[0] = 1.0; TransVar_id[1] = 0.0;
    CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config);
    val_Jacobian_i[0][0] = val_residuald[0];
    val_Jacobian_i[1][0] = val_residuald[1];
    
    TransVar_id[0] = 0.0; TransVar_id[1] = 1.0;
    CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(TransVar_i, TransVar_id, val_residual, val_residuald, config);
    val_Jacobian_i[0][1] = val_residuald[0];
    val_Jacobian_i[1][1] = val_residuald[1];
    
    //SU2_CPP2C COMMENT END
  }
  //SU2_CPP2C END CSourcePieceWise_TransLM::ComputeResidual_TransLM
}


void CSourcePieceWise_TransLM::CSourcePieceWise_TransLM__ComputeResidual_TransLM_d(su2double *TransVar_i, su2double *TransVar_id, su2double *val_residual, su2double *val_residuald, CConfig *config)
{
  su2double rey_tc, flen, re_v, strain, f_onset1, f_onset2, f_onset3, f_onset,
  f_turb, tu;
  su2double rey_tcd, f_onset1d, f_onset2d, f_onsetd;
  su2double prod, des;
  su2double prodd, desd;
  su2double f_lambda, re_theta = 0.0, re_theta_lim, r_t;
  su2double Velocity_Mag = 0.0, du_ds, theta, lambda, time_scale,
  var1, f_theta;
  su2double var1d, f_thetad;
  su2double dU_dx, dU_dy, dU_dz = 0.0;
  su2double result1;
  su2double result1d;
  su2double arg1;
  su2double arg1d;
  su2double result2;
  su2double x2;
  su2double x1;
  su2double x1d;
  su2double y1;
  su2double y1d;
  val_residuald[0] = 0.0;
  val_residual[0] = 0.0;
  val_residuald[1] = 0.0;
  val_residual[1] = 0.0;
  /* -- These lines included just so Tapenade doesn't complain --*/
//  rey = 0.0;
//  mach = 0.0;
//  tu = 0.0;
//  rey  = config->GetReynolds();
//  mach = config->GetMach();
  tu   = config->GetTurbulenceIntensity_FreeStream();
  /*--- Compute vorticity and strain (TODO: Update for 3D) ---*/
  Vorticity = fabs(PrimVar_Grad_i[1][1] - PrimVar_Grad_i[2][0]);
  /*-- Strain = sqrt(2*Sij*Sij) --*/
  result1 = pow(PrimVar_Grad_i[1][1] + PrimVar_Grad_i[2][0], 2);
  arg1 = 2.*(PrimVar_Grad_i[1][0]*PrimVar_Grad_i[1][0]+0.5*result1+
             PrimVar_Grad_i[2][1]*PrimVar_Grad_i[2][1]);
  strain = sqrt(arg1);
  /*-- Note: no incompressible for now! --*/
  if (dist_i > 0.0) {
    /*-- Intermittency eq.: --*/
    // Only operate away from wall
    result1 = pow(tu, 3);
    result2 = pow(tu, 2);
    rey_tcd = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_id[1];
    rey_tc = (4.45*result1-5.7*result2+1.37*tu+0.585)*TransVar_i[1];
    result1 = pow(tu, 2);
    flen = 0.171*result1 - 0.0083*tu + 0.0306;
    result1 = pow(dist_i, 2.);
    re_v = U_i[0]*result1/Laminar_Viscosity_i*strain;
    /*-- f_onset controls transition onset location --*/
    // Vorticity Reynolds number
    r_t = Eddy_Viscosity_i/Laminar_Viscosity_i;
    f_onset1d = -(re_v*2.193*rey_tcd/(2.193*rey_tc*(2.193*rey_tc)));
    f_onset1 = re_v/(2.193*rey_tc);
    y1 = pow(f_onset1, 4.);
    y1d = 4.*f_onset1d*pow(f_onset1, 3);
    if (f_onset1 < y1) {
      x1d = y1d;
      x1 = y1;
    } else {
      x1d = f_onset1d;
      x1 = f_onset1;
    }
    if (x1 > 2.) {
      f_onset2 = 2.;
      f_onset2d = 0.0;
    } else {
      f_onset2d = x1d;
      f_onset2 = x1;
    }
    result1 = pow(0.4*r_t, 3);
    x2 = 1. - result1;
    if (x2 < 0.)
      f_onset3 = 0.;
    else
      f_onset3 = x2;
    if (f_onset2 - f_onset3 < 0.) {
      f_onset = 0.;
      f_onsetd = 0.0;
    } else {
      f_onsetd = f_onset2d;
      f_onset = f_onset2 - f_onset3;
    }
    result1 = pow(0.25*r_t, 4);
    f_turb = exp(-result1);
    // Medida eq. 10
    arg1d = f_onsetd*TransVar_i[0] + f_onset*TransVar_id[0];
    arg1 = f_onset*TransVar_i[0];
    //result1d = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
    if (arg1 == 0) { result1d = 0.0; } else result1d = arg1d/(2.0*sqrt(arg1));
    result1 = sqrt(arg1);
    prodd = flen*c_a1*U_i[0]*strain*result1d;
    prod = flen*c_a1*U_i[0]*strain*result1;
    prodd = prodd*(1.-c_e1*TransVar_i[0]) - prod*c_e1*TransVar_id[0];
    prod = prod*(1.-c_e1*TransVar_i[0]);
    desd = c_a2*U_i[0]*Vorticity*f_turb*TransVar_id[0];
    des = c_a2*U_i[0]*Vorticity*TransVar_i[0]*f_turb;
    desd = desd*(c_e2*TransVar_i[0]-1.) + des*c_e2*TransVar_id[0];
    des = des*(c_e2*TransVar_i[0]-1.);
    val_residuald[0] = prodd - desd;
    val_residual[0] = prod - des;
    /*-- REtheta eq: --*/
    if (nDim == 2) {
      arg1 = U_i[1]*U_i[1] + U_i[2]*U_i[2];
      result1 = sqrt(arg1);
      Velocity_Mag = result1/U_i[0];
    } else
      if (nDim == 3) {
        arg1 = U_i[1]*U_i[1] + U_i[2]*U_i[2] + U_i[3]*U_i[3];
        result1 = sqrt(arg1);
        Velocity_Mag = result1/U_i[0];
      }
    /*-- Gradient of velocity magnitude ---*/
    dU_dx = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][0]+2*U_i[2
                                                                         ]/U_i[0]*PrimVar_Grad_i[2][0]);
    if (nDim == 3)
      dU_dx += 0.5*Velocity_Mag*(2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][0]);
    dU_dy = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][1]+2*U_i[2
                                                                         ]/U_i[0]*PrimVar_Grad_i[2][1]);
    if (nDim == 3)
      dU_dy += 0.5*Velocity_Mag*(2*U_i[3]/U_i[0]*PrimVar_Grad_i[3][1]);
    if (nDim == 3)
      dU_dz = 0.5*Velocity_Mag*(2*U_i[1]/U_i[0]*PrimVar_Grad_i[1][2]+2*
                                U_i[2]/U_i[0]*PrimVar_Grad_i[2][2]+2*U_i[3]/U_i[0]*
                                PrimVar_Grad_i[3][2]);
    du_ds = U_i[1]/(U_i[0]*Velocity_Mag)*dU_dx + U_i[2]/(U_i[0]*
                                                         Velocity_Mag)*dU_dy;
    // Streamwise velocity derivative
    if (nDim == 3)
      du_ds += U_i[3]/(U_i[0]*Velocity_Mag)*dU_dz;
    re_theta_lim = 20.;
    /*-- Fixed-point iterations to solve REth correlation --*/
    f_lambda = 1.;
    {
      su2double x3;
      for (int iter = 0; iter < 10; ++iter) {
        if (tu <= 1.3)
          re_theta = f_lambda*(1173.51-589.428*tu+0.2196/(tu*tu));
        else {
          result1 = pow(tu - 0.5658, -0.671);
          re_theta = 331.5*f_lambda*result1;
        }
        if (re_theta < re_theta_lim) re_theta = re_theta_lim;
        theta = re_theta*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag);
        lambda = U_i[0]*theta*theta*du_ds/Laminar_Viscosity_i;
        if (-0.1 < lambda)
          x3 = lambda;
        else
          x3 = -0.1;
        if (x3 > 0.1)
          lambda = 0.1;
        else
          lambda = x3;
        if (lambda <= 0.0) {
          result1 = pow(2./3*tu, 1.5);
          f_lambda = 1. - (-12.986*lambda-123.66*lambda*lambda-405.689
                           *lambda*lambda*lambda)*exp(-result1);
        } else
          f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-2.*tu);
      }
    }
    /*-- Calculate blending function f_theta --*/
    time_scale = 500.0*Laminar_Viscosity_i/(U_i[0]*Velocity_Mag*
                                            Velocity_Mag);
    // Deactivated the f_wake parameter...
    //theta_bl   = TransVar_i[1]*Laminar_Viscosity_i / (U_i[0]*Velocity_Mag);
    //delta_bl   = 7.5*theta_bl;
    //delta      = 50.0*Vorticity*dist_i/Velocity_Mag*delta_bl + 1e-20;
    //
    //f_wake = 1.;
    var1d = TransVar_id[0]/(1.0-1./c_e2);
    var1 = (TransVar_i[0]-1./c_e2)/(1.0-1./c_e2);    
    result1 = pow(var1, 2.0);
    result1d = 2.0*var1d*pow(var1, 1.0);
    var1d = -result1d;
    var1 = 1. - result1;
    if (var1 > 1.0) {
      f_theta = 1.0;
      f_thetad = 0.0;
    } else {
      f_thetad = var1d;
      f_theta = var1;
    }
    val_residuald[1] = c_theta*U_i[0]*(-(f_thetad*(re_theta-TransVar_i[1])
                                         )-(1.-f_theta)*TransVar_id[1])/time_scale;
    val_residual[1] = c_theta*U_i[0]/time_scale*(1.-f_theta)*(re_theta-
                                                              TransVar_i[1]);
  } else
    *val_residuald = 0.0;
  
}
