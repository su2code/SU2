/*!
 * \file turb_sources.hpp
 * \brief Numerics classes for integration of source terms in turbulence problems.
 * \version 7.4.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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

#pragma once
#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"
#include "../../scalar/scalar_sources.hpp"


/*!
 * \class CSourcePieceWise_TranLM
 * \brief Class for integrating the source terms of the Menter SST turbulence model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  su2double g_eff_i,
  g_eff_j,
  g_sep_i,
  g_sep_j;
  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double sigmat = 2.0;

  su2double Vorticity;
  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4]; /// Static storage for the Jacobian (which needs to be pointer for return type).

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] constants - SST model constants.
   * \param[in] val_kine_Inf - Freestream k, for SST with sustaining terms.
   * \param[in] val_omega_Inf - Freestream w, for SST with sustaining terms.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransLM(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config) 
      : CNumerics(val_nDim, 2, config),
        idx(val_nDim, config->GetnSpecies()) {

    
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
    Jacobian_i[0] = Jacobian_Buffer;
    Jacobian_i[1] = Jacobian_Buffer + 2;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {
  /*--- ScalarVar[0] = k, ScalarVar[0] = w, TransVar[0] = gamma, and TransVar[0] = ReThetaT ---*/
  /*--- dU/dx = PrimVar_Grad[1][0] ---*/
  AD::StartPreacc();
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(ScalarVar_i, nVar);
  AD::SetPreaccIn(ScalarVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(TransVar_i, nVar);
  AD::SetPreaccIn(TransVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+idx.Velocity(), nDim);
  AD::SetPreaccIn(Vorticity_i, 3);
  

  unsigned short iDim;
  
  su2double Corr_Rec = 0.0, Corr_F_length = 0.0, Corr_Ret = 0.0;      
  su2double dU_dx, dU_dy, dU_dz = 0.0;  
  

  su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);
  
  su2double diverg = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) diverg += PrimVar_Grad_i[iDim+idx.Velocity()][iDim];

  su2double vel_u = 0.0, vel_v = 0.0, vel_w = 0.0;
  if(nDim ==2){
      vel_u = V_i[idx.Velocity()];
      vel_v = V_i[1+idx.Velocity()];
  }
  else if(nDim ==3){
      vel_u = V_i[idx.Velocity()];
      vel_v = V_i[1+idx.Velocity()];
      vel_w = V_i[2+idx.Velocity()];
  }

  su2double Velocity_Mag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);

  AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

  Density_i = V_i[idx.Density()];
  Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
  Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

  su2double Tu = 100.0*sqrt( 2.0 * ScalarVar_i[0] / 3.0 ) / Velocity_Mag;

  Residual[0] = 0.0;       Residual[1] = 0.0;
  Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
  Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;
  
  if (dist_i > 1e-10) {

    /*--- Corr_RetC correlation*/

    if(TransVar_i[1] <= 1870){
     Corr_Rec = TransVar_i[1] - (396.035e-02 + (-120.656e-04)*TransVar_i[1] + (868.230e-06)*pow(TransVar_i[1], 2.) 
                                +(-696.506e-09)*pow(TransVar_i[1], 3.) + (174.105e-12)*pow(TransVar_i[1], 4.));
    } else {
     Corr_Rec = TransVar_i[1] - (  593.11 + (TransVar_i[1] - 1870.0) * 0.482);
    }

    /*--- F_length correlation*/

    if(TransVar_i[1] < 400){
      Corr_F_length = 398.189e-01 + (-119.270e-04)*TransVar_i[1] + (-132.567e-06) * pow(TransVar_i[1], 2.);
    } 
    else if(TransVar_i[1] >= 400 && TransVar_i[1] < 596) {
      Corr_F_length = 263.404 + (-123.939e-02) * TransVar_i[1] + (194.548e-5) * pow(TransVar_i[1], 2.) + (-101.695e-08) * pow(TransVar_i[1], 3.);
    } 
    else if(TransVar_i[1] >= 596 && TransVar_i[1] < 1200) {
      Corr_F_length = 0.5 - (TransVar_i[1] - 596.0) * 3.0e-04; 
    } 
    else {
      Corr_F_length = 0.3188;
    }    

    /*--- F_length ---*/
    su2double r_omega = Density_i*dist_i*dist_i*ScalarVar_i[1]/Laminar_Viscosity_i;
    su2double f_sub = exp(-pow(r_omega/200.0, 2));
    su2double F_length = Corr_F_length *(1.-f_sub) + 40.0*f_sub;

    /*--- F_onset ---*/
    su2double R_t = Density_i*ScalarVar_i[0]/ Laminar_Viscosity_i/ ScalarVar_i[1];
    su2double Re_v = Density_i*dist_i*dist_i*StrainMag_i/Laminar_Viscosity_i;
    su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
    su2double F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
    su2double F_onset3 = max(1.0 - pow(R_t/2.5, 3.0), 0.0);
    su2double F_onset  = max(F_onset2 - F_onset3, 0.0);

    /*-- Gradient of velocity magnitude ---*/

    dU_dx = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][0]
                              +2.*vel_v*PrimVar_Grad_i[2][0]);
    if (nDim==3)
      dU_dx += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][0]);

    dU_dy = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][1]
                              +2.*vel_v*PrimVar_Grad_i[2][1]);
    if (nDim==3)
      dU_dy += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][1]);

    if (nDim==3)
      dU_dz = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][2]
                                +2.*vel_v*PrimVar_Grad_i[2][2]
                                +2.*vel_w*PrimVar_Grad_i[3][2]);

    su2double du_ds = vel_u/Velocity_Mag*dU_dx + vel_v/Velocity_Mag*dU_dy;
    if (nDim==3)
      du_ds += vel_w/Velocity_Mag * dU_dz;
    
    /*-- Calculate blending function f_theta --*/

    su2double time_scale = 500.0*Laminar_Viscosity_i/Density_i/ Velocity_Mag /Velocity_Mag;
    su2double theta_bl   = TransVar_i[1]*Laminar_Viscosity_i / Density_i /Velocity_Mag;
    su2double delta_bl   = 7.5*theta_bl;
    su2double delta      = 50.0*VorticityMag*dist_i/Velocity_Mag*delta_bl + 1e-20;
    
    su2double re_omega = Density_i*ScalarVar_i[1]*dist_i*dist_i/Laminar_Viscosity_i;
    su2double f_wake = exp(-pow(re_omega/(1.0e+05),2));
    
    su2double var1 = (TransVar_i[0]-1.0/c_e2)/(1.0-1.0/c_e2);
    su2double var2 = 1.0 - pow(var1,2.0);
    su2double f_theta = min(max(f_wake*exp(-pow(dist_i/delta, 4)), var2), 1.0); 
    su2double f_turb = exp(-pow(R_t/4, 4));

    /*--- Corr_RetT correlation*/
    su2double Corr_Ret_lim = 20.0;
    su2double f_lambda = 1.0;
    
    su2double Retheta_Error = 200.0 , Retheta_old = 0.0;
    su2double lambda = 0.0;
    Corr_Ret = 20.0; 
    Tu = max(Tu,0.027);  

    for (int iter=0; iter<100 ; iter++) {
      
      su2double theta  = Corr_Ret * Laminar_Viscosity_i / Density_i/ Velocity_Mag;
      lambda = Density_i*theta*theta/ Laminar_Viscosity_i*du_ds;
      lambda = min(max(-0.1, lambda), 0.1);

      if (lambda<=0.0) {
        f_lambda = 1. - (-12.986*lambda - 123.66*lambda*lambda -
                         405.689*lambda*lambda*lambda)*exp(-pow(Tu/1.5, 1.5));
      } else {
        f_lambda = 1. + 0.275*(1.-exp(-35.*lambda))*exp(-Tu/0.5);
      }

      if (Tu <= 1.3) {
        Corr_Ret = f_lambda * (1173.51 - 589.428*Tu + 0.2196/Tu/Tu);
      } else {
        Corr_Ret = 331.5 * f_lambda*pow(Tu - 0.5658,-0.671);
      }
      Corr_Ret = max(Corr_Ret, Corr_Ret_lim);

      Retheta_Error = fabs(Retheta_old -Corr_Ret)/Retheta_old;

      if(Retheta_Error < 0.0000001){
        iter = 101;
      }

      Retheta_old = Corr_Ret;
    }
    
    /*-- production term of Intermeittency(Gamma) --*/
    su2double Pg = F_length*c_a1*Density_i*StrainMag_i*sqrt(F_onset*TransVar_i[0])*(1.0 - c_e1 * TransVar_i[0]);

    /*-- destruction term of Intermeittency(Gamma) --*/
    su2double Dg = c_a2*Density_i*VorticityMag*TransVar_i[0]*f_turb*(c_e2*TransVar_i[0] - 1.0);

    /*-- production term of ReThetaT --*/
    su2double Pthetat = c_theta*Density_i/time_scale * (Corr_Ret-TransVar_i[1])  * (1.0-f_theta);

    /*--- Source ---*/
    Residual[0] += (Pg - Dg)*Volume;
    Residual[1] += Pthetat*Volume; 

    /*--- Implicit part ---*/   
    Jacobian_i[0][0] = (F_length*c_a1*StrainMag_i*sqrt(F_onset)*(0.5*pow(TransVar_i[0], -0.5) -1.5*c_e1*pow(TransVar_i[0], 0.5))
                       - c_a2 * VorticityMag*f_turb*(2.0*c_e2*TransVar_i[0]-1.0) )*Volume;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -c_theta/time_scale*(1.0-f_theta)*Volume;    
  }  
  
  AD::SetPreaccOut(Residual, nVar);  
  AD::EndPreacc();

  return ResidualType<>(Residual, Jacobian_i, nullptr);
  }
};

