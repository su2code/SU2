/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
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
 * \brief Class for integrating the source terms of the LM transition model equations.
 * \ingroup SourceDiscr
 * \author S. Kang.
 */
template <class FlowIndices>
class CSourcePieceWise_TransLM final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  /*--- LM Closure constants ---*/
  const su2double c_e1 = 1.0;
  const su2double c_a1 = 2.0;
  const su2double c_e2 = 50.0;
  const su2double c_a2 = 0.06;
  const su2double sigmaf = 1.0;
  const su2double s1 = 2.0;
  const su2double c_theta = 0.03;
  const su2double sigmat = 2.0;

  su2double Residual[2];
  su2double* Jacobian_i[2];
  su2double Jacobian_Buffer[4];// Static storage for the Jacobian (which needs to be pointer for return type).

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
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

    su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                Vorticity_i[1]*Vorticity_i[1] +
                                Vorticity_i[2]*Vorticity_i[2]);
  
    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1+idx.Velocity()];
    const su2double vel_w = (nDim ==3) ? V_i[2+idx.Velocity()] : 0.0;

    const su2double Velocity_Mag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);

    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.LaminarViscosity()], V_i[idx.EddyViscosity()]);

    Density_i = V_i[idx.Density()];
    Laminar_Viscosity_i = V_i[idx.LaminarViscosity()];
    Eddy_Viscosity_i = V_i[idx.EddyViscosity()];

    Residual[0] = 0.0;       Residual[1] = 0.0;
    Jacobian_i[0][0] = 0.0;  Jacobian_i[0][1] = 0.0;
    Jacobian_i[1][0] = 0.0;  Jacobian_i[1][1] = 0.0;
  
    if (dist_i > 1e-10) {

      /*--- Corr_RetC correlation*/
      su2double Corr_Rec = 0.0;
      if(TransVar_i[1] <= 1870){
        Corr_Rec = TransVar_i[1] - (396.035e-02 + (-120.656e-04)*TransVar_i[1] + (868.230e-06)*pow(TransVar_i[1], 2.) 
                                +(-696.506e-09)*pow(TransVar_i[1], 3.) + (174.105e-12)*pow(TransVar_i[1], 4.));
      } else {
        Corr_Rec = TransVar_i[1] - (  593.11 + (TransVar_i[1] - 1870.0) * 0.482);
      }

      /*--- F_length correlation*/
      su2double Corr_F_length = 0.0;
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
      const su2double r_omega = Density_i*dist_i*dist_i*ScalarVar_i[1]/Laminar_Viscosity_i;
      const su2double f_sub = exp(-pow(r_omega/200.0, 2));
      const su2double F_length = Corr_F_length *(1.-f_sub) + 40.0*f_sub;

      /*--- F_onset ---*/
      const su2double R_t = Density_i*ScalarVar_i[0]/ Laminar_Viscosity_i/ ScalarVar_i[1];
      const su2double Re_v = Density_i*dist_i*dist_i*StrainMag_i/Laminar_Viscosity_i;
      const su2double F_onset1 = Re_v / (2.193 * Corr_Rec);
      const su2double F_onset2 = min(max(F_onset1, pow(F_onset1, 4.0)), 2.0);
      const su2double F_onset3 = max(1.0 - pow(R_t/2.5, 3.0), 0.0);
      const su2double F_onset = max(F_onset2 - F_onset3, 0.0);

      /*-- Gradient of velocity magnitude ---*/

      su2double dU_dx = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][0]
                              +2.*vel_v*PrimVar_Grad_i[2][0]);
      if (nDim==3)
        dU_dx += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][0]);

      su2double dU_dy = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][1]
                              +2.*vel_v*PrimVar_Grad_i[2][1]);
      if (nDim==3)
        dU_dy += 0.5/Velocity_Mag*( 2.*vel_w*PrimVar_Grad_i[3][1]);

      su2double dU_dz = 0.0;
      if (nDim==3)
        dU_dz = 0.5/Velocity_Mag*( 2.*vel_u*PrimVar_Grad_i[1][2]
                                +2.*vel_v*PrimVar_Grad_i[2][2]
                                +2.*vel_w*PrimVar_Grad_i[3][2]);

      su2double du_ds = vel_u/Velocity_Mag*dU_dx + vel_v/Velocity_Mag*dU_dy;
      if (nDim==3)
        du_ds += vel_w/Velocity_Mag * dU_dz;
    
      /*-- Calculate blending function f_theta --*/

      const su2double time_scale = 500.0*Laminar_Viscosity_i/Density_i/ Velocity_Mag /Velocity_Mag;
      const su2double theta_bl = TransVar_i[1]*Laminar_Viscosity_i / Density_i /Velocity_Mag;
      const su2double delta_bl = 7.5*theta_bl;
      const su2double delta = 50.0*VorticityMag*dist_i/Velocity_Mag*delta_bl + 1e-20;
    
      const su2double re_omega = Density_i*ScalarVar_i[1]*dist_i*dist_i/Laminar_Viscosity_i;
      const su2double f_wake = exp(-pow(re_omega/(1.0e+05),2));
    
      const su2double var1 = (TransVar_i[0]-1.0/c_e2)/(1.0-1.0/c_e2);
      const su2double var2 = 1.0 - pow(var1,2.0);
      const su2double f_theta = min(max(f_wake*exp(-pow(dist_i/delta, 4)), var2), 1.0);
      const su2double f_turb = exp(-pow(R_t/4, 4));

      /*--- Corr_RetT correlation*/
      su2double Corr_Ret_lim = 20.0;
      su2double f_lambda = 1.0;
    
      su2double Retheta_Error = 200.0 , Retheta_old = 0.0;
      su2double lambda = 0.0;
      su2double Corr_Ret = 20.0;
      const su2double Tu = max(100.0*sqrt( 2.0 * ScalarVar_i[0] / 3.0 ) / Velocity_Mag,0.027);
    
      for (int iter=0; iter<100 ; iter++) {
      
        su2double theta = Corr_Ret * Laminar_Viscosity_i / Density_i/ Velocity_Mag;
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
    const su2double Pg = F_length*c_a1*Density_i*StrainMag_i*sqrt(F_onset*TransVar_i[0])*(1.0 - c_e1 * TransVar_i[0]);

    /*-- destruction term of Intermeittency(Gamma) --*/
    const su2double Dg = c_a2*Density_i*VorticityMag*TransVar_i[0]*f_turb*(c_e2*TransVar_i[0] - 1.0);

    /*-- production term of ReThetaT --*/
    const su2double Pthetat = c_theta*Density_i/time_scale * (Corr_Ret-TransVar_i[1])  * (1.0-f_theta);

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

/*!
 * \class CSourcePieceWise_TranEN
 * \brief Class for integrating the source terms of the e^N transition model equations.
 * \ingroup SourceDiscr
 * \author R. Roos
 */
template <class FlowIndices>
class CSourcePieceWise_TransEN final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  su2double Residual, *Jacobian_i;
  su2double Jacobian_Buffer; /*!< \brief Static storage for the Jacobian (which needs to be pointer for return type). */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TransEN(unsigned short val_nDim, unsigned short val_nVar, const CConfig* config)
      : CNumerics(val_nDim, 1, config),
        idx(val_nDim, config->GetnSpecies()) {
    
    /*--- "Allocate" the Jacobian using the static buffer. ---*/
	Jacobian_i = &Jacobian_Buffer;
  }

  /*!
   * \brief Residual for source term integration.
   * \param[in] config - Definition of the particular problem.
   * \return A lightweight const-view (read-only) of the residual/flux and Jacobians.
   */
  ResidualType<> ComputeResidual(const CConfig* config) override {

    AD::StartPreacc();
    AD::SetPreaccIn(V_i[idx.Density()], V_i[idx.Pressure()], V_i[idx.LaminarViscosity()], StrainMag_i, ScalarVar_i[0], Volume, dist_i);
    AD::SetPreaccIn(&V_i[idx.Velocity()], nDim);
    AD::SetPreaccIn(Vorticity_i, 3);
    AD::SetPreaccIn(PrimVar_Grad_i + idx.Velocity(), nDim, nDim);
    AD::SetPreaccIn(ScalarVar_Grad_i[0], nDim);

    su2double rho 	= V_i[idx.Density()];
    su2double p 	= V_i[idx.Pressure()];
    su2double muLam = V_i[idx.LaminarViscosity()];

    const su2double VorticityMag = GeometryToolbox::Norm(3, Vorticity_i);

    const su2double vel_u = V_i[idx.Velocity()];
    const su2double vel_v = V_i[1+idx.Velocity()];
    const su2double vel_w = (nDim ==3) ? V_i[2+idx.Velocity()] : 0.0;

    const su2double vel_mag = sqrt(vel_u*vel_u + vel_v*vel_v + vel_w*vel_w);

    su2double rhoInf			= config->GetDensity_FreeStreamND();
    su2double pInf 				= config->GetPressure_FreeStreamND();
    const su2double *velInf 	= config->GetVelocity_FreeStreamND();

    su2double velInf2 = 0.0;
    for(unsigned short iDim = 0; iDim < nDim; ++iDim) {
    	velInf2 += velInf[iDim]*velInf[iDim];
    }

    su2double Gamma = config->GetGamma();

    Residual = 0.0;
    Jacobian_i[0] = 0.0;

    if (dist_i > 1e-10) {

      su2double u_e;
      if (config->GetKind_Regime() == ENUM_REGIME::COMPRESSIBLE) {

    	/*--- Estimate of the equivalent flow velocity at the edge of the boundary layer based on compressible Bernoulli's equation ---*/
		const su2double G_over_Gminus_one = Gamma/(Gamma-1);

		const su2double rho_e = pow(((pow(rhoInf,Gamma)/pInf)*p),(1/Gamma));

		u_e = sqrt(2*(G_over_Gminus_one*(pInf/rhoInf) + (velInf2/2) - G_over_Gminus_one*(p/rho_e)));

      } else {

    	/*--- Inviscid edge velocity based on incompressible Bernoulli's equation---*/
    	u_e = sqrt((rhoInf*velInf2 + 2*(p-pInf))/rho);

      }

      /*--- Local pressure-gradient parameter for the boundary layer shape factor. Minimum value of 0.328 for stability ---*/
      const su2double H_L 		= max(((StrainMag_i*dist_i)/u_e),0.328);

      /*--- Integral shape factor ---*/
      const su2double H_12 		= 13.9766*pow(H_L,4) - 22.9166*pow(H_L,3) + 13.7227*pow(H_L,2) - 1.0023*H_L + 1.6778;

      /*--- F growth parameters ---*/
      const su2double DH_12 	= (0.0616*pow(H_12,2) + 0.2339*H_12 + 3.4298)/
    		  	  	  	  	  	     (0.0047*pow(H_12,3) - 0.1056*pow(H_12,2) + 0.9350*H_12 - 1.2071);

      const su2double lH_12 	= (6.54*H_12 - 14.07)/pow(H_12,2);
      const su2double mH_12 	= (0.058*(pow((H_12 - 4),2)/(H_12 - 1)) - 0.068)*(1/lH_12);

      const su2double F_growth 	= DH_12*((1 + mH_12)*lH_12)/2;

      /*--- F crit parameters ---*/
      const su2double Re_y 		= (rho*vel_mag*dist_i)/muLam;
      const su2double k_y 		= -0.00315*pow(H_12,3) + 0.0986*pow(H_12,2) - 0.242*H_12 + 3.739;
      const su2double Re_d2_0 	= pow(10,(0.7*tanh((14/(H_12 - 1)) - 9.24) + 2.492/pow((H_12 - 1),0.43) + 0.62));
      const su2double Re_y_0	= k_y * Re_d2_0;

      short int F_crit = 0;
      if (Re_y < Re_y_0){
    	F_crit = 0;
      } else {
    	F_crit = 1;
      }

      /*--- Source term expresses stream wise growth of Tollmien_schlichting instabilities ---*/
      const su2double dn_over_dRe_d2 = 0.028*(H_12 - 1) - 0.0345*exp(-pow((3.87/(H_12 - 1) - 2.52),2));

      /*--- Production term ---*/
      const su2double P_amplification = rho*VorticityMag*F_crit*F_growth*dn_over_dRe_d2;

      /*--- Source ---*/
      Residual = P_amplification * Volume;

      /*--- Implicit part ---*/
	  Jacobian_i[0] = (rho*VorticityMag*F_crit*F_growth) * Volume;
	  
    }

    AD::SetPreaccOut(Residual);
    AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
  }
};
