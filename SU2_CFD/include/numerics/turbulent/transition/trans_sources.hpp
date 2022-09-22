/*!
 * \file trans_sources.hpp
 * \brief Numerics classes for integration of source terms in transition problems.
 * \author R. Roos
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
 * \class CSourcePieceWise_TranEN
 * \brief Class for integrating the source terms of the e^N transition model equations.
 * \ingroup SourceDiscr
 * \author R. Roos
 */
template <class FlowIndices>
class CSourcePieceWise_TransEN final : public CNumerics {
 private:
  const FlowIndices idx; /*!< \brief Object to manage the access to the flow primitives. */

  /*su2double g_eff_i,
  g_eff_j,
  g_sep_i,
  g_sep_j;*/

  su2double Vorticity;
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

    su2double VorticityMag = sqrt(Vorticity_i[0]*Vorticity_i[0] +
                                  Vorticity_i[1]*Vorticity_i[1] +
                                  Vorticity_i[2]*Vorticity_i[2]);

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

      const su2double rho_e 	= pow(((pow(rhoInf,Gamma)/pInf)*p),(1/Gamma));

      /*--- Estimate of the flow velocity at the edge of the boundary layer ---*/
      const su2double G_over_Gminus_one = Gamma/(Gamma-1);

      const su2double u_e 		= pow(2*(G_over_Gminus_one*(pInf/rhoInf) + (velInf2/2) - G_over_Gminus_one*(p/rho_e)),0.5);

      /*--- Local pressure-gradient parameter for the boundary layer shape factor ---*/
      const su2double H_L 		= (StrainMag_i*dist_i)/u_e;

      /*--- Integral shape factor ---*/
      const su2double H_12 		= 13.9766*pow(H_L,4) - 22.9166*pow(H_L,3) + 13.7227*pow(H_L,2) - 1.0023*H_L + 1.6778;

      /*--- F growth parameters ---*/
      const su2double DH_12 	= (0.0616*pow(H_12,2) + 0.2339*H_12 + 3.4298)/
    		  	  	  	  	  	     (0.0047*pow(H_12,3) - 0.1056*pow(H_12,2) + 0.9350*H_12 - 1.2071);

      const su2double lH_12 	= (6.54*H_12 - 14.07)/pow(H_12,2);
      const su2double mH_12 	= (0.058*(pow((H_12 - 4),2)/(H_12 - 1)) - 0.068)*(1/lH_12);

      const su2double F_growth 	= DH_12*((1 + mH_12)*lH_12)/2;

      /*--- F critical parameters ---*/
      const su2double Re_y 		= (rho*vel_mag*dist_i)/muLam;
      const su2double k_y 		= -0.00315*pow(H_12,3) + 0.0986*pow(H_12,2) - 0.242*H_12 + 3.739;
      const su2double Re_d2_0 	= pow(10,(0.7*tanh((14/(H_12 - 1)) - 9.24) + 2.492/pow((H_12 - 1),0.43) + 0.62));
      const su2double Re_y_0	= k_y * Re_d2_0;

      short F_crit = 0;
      if (Re_y < Re_y_0){
    	F_crit = 0;
      } else {
    	F_crit = 1;
      }

      /*--- Source term expresses streamwise growth of Tollmien_schlichting instabilities ---*/
      const su2double dn_over_dRe_d2 = 0.028*(H_12 - 1) - 0.0345*exp(-pow((3.87/(H_12 - 1) - 2.52),2));

      /*--- Production term ---*/
      const su2double P_amplification = rho*VorticityMag*F_crit*F_growth*dn_over_dRe_d2;

      /*--- Source ---*/
      Residual = P_amplification * Volume;

      /*--- Implicit part ---*/
      Jacobian_i[0] *= Volume;
    }
  
    AD::SetPreaccOut(Residual);
    AD::EndPreacc();

    return ResidualType<>(&Residual, &Jacobian_i, nullptr);
  }
};

