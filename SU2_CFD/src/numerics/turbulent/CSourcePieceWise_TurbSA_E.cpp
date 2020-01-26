/*!
 * \file CSourcePieceWise_TurbSA_E.cpp
 * \brief Implementation of numerics class CSourcePieceWise_TurbSA_E.
 * \author F. Palacios, T. Economon
 * \version 7.0.0 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../include/numerics/turbulent/CSourcePieceWise_TurbSA_E.hpp"

CSourcePieceWise_TurbSA_E::CSourcePieceWise_TurbSA_E(unsigned short val_nDim, unsigned short val_nVar,
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

CSourcePieceWise_TurbSA_E::~CSourcePieceWise_TurbSA_E(void) { }

void CSourcePieceWise_TurbSA_E::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

    //  AD::StartPreacc();
    //  AD::SetPreaccIn(V_i, nDim+6);
    //  AD::SetPreaccIn(Vorticity_i, nDim);
    //  AD::SetPreaccIn(StrainMag_i);
    //  AD::SetPreaccIn(TurbVar_i[0]);
    //  AD::SetPreaccIn(TurbVar_Grad_i[0], nDim);
    //  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);

    if (incompressible) {
      Density_i = V_i[nDim+2];
      Laminar_Viscosity_i = V_i[nDim+4];
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


    /*
     From NASA Turbulence model site. http://turbmodels.larc.nasa.gov/spalart.html
     This form was developed primarily to improve the near-wall numerical behavior of the model (i.e., the goal was to improve the convergence behavior). The reference is:
     Edwards, J. R. and Chandra, S. "Comparison of Eddy Viscosity-Transport Turbulence Models for Three-Dimensional, Shock-Separated Flowfields," AIAA Journal, Vol. 34, No. 4, 1996, pp. 756-763.
     In this modificaton Omega is replaced by Strain Rate
     */

    /*--- Evaluate Omega, here Omega is the Strain Rate ---*/

    Sbar = 0.0;
    for(iDim=0;iDim<nDim;++iDim){
        for(jDim=0;jDim<nDim;++jDim){
            Sbar+= (PrimVar_Grad_i[1+iDim][jDim]+PrimVar_Grad_i[1+jDim][iDim])*(PrimVar_Grad_i[1+iDim][jDim]);}}
    for(iDim=0;iDim<nDim;++iDim){
        Sbar-= ((2.0/3.0)*pow(PrimVar_Grad_i[1+iDim][iDim],2.0));}

    Omega= sqrt(max(Sbar,0.0));

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

        //Shat = S + TurbVar_i[0]*fv2*inv_k2_d2;
        Shat = max(S*((1.0/max(Ji,1.0e-16))+fv1),1.0e-16);

        Shat = max(Shat, 1.0e-10);
        inv_Shat = 1.0/Shat;

        /*--- Production term ---*/;

        Production = cb1*Shat*TurbVar_i[0]*Volume;

        /*--- Destruction term ---*/

        r = min(TurbVar_i[0]*inv_Shat*inv_k2_d2,10.0);
        r=tanh(r)/tanh(1.0);

        g = r + cw2*(pow(r,6.0)-r);
        g_6 = pow(g,6.0);
        glim = pow((1.0+cw3_6)/(g_6+cw3_6),1.0/6.0);
        fw = g*glim;

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
        else dShat = -S*pow(Ji,-2.0)/nu + S*dfv1;
        val_Jacobian_i[0][0] += cb1*(TurbVar_i[0]*dShat+Shat)*Volume;

        /*--- Implicit part, destruction term ---*/

        dr = (Shat-TurbVar_i[0]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
        dr=(1-pow(tanh(r),2.0))*(dr)/tanh(1.0);
        dg = dr*(1.+cw2*(6.0*pow(r,5.0)-1.0));
        dfw = dg*glim*(1.-g_6/(g_6+cw3_6));
        val_Jacobian_i[0][0] -= cw1*(dfw*TurbVar_i[0] + 2.0*fw)*TurbVar_i[0]/dist_i_2*Volume;

    }

    //  AD::SetPreaccOut(val_residual[0]);
    //  AD::EndPreacc();

}
