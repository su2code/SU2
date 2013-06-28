/*!
 * \file numerics_source_ad.cpp
 * \brief This file contains Automatically Differentiated versions
 * of appropriate source terms. These routines are produced
 * semi-automatically using python, Tapenade and some minor requirement
 * to add in small bits of code/comments
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.3
 *
 * Stanford University Unstructured (SU2) Code
 * Copyright (C) 2012 Aerospace Design Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"

void CSourcePieceWise_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourcePieceWise_Flow__SetResidual



//SU2_DIFF END CSourcePieceWise_Flow__SetResidual

}

void CSourcePieceWise_AdjDiscTurb::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourcePieceWise_Turb__SetResidual



//SU2_DIFF END CSourcePieceWise_Turb__SetResidual

}

void CSourcePieceWise_AdjDiscTurbSA::SetDirectResidual_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourcePieceWise_TurbSA__SetResidual

    double result1;
    double result1d;
    double x2;
    double x2d;
    double x1;
    double x1d;
    Density_id = U_id[0];
    Density_i = U_i[0];
    val_residuald[0] = 0.0;
    val_residual[0] = 0.0;
    /*--- Computation of divergence of velocity and vorticity ---*/
    DivVelocity = 0;
    for (iDim = 0; iDim < nDim; ++iDim)
        DivVelocity += PrimVar_Grad_i[iDim + 1][iDim];
    Vorticityd = (PrimVar_Grad_id[2][0]-PrimVar_Grad_id[1][1])*(PrimVar_Grad_i
        [2][0]-PrimVar_Grad_i[1][1]) + (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1]
        [1])*(PrimVar_Grad_id[2][0]-PrimVar_Grad_id[1][1]);
    Vorticity = (PrimVar_Grad_i[2][0]-PrimVar_Grad_i[1][1])*(PrimVar_Grad_i[2]
        [0]-PrimVar_Grad_i[1][1]);
    if (nDim == 3) {
        Vorticityd = Vorticityd + (PrimVar_Grad_id[3][1]-PrimVar_Grad_id[2][2]
            )*(PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) + (PrimVar_Grad_i[3]
            [1]-PrimVar_Grad_i[2][2])*(PrimVar_Grad_id[3][1]-PrimVar_Grad_id[2
            ][2]) + (PrimVar_Grad_id[1][2]-PrimVar_Grad_id[3][0])*(
            PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]) + (PrimVar_Grad_i[1][2]
            -PrimVar_Grad_i[3][0])*(PrimVar_Grad_id[1][2]-PrimVar_Grad_id[3][0
            ]);
        Vorticity += (PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2])*(
        PrimVar_Grad_i[3][1]-PrimVar_Grad_i[2][2]) + (PrimVar_Grad_i[1][2]-
        PrimVar_Grad_i[3][0])*(PrimVar_Grad_i[1][2]-PrimVar_Grad_i[3][0]);
    }
    Vorticityd = (Vorticity == 0.0 ? 0.0 : Vorticityd/(2.0*sqrt(Vorticity)));
    Vorticity = sqrt(Vorticity);
    if (dist_i > 0.0) {
        /*--- Production term ---*/
        dist_i_2 = dist_i*dist_i;
        nud = (Laminar_Viscosity_id*Density_i-Laminar_Viscosity_i*Density_id)/
            (Density_i*Density_i);
        nu = Laminar_Viscosity_i/Density_i;
        Jid = (TurbVar_id[0]*nu-TurbVar_i[0]*nud)/(nu*nu);
        Ji = TurbVar_i[0]/nu;
        Ji_2d = Jid*Ji + Ji*Jid;
        Ji_2 = Ji*Ji;
        Ji_3d = Ji_2d*Ji + Ji_2*Jid;
        Ji_3 = Ji_2*Ji;
        fv1d = (Ji_3d*(Ji_3+cv1_3)-Ji_3*Ji_3d)/((Ji_3+cv1_3)*(Ji_3+cv1_3));
        fv1 = Ji_3/(Ji_3+cv1_3);
        fv2d = -((Jid*(1.0+Ji*fv1)-Ji*(Jid*fv1+Ji*fv1d))/((1.0+Ji*fv1)*(1.0+Ji
            *fv1)));
        fv2 = 1.0 - Ji/(1.0+Ji*fv1);
        Omegad = Vorticityd;
        Omega = Vorticity;
        x1d = Omegad + (TurbVar_id[0]*fv2+TurbVar_i[0]*fv2d)/(k2*dist_i_2);
        x1 = Omega + TurbVar_i[0]*fv2/(k2*dist_i_2);
        if (x1 < TURB_EPS) {
            Shat = TURB_EPS;
            Shatd = 0.0;
        } else {
            Shatd = x1d;
            Shat = x1;
        }
        val_residuald[0] = cb1*Volume*(Shatd*TurbVar_i[0]+Shat*TurbVar_id[0]);
        val_residual[0] += cb1*Shat*TurbVar_i[0]*Volume;
        x2d = (TurbVar_id[0]*Shat*k2*dist_i_2-TurbVar_i[0]*k2*dist_i_2*Shatd)/
            (Shat*k2*dist_i_2*(Shat*k2*dist_i_2));
        x2 = TurbVar_i[0]/(Shat*k2*dist_i_2);
        if (x2 > 10.) {
            r = 10.;
            rd = 0.0;
        } else {
            rd = x2d;
            r = x2;
        }
        result1d = pow_d(r, rd, 6., &result1);
        gd = rd + cw2*(result1d-rd);
        g = r + cw2*(result1-r);
        g_6d = pow_d(g, gd, 6., &g_6);
        glimd = pow_d((1.0+cw3_6)/(g_6+cw3_6), -((1.0+cw3_6)*g_6d/((g_6+cw3_6)
            *(g_6+cw3_6))), 1.0/6.0, &glim);
        fwd = gd*glim + g*glimd;
        fw = g*glim;
        val_residuald[0] = val_residuald[0] - Volume*cw1*((fwd*TurbVar_i[0]+fw
            *TurbVar_id[0])*TurbVar_i[0]+fw*TurbVar_i[0]*TurbVar_id[0])/
            dist_i_2;
        val_residual[0] -= cw1*fw*TurbVar_i[0]*TurbVar_i[0]/dist_i_2*Volume;
        /*--- Diffusion term ---*/
        norm2_Grad = 0.0;
        norm2_Gradd = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            norm2_Gradd = norm2_Gradd + TurbVar_Grad_id[0][iDim]*
                TurbVar_Grad_i[0][iDim] + TurbVar_Grad_i[0][iDim]*
                TurbVar_Grad_id[0][iDim];
            norm2_Grad += TurbVar_Grad_i[0][iDim]*TurbVar_Grad_i[0][iDim];
        }
        val_residuald[0] = val_residuald[0] + cb2*Volume*norm2_Gradd/sigma;
        val_residual[0] += cb2/sigma*norm2_Grad*Volume;
    } else
        *val_residuald = 0.0;


//SU2_DIFF END CSourcePieceWise_TurbSA__SetResidual

}

void CSourcePieceWise_AdjDiscElec::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourcePieceWise_Elec__SetResidual



//SU2_DIFF END CSourcePieceWise_Elec__SetResidual

}

void CSourcePieceWise_AdjDiscLevelSet::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourcePieceWise_LevelSet__SetResidual



//SU2_DIFF END CSourcePieceWise_LevelSet__SetResidual

}

void CSourceConservative_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourceConservative_Flow__SetResidual



//SU2_DIFF END CSourceConservative_Flow__SetResidual

}

void CSourceConservative_AdjDiscTurb::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourceConservative_Turb__SetResidual



//SU2_DIFF END CSourceConservative_Turb__SetResidual

}

void CSourceRotationalFrame_AdjDiscFlow::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CSourceRotationalFrame_Flow__SetResidual



//SU2_DIFF END CSourceRotationalFrame_Flow__SetResidual

}

void CSourcePieceWise_Plasma::SetResidual_Axisymmetric_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//
	
//SU2_DIFF START CSourcePieceWise_Plasma__SetResidual_Axisymmetric

    double yinv, Gamma, Density_i, Energy_i, Energy_vib_i, Energy_el_i, 
    Pressure_i, Enthalpy_i, Enthalpy_formation_i, Velocity_i, sq_vel;
    double Density_id, Energy_id, Energy_vib_id, Pressure_id, Enthalpy_id, 
    Velocity_id, sq_veld;
    unsigned short int iDim, iSpecies, loc;

    if (Coord_i[1] > 0.0) {
        yinv = 1.0/Coord_i[1];
        *val_residuald = 0.0;
    } else {
        yinv = 0.0;
        *val_residuald = 0.0;
    }
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            loc = (nDim+3)*iSpecies;
        else
            loc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        sq_vel = 0.0;
        sq_veld = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            Velocity_id = (U_id[loc+iDim+1]*U_i[loc+0]-U_i[loc+iDim+1]*U_id[
                loc+0])/(U_i[loc+0]*U_i[loc+0]);
            Velocity_i = U_i[loc+iDim+1]/U_i[loc+0];
            sq_veld = sq_veld + Velocity_id*Velocity_i + Velocity_i*
                Velocity_id;
            sq_vel += Velocity_i*Velocity_i;
        }
        Gamma = config->GetSpecies_Gamma(iSpecies);
        Enthalpy_formation_i = config->GetEnthalpy_Formation(iSpecies);
        Density_id = U_id[loc + 0];
        Density_i = U_i[loc + 0];
        Energy_id = (U_id[loc+nDim+1]*Density_i-U_i[loc+nDim+1]*
            Density_id)/(Density_i*Density_i);
        Energy_i = U_i[loc+nDim+1]/Density_i;
        Energy_vib_i = 0.0;
        if (iSpecies < nDiatomics) {
            Energy_vib_id = (U_id[loc+nDim+2]*Density_i-U_i[loc+nDim+2]*
                Density_id)/(Density_i*Density_i);
            Energy_vib_i = U_i[loc+nDim+2]/Density_i;
        } else
            Energy_vib_id = 0.0;
        Energy_el_i = 0.0;
        Pressure_id = (Gamma-1.0)*(Density_id*(Energy_i-1.0/2.0*sq_vel-
            Enthalpy_formation_i-Energy_vib_i-Energy_el_i)+Density_i*(
            Energy_id-sq_veld/2.0-Energy_vib_id));
        Pressure_i = (Gamma-1.0)*Density_i*(Energy_i-1.0/2.0*sq_vel-
            Enthalpy_formation_i-Energy_vib_i-Energy_el_i);
        Enthalpy_id = ((U_id[loc+nDim+1]+Pressure_id)*Density_i-(U_i[loc+
            nDim+1]+Pressure_i)*Density_id)/(Density_i*Density_i);
        Enthalpy_i = (U_i[loc+nDim+1]+Pressure_i)/Density_i;
        val_residuald[loc + 0] = yinv*Volume*U_id[loc+2];
        val_residual[loc + 0] = yinv*Volume*U_i[loc+2];
        val_residuald[loc + 1] = (yinv*Volume*(U_id[loc+1]*U_i[loc+2]+U_i[loc+
            1]*U_id[loc+2])*U_i[loc+0]-yinv*Volume*U_i[loc+1]*U_i[loc+2]*U_id[
            loc+0])/(U_i[loc+0]*U_i[loc+0]);
        val_residual[loc + 1] = yinv*Volume*U_i[loc+1]*U_i[loc+2]/U_i[loc+0];
        val_residuald[loc + 2] = (yinv*Volume*(U_id[loc+2]*U_i[loc+2]+U_i[loc+
            2]*U_id[loc+2])*U_i[loc+0]-yinv*Volume*(U_i[loc+2]*U_i[loc+2])*
            U_id[loc+0])/(U_i[loc+0]*U_i[loc+0]);
        val_residual[loc + 2] = yinv*Volume*U_i[loc+2]*U_i[loc+2]/U_i[loc+0];
        val_residuald[loc + 3] = yinv*Volume*(Enthalpy_id*U_i[loc+2]+
            Enthalpy_i*U_id[loc+2]);
        val_residual[loc + 3] = yinv*Volume*Enthalpy_i*U_i[loc+2];
        if (iSpecies < nDiatomics) {
            val_residuald[loc + 4] = (yinv*Volume*(U_id[loc+4]*U_i[loc+2]+U_i[
                loc+4]*U_id[loc+2])*U_i[loc+0]-yinv*Volume*U_i[loc+4]*U_i[loc+
                2]*U_id[loc+0])/(U_i[loc+0]*U_i[loc+0]);
            val_residual[loc + 4] = yinv*Volume*U_i[loc+4]*U_i[loc+2]/U_i[loc+
                0];
        }
    }


//SU2_DIFF END CSourcePieceWise_Plasma__SetResidual_Axisymmetric
}  



void CSourcePieceWise_Plasma::SetResidual_Chemistry_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//
	
//SU2_DIFF START CSourcePieceWise_Plasma__SetResidual_Chemistry

    double T_min;
    // Minimum temperature for the modified temperature calculations.
    double epsilon;
    // Parameter for the modified temperature calculations.
    unsigned short int iSpecies, jSpecies, iReaction, iVar, iDim, ii;
    unsigned short int iLoc, jLoc;
    unsigned short int counterFwd, counterBkw;
    double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, 
    Gas_constant;
    double Energy_vibd, Vel2d;
    float arg1;
    float arg1d;
    double arg10;
    double arg10d;
    float result1;
    float result1d;
    double arg2;
    double arg2d;
    int ii1;
    T_min = 800;
    epsilon = 80;
  
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        Temp_tr_id[ii1] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        Energy_vib = 0.0;
        Energy_el = 0.0;
        Vel2 = 0.0;
        Vel2d = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            Vel2d = Vel2d + (((U_id[iLoc+iDim+1]*U_i[iLoc+0]-U_i[iLoc+iDim+1]*
                U_id[iLoc+0])*U_i[iLoc+iDim+1]/(U_i[iLoc+0]*U_i[iLoc+0])+U_i[
                iLoc+iDim+1]*U_id[iLoc+iDim+1]/U_i[iLoc+0])*U_i[iLoc+0]-U_i[
                iLoc+iDim+1]*U_i[iLoc+iDim+1]*U_id[iLoc+0]/U_i[iLoc+0])/(U_i[
                iLoc+0]*U_i[iLoc+0]);
            Vel2 += U_i[iLoc+iDim+1]/U_i[iLoc+0]*U_i[iLoc+iDim+1]/U_i[iLoc+0];
        }
        if (iSpecies < nDiatomics) {
            Energy_vibd = (U_id[iLoc+nDim+2]*U_i[iLoc+0]-U_i[iLoc+nDim+2
                ]*U_id[iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]);
            Energy_vib = U_i[iLoc+nDim+2]/U_i[iLoc+0];
        } else
            Energy_vibd = 0.0;
            Gamma = config->GetSpecies_Gamma(iSpecies);
        		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
        		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
        Temp_tr_id[iSpecies] = (Gamma-1.0)*((U_id[iLoc+nDim+1]*U_i[iLoc]-
            U_i[iLoc+nDim+1]*U_id[iLoc])/(U_i[iLoc]*U_i[iLoc])-0.5*Vel2d-
            Energy_vibd)/Gas_constant;
        Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant*(U_i[iLoc+nDim+1]/
            U_i[iLoc]-0.5*Vel2-Enthalpy_formation-Energy_vib-Energy_el);
    }
    /*--- Initialize all components of the residual and Jacobian to zero ---*/
    for (iVar = 0; iVar < nVar; ++iVar) {
        val_residuald[iVar] = 0.0;
        val_residual[iVar] = 0.0;
    }
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        bkwRxnd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        Keqd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        fwdRxnd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        ReactionRateFwdd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        ReactionRateBkwd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        T_rxnbd[ii1] = 0.0;
    for (ii1 = 0; ii1 < nReactions; ++ii1)
        T_rxnfd[ii1] = 0.0;
    for (iReaction = 0; iReaction < nReactions; ++iReaction) {
        /*--- Calculate the rate-controlling temperatures ---*/
        //NOTE:  This implementation takes the geometric mean of the TR temps of all particpants.
        //       This is NOT consistent with the Park chemical models, and is only a first-cut approach.
        T_rxnfd[iReaction] = 0.0;
        T_rxnf[iReaction] = 1.0;
        T_rxnbd[iReaction] = 0.0;
        T_rxnb[iReaction] = 1.0;
        counterFwd = 0;
        counterBkw = 0;
        for (ii = 0; ii < 3; ++ii) {
            iSpecies = Reactions[iReaction][0][ii];
            jSpecies = Reactions[iReaction][1][ii];
            /*--- Reactants ---*/
            if (iSpecies != nSpecies) {
                T_rxnfd[iReaction] = T_rxnfd[iReaction]*Temp_tr_i[iSpecies] + 
                    T_rxnf[iReaction]*Temp_tr_id[iSpecies];
                T_rxnf[iReaction] *= Temp_tr_i[iSpecies];
                counterFwd++;
            }
            /*--- Products ---*/
            if (jSpecies != nSpecies) {
                T_rxnbd[iReaction] = T_rxnbd[iReaction]*Temp_tr_i[jSpecies] + 
                    T_rxnb[iReaction]*Temp_tr_id[jSpecies];
                T_rxnb[iReaction] *= Temp_tr_i[jSpecies];
                counterBkw++;
            }
        }
        arg1d = T_rxnfd[iReaction]/(counterFwd*T_rxnf[iReaction]);
        arg1 = 1.0/counterFwd*log(T_rxnf[iReaction]);
        T_rxnfd[iReaction] = arg1d*exp(arg1);
        T_rxnf[iReaction] = exp(arg1);
        arg1d = T_rxnbd[iReaction]/(counterBkw*T_rxnb[iReaction]);
        arg1 = 1.0/counterBkw*log(T_rxnb[iReaction]);
        T_rxnbd[iReaction] = arg1d*exp(arg1);
        T_rxnb[iReaction] = exp(arg1);
        /*--- Apply a modified temperature to ease the stiffness at low temperatures ---
        */
        arg10d = T_rxnfd[iReaction]*(T_rxnf[iReaction]-T_min) + (T_rxnf[
            iReaction]-T_min)*T_rxnfd[iReaction];
        arg10 = (T_rxnf[iReaction]-T_min)*(T_rxnf[iReaction]-T_min) + epsilon*
            epsilon;
        result1d = (arg10 == 0.0 ? 0.0 : arg10d/(2.0*sqrt(arg10)));
        result1 = sqrt(arg10);
        T_rxnfd[iReaction] = 0.5*(T_rxnfd[iReaction]+result1d);
        T_rxnf[iReaction] = 0.5*(T_rxnf[iReaction]+T_min+result1);
        arg10d = T_rxnbd[iReaction]*(T_rxnb[iReaction]-T_min) + (T_rxnb[
            iReaction]-T_min)*T_rxnbd[iReaction];
        arg10 = (T_rxnb[iReaction]-T_min)*(T_rxnb[iReaction]-T_min) + epsilon*
            epsilon;
        result1d = (arg10 == 0.0 ? 0.0 : arg10d/(2.0*sqrt(arg10)));
        result1 = sqrt(arg10);
        T_rxnbd[iReaction] = 0.5*(T_rxnbd[iReaction]+result1d);
        T_rxnb[iReaction] = 0.5*(T_rxnb[iReaction]+T_min+result1);
        /*--- Calculate equilibrium extent of reaction ---*/
        		GetEq_Rxn_Coefficients(EqRxnConstants, config);
        //NOTE: Scalabrin implementation
        arg10d = EqRxnConstants[iReaction][0]*T_rxnbd[iReaction]/10000.0 - 
            EqRxnConstants[iReaction][2]*T_rxnbd[iReaction]/T_rxnb[iReaction] 
            - EqRxnConstants[iReaction][3]*10000.0*T_rxnbd[iReaction]/(T_rxnb[
            iReaction]*T_rxnb[iReaction]) - EqRxnConstants[iReaction][4]*2*(
            10000.0*10000.0)*T_rxnbd[iReaction]/pow(T_rxnb[iReaction],3);
        arg10 = EqRxnConstants[iReaction][0]*(T_rxnb[iReaction]/10000.0) + 
            EqRxnConstants[iReaction][1] + EqRxnConstants[iReaction][2]*log(
            10000.0/T_rxnb[iReaction]) + EqRxnConstants[iReaction][3]*(10000.0
            /T_rxnb[iReaction]) + EqRxnConstants[iReaction][4]*(10000.0/T_rxnb
            [iReaction])*(10000.0/T_rxnb[iReaction]);
        Keqd[iReaction] = arg10d*exp(arg10);
        Keq[iReaction] = exp(arg10);
        /*--- Calculate reaction rate coefficients ---*/
        arg10d = eta[iReaction]*T_rxnfd[iReaction]/T_rxnf[iReaction];
        arg10 = eta[iReaction]*log(T_rxnf[iReaction]);
        arg2d = -(-(theta[iReaction]*T_rxnfd[iReaction])/(T_rxnf[iReaction]*
            T_rxnf[iReaction]));
        arg2 = -theta[iReaction]/T_rxnf[iReaction];
        ReactionRateFwdd[iReaction] = Cf[iReaction]*(arg10d*exp(arg10)*exp(
            arg2)+exp(arg10)*arg2d*exp(arg2));
        ReactionRateFwd[iReaction] = Cf[iReaction]*exp(arg10)*exp(arg2);
        arg10d = eta[iReaction]*T_rxnbd[iReaction]/T_rxnb[iReaction];
        arg10 = eta[iReaction]*log(T_rxnb[iReaction]);
        arg2d = -(-(theta[iReaction]*T_rxnbd[iReaction])/(T_rxnb[iReaction]*
            T_rxnb[iReaction]));
        arg2 = -theta[iReaction]/T_rxnb[iReaction];
        ReactionRateBkwd[iReaction] = (Cf[iReaction]*(arg10d*exp(arg10)*exp(
            arg2)+exp(arg10)*arg2d*exp(arg2))*Keq[iReaction]-Cf[iReaction]*exp
            (arg10)*exp(arg2)*Keqd[iReaction])/(Keq[iReaction]*Keq[iReaction])
        ;
        ReactionRateBkw[iReaction] = Cf[iReaction]*exp(arg10)*exp(arg2)/Keq[
            iReaction];
        //		ReactionRateBkw[iReaction] = ReactionRateFwd[iReaction] / Keq[iReaction];
        fwdRxnd[iReaction] = 0.0;
        fwdRxn[iReaction] = 1.0;
        bkwRxnd[iReaction] = 0.0;
        bkwRxn[iReaction] = 1.0;
        for (ii = 0; ii < 3; ++ii) {
            /*--- Reactants ---*/
            iSpecies = Reactions[iReaction][0][ii];
            if (iSpecies != nSpecies) {
                if (iSpecies < nDiatomics)
                    iLoc = (nDim+3)*iSpecies;
                else
                    iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics
                        );
                fwdRxnd[iReaction] = fwdRxnd[iReaction]*0.001*U_i[iLoc+0]/
                    Molar_Mass[iSpecies] + fwdRxn[iReaction]*0.001*U_id[iLoc+0
                    ]/Molar_Mass[iSpecies];
                fwdRxn[iReaction] *= 0.001*U_i[iLoc+0]/Molar_Mass[iSpecies];
            }
            /*--- Products ---*/
            jSpecies = Reactions[iReaction][1][ii];
            if (jSpecies != nSpecies) {
                if (jSpecies < nDiatomics)
                    jLoc = (nDim+3)*jSpecies;
                else
                    jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics
                        );
                bkwRxnd[iReaction] = bkwRxnd[iReaction]*0.001*U_i[jLoc+0]/
                    Molar_Mass[jSpecies] + bkwRxn[iReaction]*0.001*U_id[jLoc+0
                    ]/Molar_Mass[jSpecies];
                bkwRxn[iReaction] *= 0.001*U_i[jLoc+0]/Molar_Mass[jSpecies];
            }
        }
        fwdRxnd[iReaction] = 1000.0*(ReactionRateFwdd[iReaction]*fwdRxn[
            iReaction]+ReactionRateFwd[iReaction]*fwdRxnd[iReaction]);
        fwdRxn[iReaction] = 1000.0*ReactionRateFwd[iReaction]*fwdRxn[iReaction
            ];
        bkwRxnd[iReaction] = 1000.0*(ReactionRateBkwd[iReaction]*bkwRxn[
            iReaction]+ReactionRateBkw[iReaction]*bkwRxnd[iReaction]);
        bkwRxn[iReaction] = 1000.0*ReactionRateBkw[iReaction]*bkwRxn[iReaction
            ];
    }
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        w_dotd[iSpecies] = 0.0;
        w_dot[iSpecies] = 0.0;
    }
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        w_dotd[ii1] = 0.0;
    for (iReaction = 0; iReaction < nReactions; ++iReaction)
        for (ii = 0; ii < 3; ++ii) {
            /*--- Products ---*/
            iSpecies = Reactions[iReaction][1][ii];
            if (iSpecies != nSpecies) {
                w_dotd[iSpecies] = w_dotd[iSpecies] + Molar_Mass[iSpecies]*(
                    fwdRxnd[iReaction]-bkwRxnd[iReaction]);
                w_dot[iSpecies] += Molar_Mass[iSpecies]*(fwdRxn[iReaction]-
                bkwRxn[iReaction]);
            }
            /*--- Reactants ---*/
            iSpecies = Reactions[iReaction][0][ii];
            if (iSpecies != nSpecies) {
                w_dotd[iSpecies] = w_dotd[iSpecies] - Molar_Mass[iSpecies]*(
                    fwdRxnd[iReaction]-bkwRxnd[iReaction]);
                w_dot[iSpecies] -= Molar_Mass[iSpecies]*(fwdRxn[iReaction]-
                bkwRxn[iReaction]);
            }
        }
    *val_residuald = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        val_residuald[iLoc] = Volume*w_dotd[iSpecies];
        val_residual[iLoc] = w_dot[iSpecies]*Volume;
        for (iDim = 0; iDim < nDim; ++iDim) {
            val_residuald[iLoc + 1 + iDim] = Volume*((w_dotd[iSpecies]*U_i[
                iLoc+1+iDim]+w_dot[iSpecies]*U_id[iLoc+1+iDim])*U_i[iLoc]-
                w_dot[iSpecies]*U_i[iLoc+1+iDim]*U_id[iLoc])/(U_i[iLoc]*U_i[
                iLoc]);
            val_residual[iLoc + 1 + iDim] = w_dot[iSpecies]*U_i[iLoc+1+iDim]/
                U_i[iLoc]*Volume;
        }
        val_residuald[iLoc + 1 + nDim] = Volume*((w_dotd[iSpecies]*U_i[iLoc
            +1+nDim]+w_dot[iSpecies]*U_id[iLoc+1+nDim])*U_i[iLoc]-w_dot[
            iSpecies]*U_i[iLoc+1+nDim]*U_id[iLoc])/(U_i[iLoc]*U_i[iLoc]);
        val_residual[iLoc + 1 + nDim] = w_dot[iSpecies]*U_i[iLoc+1+nDim]
            /U_i[iLoc]*Volume;
        if (iSpecies < nDiatomics) {
            val_residuald[iLoc + 2 + nDim] = Volume*((w_dotd[iSpecies]*U_i[
                iLoc+2+nDim]+w_dot[iSpecies]*U_id[iLoc+2+nDim])*U_i[iLoc
                ]-w_dot[iSpecies]*U_i[iLoc+2+nDim]*U_id[iLoc])/(U_i[iLoc]*
                U_i[iLoc]);
            val_residual[iLoc + 2 + nDim] = w_dot[iSpecies]*U_i[iLoc+2+
                nDim]/U_i[iLoc]*Volume;
        }
    }


//SU2_DIFF END CSourcePieceWise_Plasma__SetResidual_Chemistry
}

void CSourcePieceWise_Plasma::SetResidual_MomentumExch_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//
	
//SU2_DIFF START CSourcePieceWise_Plasma__SetResidual_MomentumExch

    unsigned short int iDim, iSpecies, jSpecies, iVar, iLoc, jLoc;
    double collisionFreq, collisionArea, velocity_Species_i, 
    velocity_Species_j, collisionVelocity;
    double collisionFreqd, collisionAread, velocity_Species_id, 
    velocity_Species_jd, collisionVelocityd;
    double T_control;
    double T_controld;
    double radius_electronIonCollision, electron_temperature;
    double radius_electronIonCollisiond, electron_temperatured;
    double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, 
    Gas_constant;
    double Energy_vibd, Vel2d;
    double arg1;
    double arg1d;
    int ii2;
    int ii1;

    /*--- Initialization ---*/
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies)
        for (iDim = 0; iDim < nDim; ++iDim) {
            Pd[iSpecies][iDim] = 0.0;
            P[iSpecies][iDim] = 0.0;
        }
    for (iVar = 0; iVar < nVar; ++iVar) {
        val_residuald[iVar] = 0.0;
        val_residual[iVar] = 0.0;
    }
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        Temp_tr_id[ii1] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        Energy_vib = 0.0;
        Energy_el = 0.0;
        Vel2 = 0.0;
        Vel2d = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            Vel2d = Vel2d + (((U_id[iLoc+iDim+1]*U_i[iLoc+0]-U_i[iLoc+iDim+1]*
                U_id[iLoc+0])*U_i[iLoc+iDim+1]/(U_i[iLoc+0]*U_i[iLoc+0])+U_i[
                iLoc+iDim+1]*U_id[iLoc+iDim+1]/U_i[iLoc+0])*U_i[iLoc+0]-U_i[
                iLoc+iDim+1]*U_i[iLoc+iDim+1]*U_id[iLoc+0]/U_i[iLoc+0])/(U_i[
                iLoc+0]*U_i[iLoc+0]);
            Vel2 += U_i[iLoc+iDim+1]/U_i[iLoc+0]*U_i[iLoc+iDim+1]/U_i[iLoc+0];
        }
        if (iSpecies < nDiatomics) {
            Energy_vibd = (U_id[iLoc+nDim+2]*U_i[iLoc+0]-U_i[iLoc+nDim+2
                ]*U_id[iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]);
            Energy_vib = U_i[iLoc+nDim+2]/U_i[iLoc+0];
        } else
            Energy_vibd = 0.0;
        		Gamma = config->GetSpecies_Gamma(iSpecies);
        		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
        		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
        Temp_tr_id[iSpecies] = (Gamma-1.0)*((U_id[iLoc+nDim+1]*U_i[iLoc]-
            U_i[iLoc+nDim+1]*U_id[iLoc])/(U_i[iLoc]*U_i[iLoc])-0.5*Vel2d-
            Energy_vibd)/Gas_constant;
        Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant*(U_i[iLoc+nDim+1]/
            U_i[iLoc]-0.5*Vel2-Enthalpy_formation-Energy_vib-Energy_el);
    }
    electron_temperatured = Temp_tr_id[nSpecies - 1];
    electron_temperature = Temp_tr_i[nSpecies - 1];
    *val_residuald = 0.0;
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        for (ii2 = 0; ii2 < nDim; ++ii2)
            Pd[ii1][ii2] = 0.0;
    /*--- Solve for momentum exchange between species due to collisions ---*/
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        for (jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
            if (jSpecies < nDiatomics)
                jLoc = (nDim+3)*jSpecies;
            else
                jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
            if (iSpecies != jSpecies) {
                /*				collisionArea = PI_NUMBER * ((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0)

																																	* ((Molecular_Diameter[iSpecies] + Molecular_Diameter[jSpecies])/2.0); 
                */
                arg1d = Temp_tr_id[iSpecies]*Temp_tr_i[jSpecies] + Temp_tr_i[
                    iSpecies]*Temp_tr_id[jSpecies];
                arg1 = Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies];
                T_controld = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
                T_control = sqrt(arg1);
                /*collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 

                      * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 

                            + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
                */
                arg1d = (Omega11[iSpecies][jSpecies][0]*(T_controld*log(
                    T_control)/T_control+log(T_control)*T_controld/T_control)+
                    Omega11[iSpecies][jSpecies][1]*T_controld/T_control)*log(
                    T_control) + (Omega11[iSpecies][jSpecies][0]*log(T_control
                    )*log(T_control)+Omega11[iSpecies][jSpecies][1]*log(
                    T_control)+Omega11[iSpecies][jSpecies][2])*T_controld/
                    T_control;
                arg1 = (Omega11[iSpecies][jSpecies][0]*log(T_control)*log(
                    T_control)+Omega11[iSpecies][jSpecies][1]*log(T_control)+
                    Omega11[iSpecies][jSpecies][2])*log(T_control);
                collisionAread = 1E-20*Omega11[iSpecies][jSpecies][3]*arg1d*
                    exp(arg1);
                collisionArea = 1E-20*Omega11[iSpecies][jSpecies][3]*exp(arg1);
                /*----- An if-condition for the case when a positive charge collides with an electron -----
                */
                if ((ChargeNumber[iSpecies] == 1 && jSpecies == nSpecies - 1) || 
                (iSpecies == nSpecies - 1 && ChargeNumber[jSpecies] == 1)) {
                    radius_electronIonCollisiond = -(ELECTRON_CHARGE*ELECTRON_CHARGE*32.0*
                        FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperatured/(32.0*FREE_PERMITTIVITY*
                        BOLTZMANN_CONSTANT*electron_temperature*(32.0*FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*
                        electron_temperature)));
                    radius_electronIonCollision = ELECTRON_CHARGE*ELECTRON_CHARGE/(32.0*
                        FREE_PERMITTIVITY*BOLTZMANN_CONSTANT*electron_temperature);
                    collisionAread = PI_NUMBER*(radius_electronIonCollisiond*
                        radius_electronIonCollision+
                        radius_electronIonCollision*
                        radius_electronIonCollisiond);
                    collisionArea = PI_NUMBER*radius_electronIonCollision*
                        radius_electronIonCollision;
                }
                arg1d = 8.0*BOLTZMANN_CONSTANT*Temp_tr_id[iSpecies]/(PI_NUMBER*
                    Molecular_Mass[iSpecies]);
                arg1 = 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass
                    [iSpecies]);
                velocity_Species_id = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(
                    arg1)));
                velocity_Species_i = sqrt(arg1);
                arg1d = 8.0*BOLTZMANN_CONSTANT*Temp_tr_id[jSpecies]/(PI_NUMBER*
                    Molecular_Mass[jSpecies]);
                arg1 = 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass
                    [jSpecies]);
                velocity_Species_jd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(
                    arg1)));
                velocity_Species_j = sqrt(arg1);
                arg1d = velocity_Species_id*velocity_Species_i + 
                    velocity_Species_i*velocity_Species_id + 
                    velocity_Species_jd*velocity_Species_j + 
                    velocity_Species_j*velocity_Species_jd;
                arg1 = velocity_Species_i*velocity_Species_i + 
                    velocity_Species_j*velocity_Species_j;
                collisionVelocityd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1
                    )));
                collisionVelocity = sqrt(arg1);
                collisionFreqd = (U_id[jLoc+0]*collisionArea*collisionVelocity
                    +U_i[jLoc+0]*(collisionAread*collisionVelocity+
                    collisionArea*collisionVelocityd))/(Molecular_Mass[
                    iSpecies]+Molecular_Mass[jSpecies]);
                collisionFreq = U_i[jLoc+0]*collisionArea*collisionVelocity/(
                    Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
                for (iDim = 0; iDim < nDim; ++iDim) {
                    Pd[iSpecies][iDim] = Pd[iSpecies][iDim] + (U_id[iLoc+0]*
                        collisionFreq+U_i[iLoc+0]*collisionFreqd)*(U_i[jLoc+1+
                        iDim]/U_i[jLoc+0]-U_i[iLoc+1+iDim]/U_i[iLoc+0]) + U_i[
                        iLoc+0]*collisionFreq*((U_id[jLoc+1+iDim]*U_i[jLoc+0]-
                        U_i[jLoc+1+iDim]*U_id[jLoc+0])/(U_i[jLoc+0]*U_i[jLoc+0
                        ])-(U_id[iLoc+1+iDim]*U_i[iLoc+0]-U_i[iLoc+1+iDim]*
                        U_id[iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]));
                    P[iSpecies][iDim] += U_i[iLoc+0]*collisionFreq*(U_i[jLoc+1
                    +iDim]/U_i[jLoc+0]-U_i[iLoc+1+iDim]/U_i[iLoc+0]);
                }
            }
        }
        for (iDim = 0; iDim < nDim; ++iDim) {
            val_residuald[iLoc + 1 + iDim] = Volume*Pd[iSpecies][iDim];
            val_residual[iLoc + 1 + iDim] = P[iSpecies][iDim]*Volume;
            val_residuald[iLoc + nDim + 1] = val_residuald[iLoc + nDim +
                1] + Volume*((Pd[iSpecies][iDim]*U_i[iLoc+1+iDim]+P[iSpecies][
                iDim]*U_id[iLoc+1+iDim])*U_i[iLoc+0]-P[iSpecies][iDim]*U_i[
                iLoc+1+iDim]*U_id[iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]);
            val_residual[iLoc + nDim + 1] += P[iSpecies][iDim]*U_i[iLoc+1+
            iDim]/U_i[iLoc+0]*Volume;
        }
    }


//SU2_DIFF END CSourcePieceWise_Plasma__SetResidual_MomentumExch	
}

void CSourcePieceWise_Plasma::SetResidual_EnergyExch_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//
	
//SU2_DIFF START CSourcePieceWise_Plasma__SetResidual_EnergyExch

    unsigned short int iDim, iSpecies, jSpecies, iLoc, jLoc, iVar;
    double collisionFreq, collisionArea, velocity_Species_i, 
    velocity_Species_j, collisionVelocity, vel_dot_prod;
    double collisionFreqd, collisionAread, velocity_Species_id, 
    velocity_Species_jd, collisionVelocityd;
    double coefficient;
    double coefficientd;
    double T_control;
    double T_controld;
    double Energy_vib, Energy_el, Vel2, Gamma, Enthalpy_formation, 
    Gas_constant;
    double Energy_vibd, Vel2d;
    double arg1;
    double arg1d;
    float arg10;
    float arg10d;
    double result1;
    double result1d;
    float arg2;
    float arg2d;
    float arg3;
    double arg4;
    double arg4d;
    double Q_tvd[nSpecies];
    int ii1;

    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        CharVibTempd[iSpecies] = 0.0;
    }
  
    /*--- Energy transport via elastic collisions ---
Comment: Two terms, both from kinetic theory.  The first term accounts for changes in species temperature from

	 elastic encounters between particles, the second term accounts for frictional heating between species
    */
    //Comment: From Lee 1985, and originally from Engineering Magnetohydrodynamics by Sutton (1965)
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        Q_elasticd[iSpecies] = 0.0;
        Q_elastic[iSpecies] = 0.0;
        Q_tvd[iSpecies] = 0.0;
        Q_tv[iSpecies] = 0.0;
    }
    for (iVar = 0; iVar < nVar; ++iVar) {
        val_residuald[iVar] = 0.0;
        val_residual[iVar] = 0.0;
    }
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        SpeciesPressure_id[ii1] = 0.0;
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        Temp_tr_id[ii1] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        Energy_vib = 0.0;
        Energy_el = 0.0;
        Vel2 = 0.0;
        Vel2d = 0.0;
        for (iDim = 0; iDim < nDim; ++iDim) {
            Vel2d = Vel2d + (((U_id[iLoc+iDim+1]*U_i[iLoc+0]-U_i[iLoc+iDim+1]*
                U_id[iLoc+0])*U_i[iLoc+iDim+1]/(U_i[iLoc+0]*U_i[iLoc+0])+U_i[
                iLoc+iDim+1]*U_id[iLoc+iDim+1]/U_i[iLoc+0])*U_i[iLoc+0]-U_i[
                iLoc+iDim+1]*U_i[iLoc+iDim+1]*U_id[iLoc+0]/U_i[iLoc+0])/(U_i[
                iLoc+0]*U_i[iLoc+0]);
            Vel2 += U_i[iLoc+iDim+1]/U_i[iLoc+0]*U_i[iLoc+iDim+1]/U_i[iLoc+0];
        }
        if (iSpecies < nDiatomics) {
            Energy_vibd = (U_id[iLoc+nDim+2]*U_i[iLoc+0]-U_i[iLoc+nDim+2
                ]*U_id[iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]);
            Energy_vib = U_i[iLoc+nDim+2]/U_i[iLoc+0];
        } else
            Energy_vibd = 0.0;
        		Gamma = config->GetSpecies_Gamma(iSpecies);
        		Enthalpy_formation = config->GetEnthalpy_Formation(iSpecies);
        		Gas_constant = config->GetSpecies_Gas_Constant(iSpecies);    
        Temp_tr_id[iSpecies] = (Gamma-1.0)*((U_id[iLoc+nDim+1]*U_i[iLoc]-
            U_i[iLoc+nDim+1]*U_id[iLoc])/(U_i[iLoc]*U_i[iLoc])-0.5*Vel2d-
            Energy_vibd)/Gas_constant;
        Temp_tr_i[iSpecies] = (Gamma-1.0)/Gas_constant*(U_i[iLoc+nDim+1]/
            U_i[iLoc]-0.5*Vel2-Enthalpy_formation-Energy_vib-Energy_el);
        SpeciesPressure_id[iSpecies] = (Gamma-1.0)*(U_id[iLoc]*(U_i[iLoc+
            nDim+1]/U_i[iLoc]-0.5*Vel2-Enthalpy_formation-Energy_vib-
            Energy_el)+U_i[iLoc]*((U_id[iLoc+nDim+1]*U_i[iLoc]-U_i[iLoc+
            nDim+1]*U_id[iLoc])/(U_i[iLoc]*U_i[iLoc])-0.5*Vel2d-Energy_vibd
            ));
        SpeciesPressure_i[iSpecies] = (Gamma-1.0)*U_i[iLoc]*(U_i[iLoc+nDim+
            1]/U_i[iLoc]-0.5*Vel2-Enthalpy_formation-Energy_vib-Energy_el);
    }
    *val_residuald = 0.0;
    for (ii1 = 0; ii1 < nSpecies; ++ii1)
        Q_elasticd[ii1] = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
        if (iSpecies < nDiatomics)
            iLoc = (nDim+3)*iSpecies;
        else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        for (jSpecies = 0; jSpecies < nSpecies; ++jSpecies)
            if (jSpecies != iSpecies) {
                if (jSpecies < nDiatomics)
                    jLoc = (nDim+3)*jSpecies;
                else
                    jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics
                        );
                /*				collisionArea = PI_NUMBER * ((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0)

																																	* ((Molecular_Diameter[iSpecies]+Molecular_Diameter[jSpecies])/2.0);
                */
                arg1d = Temp_tr_id[iSpecies]*Temp_tr_i[jSpecies] + Temp_tr_i[
                    iSpecies]*Temp_tr_id[jSpecies];
                arg1 = Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies];
                T_controld = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
                T_control = sqrt(arg1);
                /*        collisionArea = 1E-20 * Omega11[iSpecies][jSpecies][3] 

                        * pow(T_control, Omega11[iSpecies][jSpecies][0]*log(T_control)*log(T_control) 

                              + Omega11[iSpecies][jSpecies][1]*log(T_control) + Omega11[iSpecies][jSpecies][2]);
                */
                arg1d = (Omega11[iSpecies][jSpecies][0]*(T_controld*log(
                    T_control)/T_control+log(T_control)*T_controld/T_control)+
                    Omega11[iSpecies][jSpecies][1]*T_controld/T_control)*log(
                    T_control) + (Omega11[iSpecies][jSpecies][0]*log(T_control
                    )*log(T_control)+Omega11[iSpecies][jSpecies][1]*log(
                    T_control)+Omega11[iSpecies][jSpecies][2])*T_controld/
                    T_control;
                arg1 = (Omega11[iSpecies][jSpecies][0]*log(T_control)*log(
                    T_control)+Omega11[iSpecies][jSpecies][1]*log(T_control)+
                    Omega11[iSpecies][jSpecies][2])*log(T_control);
                collisionAread = 1E-20*Omega11[iSpecies][jSpecies][3]*arg1d*
                    exp(arg1);
                collisionArea = 1E-20*Omega11[iSpecies][jSpecies][3]*exp(arg1)
                ;
                arg1d = 8.0*BOLTZMANN_CONSTANT*Temp_tr_id[iSpecies]/(PI_NUMBER*
                    Molecular_Mass[iSpecies]);
                arg1 = 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molecular_Mass
                    [iSpecies]);
                velocity_Species_id = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(
                    arg1)));
                velocity_Species_i = sqrt(arg1);
                arg1d = 8.0*BOLTZMANN_CONSTANT*Temp_tr_id[jSpecies]/(PI_NUMBER*
                    Molecular_Mass[jSpecies]);
                arg1 = 8.0*BOLTZMANN_CONSTANT*Temp_tr_i[jSpecies]/(PI_NUMBER*Molecular_Mass
                    [jSpecies]);
                velocity_Species_jd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(
                    arg1)));
                velocity_Species_j = sqrt(arg1);
                arg1d = velocity_Species_id*velocity_Species_i + 
                    velocity_Species_i*velocity_Species_id + 
                    velocity_Species_jd*velocity_Species_j + 
                    velocity_Species_j*velocity_Species_jd;
                arg1 = velocity_Species_i*velocity_Species_i + 
                    velocity_Species_j*velocity_Species_j;
                collisionVelocityd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1
                    )));
                collisionVelocity = sqrt(arg1);
                collisionFreqd = (U_id[jLoc+0]*collisionArea*collisionVelocity
                    +U_i[jLoc+0]*(collisionAread*collisionVelocity+
                    collisionArea*collisionVelocityd))/(Molecular_Mass[
                    iSpecies]+Molecular_Mass[jSpecies]);
                collisionFreq = U_i[jLoc+0]*collisionArea*collisionVelocity/(
                    Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies]);
                vel_dot_prod = 0.0;
                for (iDim = 0; iDim < nDim; ++iDim)
                    vel_dot_prod += (U_i[iLoc+1+iDim]-U_i[jLoc+1+iDim])*U_i[
                    iLoc+1+iDim];
                /* Exchange from Lee and Sutton, heavy particles */
                coefficientd = BOLTZMANN_CONSTANT*3.0*(2.0*U_id[iLoc+0]*collisionFreq/(
                    Molecular_Mass[iSpecies]+Molecular_Mass[jSpecies])+2.0*U_i
                    [iLoc+0]*collisionFreqd/(Molecular_Mass[iSpecies]+
                    Molecular_Mass[jSpecies]))/2.0;
                coefficient = 2.0*U_i[iLoc+0]/(Molecular_Mass[iSpecies]+
                    Molecular_Mass[jSpecies])*collisionFreq*3.0/2.0*BOLTZMANN_CONSTANT;
                Q_elasticd[iSpecies] = Q_elasticd[iSpecies] + coefficientd*(
                    Temp_tr_i[jSpecies]-Temp_tr_i[iSpecies]) + coefficient*(
                    Temp_tr_id[jSpecies]-Temp_tr_id[iSpecies]);
                Q_elastic[iSpecies] += coefficient*(Temp_tr_i[jSpecies]-
                Temp_tr_i[iSpecies]);
            }
        val_residuald[iLoc + 1 + nDim] = Volume*Q_elasticd[iSpecies];
        val_residual[iLoc + 1 + nDim] = Q_elastic[iSpecies]*Volume;
    }
    /*--- Translational-rotational & vibrational energy exchange via inelastic collisions ---
    */
    //Comment: Landau-Teller formulation
    //Comment: Based on Scalabrin and Candler.  May need to re-visit with much more detail and derive for my formulation
    double tau_sr, tau_ps, LimitingXSection, AvgMolecularSpeed, ReducedMass, 
    A_sr, estar_vs, e_vs, q_tr_vs;
    double tau_srd, tau_psd, LimitingXSectiond, AvgMolecularSpeedd, estar_vsd
    , e_vsd, q_tr_vsd;
    double MixturePressure, MixtureNumDensity;
    double MixtureNumDensityd;
    /*--- Calculate mixture quantities ---*/
    MixturePressure = 0.0;
    MixtureNumDensity = 0.0;
    MixtureNumDensityd = 0.0;
    for (iSpecies = 0; iSpecies < nSpecies; ++iSpecies) {
      if (iSpecies < nDiatomics) {
            iLoc = (nDim+3)*iSpecies;
            CharVibTemp[iSpecies] = config->GetCharVibTemp(iSpecies);
      } else
            iLoc = (nDim+3)*nDiatomics + (nDim+2)*(iSpecies-nDiatomics);
        MixturePressure += SpeciesPressure_i[iSpecies];
        MixtureNumDensityd = MixtureNumDensityd + U_id[iLoc+0]/Molecular_Mass[
            iSpecies];
        MixtureNumDensity += U_i[iLoc+0]/Molecular_Mass[iSpecies];
    }
    for (iSpecies = 0; iSpecies < nDiatomics; ++iSpecies) {
        iLoc = (nDim+3)*iSpecies;
        Q_tvd[iSpecies] = 0.0;
        Q_tv[iSpecies] = 0.0;
        for (jSpecies = 0; jSpecies < nSpecies; ++jSpecies) {
            if (jSpecies < nDiatomics)
                jLoc = (nDim+3)*jSpecies;
            else
                jLoc = (nDim+3)*nDiatomics + (nDim+2)*(jSpecies-nDiatomics);
            /*--- Calculate Landau-Teller relaxation time ---*/
            ReducedMass = Molar_Mass[iSpecies]*Molar_Mass[jSpecies]/(
                Molar_Mass[iSpecies]+Molar_Mass[jSpecies]);
            //			LimitingXSection = 1E-20*pow(50000/Temp_tr_i[iSpecies], 2.0);
            arg10d = -(2.0*Temp_tr_id[iSpecies]/Temp_tr_i[iSpecies]);
            arg10 = 2.0*log(50000/Temp_tr_i[iSpecies]);
            LimitingXSectiond = 1E-20*arg10d*exp(arg10);
            LimitingXSection = 1E-20*exp(arg10);
            arg1d = 8.0*UNIVERSAL_GAS_CONSTANT*Temp_tr_id[iSpecies]/(PI_NUMBER*Molar_Mass[
                iSpecies]);
            arg1 = 8.0*UNIVERSAL_GAS_CONSTANT*Temp_tr_i[iSpecies]/(PI_NUMBER*Molar_Mass[
                iSpecies]);
            AvgMolecularSpeedd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
            AvgMolecularSpeed = sqrt(arg1);
            //			A_sr = 1.16 * 1E-3 * sqrt(ReducedMass) * pow(CharVibTemp[iSpecies], 4.0/3.0);
            result1 = sqrt(ReducedMass);
            arg10 = 4.0/3.0*log(CharVibTemp[iSpecies]);
            A_sr = 1.16*1E-3*result1*exp(arg10);
            /*			tau_sr =   101325/(SpeciesPressure_i[iSpecies]+SpeciesPressure_i[jSpecies])

																			* exp(A_sr*(pow(sqrt(Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies]),-1.0/3.0) - 0.015*pow(ReducedMass,0.25)) - 18.42);
            */
            arg1d = Temp_tr_id[iSpecies]*Temp_tr_i[jSpecies] + Temp_tr_i[
                iSpecies]*Temp_tr_id[jSpecies];
            arg1 = Temp_tr_i[iSpecies]*Temp_tr_i[jSpecies];
            result1d = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
            result1 = sqrt(arg1);
            arg2d = -(result1d/(3.0*result1));
            arg2 = -1.0/3.0*log(result1);
            arg3 = 0.25*log(ReducedMass);
            arg4d = A_sr*arg2d*exp(arg2);
            arg4 = A_sr*(exp(arg2)-0.015*exp(arg3)) - 18.42;
            tau_srd = 101325*arg4d*exp(arg4)/(SpeciesPressure_i[iSpecies]+
                SpeciesPressure_i[jSpecies]) - 101325*(SpeciesPressure_id[
                iSpecies]+SpeciesPressure_id[jSpecies])*exp(arg4)/((
                SpeciesPressure_i[iSpecies]+SpeciesPressure_i[jSpecies])*(
                SpeciesPressure_i[iSpecies]+SpeciesPressure_i[jSpecies]));
            tau_sr = 101325/(SpeciesPressure_i[iSpecies]+SpeciesPressure_i[
                jSpecies])*exp(arg4);
            tau_psd = -(((LimitingXSectiond*AvgMolecularSpeed+LimitingXSection
                *AvgMolecularSpeedd)*MixtureNumDensity+LimitingXSection*
                AvgMolecularSpeed*MixtureNumDensityd)/(LimitingXSection*
                AvgMolecularSpeed*MixtureNumDensity*(LimitingXSection*
                AvgMolecularSpeed*MixtureNumDensity)));
            tau_ps = 1.0/(LimitingXSection*AvgMolecularSpeed*MixtureNumDensity
                );
            arg1d = -(CharVibTemp[iSpecies]*Temp_tr_id[jSpecies]/(Temp_tr_i[
                jSpecies]*Temp_tr_i[jSpecies]));
            arg1 = CharVibTemp[iSpecies]/Temp_tr_i[jSpecies];
            estar_vsd = -(UNIVERSAL_GAS_CONSTANT*CharVibTemp[iSpecies]*arg1d*exp(arg1)/
                Molar_Mass[iSpecies]/((exp(arg1)-1.0)*(exp(arg1)-1.0)));
            estar_vs = UNIVERSAL_GAS_CONSTANT/Molar_Mass[iSpecies]*CharVibTemp[iSpecies]/(exp
                (arg1)-1.0);
            e_vsd = (U_id[iLoc+nDim+2]*U_i[iLoc+0]-U_i[iLoc+nDim+2]*U_id
                [iLoc+0])/(U_i[iLoc+0]*U_i[iLoc+0]);
            e_vs = U_i[iLoc+nDim+2]/U_i[iLoc+0];
            /*--- Energy transferred per-molecule from species r to species s ---
            */
            q_tr_vsd = (Molecular_Mass[iSpecies]*(estar_vsd-e_vsd)*(tau_sr+
                tau_ps)-Molecular_Mass[iSpecies]*(estar_vs-e_vs)*(tau_srd+
                tau_psd))/((tau_sr+tau_ps)*(tau_sr+tau_ps));
            q_tr_vs = Molecular_Mass[iSpecies]*(estar_vs-e_vs)/(tau_sr+tau_ps)
            ;
            /*--- Convert to energy per volume for r and s species and multiply by volume for residual ---
            */
            //			Q_tv[iSpecies] += U_i[iLoc+0] * (estar_vs - e_vs)/(tau_sr + tau_ps);
            val_residuald[iLoc + nDim + 2] = val_residuald[iLoc + nDim +
                2] + Volume*(q_tr_vsd*U_i[iLoc+0]+q_tr_vs*U_id[iLoc+0])/
                Molecular_Mass[iSpecies];
            val_residual[iLoc + nDim + 2] += q_tr_vs*U_i[iLoc+0]/
            Molecular_Mass[iSpecies]*Volume;
            val_residuald[iLoc + nDim + 1] = val_residuald[iLoc + nDim +
                1] + Volume*(q_tr_vsd*U_i[iLoc+0]+q_tr_vs*U_id[iLoc+0])/
                Molecular_Mass[iSpecies];
            val_residual[iLoc + nDim + 1] += q_tr_vs*U_i[iLoc+0]/
            Molecular_Mass[iSpecies]*Volume;
            val_residuald[jLoc + nDim + 1] = val_residuald[jLoc + nDim +
                1] - Volume*(q_tr_vsd*U_i[jLoc+0]+q_tr_vs*U_id[jLoc+0])/
                Molecular_Mass[jSpecies];
            val_residual[jLoc + nDim + 1] -= q_tr_vs*U_i[jLoc+0]/
            Molecular_Mass[jSpecies]*Volume;
        }
    }


//SU2_DIFF END CSourcePieceWise_Plasma__SetResidual_EnergyExch
	
}
