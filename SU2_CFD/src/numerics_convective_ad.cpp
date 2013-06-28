/*!
 * \file numerics_convective_ad.cpp
 * \brief This file contains Automatically Differentiated versions
 * of appropriate convective terms. These routines are produced
 * semi-automatically using python, Tapenade and some minor requirement
 * to add in small bits of code/comments
 * \author Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
 * \version 2.0.4
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
#include "../include/math_ad.hpp"



void CUpwRoe_AdjDiscFlow::SetDirectResidual_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwRoe_Flow__SetResidual

    double arg1;
    double arg1d;
    double tmpresult;
    int ii2;
    int ii1;
    Area = 0;
    for (iDim = 0; iDim < nDim; ++iDim)
        Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);
    /*-- Unit Normal ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
        UnitaryNormald[iDim] = 0.0;
        UnitaryNormal[iDim] = Normal[iDim]/Area;
    }
    /*--- Conserved variables at point i,

        Need to recompute SoundSpeed / Pressure / Enthalpy in 

        case of 2nd order reconstruction ---*/
    Density_id = U_id[0];
    Density_i = U_i[0];
    sq_vel = 0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Velocity_id[ii1] = 0.0;
    sq_veld = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Velocity_id[iDim] = (U_id[iDim+1]*Density_i-U_i[iDim+1]*Density_id)/(
            Density_i*Density_i);
        Velocity_i[iDim] = U_i[iDim+1]/Density_i;
        sq_veld = sq_veld + Velocity_id[iDim]*Velocity_i[iDim] + Velocity_i[
            iDim]*Velocity_id[iDim];
        sq_vel += Velocity_i[iDim]*Velocity_i[iDim];
    }
    Energy_id = (U_id[nDim+1]*Density_i-U_i[nDim+1]*Density_id)/(Density_i*
        Density_i);
    Energy_i = U_i[nDim+1]/Density_i;
    arg1d = Gamma*Gamma_Minus_One*(Energy_id-0.5*sq_veld);
    arg1 = Gamma*Gamma_Minus_One*(Energy_i-0.5*sq_vel);
    SoundSpeed_id = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
    SoundSpeed_i = sqrt(arg1);
    Pressure_id = ((SoundSpeed_id*SoundSpeed_i+SoundSpeed_i*SoundSpeed_id)*
        Density_i+SoundSpeed_i*SoundSpeed_i*Density_id)/Gamma;
    Pressure_i = SoundSpeed_i*SoundSpeed_i*Density_i/Gamma;
    Enthalpy_id = ((U_id[nDim+1]+Pressure_id)*Density_i-(U_i[nDim+1]+
        Pressure_i)*Density_id)/(Density_i*Density_i);
    Enthalpy_i = (U_i[nDim+1]+Pressure_i)/Density_i;
    /*--- Conserved variables at point j,

        Need to recompute SoundSpeed / Pressure / Enthalpy in 

        case of 2nd order reconstruction ---*/
    Density_jd = U_jd[0];
    Density_j = U_j[0];
    sq_vel = 0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Velocity_jd[ii1] = 0.0;
    sq_veld = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Velocity_jd[iDim] = (U_jd[iDim+1]*Density_j-U_j[iDim+1]*Density_jd)/(
            Density_j*Density_j);
        Velocity_j[iDim] = U_j[iDim+1]/Density_j;
        sq_veld = sq_veld + Velocity_jd[iDim]*Velocity_j[iDim] + Velocity_j[
            iDim]*Velocity_jd[iDim];
        sq_vel += Velocity_j[iDim]*Velocity_j[iDim];
    }
    Energy_jd = (U_jd[nDim+1]*Density_j-U_j[nDim+1]*Density_jd)/(Density_j*
        Density_j);
    Energy_j = U_j[nDim+1]/Density_j;
    arg1d = Gamma*Gamma_Minus_One*(Energy_jd-0.5*sq_veld);
    arg1 = Gamma*Gamma_Minus_One*(Energy_j-0.5*sq_vel);
    SoundSpeed_jd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
    SoundSpeed_j = sqrt(arg1);
    Pressure_jd = ((SoundSpeed_jd*SoundSpeed_j+SoundSpeed_j*SoundSpeed_jd)*
        Density_j+SoundSpeed_j*SoundSpeed_j*Density_jd)/Gamma;
    Pressure_j = SoundSpeed_j*SoundSpeed_j*Density_j/Gamma;
    Enthalpy_jd = ((U_jd[nDim+1]+Pressure_jd)*Density_j-(U_j[nDim+1]+
        Pressure_j)*Density_jd)/(Density_j*Density_j);
    Enthalpy_j = (U_j[nDim+1]+Pressure_j)/Density_j;
    /*--- Roe-averaged variables at interface between i & j ---*/
    Rd = (Density_j/Density_i == 0.0 ? 0.0 : (Density_jd*Density_i-Density_j*
        Density_id)/(Density_i*Density_i*2.0*sqrt(Density_j/Density_i)));
    R = sqrt(Density_j/Density_i);
    RoeDensityd = Rd*Density_i + R*Density_id;
    RoeDensity = R*Density_i;
    sq_vel = 0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        RoeVelocityd[ii1] = 0.0;
    sq_veld = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        RoeVelocityd[iDim] = ((Rd*Velocity_j[iDim]+R*Velocity_jd[iDim]+
            Velocity_id[iDim])*(R+1)-(R*Velocity_j[iDim]+Velocity_i[iDim])*Rd)
            /((R+1)*(R+1));
        RoeVelocity[iDim] = (R*Velocity_j[iDim]+Velocity_i[iDim])/(R+1);
        sq_veld = sq_veld + RoeVelocityd[iDim]*RoeVelocity[iDim] + RoeVelocity
            [iDim]*RoeVelocityd[iDim];
        sq_vel += RoeVelocity[iDim]*RoeVelocity[iDim];
    }
    RoeEnthalpyd = ((Rd*Enthalpy_j+R*Enthalpy_jd+Enthalpy_id)*(R+1)-(R*
        Enthalpy_j+Enthalpy_i)*Rd)/((R+1)*(R+1));
    RoeEnthalpy = (R*Enthalpy_j+Enthalpy_i)/(R+1);
    arg1d = (Gamma-1)*(RoeEnthalpyd-0.5*sq_veld);
    arg1 = (Gamma-1)*(RoeEnthalpy-0.5*sq_vel);
    RoeSoundSpeedd = (arg1 == 0.0 ? 0.0 : arg1d/(2.0*sqrt(arg1)));
    RoeSoundSpeed = sqrt(arg1);
    /*--- Compute Proj_flux_tensor_i ---*/
    if (nDim == 2) {
        double rhou = Density_i*Velocity_i[0];
        double rhoud = Density_id*Velocity_i[0] + Density_i*Velocity_id[0];
        double rhov = Density_i*Velocity_i[1];
        double rhovd = Density_id*Velocity_i[1] + Density_i*Velocity_id[1];
        for (ii1 = 0; ii1 < nVar; ++ii1)
            Proj_flux_tensor_id[ii1] = 0.0;
        Proj_flux_tensor_id[0] = Normal[0]*rhoud;
        Proj_flux_tensor_i[0] = rhou*Normal[0];
        Proj_flux_tensor_id[1] = Normal[0]*(rhoud*Velocity_i[0]+rhou*
            Velocity_id[0]+Pressure_id);
        Proj_flux_tensor_i[1] = (rhou*Velocity_i[0]+Pressure_i)*Normal[0];
        Proj_flux_tensor_id[2] = Normal[0]*(rhoud*Velocity_i[1]+rhou*
            Velocity_id[1]);
        Proj_flux_tensor_i[2] = rhou*Velocity_i[1]*Normal[0];
        Proj_flux_tensor_id[3] = Normal[0]*(rhoud*Enthalpy_i+rhou*Enthalpy_id)
        ;
        Proj_flux_tensor_i[3] = rhou*Enthalpy_i*Normal[0];
        Proj_flux_tensor_id[0] = Proj_flux_tensor_id[0] + Normal[1]*rhovd;
        Proj_flux_tensor_i[0] += rhov*Normal[1];
        Proj_flux_tensor_id[1] = Proj_flux_tensor_id[1] + Normal[1]*(rhovd*
            Velocity_i[0]+rhov*Velocity_id[0]);
        Proj_flux_tensor_i[1] += rhov*Velocity_i[0]*Normal[1];
        Proj_flux_tensor_id[2] = Proj_flux_tensor_id[2] + Normal[1]*(rhovd*
            Velocity_i[1]+rhov*Velocity_id[1]+Pressure_id);
        Proj_flux_tensor_i[2] += (rhov*Velocity_i[1]+Pressure_i)*Normal[1];
        Proj_flux_tensor_id[3] = Proj_flux_tensor_id[3] + Normal[1]*(rhovd*
            Enthalpy_i+rhov*Enthalpy_id);
        Proj_flux_tensor_i[3] += rhov*Enthalpy_i*Normal[1];
    } else {
        double rhou = Density_i*Velocity_i[0];
        double rhoud = Density_id*Velocity_i[0] + Density_i*Velocity_id[0];
        double rhov = Density_i*Velocity_i[1];
        double rhovd = Density_id*Velocity_i[1] + Density_i*Velocity_id[1];
        double rhow = Density_i*Velocity_i[2];
        double rhowd = Density_id*Velocity_i[2] + Density_i*Velocity_id[2];
        for (ii1 = 0; ii1 < nVar; ++ii1)
            Proj_flux_tensor_id[ii1] = 0.0;
        Proj_flux_tensor_id[0] = Normal[0]*rhoud;
        Proj_flux_tensor_i[0] = rhou*Normal[0];
        Proj_flux_tensor_id[1] = Normal[0]*(rhoud*Velocity_i[0]+rhou*
            Velocity_id[0]+Pressure_id);
        Proj_flux_tensor_i[1] = (rhou*Velocity_i[0]+Pressure_i)*Normal[0];
        Proj_flux_tensor_id[2] = Normal[0]*(rhoud*Velocity_i[1]+rhou*
            Velocity_id[1]);
        Proj_flux_tensor_i[2] = rhou*Velocity_i[1]*Normal[0];
        Proj_flux_tensor_id[3] = Normal[0]*(rhoud*Velocity_i[2]+rhou*
            Velocity_id[2]);
        Proj_flux_tensor_i[3] = rhou*Velocity_i[2]*Normal[0];
        Proj_flux_tensor_id[4] = Normal[0]*(rhoud*Enthalpy_i+rhou*Enthalpy_id)
        ;
        Proj_flux_tensor_i[4] = rhou*Enthalpy_i*Normal[0];
        Proj_flux_tensor_id[0] = Proj_flux_tensor_id[0] + Normal[1]*rhovd;
        Proj_flux_tensor_i[0] += rhov*Normal[1];
        Proj_flux_tensor_id[1] = Proj_flux_tensor_id[1] + Normal[1]*(rhovd*
            Velocity_i[0]+rhov*Velocity_id[0]);
        Proj_flux_tensor_i[1] += rhov*Velocity_i[0]*Normal[1];
        Proj_flux_tensor_id[2] = Proj_flux_tensor_id[2] + Normal[1]*(rhovd*
            Velocity_i[1]+rhov*Velocity_id[1]+Pressure_id);
        Proj_flux_tensor_i[2] += (rhov*Velocity_i[1]+Pressure_i)*Normal[1];
        Proj_flux_tensor_id[3] = Proj_flux_tensor_id[3] + Normal[1]*(rhovd*
            Velocity_i[2]+rhov*Velocity_id[2]);
        Proj_flux_tensor_i[3] += rhov*Velocity_i[2]*Normal[1];
        Proj_flux_tensor_id[4] = Proj_flux_tensor_id[4] + Normal[1]*(rhovd*
            Enthalpy_i+rhov*Enthalpy_id);
        Proj_flux_tensor_i[4] += rhov*Enthalpy_i*Normal[1];
        Proj_flux_tensor_id[0] = Proj_flux_tensor_id[0] + Normal[2]*rhowd;
        Proj_flux_tensor_i[0] += rhow*Normal[2];
        Proj_flux_tensor_id[1] = Proj_flux_tensor_id[1] + Normal[2]*(rhowd*
            Velocity_i[0]+rhow*Velocity_id[0]);
        Proj_flux_tensor_i[1] += rhow*Velocity_i[0]*Normal[2];
        Proj_flux_tensor_id[2] = Proj_flux_tensor_id[2] + Normal[2]*(rhowd*
            Velocity_i[1]+rhow*Velocity_id[1]);
        Proj_flux_tensor_i[2] += rhow*Velocity_i[1]*Normal[2];
        Proj_flux_tensor_id[3] = Proj_flux_tensor_id[3] + Normal[2]*(rhowd*
            Velocity_i[2]+rhow*Velocity_id[2]+Pressure_id);
        Proj_flux_tensor_i[3] += (rhow*Velocity_i[2]+Pressure_i)*Normal[2];
        Proj_flux_tensor_id[4] = Proj_flux_tensor_id[4] + Normal[2]*(rhowd*
            Enthalpy_i+rhow*Enthalpy_id);
        Proj_flux_tensor_i[4] += rhow*Enthalpy_i*Normal[2];
    }
    /*--- Compute Proj_flux_tensor_j ---*/
    if (nDim == 2) {
        double rhou = Density_j*Velocity_j[0];
        double rhoud = Density_jd*Velocity_j[0] + Density_j*Velocity_jd[0];
        double rhov = Density_j*Velocity_j[1];
        double rhovd = Density_jd*Velocity_j[1] + Density_j*Velocity_jd[1];
        for (ii1 = 0; ii1 < nVar; ++ii1)
            Proj_flux_tensor_jd[ii1] = 0.0;
        Proj_flux_tensor_jd[0] = Normal[0]*rhoud;
        Proj_flux_tensor_j[0] = rhou*Normal[0];
        Proj_flux_tensor_jd[1] = Normal[0]*(rhoud*Velocity_j[0]+rhou*
            Velocity_jd[0]+Pressure_jd);
        Proj_flux_tensor_j[1] = (rhou*Velocity_j[0]+Pressure_j)*Normal[0];
        Proj_flux_tensor_jd[2] = Normal[0]*(rhoud*Velocity_j[1]+rhou*
            Velocity_jd[1]);
        Proj_flux_tensor_j[2] = rhou*Velocity_j[1]*Normal[0];
        Proj_flux_tensor_jd[3] = Normal[0]*(rhoud*Enthalpy_j+rhou*Enthalpy_jd)
        ;
        Proj_flux_tensor_j[3] = rhou*Enthalpy_j*Normal[0];
        Proj_flux_tensor_jd[0] = Proj_flux_tensor_jd[0] + Normal[1]*rhovd;
        Proj_flux_tensor_j[0] += rhov*Normal[1];
        Proj_flux_tensor_jd[1] = Proj_flux_tensor_jd[1] + Normal[1]*(rhovd*
            Velocity_j[0]+rhov*Velocity_jd[0]);
        Proj_flux_tensor_j[1] += rhov*Velocity_j[0]*Normal[1];
        Proj_flux_tensor_jd[2] = Proj_flux_tensor_jd[2] + Normal[1]*(rhovd*
            Velocity_j[1]+rhov*Velocity_jd[1]+Pressure_jd);
        Proj_flux_tensor_j[2] += (rhov*Velocity_j[1]+Pressure_j)*Normal[1];
        Proj_flux_tensor_jd[3] = Proj_flux_tensor_jd[3] + Normal[1]*(rhovd*
            Enthalpy_j+rhov*Enthalpy_jd);
        Proj_flux_tensor_j[3] += rhov*Enthalpy_j*Normal[1];
    } else {
        double rhou = Density_j*Velocity_j[0];
        double rhoud = Density_jd*Velocity_j[0] + Density_j*Velocity_jd[0];
        double rhov = Density_j*Velocity_j[1];
        double rhovd = Density_jd*Velocity_j[1] + Density_j*Velocity_jd[1];
        double rhow = Density_j*Velocity_j[2];
        double rhowd = Density_jd*Velocity_j[2] + Density_j*Velocity_jd[2];
        for (ii1 = 0; ii1 < nVar; ++ii1)
            Proj_flux_tensor_jd[ii1] = 0.0;
        Proj_flux_tensor_jd[0] = Normal[0]*rhoud;
        Proj_flux_tensor_j[0] = rhou*Normal[0];
        Proj_flux_tensor_jd[1] = Normal[0]*(rhoud*Velocity_j[0]+rhou*
            Velocity_jd[0]+Pressure_jd);
        Proj_flux_tensor_j[1] = (rhou*Velocity_j[0]+Pressure_j)*Normal[0];
        Proj_flux_tensor_jd[2] = Normal[0]*(rhoud*Velocity_j[1]+rhou*
            Velocity_jd[1]);
        Proj_flux_tensor_j[2] = rhou*Velocity_j[1]*Normal[0];
        Proj_flux_tensor_jd[3] = Normal[0]*(rhoud*Velocity_j[2]+rhou*
            Velocity_jd[2]);
        Proj_flux_tensor_j[3] = rhou*Velocity_j[2]*Normal[0];
        Proj_flux_tensor_jd[4] = Normal[0]*(rhoud*Enthalpy_j+rhou*Enthalpy_jd)
        ;
        Proj_flux_tensor_j[4] = rhou*Enthalpy_j*Normal[0];
        Proj_flux_tensor_jd[0] = Proj_flux_tensor_jd[0] + Normal[1]*rhovd;
        Proj_flux_tensor_j[0] += rhov*Normal[1];
        Proj_flux_tensor_jd[1] = Proj_flux_tensor_jd[1] + Normal[1]*(rhovd*
            Velocity_j[0]+rhov*Velocity_jd[0]);
        Proj_flux_tensor_j[1] += rhov*Velocity_j[0]*Normal[1];
        Proj_flux_tensor_jd[2] = Proj_flux_tensor_jd[2] + Normal[1]*(rhovd*
            Velocity_j[1]+rhov*Velocity_jd[1]+Pressure_jd);
        Proj_flux_tensor_j[2] += (rhov*Velocity_j[1]+Pressure_j)*Normal[1];
        Proj_flux_tensor_jd[3] = Proj_flux_tensor_jd[3] + Normal[1]*(rhovd*
            Velocity_j[2]+rhov*Velocity_jd[2]);
        Proj_flux_tensor_j[3] += rhov*Velocity_j[2]*Normal[1];
        Proj_flux_tensor_jd[4] = Proj_flux_tensor_jd[4] + Normal[1]*(rhovd*
            Enthalpy_j+rhov*Enthalpy_jd);
        Proj_flux_tensor_j[4] += rhov*Enthalpy_j*Normal[1];
        Proj_flux_tensor_jd[0] = Proj_flux_tensor_jd[0] + Normal[2]*rhowd;
        Proj_flux_tensor_j[0] += rhow*Normal[2];
        Proj_flux_tensor_jd[1] = Proj_flux_tensor_jd[1] + Normal[2]*(rhowd*
            Velocity_j[0]+rhow*Velocity_jd[0]);
        Proj_flux_tensor_j[1] += rhow*Velocity_j[0]*Normal[2];
        Proj_flux_tensor_jd[2] = Proj_flux_tensor_jd[2] + Normal[2]*(rhowd*
            Velocity_j[1]+rhow*Velocity_jd[1]);
        Proj_flux_tensor_j[2] += rhow*Velocity_j[1]*Normal[2];
        Proj_flux_tensor_jd[3] = Proj_flux_tensor_jd[3] + Normal[2]*(rhowd*
            Velocity_j[2]+rhow*Velocity_jd[2]+Pressure_jd);
        Proj_flux_tensor_j[3] += (rhow*Velocity_j[2]+Pressure_j)*Normal[2];
        Proj_flux_tensor_jd[4] = Proj_flux_tensor_jd[4] + Normal[2]*(rhowd*
            Enthalpy_j+rhow*Enthalpy_jd);
        Proj_flux_tensor_j[4] += rhow*Enthalpy_j*Normal[2];
    }
    /*--- Compute P and Lambda (do it with the Normal) ---*/
    double sqvel, rhooc, rhoxc, c2;
    double sqveld, rhoocd, rhoxcd;
    rhoocd = (RoeDensityd*RoeSoundSpeed-RoeDensity*RoeSoundSpeedd)/(
        RoeSoundSpeed*RoeSoundSpeed);
    rhooc = RoeDensity/RoeSoundSpeed;
    rhoxcd = RoeDensityd*RoeSoundSpeed + RoeDensity*RoeSoundSpeedd;
    rhoxc = RoeDensity*RoeSoundSpeed;
    c2 = RoeSoundSpeed*RoeSoundSpeed;
//    rhooc, rhoxc, c2;
    if (nDim == 2) {
        sqveld = RoeVelocityd[0]*RoeVelocity[0] + RoeVelocity[0]*RoeVelocityd[
            0] + RoeVelocityd[1]*RoeVelocity[1] + RoeVelocity[1]*RoeVelocityd[
            1];
        sqvel = RoeVelocity[0]*RoeVelocity[0] + RoeVelocity[1]*RoeVelocity[1];
        P_Tensord[0][0] = 0.0;
        P_Tensor[0][0] = 1.0;
        P_Tensord[0][1] = 0.0;
        P_Tensor[0][1] = 0.0;
        for (ii1 = 0; ii1 < nVar; ++ii1)
            for (ii2 = 0; ii2 < nVar; ++ii2)
                P_Tensord[ii1][ii2] = 0.0;
        P_Tensord[0][2] = 0.5*rhoocd;
        P_Tensor[0][2] = 0.5*rhooc;
        P_Tensord[0][3] = 0.5*rhoocd;
        P_Tensor[0][3] = 0.5*rhooc;
        P_Tensord[1][0] = RoeVelocityd[0];
        P_Tensor[1][0] = RoeVelocity[0];
        P_Tensord[1][1] = UnitaryNormal[1]*RoeDensityd;
        P_Tensor[1][1] = RoeDensity*UnitaryNormal[1];
        P_Tensord[1][2] = 0.5*(RoeVelocityd[0]*rhooc+RoeVelocity[0]*rhoocd+
            UnitaryNormal[0]*RoeDensityd);
        P_Tensor[1][2] = 0.5*(RoeVelocity[0]*rhooc+UnitaryNormal[0]*RoeDensity
            );
        P_Tensord[1][3] = 0.5*(RoeVelocityd[0]*rhooc+RoeVelocity[0]*rhoocd-
            UnitaryNormal[0]*RoeDensityd);
        P_Tensor[1][3] = 0.5*(RoeVelocity[0]*rhooc-UnitaryNormal[0]*RoeDensity
            );
        P_Tensord[2][0] = RoeVelocityd[1];
        P_Tensor[2][0] = RoeVelocity[1];
        P_Tensord[2][1] = -(UnitaryNormal[0]*RoeDensityd);
        P_Tensor[2][1] = -RoeDensity*UnitaryNormal[0];
        P_Tensord[2][2] = 0.5*(RoeVelocityd[1]*rhooc+RoeVelocity[1]*rhoocd+
            UnitaryNormal[1]*RoeDensityd);
        P_Tensor[2][2] = 0.5*(RoeVelocity[1]*rhooc+UnitaryNormal[1]*RoeDensity
            );
        P_Tensord[2][3] = 0.5*(RoeVelocityd[1]*rhooc+RoeVelocity[1]*rhoocd-
            UnitaryNormal[1]*RoeDensityd);
        P_Tensor[2][3] = 0.5*(RoeVelocity[1]*rhooc-UnitaryNormal[1]*RoeDensity
            );
        P_Tensord[3][0] = 0.5*sqveld;
        P_Tensor[3][0] = 0.5*sqvel;
        P_Tensord[3][1] = UnitaryNormal[1]*(RoeDensityd*RoeVelocity[0]+
            RoeDensity*RoeVelocityd[0]) - UnitaryNormal[0]*(RoeDensityd*
            RoeVelocity[1]+RoeDensity*RoeVelocityd[1]);
        P_Tensor[3][1] = RoeDensity*RoeVelocity[0]*UnitaryNormal[1] - 
            RoeDensity*RoeVelocity[1]*UnitaryNormal[0];
        P_Tensord[3][2] = 0.5*(0.5*(sqveld*rhooc+sqvel*rhoocd)+UnitaryNormal[0
            ]*(RoeDensityd*RoeVelocity[0]+RoeDensity*RoeVelocityd[0])+
            UnitaryNormal[1]*(RoeDensityd*RoeVelocity[1]+RoeDensity*
            RoeVelocityd[1])+rhoxcd/Gamma_Minus_One);
        P_Tensor[3][2] = 0.5*(0.5*sqvel*rhooc+RoeDensity*RoeVelocity[0]*
            UnitaryNormal[0]+RoeDensity*RoeVelocity[1]*UnitaryNormal[1]+rhoxc/
            Gamma_Minus_One);
        P_Tensord[3][3] = 0.5*(0.5*(sqveld*rhooc+sqvel*rhoocd)-UnitaryNormal[0
            ]*(RoeDensityd*RoeVelocity[0]+RoeDensity*RoeVelocityd[0])-
            UnitaryNormal[1]*(RoeDensityd*RoeVelocity[1]+RoeDensity*
            RoeVelocityd[1])+rhoxcd/Gamma_Minus_One);
        P_Tensor[3][3] = 0.5*(0.5*sqvel*rhooc-RoeDensity*RoeVelocity[0]*
            UnitaryNormal[0]-RoeDensity*RoeVelocity[1]*UnitaryNormal[1]+rhoxc/
            Gamma_Minus_One);
    } else {
        sqveld = RoeVelocityd[0]*RoeVelocity[0] + RoeVelocity[0]*RoeVelocityd[
            0] + RoeVelocityd[1]*RoeVelocity[1] + RoeVelocity[1]*RoeVelocityd[
            1] + RoeVelocityd[2]*RoeVelocity[2] + RoeVelocity[2]*RoeVelocityd[
            2];
        sqvel = RoeVelocity[0]*RoeVelocity[0] + RoeVelocity[1]*RoeVelocity[1] 
            + RoeVelocity[2]*RoeVelocity[2];
        P_Tensord[0][0] = 0.0;
        P_Tensor[0][0] = UnitaryNormal[0];
        P_Tensord[0][1] = 0.0;
        P_Tensor[0][1] = UnitaryNormal[1];
        P_Tensord[0][2] = 0.0;
        P_Tensor[0][2] = UnitaryNormal[2];
        for (ii1 = 0; ii1 < nVar; ++ii1)
            for (ii2 = 0; ii2 < nVar; ++ii2)
                P_Tensord[ii1][ii2] = 0.0;
        P_Tensord[0][3] = 0.5*rhoocd;
        P_Tensor[0][3] = 0.5*rhooc;
        P_Tensord[0][4] = 0.5*rhoocd;
        P_Tensor[0][4] = 0.5*rhooc;
        P_Tensord[1][0] = UnitaryNormal[0]*RoeVelocityd[0];
        P_Tensor[1][0] = RoeVelocity[0]*UnitaryNormal[0];
        P_Tensord[1][1] = UnitaryNormal[1]*RoeVelocityd[0] - UnitaryNormal[2]*
            RoeDensityd;
        P_Tensor[1][1] = RoeVelocity[0]*UnitaryNormal[1] - RoeDensity*
            UnitaryNormal[2];
        P_Tensord[1][2] = UnitaryNormal[2]*RoeVelocityd[0] + UnitaryNormal[1]*
            RoeDensityd;
        P_Tensor[1][2] = RoeVelocity[0]*UnitaryNormal[2] + RoeDensity*
            UnitaryNormal[1];
        P_Tensord[1][3] = 0.5*(RoeVelocityd[0]*rhooc+RoeVelocity[0]*rhoocd+
            UnitaryNormal[0]*RoeDensityd);
        P_Tensor[1][3] = 0.5*(RoeVelocity[0]*rhooc+RoeDensity*UnitaryNormal[0]
            );
        P_Tensord[1][4] = 0.5*(RoeVelocityd[0]*rhooc+RoeVelocity[0]*rhoocd-
            UnitaryNormal[0]*RoeDensityd);
        P_Tensor[1][4] = 0.5*(RoeVelocity[0]*rhooc-RoeDensity*UnitaryNormal[0]
            );
        P_Tensord[2][0] = UnitaryNormal[0]*RoeVelocityd[1] + UnitaryNormal[2]*
            RoeDensityd;
        P_Tensor[2][0] = RoeVelocity[1]*UnitaryNormal[0] + RoeDensity*
            UnitaryNormal[2];
        P_Tensord[2][1] = UnitaryNormal[1]*RoeVelocityd[1];
        P_Tensor[2][1] = RoeVelocity[1]*UnitaryNormal[1];
        P_Tensord[2][2] = UnitaryNormal[2]*RoeVelocityd[1] - UnitaryNormal[0]*
            RoeDensityd;
        P_Tensor[2][2] = RoeVelocity[1]*UnitaryNormal[2] - RoeDensity*
            UnitaryNormal[0];
        P_Tensord[2][3] = 0.5*(RoeVelocityd[1]*rhooc+RoeVelocity[1]*rhoocd+
            UnitaryNormal[1]*RoeDensityd);
        P_Tensor[2][3] = 0.5*(RoeVelocity[1]*rhooc+RoeDensity*UnitaryNormal[1]
            );
        P_Tensord[2][4] = 0.5*(RoeVelocityd[1]*rhooc+RoeVelocity[1]*rhoocd-
            UnitaryNormal[1]*RoeDensityd);
        P_Tensor[2][4] = 0.5*(RoeVelocity[1]*rhooc-RoeDensity*UnitaryNormal[1]
            );
        P_Tensord[3][0] = UnitaryNormal[0]*RoeVelocityd[2] - UnitaryNormal[1]*
            RoeDensityd;
        P_Tensor[3][0] = RoeVelocity[2]*UnitaryNormal[0] - RoeDensity*
            UnitaryNormal[1];
        P_Tensord[3][1] = UnitaryNormal[1]*RoeVelocityd[2] + UnitaryNormal[0]*
            RoeDensityd;
        P_Tensor[3][1] = RoeVelocity[2]*UnitaryNormal[1] + RoeDensity*
            UnitaryNormal[0];
        P_Tensord[3][2] = UnitaryNormal[2]*RoeVelocityd[2];
        P_Tensor[3][2] = RoeVelocity[2]*UnitaryNormal[2];
        P_Tensord[3][3] = 0.5*(RoeVelocityd[2]*rhooc+RoeVelocity[2]*rhoocd+
            UnitaryNormal[2]*RoeDensityd);
        P_Tensor[3][3] = 0.5*(RoeVelocity[2]*rhooc+RoeDensity*UnitaryNormal[2]
            );
        P_Tensord[3][4] = 0.5*(RoeVelocityd[2]*rhooc+RoeVelocity[2]*rhoocd-
            UnitaryNormal[2]*RoeDensityd);
        P_Tensor[3][4] = 0.5*(RoeVelocity[2]*rhooc-RoeDensity*UnitaryNormal[2]
            );
        P_Tensord[4][0] = 0.5*UnitaryNormal[0]*sqveld + UnitaryNormal[2]*(
            RoeDensityd*RoeVelocity[1]+RoeDensity*RoeVelocityd[1]) - 
            UnitaryNormal[1]*(RoeDensityd*RoeVelocity[2]+RoeDensity*
            RoeVelocityd[2]);
        P_Tensor[4][0] = 0.5*sqvel*UnitaryNormal[0] + RoeDensity*RoeVelocity[1
            ]*UnitaryNormal[2] - RoeDensity*RoeVelocity[2]*UnitaryNormal[1];
        P_Tensord[4][1] = 0.5*UnitaryNormal[1]*sqveld - UnitaryNormal[2]*(
            RoeDensityd*RoeVelocity[0]+RoeDensity*RoeVelocityd[0]) + 
            UnitaryNormal[0]*(RoeDensityd*RoeVelocity[2]+RoeDensity*
            RoeVelocityd[2]);
        P_Tensor[4][1] = 0.5*sqvel*UnitaryNormal[1] - RoeDensity*RoeVelocity[0
            ]*UnitaryNormal[2] + RoeDensity*RoeVelocity[2]*UnitaryNormal[0];
        P_Tensord[4][2] = 0.5*UnitaryNormal[2]*sqveld + UnitaryNormal[1]*(
            RoeDensityd*RoeVelocity[0]+RoeDensity*RoeVelocityd[0]) - 
            UnitaryNormal[0]*(RoeDensityd*RoeVelocity[1]+RoeDensity*
            RoeVelocityd[1]);
        P_Tensor[4][2] = 0.5*sqvel*UnitaryNormal[2] + RoeDensity*RoeVelocity[0
            ]*UnitaryNormal[1] - RoeDensity*RoeVelocity[1]*UnitaryNormal[0];
        P_Tensord[4][3] = 0.5*(0.5*(sqveld*rhooc+sqvel*rhoocd)+RoeDensityd*(
            RoeVelocity[0]*UnitaryNormal[0]+RoeVelocity[1]*UnitaryNormal[1]+
            RoeVelocity[2]*UnitaryNormal[2])+RoeDensity*(UnitaryNormal[0]*
            RoeVelocityd[0]+UnitaryNormal[1]*RoeVelocityd[1]+UnitaryNormal[2]*
            RoeVelocityd[2])+rhoxcd/Gamma_Minus_One);
        P_Tensor[4][3] = 0.5*(0.5*sqvel*rhooc+RoeDensity*(RoeVelocity[0]*
            UnitaryNormal[0]+RoeVelocity[1]*UnitaryNormal[1]+RoeVelocity[2]*
            UnitaryNormal[2])+rhoxc/Gamma_Minus_One);
        P_Tensord[4][4] = 0.5*(0.5*(sqveld*rhooc+sqvel*rhoocd)-RoeDensityd*(
            RoeVelocity[0]*UnitaryNormal[0]+RoeVelocity[1]*UnitaryNormal[1]+
            RoeVelocity[2]*UnitaryNormal[2])-RoeDensity*(UnitaryNormal[0]*
            RoeVelocityd[0]+UnitaryNormal[1]*RoeVelocityd[1]+UnitaryNormal[2]*
            RoeVelocityd[2])+rhoxcd/Gamma_Minus_One);
        P_Tensor[4][4] = 0.5*(0.5*sqvel*rhooc-RoeDensity*(RoeVelocity[0]*
            UnitaryNormal[0]+RoeVelocity[1]*UnitaryNormal[1]+RoeVelocity[2]*
            UnitaryNormal[2])+rhoxc/Gamma_Minus_One);
    }
    ProjVelocity = 0.0;
    ProjVelocity_i = 0.0;
    ProjVelocity_j = 0.0;
    ProjVelocityd = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        ProjVelocityd = ProjVelocityd + UnitaryNormal[iDim]*RoeVelocityd[iDim]
        ;
        ProjVelocity += RoeVelocity[iDim]*UnitaryNormal[iDim];
        ProjVelocity_i += Velocity_i[iDim]*UnitaryNormal[iDim];
        ProjVelocity_j += Velocity_j[iDim]*UnitaryNormal[iDim];
    }
    for (ii1 = 0; ii1 < nVar; ++ii1)
        Lambdad[ii1] = 0.0;
    /*--- Flow eigenvalues and Entropy correctors ---*/
    for (iDim = 0; iDim < nDim; ++iDim) {
        Lambdad[iDim] = ProjVelocityd;
        Lambda[iDim] = ProjVelocity;
    }
    Lambdad[nVar - 2] = ProjVelocityd + RoeSoundSpeedd;
    Lambda[nVar - 2] = ProjVelocity + RoeSoundSpeed;
    Lambdad[nVar - 1] = ProjVelocityd - RoeSoundSpeedd;
    Lambda[nVar - 1] = ProjVelocity - RoeSoundSpeed;
    /*--- Entropy correction ---
	for (iVar = 0; iVar < nVar; iVar++)

		if ( fabs(Lambda[iVar]) < Epsilon[iVar] )

			Lambda[iVar] = (Lambda[iVar]*Lambda[iVar] + Epsilon[iVar]*Epsilon[iVar])/(2.0*Epsilon[iVar]);

		else

			Lambda[iVar] = fabs(Lambda[iVar]);*/
    for (iVar = 0; iVar < nVar; ++iVar) {
        Lambdad[iVar] = fabs_d(Lambda[iVar], Lambdad[iVar], &tmpresult);
        Lambda[iVar] = tmpresult;
    }
    /*--- Compute wave amplitudes (characteristics) ---*/
    proj_delta_vel = 0.0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        delta_veld[ii1] = 0.0;
    proj_delta_veld = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        delta_veld[iDim] = Velocity_jd[iDim] - Velocity_id[iDim];
        delta_vel[iDim] = Velocity_j[iDim] - Velocity_i[iDim];
        proj_delta_veld = proj_delta_veld + Normal[iDim]*delta_veld[iDim];
        proj_delta_vel += delta_vel[iDim]*Normal[iDim];
    }
    delta_pd = Pressure_jd - Pressure_id;
    delta_p = Pressure_j - Pressure_i;
    delta_rhod = Density_jd - Density_id;
    delta_rho = Density_j - Density_i;
    proj_delta_veld = proj_delta_veld/Area;
    proj_delta_vel = proj_delta_vel/Area;
    if (nDim == 2) {
        for (ii1 = 0; ii1 < nVar; ++ii1)
            delta_waved[ii1] = 0.0;
        delta_waved[0] = delta_rhod - (delta_pd*(RoeSoundSpeed*RoeSoundSpeed)-
            delta_p*(RoeSoundSpeedd*RoeSoundSpeed+RoeSoundSpeed*RoeSoundSpeedd
            ))/(RoeSoundSpeed*RoeSoundSpeed*(RoeSoundSpeed*RoeSoundSpeed));
        delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
        delta_waved[1] = UnitaryNormal[1]*delta_veld[0] - UnitaryNormal[0]*
            delta_veld[1];
        delta_wave[1] = UnitaryNormal[1]*delta_vel[0] - UnitaryNormal[0]*
            delta_vel[1];
        delta_waved[2] = proj_delta_veld + (delta_pd*RoeDensity*RoeSoundSpeed-
            delta_p*(RoeDensityd*RoeSoundSpeed+RoeDensity*RoeSoundSpeedd))/(
            RoeDensity*RoeSoundSpeed*(RoeDensity*RoeSoundSpeed));
        delta_wave[2] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
        delta_waved[3] = (delta_pd*RoeDensity*RoeSoundSpeed-delta_p*(
            RoeDensityd*RoeSoundSpeed+RoeDensity*RoeSoundSpeedd))/(RoeDensity*
            RoeSoundSpeed*(RoeDensity*RoeSoundSpeed)) - proj_delta_veld;
        delta_wave[3] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
        //*val_residuald = 0.0;
    } else {
        for (ii1 = 0; ii1 < nVar; ++ii1)
            delta_waved[ii1] = 0.0;
        // nDim == 3
        delta_waved[0] = delta_rhod - (delta_pd*(RoeSoundSpeed*RoeSoundSpeed)-
            delta_p*(RoeSoundSpeedd*RoeSoundSpeed+RoeSoundSpeed*RoeSoundSpeedd
            ))/(RoeSoundSpeed*RoeSoundSpeed*(RoeSoundSpeed*RoeSoundSpeed));
        delta_wave[0] = delta_rho - delta_p/(RoeSoundSpeed*RoeSoundSpeed);
        delta_waved[1] = UnitaryNormal[0]*delta_veld[2] - UnitaryNormal[2]*
            delta_veld[0];
        delta_wave[1] = UnitaryNormal[0]*delta_vel[2] - UnitaryNormal[2]*
            delta_vel[0];
        delta_waved[2] = UnitaryNormal[1]*delta_veld[0] - UnitaryNormal[0]*
            delta_veld[1];
        delta_wave[2] = UnitaryNormal[1]*delta_vel[0] - UnitaryNormal[0]*
            delta_vel[1];
        delta_waved[3] = proj_delta_veld + (delta_pd*RoeDensity*RoeSoundSpeed-
            delta_p*(RoeDensityd*RoeSoundSpeed+RoeDensity*RoeSoundSpeedd))/(
            RoeDensity*RoeSoundSpeed*(RoeDensity*RoeSoundSpeed));
        delta_wave[3] = proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
        delta_waved[4] = (delta_pd*RoeDensity*RoeSoundSpeed-delta_p*(
            RoeDensityd*RoeSoundSpeed+RoeDensity*RoeSoundSpeedd))/(RoeDensity*
            RoeSoundSpeed*(RoeDensity*RoeSoundSpeed)) - proj_delta_veld;
        delta_wave[4] = -proj_delta_vel + delta_p/(RoeDensity*RoeSoundSpeed);
        //*val_residuald = 0.0;
    }
    /*--- Roe's Flux approximation ---*/
    for (iVar = 0; iVar < nVar; ++iVar) {
        val_residuald[iVar] = 0.5*(Proj_flux_tensor_id[iVar]+
            Proj_flux_tensor_jd[iVar]);
        val_residual[iVar] = 0.5*(Proj_flux_tensor_i[iVar]+Proj_flux_tensor_j[
            iVar]);
        for (jVar = 0; jVar < nVar; ++jVar) {
            val_residuald[iVar] = val_residuald[iVar] - 0.5*Area*((Lambdad[
                jVar]*delta_wave[jVar]+Lambda[jVar]*delta_waved[jVar])*
                P_Tensor[iVar][jVar]+Lambda[jVar]*delta_wave[jVar]*P_Tensord[
                iVar][jVar]);
            val_residual[iVar] -= 0.5*Lambda[jVar]*delta_wave[jVar]*P_Tensor[
            iVar][jVar]*Area;
        }
    }


//SU2_DIFF END CUpwRoe_Flow__SetResidual

}

void CUpwSca_AdjDiscTurbSA::SetDirectResidual_ad(double *val_residual, double *val_residuald, CConfig *config) {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwSca_TurbSA__SetResidual

    double result1;
    double result1d;
    int ii1;
    Density_id = U_id[0];
    Density_i = U_i[0];
    Density_jd = U_jd[0];
    Density_j = U_j[0];
    q_ij = 0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Velocity_id[ii1] = 0.0;
    for (ii1 = 0; ii1 < nDim; ++ii1)
        Velocity_jd[ii1] = 0.0;
    q_ijd = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Velocity_id[iDim] = (U_id[iDim+1]*Density_i-U_i[iDim+1]*Density_id)/(
            Density_i*Density_i);
        Velocity_i[iDim] = U_i[iDim+1]/Density_i;
        Velocity_jd[iDim] = (U_jd[iDim+1]*Density_j-U_j[iDim+1]*Density_jd)/(
            Density_j*Density_j);
        Velocity_j[iDim] = U_j[iDim+1]/Density_j;
        q_ijd = q_ijd + 0.5*Normal[iDim]*(Velocity_id[iDim]+Velocity_jd[iDim])
        ;
        q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
    result1d = fabs_d(q_ij, q_ijd, &result1);
    a0d = 0.5*(q_ijd+result1d);
    a0 = 0.5*(q_ij+result1);
    result1d = fabs_d(q_ij, q_ijd, &result1);
    a1d = 0.5*(q_ijd-result1d);
    a1 = 0.5*(q_ij-result1);
    val_residuald[0] = a0d*TurbVar_i[0] + a0*TurbVar_id[0] + a1d*TurbVar_j[0] 
        + a1*TurbVar_jd[0];
    val_residual[0] = a0*TurbVar_i[0] + a1*TurbVar_j[0];


//SU2_DIFF END CUpwSca_TurbSA__SetResidual

}

void CUpwRoeArtComp_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwRoeArtComp_Flow__SetResidual



//SU2_DIFF END CUpwRoeArtComp_Flow__SetResidual

}

void CUpwLin_AdjDiscLevelSet::SetDirectResidual_ad ()  {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwLin_LevelSet__SetResidual



//SU2_DIFF END CUpwLin_LevelSet__SetResidual

}

void CUpwSca_AdjDiscTurb::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwSca_Turb__SetResidual



//SU2_DIFF END CUpwSca_Turb__SetResidual

}

void CCentJST_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CCentJST_Flow__SetResidual



//SU2_DIFF END CCentJST_Flow__SetResidual

}

void CCentJSTArtComp_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CCentJSTArtComp_Flow__SetResidual



//SU2_DIFF END CCentJSTArtComp_Flow__SetResidual

}

void CCentLax_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CCentLax_Flow__SetResidual



//SU2_DIFF END CCentLax_Flow__SetResidual

}

void CCentLaxArtComp_AdjDiscFlow::SetDirectResidual_ad () {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CCentLaxArtComp_Flow__SetResidual



//SU2_DIFF END CCentLaxArtComp_Flow__SetResidual

}

void CUpwRoe_AdjDiscPlasmaDiatomic::SetDirectResidual_ad() {
//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CUpwRoe_PlasmaDiatomic__SetResidual



//SU2_DIFF END CUpwRoe_PlasmaDiatomic__SetResidual

}
