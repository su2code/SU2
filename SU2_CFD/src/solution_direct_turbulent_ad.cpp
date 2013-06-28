/*!
 * \file solution_direct_turbulent_ad.cpp
 * \brief Main subrotuines for solving direct problems (Euler, Navier-Stokes, etc.).
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

#include "../include/solution_structure.hpp"

void CAdjTurbSolution::CalcGradient_GG_ad(double *val_U_i, double *val_U_id, double **val_U_js, double **val_U_jsd,
		unsigned short nNeigh, unsigned short numVar, double **Normals, double **grad_U_i, double **grad_U_id, CConfig *config,
		CGeometry *geometry, unsigned long iPoint){

	double Volume;
	Volume = geometry->node[iPoint]->GetVolume();

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CTurbSolution__CalcGradient_GG

    unsigned short int iVar, iNeigh, jVar, kDim;
    /*--- Set Gradient to Zero ---*/
    for (iVar = 0; iVar < numVar; ++iVar)
        for (kDim = 0; kDim < nDim; ++kDim) {
            grad_U_id[iVar][kDim] = 0.0;
            grad_U_i[iVar][kDim] = 0.0;
        }
    **grad_U_id = 0.0;
    for (iNeigh = 0; iNeigh < nNeigh; ++iNeigh)
        for (jVar = 0; jVar < numVar; ++jVar)
            for (kDim = 0; kDim < nDim; ++kDim) {
                grad_U_id[jVar][kDim] = grad_U_id[jVar][kDim] + 0.5*Normals[
                    iNeigh][kDim]*(val_U_id[jVar]+val_U_jsd[iNeigh][jVar])/
                    Volume;
                grad_U_i[jVar][kDim] += 0.5*(val_U_i[jVar]+val_U_js[iNeigh][
                jVar])*Normals[iNeigh][kDim]/Volume;
            }


//SU2_DIFF END CTurbSolution__CalcGradient_GG
    
//    for (unsigned short jPos = 0; jPos < numVar; jPos++){
//        cout << "@i: " << val_U_i[jPos] << endl;;
//        for (unsigned short jNode = 0; jNode < nNeigh; jNode++){
//            cout << "@ks: " << val_U_js[jNode][jPos] << endl;
//        }
//    }
    

}

void CAdjTurbSolution::CalcGradient_LS_ad(double *val_U_i, double *val_U_id, double **val_U_js, double **val_U_jsd,
		unsigned short nNeigh, double *coords_i, double **coords_js, double **grad_U_i, double **grad_U_id,
		CConfig *config){

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CTurbSolution__CalcGradient_LS




//SU2_DIFF END CTurbSolution__CalcGradient_LS

}

void CAdjTurbSolution::CalcPrimVar_Compressible_ad(double *val_Vars, double *val_Varsd, double Gamma, double Gas_Constant, unsigned short numVar,
		double turb_ke, double* Primitive, double* Primitived, CConfig *config){

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CTurbSolution__CalcPrimVar_Compressible

    unsigned short int iDim;
    double Velocity2;
    double Velocity2d;
    //SetVelocity2();								 // Compute the modulus of the velocity.
    Velocity2 = 0.0;
    Velocity2d = 0.0;
    for (iDim = 0; iDim < nDim; ++iDim) {
        Velocity2d = Velocity2d + ((val_Varsd[iDim+1]*val_Vars[iDim+1]+
            val_Vars[iDim+1]*val_Varsd[iDim+1])*(val_Vars[0]*val_Vars[0])-
            val_Vars[iDim+1]*val_Vars[iDim+1]*(val_Varsd[0]*val_Vars[0]+
            val_Vars[0]*val_Varsd[0]))/(val_Vars[0]*val_Vars[0]*(val_Vars[0]*
            val_Vars[0]));
        Velocity2 += val_Vars[iDim+1]*val_Vars[iDim+1]/(val_Vars[0]*val_Vars[0
        ]);
    }
    //SetPressure(Gamma, turb_ke);   // Requires Velocity2 computation.
    Primitived[nDim + 1] = (Gamma-1.0)*(val_Varsd[0]*(val_Vars[numVar-1]/
        val_Vars[0]-0.5*Velocity2-turb_ke)+val_Vars[0]*((val_Varsd[numVar-1]*
        val_Vars[0]-val_Vars[numVar-1]*val_Varsd[0])/(val_Vars[0]*val_Vars[0]
        )-0.5*Velocity2d));
    Primitive[nDim + 1] = (Gamma-1.0)*val_Vars[0]*(val_Vars[numVar-1]/
        val_Vars[0]-0.5*Velocity2-turb_ke);
    //note: if turb_ke used, need to include its sensitivity...
    //SetTemperature(Gas_Constant);  // Requires pressure computation.
    Primitived[0] = (Primitived[nDim+1]*Gas_Constant*val_Vars[0]-Primitive[
        nDim+1]*Gas_Constant*val_Varsd[0])/(Gas_Constant*val_Vars[0]*(
        Gas_Constant*val_Vars[0]));
    Primitive[0] = Primitive[nDim+1]/(Gas_Constant*val_Vars[0]);
    //SetLaminarViscosity();				 // Requires temperature computation.
    //CalcLaminarViscosity( Solution, LaminarViscosity, config);
    for (iDim = 0; iDim < nDim; ++iDim) {
        Primitived[iDim + 1] = (val_Varsd[iDim+1]*val_Vars[0]-val_Vars[iDim+1]
            *val_Varsd[0])/(val_Vars[0]*val_Vars[0]);
        Primitive[iDim + 1] = val_Vars[iDim+1]/val_Vars[0];
    }


//SU2_DIFF END CTurbSolution__CalcPrimVar_Compressible

}

void CAdjTurbSolution::CalcLaminarViscosity_ad(double *val_U_i, double *val_U_id,
		double *val_laminar_viscosity_i, double *val_laminar_viscosity_id, CConfig *config){

	double Gas_Constant = config->GetGas_Constant();
	double Temperature_Ref = config->GetTemperature_Ref();
	double Viscosity_Ref = config->GetViscosity_Ref();

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CTurbSolution__CalcLaminarViscosity

    double Temperature, Temperature_Dim;
    double Temperatured, Temperature_Dimd;
    double Density, Pressure;
    double Densityd, Pressured;
    double result1;
    double result1d;
    Densityd = val_U_id[0];
    Density = val_U_i[0];
    // Not incompressible
    if (nDim == 2) {
        Pressured = Gamma_Minus_One*(val_U_id[3]-((val_U_id[1]*val_U_i[1]+
            val_U_i[1]*val_U_id[1]+val_U_id[2]*val_U_i[2]+val_U_i[2]*val_U_id[
            2])*2.0*val_U_i[0]-(val_U_i[1]*val_U_i[1]+val_U_i[2]*val_U_i[2])*
            2.0*val_U_id[0])/(2.0*val_U_i[0]*(2.0*val_U_i[0])));
        Pressure = Gamma_Minus_One*(val_U_i[3]-(val_U_i[1]*val_U_i[1]+val_U_i[
            2]*val_U_i[2])/(2.0*val_U_i[0]));
    } else {
        Pressured = Gamma_Minus_One*(val_U_id[3]-((val_U_id[1]*val_U_i[1]+
            val_U_i[1]*val_U_id[1]+val_U_id[2]*val_U_i[2]+val_U_i[2]*val_U_id[
            2]+val_U_id[3]*val_U_i[3]+val_U_i[3]*val_U_id[3])*2.0*val_U_i[0]-(
            val_U_i[1]*val_U_i[1]+val_U_i[2]*val_U_i[2]+val_U_i[3]*val_U_i[3])
            *2.0*val_U_id[0])/(2.0*val_U_i[0]*(2.0*val_U_i[0])));
        Pressure = Gamma_Minus_One*(val_U_i[3]-(val_U_i[1]*val_U_i[1]+val_U_i[
            2]*val_U_i[2]+val_U_i[3]*val_U_i[3])/(2.0*val_U_i[0]));
    }
    Temperatured = (Pressured*Gas_Constant*Density-Pressure*Gas_Constant*
        Densityd)/(Gas_Constant*Density*(Gas_Constant*Density));
    Temperature = Pressure/(Gas_Constant*Density);
    /*--- Calculate viscosity from a non-dim. Sutherland's Law ---*/
    Temperature_Dimd = Temperature_Ref*Temperatured;
    Temperature_Dim = Temperature*Temperature_Ref;
    result1d = pow_d(Temperature_Dim/300.0, Temperature_Dimd/300.0, 3.0/2.0, &
        result1);
    val_laminar_viscosity_id[0] = 1.853E-5*((300.0+110.3)*result1d*(
        Temperature_Dim+110.3)-result1*(300.0+110.3)*Temperature_Dimd)/((
        Temperature_Dim+110.3)*(Temperature_Dim+110.3));
    val_laminar_viscosity_i[0] = 1.853E-5*(result1*(300.0+110.3)/(
        Temperature_Dim+110.3));
    val_laminar_viscosity_id[0] = val_laminar_viscosity_id[0]/Viscosity_Ref;
    val_laminar_viscosity_i[0] = val_laminar_viscosity_i[0]/Viscosity_Ref;


//SU2_DIFF END CTurbSolution__CalcLaminarViscosity

}
// CORRECT PLACE?**********
void CAdjTurbSolution::CalcEddyViscosity_ad(double *val_FlowVars, double *val_FlowVarsd,
		double val_laminar_viscosity, double val_laminar_viscosityd,
		double *val_TurbVar, double *val_TurbVard, double *val_eddy_viscosity, double *val_eddy_viscosityd,
		CConfig *config){

//************************************************//
// Please do not delete //SU2_DIFF comment lines  //
//************************************************//

//SU2_DIFF START CTurbSASolution__CalcEddyViscosity

    double Ji, Ji_3, fv1;
    double Jid, Ji_3d, fv1d;
    double rho, mu, nu, nu_hat;
    double rhod, mud, nud, nu_hatd;
    double cv1_3 = 7.1*7.1*7.1;
    //	if (incompressible) {
    //		rho = ??;
    //		mu  = val_laminar_viscosity;
    //	}
    //	else {
    rhod = val_FlowVarsd[0];
    rho = val_FlowVars[0];
    mud = val_laminar_viscosityd;
    mu = val_laminar_viscosity;
    //	}
    nud = (mud*rho-mu*rhod)/(rho*rho);
    nu = mu/rho;
    nu_hatd = val_TurbVard[0];
    nu_hat = val_TurbVar[0];
    Jid = (nu_hatd*nu-nu_hat*nud)/(nu*nu);
    Ji = nu_hat/nu;
    Ji_3d = (Jid*Ji+Ji*Jid)*Ji + Ji*Ji*Jid;
    Ji_3 = Ji*Ji*Ji;
    fv1d = (Ji_3d*(Ji_3+cv1_3)-Ji_3*Ji_3d)/((Ji_3+cv1_3)*(Ji_3+cv1_3));
    fv1 = Ji_3/(Ji_3+cv1_3);
    val_eddy_viscosityd[0] = (rhod*fv1+rho*fv1d)*nu_hat + rho*fv1*nu_hatd;
    val_eddy_viscosity[0] = rho*fv1*nu_hat;


//SU2_DIFF END CTurbSASolution__CalcEddyViscosity

}
