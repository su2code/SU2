/*!
 * fluid_model_flp.cpp
 * \brief Source of the ideal gas model.
 * \author: S.Vitale, G.Gori, M.Pini, A.Guardone, P.Colonna, T.P. van der Stelt
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

#include "../include/fluid_model.hpp"
#include "fluidprop.h"

#define LEN_FLUIDNAMES 256
#define LEN_COMPONENTS 32


CFluidProp::CFluidProp() : CFluidModel() {

}


CFluidProp::CFluidProp( string thermolib, int ncomp, string* comp, double* conc, double T_ref, double P_ref, double rho_ref ) : CFluidModel() {

	string  libraryName[10];

	libraryName[0] = "freeStanMix";
	libraryName[1] = "IF97";
	libraryName[2] = "GasMix";
	libraryName[3] = "PCP-SAFT";
	libraryName[4] = "RefProp";
	libraryName[5] = "ScalingLaws";
	libraryName[6] = "StanMix";
	libraryName[7] = "TPSI";
	libraryName[8] = "vThermo";
	libraryName[9] = "LuTEoS";

	cout << " " << "\n";

	init_fluidprop();
	
        fluidprop_setunits( "SI", " ", " ", " ");
        //fluidprop_setrefstate( T_ref, P_ref); //, rho_ref);

	cout << "------- CHECK FLUIDPROP LIBRARIES -------" << "\n";

	int i;
	for ( i=0; i<=9; i++ )
	  {
	  int  version[4];
	  if ( fluidprop_getversion( libraryName[i].c_str(), version ) )
	  {
	    printf( "%s version %d.%d.%d.%d available\n", libraryName[i].c_str(), version[0], version[1], version[2], version[3] );
	  }
	  else
	  {
	    printf( "%s library not detected\n", libraryName[i].c_str());
	    throw(-1);
	  }

	}

	ThermoLib = thermolib;
	nComp = ncomp;
	Comp= comp;
	Conc = conc;

	printf("Selected Library       : %s\n", ThermoLib.c_str() );

	printf("Selected Components    : ");
	for( int i = 0; i < nComp; i++)
	   printf("%s, ", Comp[i].c_str());
	printf("\n");
    
	printf("Selected Concentrations: ");
	for( int i = 0; i < nComp; i++)
	   printf("%f, ", Conc[i]);
	printf("\n");
    
	char LocalComp[20][LEN_COMPONENTS];
	double LocalConc[20];
	for( int i = 0; i < nComp; i++)
	{
	   strcpy( LocalComp[i], Comp[i].c_str());
	   LocalConc[i] = Conc[i];
	}
	fluidprop_setfluid( ThermoLib.c_str(), nComp, LocalComp[0], LEN_COMPONENTS, LocalConc );

	cout << "-----------------------------------------" << "\n";
	cout << " " << "\n";

	// Store the remaining referenced properties
	h_ref = e_ref;
	s_ref = e_ref/T_ref;
	dPdrho_e_ref = e_ref;
	dPde_rho_ref = rho_ref;
	dTdrho_e_ref = rho_ref/T_ref;
	dTde_rho_ref = e_ref/T_ref;
 }


CFluidProp::~CFluidProp(void) {

}


const bool single_phase = true;


void CFluidProp::SetTDState_rhoe (double rho, double e ){

	const char* pair;
        if (single_phase)
           pair = "vu_1ph";
        else
           pair = "vu";       

	//rho = rho*rho_ref;
	//e = e*e_ref;
	Density = rho;
	StaticEnergy = e;
	double v = 1.0/rho;

        struct fluidstate_t output;
        double x[20], y[20];
	fluidprop_allprops( pair, v, e, &output, x, y); 
        Pressure = output.P;
        Temperature = output.T;
	SoundSpeed2 = pow( output.c, 2);
	Entropy = output.s;
	dPdrho_e = output.alpha;
	dPde_rho = output.beta;
	dTdrho_e = Temperature * dPde_rho * pow(v,2);
	dTde_rho = 1. / output.cv;

        /*Pressure = fluidprop_pressure ( pair, v, e );
	Temperature = fluidprop_temperature ( pair, v, e );
	SoundSpeed2 = pow( fluidprop_soundspeed ( pair, v, e ), 2);
	Entropy = fluidprop_entropy ( pair, v, e );
	dPdrho_e = fluidprop_alpha ( pair, v, e );
	dPde_rho = fluidprop_beta ( pair, v, e );
	dTdrho_e = Temperature * dPde_rho * pow(v,2);
	dTde_rho = 1. / fluidprop_heatcapv ( pair, v, e );*/

}

void CFluidProp::SetTDState_PT (double P, double T ){

	const char* pair;
        if (single_phase)
           pair = "PT_1ph";
        else
           pair = "PT";       

	//const char* pair = "PT";
	//P = P*P_ref;
	//T = T*T_ref;
	Pressure = P;
	Temperature = T;

        struct fluidstate_t output;
        double x[20], y[20];
        fluidprop_allprops( pair, P, T, &output, x, y); 
        Density = output.d;
	StaticEnergy = output.u;
	SoundSpeed2 = pow( output.c, 2);
	Entropy = output.s;
	dPdrho_e = output.alpha;
	dPde_rho = output.beta;
	dTdrho_e = Temperature * dPde_rho * pow(output.v,2);
	dTde_rho = 1. / output.cv;
	
        /*Density = fluidprop_density ( pair, P, T );
	StaticEnergy = fluidprop_intenergy ( pair, P, T );
	SoundSpeed2 = pow( fluidprop_soundspeed ( pair, P, T ),2);
	Entropy = fluidprop_entropy ( pair, P, T );
	dPdrho_e = fluidprop_alpha ( pair, P, T );
	dPde_rho = fluidprop_beta ( pair, P, T );
	dTdrho_e = T * dPde_rho / pow(Density,2);
	dTde_rho = 1. / fluidprop_heatcapv ( pair, P, T );*/

}

void CFluidProp::SetTDState_Prho (double P, double rho ){

	const char* pair;
        if (single_phase)
           pair = "Pd_1ph";
        else
           pair = "Pd";       
	//const char* pair = "Pd";
	//P = P*P_ref;
	//rho = rho*rho_ref;
	Pressure = P;
	Density = rho;

        struct fluidstate_t output;
        double x[20], y[20];
        fluidprop_allprops( pair, P, rho, &output, x, y); 
	Temperature = output.T;
	StaticEnergy = output.u;
	SoundSpeed2 = pow( output.c, 2);
	Entropy = output.s;
	dPdrho_e = output.alpha;
	dPde_rho = output.beta;
	dTdrho_e = Temperature * dPde_rho * pow(output.v,2);
	dTde_rho = 1. / output.cv;
	
	/*Temperature = fluidprop_temperature ( pair, P, rho );
	StaticEnergy = fluidprop_intenergy ( pair, P, rho );
	SoundSpeed2 = pow( fluidprop_soundspeed ( pair, P, rho ),2);
	Entropy = fluidprop_entropy ( pair, P, rho );
	dPdrho_e = fluidprop_alpha ( pair, P, rho );
	dPde_rho = fluidprop_beta ( pair, P, rho );
	dTdrho_e = Temperature * dPde_rho / pow(rho,2);
	dTde_rho = 1. / fluidprop_heatcapv ( pair, P, rho );*/

}

void CFluidProp::SetEnergy_Prho (double P, double rho ){

	const char* pair;
        if (single_phase)
           pair = "Pd_1ph";
        else
           pair = "Pd";       

	//const char* pair = "Pd";
	//P = P*P_ref;
	//rho = rho*rho_ref;
	
	StaticEnergy = fluidprop_intenergy ( pair, P, rho );

}

void CFluidProp::SetTDState_hs (double h, double s ){

	const char* pair;
        if (single_phase)
           pair = "hs_1ph";
        else
           pair = "hs";       
	
	//const char* pair = "hs";
	//h = h*h_ref;
	//s = s*s_ref;
	Entropy = s;

        struct fluidstate_t output;
        double x[20], y[20];
        fluidprop_allprops( pair, h, s, &output, x, y); 
        Pressure = output.P;
	Temperature = output.T;
        Density = output.d;
	StaticEnergy = output.u;
	SoundSpeed2 = pow( output.c, 2);
	dPdrho_e = output.alpha;
	dPde_rho = output.beta;
	dTdrho_e = Temperature * dPde_rho * pow(output.v,2);
	dTde_rho = 1. / output.cv;
	
	/*Pressure = fluidprop_pressure ( pair, h, s );
	Temperature = fluidprop_temperature ( pair, h, s );
	Density = fluidprop_density ( pair, h, s );
	StaticEnergy = h - Pressure/Density;
	SoundSpeed2 = pow( fluidprop_soundspeed ( pair, h, s ),2);
	dPdrho_e = fluidprop_alpha ( pair, h, s );
	dPde_rho = fluidprop_beta ( pair, h, s );
	dTdrho_e = Temperature * dPde_rho / pow(Density,2);
	dTde_rho = 1. / fluidprop_heatcapv ( pair, h, s );*/

}

void CFluidProp::SetTDState_rhoT (double rho, double T ){
    
	const char* pair;
        if (single_phase)
           pair = "Td_1ph";
        else
           pair = "Td";       

	//const char* pair = "Td";
	//T = T*T_ref;
	//rho = rho*rho_ref;
	Density = rho;
	Temperature = T;
    
        struct fluidstate_t output;
        double x[20], y[20];
        fluidprop_allprops( pair, T, rho, &output, x, y); 
        Pressure = output.P;
	Temperature = output.T;
	SoundSpeed2 = pow( output.c, 2);
	Entropy = output.s;
	dPdrho_e = output.alpha;
	dPde_rho = output.beta;
	dTdrho_e = Temperature * dPde_rho * pow(output.v,2);
	dTde_rho = 1. / output.cv;
	
	/*Pressure = fluidprop_pressure ( pair, T, rho );
	Temperature = fluidprop_temperature ( pair, T, rho );
	SoundSpeed2 = pow( fluidprop_soundspeed ( pair, T, rho), 2);
	Entropy = fluidprop_entropy ( pair, T, rho );
	dPdrho_e = fluidprop_alpha ( pair, T, rho );
	dPde_rho = fluidprop_beta ( pair, T, rho );
	dTdrho_e = T * dPde_rho / pow(rho,2);
	dTde_rho = 1. / fluidprop_heatcapv ( pair, T, rho );*/
    
}

void CFluidProp::SetTDState_NonDim () {

	Pressure = Pressure/P_ref;
	Temperature = Temperature/T_ref;
	Density = Density/rho_ref;
	StaticEnergy = StaticEnergy/e_ref;
	SoundSpeed2 = SoundSpeed2/pow(v_ref,2);
	dPdrho_e = dPdrho_e/dPdrho_e_ref;
	dPde_rho = dPde_rho/dPde_rho_ref;
	dTdrho_e = dTdrho_e/dTdrho_e_ref;
	dTde_rho = dTde_rho/dTde_rho_ref;

}








