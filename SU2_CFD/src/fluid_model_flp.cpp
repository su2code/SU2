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

#ifdef HAVE_FluidProp

#include "../include/fluid_model.hpp"
#include "fluidprop.h"

#define LEN_FLUIDNAMES 256
#define LEN_COMPONENTS 32

long i = 0;

CFluidProp::CFluidProp() : CFluidModel() {

}

CFluidProp::CFluidProp( string thermolib, int ncomp, string* comp, double* conc, bool HasSinglePhaseOnly, string tablename, double T_ref, double P_ref, double rho_ref, int ErrorLevel ) : CFluidModel() {

        // Copy fluid data to CFluidProp object
  	ThermoLib = thermolib;
	nComp = ncomp;
	Comp= comp;
	Conc = conc;
        SinglePhaseOnly = HasSinglePhaseOnly;
        TableName = tablename;

        this->ErrorLevel = ErrorLevel;

	// Prepare composition array's for fluidprop_setfluid 
	char LocalComp[20][LEN_COMPONENTS];
	double LocalConc[20];
	for( int i = 0; i < nComp; i++) {
	   strcpy( LocalComp[i], Comp[i].c_str());
	   LocalConc[i] = Conc[i];
	}

        // Intialize FluidProp (it opens the libraries) 
	if (!fluidprop_isinit()) init_fluidprop();

        // Define the working fluid 
	printf("SetFluid...\n");
	fluidprop_setfluid( ThermoLib.c_str(), nComp, LocalComp[0], LEN_COMPONENTS, LocalConc );

        // Table generator
        /*fluidprop_setunits( "SI", " ", " ", " ");
        string table_name      = "CO2_table";
        string table_type      = "consistent";  // or "non-consisent"
        string interpol_method = "bicubic2";
        string mesh_name       = "CO2_mesh";
        string search_alg      = "kdtree"; 
        int    phase           = 2;             // phase = 1 --> grid covers only single phase area, phase = 2 --> grid includes 2-phase VLE area
        fluidprop_settable( table_name.c_str(), table_type.c_str(), interpol_method.c_str(), mesh_name.c_str(), search_alg.c_str(), phase);
        printf( "FluidProp error message: %s\n", fluidprop_geterror());
	fluidprop_setfluid( ThermoLib.c_str(), nComp, LocalComp[0], LEN_COMPONENTS, LocalConc );*/
        //throw(-1);

	// In case of using LuT, read the table
        bool LuT = ThermoLib.substr(0,3) == "LuT";
        if (LuT) printf("UseTable...\n");
        if (LuT) fluidprop_usetable( TableName.c_str());
        //if (LuT) fluidprop_setfluid( ThermoLib.c_str(), nComp, LocalComp[0], LEN_COMPONENTS, LocalConc );
        if (LuT) printf("UseTable finished!\n");

        if (strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
           throw(-1);
        }
        else {
           
           // Set units to SI (as a result also setrefstate_nondim expects input in SI units)
           fluidprop_setunits( "SI", " ", " ", " ");
           
           // Set reference state for non-dimensionalization 
           fluidprop_setrefstate_nondim( T_ref, P_ref, 1./rho_ref);
           
           if (strcmp( fluidprop_geterror(), "No errors")) {
              printf( "FluidProp error message: %s\n", fluidprop_geterror());
              throw(-1);
           }
           else {
              int version[4];
              if (LuT)
                 fluidprop_getversion( "LuTEoS", version ); 
	      else
                 fluidprop_getversion( thermolib.c_str(), version ); 

              printf("-----------------------------------------------------\n");
              printf("FluidProp fluid specification\n");
              printf("-----------------------------------------------------\n");
              if (LuT) {
 	         printf("   Selected Library       : LuTEoS version %d.%d.%d.%d\n", version[0], version[1], version[2], version[3]);
 	         printf("   Selected Table         : %s\n", TableName.c_str());
	      }
              else
 	         printf("   Selected Library       : %s version %d.%d.%d.%d\n", thermolib.c_str(), version[0], version[1], version[2], version[3]);

	      printf("   Selected Components    : ");
	      for( int i = 0; i < nComp; i++)
	         printf("%s, ", Comp[i].c_str());
              printf("\n");
    
	      printf("   Selected Concentrations: ");
	      for( int i = 0; i < nComp; i++)
	         printf("%f, ", Conc[i]);
	      printf("\n"); 
  
              printf("   Error message          : %s\n", fluidprop_geterror());
              printf("-----------------------------------------------------\n\n");
	   }
        }
}


CFluidProp::~CFluidProp(void) {

}


void CFluidProp::SetTDState_rhoe (double rho, double e ){

	const char* pair;
        if (SinglePhaseOnly)
           pair = "du_1ph";
        else
           pair = "du";       

	Density = rho;
	StaticEnergy = e;

        struct fluidstate_su2_t output;
        if (ErrorLevel > 1) printf( "%ld.  rho = %f, u   = %f: ", i++, rho, e);
	fluidprop_allprops_su2( pair, rho, e, &output); 
        if (ErrorLevel > 1) printf( "P   = %f, T   = %f, s = %f, c = %f, cv = %f \n", 
                                    output.P, output.T, output.s, sqrt(output.c2), output.cv);

        Pressure    = output.P;
        Temperature = output.T;
	SoundSpeed2 = output.c2;
	Entropy     = output.s;
	dPdrho_e    = output.alpha;
	dPde_rho    = output.beta;
	dTdrho_e    = Temperature * dPde_rho / pow(rho,2);
	dTde_rho    = 1. / output.cv;

        if (ErrorLevel > 0 && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "rho = %f, u = %f\n",rho, e);
        }
}

void CFluidProp::SetTDState_PT (double P, double T ){

	const char* pair;
        if (SinglePhaseOnly)
           pair = "PT_1ph";
        else
           pair = "PT";       

	Pressure = P;
	Temperature = T;

        struct fluidstate_su2_t output;

        if (ErrorLevel > 1) printf( "%ld.  P =   %f, T =   %f: ", i++, P, T);
        fluidprop_allprops_su2( pair, P, T, &output); 
        if (ErrorLevel > 1) printf( "rho = %f, u   = %f, s = %f, c = %f, cv = %f \n", 
                                    output.d, output.u, output.s, sqrt(output.c2), output.cv);

        Density      = output.d;
	StaticEnergy = output.u;
	SoundSpeed2  = output.c2;
	Entropy      = output.s;
	dPdrho_e     = output.alpha;
	dPde_rho     = output.beta;
	dTdrho_e     = Temperature * dPde_rho / pow(output.d,2);
	dTde_rho     = 1. / output.cv;
	       
        if (ErrorLevel > 0  && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "P = %f, T = %f, u = %f\n", P, T, output.u);
        }
}

void CFluidProp::SetTDState_Prho (double P, double rho ){

	const char* pair;
        if (SinglePhaseOnly)
           pair = "Pd_1ph";
        else
           pair = "Pd";       

	Pressure = P;
	Density = rho;

        struct fluidstate_su2_t output;
        if (ErrorLevel > 1) printf( "%ld.  P   = %f, rho = %f: ", i++, P, rho);
        fluidprop_allprops_su2( pair, P, rho, &output); 
        if (ErrorLevel > 1) printf( "T   = %f, u   = %f, s = %f, c = %f, cv = %f \n", 
                                    output.T, output.u, output.s, sqrt(output.c2), output.cv);

	Temperature  = output.T;
	StaticEnergy = output.u;
	SoundSpeed2  = output.c2;
	Entropy      = output.s;
	dPdrho_e     = output.alpha;
	dPde_rho     = output.beta;
	dTdrho_e     = Temperature * dPde_rho / pow(rho,2);
	dTde_rho     = 1. / output.cv;
	
        if (ErrorLevel > 0 && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "P = %f, rho = %.16lf, u = %f, T = %f\n", P, rho, output.u, output.T);
        }

}

void CFluidProp::SetEnergy_Prho (double P, double rho ){

	const char* pair;
        if (SinglePhaseOnly)
           pair = "Pd_1ph";
        else
           pair = "Pd";       

        if (ErrorLevel > 1) printf( "%ld.  P =   %f, rho = %f: ", i++, P, rho);
	StaticEnergy = fluidprop_intenergy ( pair, P, rho );
        if (ErrorLevel > 1) printf( "StaticEnergy = %f\n", StaticEnergy);

        if (ErrorLevel > 0 && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "P = %f, rho = %f, u = %f\n", P, rho, StaticEnergy);
        }
}

void CFluidProp::SetTDState_hs (double h, double s ){

	const char* pair;
        if (SinglePhaseOnly)
           pair = "hs_1ph";
        else
           pair = "hs";       
	
	Entropy = s;

        struct fluidstate_su2_t output;
        if (ErrorLevel > 1) printf( "%ld.  h =   %f, s   = %f: ", i++, h, s);
        fluidprop_allprops_su2( pair, h, s, &output); 
        if (ErrorLevel > 1) printf( "P   = %f, rho = %f, T = %f, u = %f, c = %f, cv = %f \n", 
                                    output.P, output.d, output.T, output.u, sqrt(output.c2), output.cv);

        Pressure     = output.P;
	Temperature  = output.T;
        Density      = output.d;
	StaticEnergy = output.u;
	SoundSpeed2  = output.c2;
	dPdrho_e     = output.alpha;
	dPde_rho     = output.beta;
	dTdrho_e     = Temperature * dPde_rho / pow(output.d,2);
	dTde_rho     = 1. / output.cv;

        if (ErrorLevel > 0 && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "h = %f, s = %f, u = %f\n", h, s, output.u);
        }

}

void CFluidProp::SetTDState_rhoT (double rho, double T ){
    
	const char* pair;
        if (SinglePhaseOnly)
           pair = "Td_1ph";
        else
           pair = "Td";       

	Density = rho;
	Temperature = T;
    
        struct fluidstate_su2_t output;
        if (ErrorLevel > 1) printf( "%ld.  T   = %f, rho = %f: ", i++, T, rho);
        fluidprop_allprops_su2( pair, T, rho, &output); 
        if (ErrorLevel > 1) printf( "P   = %f, u   = %f, s = %f, c = %f, cv = %f \n", 
                                    output.P, output.u, output.s, sqrt(output.c2), output.cv);

        Pressure     = output.P;
	StaticEnergy = output.u;
	SoundSpeed2  = output.c2;
	Entropy      = output.s;
	dPdrho_e     = output.alpha;
	dPde_rho     = output.beta;
	dTdrho_e     = Temperature * dPde_rho / pow(rho,2);
	dTde_rho     = 1. / output.cv;
	
        if (ErrorLevel > 0 && strcmp( fluidprop_geterror(), "No errors")) {
           printf( "FluidProp error message: %s\n", fluidprop_geterror());
	   printf( "rho = %f, T = %f, u = %f\n", rho, T, output.u);
        }
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

#endif






