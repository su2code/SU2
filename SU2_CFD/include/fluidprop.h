/*============================================================================================! */
/*                                                                                            ! */
/*  FluidProp      Teus van der Stelt        Copyright � 2004-2008. All Rights Reserved.      ! */
/*  ========       Piero Colonna             Energy Technology Section, TU Delft.             ! */
/*                                                                                            ! */
/*--------------------------------------------------------------------------------------------! */
/*                                                                                            ! */
/*  Module name   : FluidProp                                                                 ! */
/*                                                                                            ! */
/*  Description   : This module loads and provides access to the FluidProp library.           ! */
/*                                                                                            ! */
/*  Creation date : 21-07-2008, Henk Seubers, VORtech                                         ! */
/*                                                                                            ! */
/*============================================================================================! */

#ifndef _NL_TUDELFT_FLUIDPROP_FLUIDPROP_H
#define _NL_TUDELFT_FLUIDPROP_FLUIDPROP_H

#include "size3264bits.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Define BOOL if it is not already defined */
#ifndef BOOL
   #ifdef __cplusplus
      typedef bool BOOL;
   #else
      typedef int BOOL;
   #endif
#endif

/* Define TRUE if it is not already defined */
#ifndef TRUE
   #define TRUE ((BOOL)1)
#endif

/* Define FALSE if it is not already defined */
#ifndef FALSE
   #define FALSE ((BOOL)0)
#endif

#define FLUIDPROP_ALLPROPS_NFIELDS 23
#define FLUIDPROP_ALLPROPSSAT_NFIELDS 33

struct fluidstate_t {
        double P;
        double T; 
        double v; 
        double d; 
        double h; 
        double s; 
        double u; 
        double q; 
        double cv; 
        double cp; 
        double c; 
        double alpha; 
        double beta; 
        double chi; 
        double fi; 
        double ksi; 
        double psi; 
        double zeta; 
        double theta;     
        double kappa; 
        double gamma; 
        double eta; 
        double lambda; 
        double sigma; 
        double d_liq;   
        double d_vap; 
        double h_liq; 
        double h_vap; 
        double T_sat;         
        double dd_liq_dP; 
        double dd_vap_dP; 
        double dh_liq_dP;    
        double dh_vap_dP; 
        double dT_sat_dP;
};

/* ==================================================! */
/*  Procedures for object initialization and cleanup ! */
/* ==================================================! */

/* Is library loaded?  */
BOOL fluidprop_isinit();

/* Load library */
BOOL init_fluidprop();

/* Unload library */
void uninit_fluidprop();

/* Create object */
void fluidprop_createobject( const char* modelname);

/* Returns version array (4 longs) */
BOOL fluidprop_getversion( const char* modelname, int* version);

/* Return error string */
char* fluidprop_geterror();

/* Return fluidprop install path */
const char* fluidprop_getinstallpath();

/* Return fluidprop data path*/
const char* fluidprop_getdatapath();

/* functions calculating thermodynamic properties */
void fluidprop_setfluid( const char* ModelName, int nComp, const char* Comp, int len_Comp, const double* Conc);

/* returns current fluid */
void fluidprop_getfluid( char* ModelName, int len_ModelName, int* nComp, char* Comp, int len_Comp, double* Conc);

/* function setting table parameters for look-up table approach */
void fluidprop_settable( const char* GenerateMode, const char* InterpolMethod, const char* StepFunc, 
			 const char* SearchAlg, int nRow, int nCol, int nSat, double Tlow, double Thigh, 
                         double StepMin, double StepMax);

/* returns all fluid names */
void fluidprop_getfluidnames( const char* LongShort, const char* ModelName, int* nFluids, char* FluidNames, int len_FluidNames);

/* returns all properties in OutputState */
void fluidprop_allpropssat( const char* InputSpec, double Input1, double Input2, struct fluidstate_t* OutputState, double* x, double* y);

/* returns all properties in OutputState except derivatives */
void fluidprop_allprops( const char* InputSpec, double Input1, double Input2, struct fluidstate_t* OutputState, double* x, double* y);

/* returns all properties needed by joe in OutputState */
void fluidprop_allpropsjoe( const char* InputSpec, double Input1, double Input2, struct fluidstatejoe_t* OutputState);

double fluidprop_solve( const char* FuncSpec, double FuncVal, const char* InputSpec, int TargetProp, double FixedVal, double MinVal, double MaxVal);

double fluidprop_pressure( const char* InputSpec, double Input1, double Input2);

double fluidprop_temperature( const char* InputSpec, double Input1, double Input2);

double fluidprop_specvolume( const char* InputSpec, double Input1, double Input2);

double fluidprop_density( const char* InputSpec, double Input1, double Input2);

double fluidprop_enthalpy( const char* InputSpec, double Input1, double Input2);

double fluidprop_entropy( const char* InputSpec, double Input1, double Input2);

double fluidprop_intenergy( const char* InputSpec, double Input1, double Input2);

double fluidprop_vaporqual( const char* InputSpec, double Input1, double Input2);

/* returns the liquid phase composition in the array pointed to by Output (size=#components in mixture) */
void fluidprop_liquidcmp( const char* InputSpec, double Input1, double Input2, double* Output);

/* returns the liquid phase composition in the array pointed to by Output (size=#components in mixture) */
void fluidprop_vaporcmp( const char* InputSpec, double Input1, double Input2, double* Output);

double fluidprop_heatcapv( const char* InputSpec, double Input1, double Input2);

double fluidprop_heatcapp( const char* InputSpec, double Input1, double Input2);

double fluidprop_soundspeed( const char* InputSpec, double Input1, double Input2);

double fluidprop_alpha( const char* InputSpec, double Input1, double Input2);

double fluidprop_beta( const char* InputSpec, double Input1, double Input2);

double fluidprop_chi( const char* InputSpec, double Input1, double Input2);

double fluidprop_fi( const char* InputSpec, double Input1, double Input2);

double fluidprop_ksi( const char* InputSpec, double Input1, double Input2);

double fluidprop_psi( const char* InputSpec, double Input1, double Input2);

double fluidprop_zeta( const char* InputSpec, double Input1, double Input2);

double fluidprop_theta( const char* InputSpec, double Input1, double Input2);

double fluidprop_kappa( const char* InputSpec, double Input1, double Input2);

double fluidprop_gamma( const char* InputSpec, double Input1, double Input2);

double fluidprop_viscosity( const char* InputSpec, double Input1, double Input2);

double fluidprop_thermcond( const char* InputSpec, double Input1, double Input2);

double fluidprop_surftens( const char* InputSpec, double Input1, double Input2);

double fluidprop_mmol();

double fluidprop_tcrit();

double fluidprop_pcrit();

double fluidprop_tmin();

double fluidprop_tmax();

void fluidprop_allinfo( double* M_mol, double* T_crit, double* P_crit, double* T_min, double* T_max);

void fluidprop_setunits( const char* UnitSet, const char* MassOrMole, const char* Properties, const char* Units);

void fluidprop_setrefstate( double T_ref, double P_ref);

/* returns P, T, v, h, s, u, alpha, beta, chi, fi, zeta, gamma in OutputState */
void fluidprop_zflow_vu(struct fluidstate_t* OutputState);

void fluidprop_saturationline( int* nPnts, double* T, double* P, double* x, double* y);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _NL_TUDELFT_FLUIDPROP_FLUIDPROP_H */
