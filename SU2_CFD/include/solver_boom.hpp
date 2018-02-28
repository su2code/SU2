#pragma once

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

#include "../../Common/include/mpi_structure.hpp"

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"
#include "../../SU2_CFD/include/SUBoom.hpp"  // For things like mergesort, AtmosISA...

using namespace std;

class CBoom_AugBurgers{
public:
  unsigned long n_prop;
  unsigned long **pointID_original;
  su2double ***Coord_original;

  /*---Flight variables---*/
  su2double flt_h;
  su2double flt_M;
  su2double flt_psi, flt_gamma, flt_mu;
  su2double flt_heading[3];

  /*---Atmosphere variables---*/
  su2double atm_g, atm_R, atm_cp;
  su2double T_inf, a_inf, p_inf, rho_inf;

  /*---Ray variables---*/
  unsigned short ray_N_phi;	// Number of rays
  su2double ray_r0;			// Initial ray distance
  su2double ray_t0;			// Initial ray time
  su2double *ray_phi;		// Ray angles
  su2double *ray_x, *ray_y, ray_z, ray_A;
  su2double ray_lambda, *ray_gamma, *ray_theta;
  bool ground_flag;     // Whether or not we've propagated to the ground

  /*---ABE variables---*/
  su2double p0, 			// Reference pressure
  			w0,				// Reference angular frequency
  			rho0,			// Ambient density (function of altitude)
  			c0,				// Ambient sound speed (function of altitude)
  			dsigma,			// Step size
        dz,       // Change in altitude
        dtau,     // Grid spacing
  			xbar,			// Shock formation distance of plane wave
  			beta,			// Coefficient of nonlinearity
  			C_nu_O2,		// Dispersion parameter for O2 (dimensionless)
  			C_nu_N2,		// Dispersion parameter for N2 (dimensionless)
  			theta_nu_O2,	// Relaxation time parameter for O2 (dimensionless)
  			theta_nu_N2,	// Relaxation time parameter for N2 (dimensionless)
  			m_nu_O2,		// Dispersion parameter for O2 (dimensionless)
  			m_nu_N2,		// Dispersion parameter for N2 (dimensionless)
  			tau_nu_O2,		// Relaxation time for O2 (dimensionless)
  			tau_nu_N2,		// Relaxation time for N2 (dimensionless)
  			alpha0_tv,		// Thermoviscous attenuation coefficient
  			delta,			// Diffusivity of sound
  			Gamma,			// Thermoviscous parameter (dimensionless)
  			mu,				// Viscosity
  			kappa;			// Thermal conduction coefficient

  /*---Sensitivity---*/
  unsigned short nDim;
  unsigned long *nPanel, nSig, *nPointID, nPanel_tot, nPointID_loc, *nPointID_proc, nPointID_tot;
  unsigned long **PointID;
  su2double ***dJdU;

  /*---Ray data class for passing into RK4---*/
  class RayData{
    public:
      su2double L;
      su2double T;
      su2double c0;
      su2double nu;
      unsigned long M;
  };

  /*---Signal class for storing pressure signal---*/
  class Signal{
    public:
      unsigned long *len;   // signal length
      su2double **x,        // x-coordinates for every azimuth
                *t,         // time for signal points at a single azimuth
                *t_prime,   // retarded time
                *tau,       // retarded time (dimensionless)
                *taud,      // distorted time
                **p_prime,  // pressure signal for every azimuth
                *P,         // pressure signal (dimensionless)
                *dP_att,    // change in pressure (Attenuation)
                *dP_rel,    // change in pressure (Relaxation)
                *dP_spr,    // change in pressure (Spreading)
                *dP_str;    // change in pressure (Stratification)
  };

  Signal signal;

  /*---Main functions for boom calculations---*/
  CBoom_AugBurgers();
  CBoom_AugBurgers(CSolver *solver, CConfig *config, CGeometry *geometry);
  ~CBoom_AugBurgers(void);

  void SearchLinear(CConfig *config, CGeometry *geometry, 
                    const su2double r0, const su2double *phi);
  void ExtractLine(CGeometry *geometry, const su2double r0, unsigned short iPhi);
  void ExtractPressure(CSolver *solver, CConfig *config, CGeometry *geometry, unsigned short iPhi);
  bool InsideElem(CGeometry *geometry, su2double r0, su2double phi, unsigned long jElem, su2double *p0, su2double *p1);
  int Intersect2D(su2double r0, su2double *Coord_i, su2double *Coord_ip1, su2double *p0, su2double *p1);
  int Intersect3D(su2double r0, su2double phi, int nCoord, su2double **Coord_i, su2double *p1);

  void AtmosISA(su2double& h0, su2double& T, su2double& a, su2double& p,
                su2double& rho, su2double& g);
  
  void PropagateSignal(unsigned short iPhi);
  void Preprocessing(unsigned short iPhi, unsigned long iIter);
  void CreateUniformGridSignal(unsigned short iPhi);
  void CreateInitialRayTube(unsigned short iPhi);
  void Nonlinearity(unsigned short iPhi);
  void Attenuation(unsigned short iPhi);
  void Relaxation(unsigned short iPhi);
  void Spreading(unsigned short iPhi);
  void Stratification(unsigned short iPhi);
  void Iterate(unsigned short iPhi);
  void WriteSensitivities();

};
