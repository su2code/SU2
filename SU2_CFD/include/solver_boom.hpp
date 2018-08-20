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
  unsigned short Kind_Boom_Cost, Kind_Step;
  unsigned long n_prop, m_pow_2;
  unsigned long **pointID_original;
  su2double ***Coord_original;
  su2double CFL_reduce, Step_size, Step_growth;
  bool AD_flag;

  /*---Flight variables---*/
  su2double flt_h;
  su2double flt_M, flt_U;
  su2double flt_psi,    // Heading angle
            flt_gamma,  // Climb angle
            flt_mu;     // Mach angle
  su2double flt_heading[3];

  /*---Atmosphere variables---*/
  su2double atm_g, atm_R, atm_cp;
  su2double T_inf, a_inf, p_inf, rho_inf;

  /*---Ray variables---*/
  unsigned short ray_N_phi;	// Number of rays
  su2double ray_r0;			// Initial ray distance
  su2double ray_t0;			// Initial ray time
  su2double *ray_c0;    // Snell's constant
  su2double *ray_phi;		// Ray angles
  su2double *ray_nu;    // Heading angle of wave normal
  su2double ray_dt, ray_dphi;  // Perturbations for ray tube
  su2double *ray_x, *ray_y, ray_z, ray_A;
  su2double ray_lambda, *ray_gamma, *ray_theta;
  bool ground_flag;     // Whether or not we've propagated to the ground

  /*---Scale variables---*/
  su2double scale_L;

  /*---ABE variables---*/
  su2double p0, 			// Reference pressure
        p_peak,   // Reference pressure for shock formation distance
        f0,       // Reference frequency
  			w0,				// Reference angular frequency
        w0_g,     // Angular frequency at the ground
  			rho0,			// Ambient density (function of altitude)
  			c0,				// Ambient sound speed (function of altitude)
        M_a,      // Acoustic Mach number
        dx_avg,   // Average spacing of signal
  			dsigma,			// Step size
        dsigma_old, // Previous step size
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

  /*---Cost function---*/
  su2double *PLdB;

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
                *dP;        // change in pressure
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
  void HumidityISO(su2double& z0, su2double& p_inf, su2double& T_inf, su2double& h);

  void PropagateSignal(unsigned short iPhi);
  void Preprocessing(unsigned short iPhi, unsigned long iIter);
  void CreateUniformGridSignal(unsigned short iPhi);
  void CreateInitialRayTube(unsigned short iPhi);
  void DetermineStepSize(unsigned short iPhi, unsigned long iIter);
  void Nonlinearity(unsigned short iPhi);
  void Attenuation(unsigned short iPhi);
  void Relaxation(unsigned short iPhi, unsigned long iIter);
  void Scaling(unsigned short iPhi);
  void PerceivedLoudness(unsigned short iPhi);
  void FFT(unsigned long m, su2double *x, su2double *y);
  void FourierTransform(unsigned short iPhi, su2double *w, su2double *p_of_w, su2double& p_dc, unsigned short n_sample);
  // void FourierTransform(unsigned short iPhi, su2double *w, su2double *p_of_w, su2double *f_min, su2double *f_max, unsigned short n_band, unsigned short n_sample_per_band);
  void MarkVII(unsigned short iPhi, su2double *SPL_band, su2double *fc, unsigned short n_band);
  void AcousticEnergy(unsigned short iPhi);
  void WriteGroundPressure(unsigned short iPhi);
  void WriteSensitivities();

};
