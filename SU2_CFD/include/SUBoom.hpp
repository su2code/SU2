#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/matrix_structure.hpp"
#include "../../Common/include/vector_structure.hpp"
#include "../../Common/include/linear_solvers_structure.hpp"
#include "../../Common/include/grid_movement_structure.hpp"
#include "../../SU2_CFD/include/solver_structure.hpp"

using namespace std;

class SUBoom{
public:
  /*---Flight variables---*/
  double flt_h;
  double flt_M;
  double flt_psi, flt_gamma, flt_mu;
  double flt_heading[3];

  /*---Atmosphere variables---*/
  double atm_g;
  unsigned int atm_noise_flag;

  double T_inf, a_inf, p_inf, rho_inf;
  double z[45001], T_of_z[45001], a_of_z[45001], p_of_z[45001], rho_of_z[45001];

  /*---Scale factors---*/
  double scale_L, scale_T, scale_p, scale_m, scale_z;
  double scale_C1, scale_C2;

  /*---Ray variables---*/
  unsigned int ray_N_phi;
  double ray_r0;
  double *ray_t0;
  double *ray_phi;
  double **ray_nu;
  double **ray_theta0;
  double **ray_c0;
  double ***theta;
  double ***x_of_z, ***y_of_z, ***t_of_z;
  //double ***a_of_z, ***rho_of_z, ***p_of_z;
  double ***dxdt, ***dydt, ***dzdt, ***tangent;
  double **ray_A;
  double **ray_C1, **ray_C2, **ray_dC1, **ray_dC2;

  /*---Tolerance variables---*/
  double tol_l, tol_m;
  double tol_dp, tol_dr, tol_dphi;
  double tol_integrate;

  /*---Boom strength---*/
  double p_rise, p_max, p_rise2;

  /*---Ray data class for passing into RK4---*/
  class RayData{
    public:
      double L;
      double T;
      double c0;
      double nu;

      double diff_m2, diff_p2;
      double scale_C1, scale_C2;
      double *t, *C1, *C2, *dC1, *dC2;
      double M;

      //void GetAtmosphericData(double& a, double& rho, double& p, double h);
  };

  /*---Signal class for storing pressure signal---*/
  class Signal{
    public:
      int original_len, M, final_M;
      double *m, *dp, *l, *fvec;
      double *x, *original_T, *original_p;
      double *final_T, *final_p;
  };

  Signal signal;

  /*---Main functions for boom calculations---*/
  SUBoom();
  SUBoom(CSolver *solver, CConfig *config, CGeometry *geometry);
/*double h0, double M, double psi, double gamma, double g,
         unsigned int noise_flag, unsigned int N_phi,
         double phi[], double dphi, double r0, double dr,
	     double m, double dp, double l);*/

  void AtmosISA(double& h0, double& T, double& a, double& p,
                double& rho, double& g);
  void ConditionAtmosphericData();
  void ScaleFactors();
  void InitialWaveNormals();
  void RayTracer();
  void RayTubeArea();
  void FindInitialRayTime();
  void ODETerms();
  void DistanceToTime();
  void CreateSignature(); // Note: p is pressure disturbance (relative)
  void PropagateSignal();

  void GetAtmosphericData(double& a, double& rho, double& p, double h);
  void Sph2Cart(double& nx, double& ny, double& nz, double az, double elev,
                  double r);
  double *rk4(double x0, int m, double y0[], double dx, RayData data,
              double *f(double x, int m, double y[], RayData data));
  double *SplineGetDerivs(double x[], double y[], int N);
  double matchr(int i, int j, double h_L, double r0);
  double *ClipLambdaZeroSegment(double fvec[], int &M);
};
double *derivs(double x, int m, double y[], SUBoom::RayData data);
double *derivsProp(double t, int m, double y[], SUBoom::RayData data);
double **WaveformToPressureSignal(double fvec[], int M, int &Msig);
double EvaluateSpline(double x, int N, double t[], double fit[], double coeffs[]);
