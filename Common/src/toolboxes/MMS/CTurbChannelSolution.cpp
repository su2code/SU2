/*!
 * \file CTurbChannelSolution.cpp
 * \brief Implementations of the member functions of CTurbChannelSolution
 * \author E.Molina
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include "../../../include/toolboxes/MMS/CTurbChannelSolution.hpp"
#include <random>

CTurbChannelSolution::CTurbChannelSolution(void) : CVerificationSolution() { }

CTurbChannelSolution::CTurbChannelSolution(unsigned short val_nDim,
                           unsigned short val_nVar,
                           unsigned short val_iMesh,
                           CConfig*       config)
  : CVerificationSolution(val_nDim, val_nVar, val_iMesh, config) {
    
  
  /*--- Perform some sanity and error checks for this solution here. ---*/

  if((config->GetUnsteady_Simulation() != TIME_STEPPING) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_1ST) &&
     (config->GetUnsteady_Simulation() != DT_STEPPING_2ND))
    SU2_MPI::Error("Unsteady mode must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

  if(config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("Compressible flow equations must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

  if((config->GetKind_Solver() != NAVIER_STOKES) &&
     (config->GetKind_Solver() != FEM_NAVIER_STOKES))
    SU2_MPI::Error("Navier Stokes equations must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

  if((config->GetKind_FluidModel() != STANDARD_AIR) &&
     (config->GetKind_FluidModel() != IDEAL_GAS))
    SU2_MPI::Error("Standard air or ideal gas must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

  if(config->GetKind_ViscosityModel() != CONSTANT_VISCOSITY)
    SU2_MPI::Error("Constant viscosity must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);
  
  if(config->GetRef_NonDim() != FREESTREAM_PRESS_EQ_ONE)
    SU2_MPI::Error("FREESTREAM_PRESS_EQ_ONE must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

  if(config->GetKind_ConductivityModel() != CONSTANT_PRANDTL)
    SU2_MPI::Error("Constant Prandtl number must be selected for the Turbulent Channel",
                   CURRENT_FUNCTION);

    
  /*--- Store specific parameters here. ---*/
    
  ReynoldsFriction = 5200.0;  // Friction Reynolds Number.
  TWall            = 273.15; // Wall Temperature
  dx               = 0.2;    // Mesh spacing in the x-direction
  dz               = 0.1;    // Mesh spacing in the z-direction
  Constant_Turb    = 10.;    // Constant to multiply the velocity fluctuation.
  CurrentTime      = config->GetDelta_UnstTime();
    
  /*--- Turbulent flow. Use the relation of Malaspinas and Sagaut,JCP 275 2014,
  to compute the Reynolds number. ---*/
    
  ReynoldsMeanVelocity = pow((8./ 0.073),(4./7.))
                        * pow(ReynoldsFriction,(8./7.));
    
  /*--- Useful coefficients  ---*/
  RGas            = config->GetGas_Constant();
  su2double Gamma = config->GetGamma();
  su2double Cv    = RGas/(Gamma - 1.0);
  su2double Prandtl_Lam = config->GetPrandtl_Lam();
  muLam = config->GetViscosity_FreeStream();
  su2double rhoDim = config->GetDensity_FreeStream();
    
  rhoRef  = config->GetDensity_Ref();
  uRef    = config->GetVelocity_Ref();
  
  ovGm1 = 1.0/(Gamma - 1.0);

  // Friction velocity
  uTau = (2.0*ReynoldsFriction*muLam)/(rhoDim*config->GetLength_Reynolds());
    
  // Friction Length
  lTau = (muLam/rhoDim) / uTau;

  // Compute the wall shear stress and the body force.
  tauWall = rhoDim*uTau*uTau;

  fBodyX = 2.0*tauWall/config->GetLength_Reynolds();

  // The velocity profile is approximated with the following fit
  // u = a0*(1-|2*y/h|^alpha), where a0 = uMean*(alpha+1)/alpha.
  // The coefficient alpha can be computed from the friction and mean
  // velocity Reynolds numbers and uMean from the Mach number.
  // Compute these coefficients below.
  su2double uMean = config->GetModVel_FreeStream();
  alpha = 2.0*ReynoldsFriction*ReynoldsFriction/ReynoldsMeanVelocity - 1.0;
  a0    = uMean*(alpha+1.0)/alpha;

  // From a simple analysis of the fully developed equations it follows
  // that the reduced temperature distribution is theta = 1-|2*y/h|^(alpha+2).
  // When the removed heat equals the work of the body force vector, the
  // temperature in the middle of the channel can be estimated. This is
  // done below.
  halfChan = 0.5*config->GetLength_Reynolds();
  su2double kCondLam = muLam*Gamma*Cv/Prandtl_Lam;
  TMiddle  = TWall + halfChan*halfChan*uMean*fBodyX/(kCondLam*(alpha+2.0));

  // Compute the fully developed RANS solution, if needed.
  ComputeFullyDevelopedRANS(config, yRANS, rhoRANS,
                              uRANS, kRANS, omegaRANS,mutRANS);
  
  // Compute the random numbers for the synthetic turbulence.
  STG_Preprocessing(PhaseMode, RandUnitVec, RandUnitNormal);
  
  
  
  /*--- Write a message that the solution is initialized for the
   Turbulent Channel test case. ---*/
  
  if ((rank == MASTER_NODE) && (val_iMesh == MESH_0)) {
    cout << endl;
    cout << "Warning: Fluid properties and solution are being " << endl;
    cout << "         initialized for the Turbulent Channel case!!!" << endl;
    cout << setprecision(9) << "Friction Reynolds Number: " << ReynoldsFriction << endl;
    cout << setprecision(9) << "Mean Flow Reynolds Number: " << ReynoldsMeanVelocity << endl;
    cout << setprecision(9) << "Body Force: " << fBodyX  << endl;
    cout << setprecision(6) << "YPlus~1.: " << lTau << endl;
    cout << setprecision(4) << "dt+: " << config->GetDelta_UnstTime() * uTau/lTau << endl;
    cout << setprecision(4) << "Total CTUs: " << config->GetTotal_UnstTime() * uMean/ (2.0* PI_NUMBER) << endl;
    cout << endl << flush;
  }
}

CTurbChannelSolution::~CTurbChannelSolution(void) { }

void CTurbChannelSolution::GetSolution(const su2double *val_coords,
                               const su2double val_t,
                               su2double       *val_solution) {

  /* The initial conditions are set for the turbulent channel flow test case.*/

  su2double val_coordsZ      = 0.0;
  if (nDim == 3) val_coordsZ = val_coords[2];

  // Determine the y-coordinate:- Remember that the channel grid is from 0 -> 2*delta
  //                            - RANS equations are only solved for the upper half of the domain;

  su2double y = fabs(val_coords[1]/halfChan - halfChan);

  // Search for the y-coordinate in yRANS.
  unsigned int i;
  for(i=1; i<yRANS.size(); ++i)
    if((y >= yRANS[i-1]) && (y <= yRANS[i])) break;

  // Determine the interpolation weights.
  const su2double wi   = (y - yRANS[i-1])/(yRANS[i] - yRANS[i-1]);
  const su2double wim1 = 1.0 - wi;

  /* Compute the primitive variables and turbulent variables. */
  su2double rho     = (wim1*rhoRANS[i-1] + wi*rhoRANS[i]);
  su2double u       = (wim1*uRANS[i-1]   + wi*uRANS[i]);
  su2double v       = 0.0;
  su2double w       = 0.0;
  su2double turb_ke = (wim1*kRANS[i-1]     + wi*kRANS[i]);
  su2double omega   = (wim1*omegaRANS[i-1] + wi*omegaRANS[i]);
  su2double mut     = (wim1*mutRANS[i-1] + wi*mutRANS[i]);
  su2double dist    = (wim1*distRANS[i-1] + wi*distRANS[i]);
  su2double dudy    = (wim1*dudyRANS[i-1] + wi*dudyRANS[i]);
  
  //Determine the turbulent length scale for the SST model.
  const su2double lengthTurb = sqrt(turb_ke) / (0.09 * omega);

  // Determine the Strain Magnitude
  su2double StrainMag = 0.0;
  StrainMag += 2.0*pow(0.5*dudy, 2.0);
  StrainMag = sqrt(2.0*StrainMag);
  
  // Kolmogorov length scale and wave number
  su2double epsilon = 2.0 * (muLam/rho) * pow(StrainMag,2.0); // TODO: Check the energy dissipation rate eq.
  su2double k_neta = 2.0 * PI_NUMBER / pow( pow((muLam/rho),3.0) / epsilon , 0.25);
  
  // Calculate the rate of strain tensor
  su2double S_ij[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  S_ij[0][1] = dudy;
  S_ij[1][0] = S_ij[0][1];
  
  // Divergence of the velocity component. It should be zero 1D.
  su2double divVel = 0.0;
  
  // Using rate of strain matrix, calculate Reynolds stress tensor
  su2double delta[3][3] = {{1.0, 0.0, 0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  su2double R_ij[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  for (unsigned short iDim = 0; iDim < 3; iDim++){
    for (unsigned short jDim = 0; jDim < 3; jDim++){
      R_ij[iDim][jDim] = TWO3 * rho * turb_ke * delta[iDim][jDim]
      - mut * (2. * S_ij[iDim][jDim] - TWO3 * divVel * delta[iDim][jDim]);
    }
  }
  
  // Cholesky decomposition of the Reynolds Stress tensor
  su2double a_ij[3][3] = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};
  a_ij[0][0] = sqrt(R_ij[0][0]);
  a_ij[1][0] = R_ij[1][0] / a_ij[0][0];
  a_ij[1][1] = sqrt(R_ij[1][1] - pow(a_ij[1][0],2.0));
  a_ij[2][0] = R_ij[2][0] / a_ij[0][0];
  a_ij[2][1] = (R_ij[2][1] - a_ij[1][0] * a_ij[2][0]) / a_ij[1][1];
  a_ij[2][2] = sqrt(R_ij[2][2] - pow(a_ij[2][0],2.0) - pow(a_ij[2][1],2.0));
 
  // Wave number that corresponds to the wavelength of the most energy containing mode.
  su2double l_e = min(2.0 * dist, 3.0 * lengthTurb);
  su2double k_e = 2.0 * PI_NUMBER / l_e;
  
  // Number of wave numbers
  int nWave = (int)WaveNumbers.size();
   
  // Calculate the local Nyquist wave number
  su2double k_cut = 2.0 * PI_NUMBER / min_lengthNyquist;

  // Initialize vector and auxiliary vector of velocity fluctuations.
  su2double VelTurb[3]    = {0.0, 0.0, 0.0};
  su2double VelAuxTurb[3] = {0.0, 0.0, 0.0};

  // Loop over the wave numbers and calculate the spectrum of the kinetic energy.
  su2double E_k_sum = 0.0;
  vector<su2double> E_k, f_neta_, f_cut_, NormalizedAmplitude_n_;
  for(int j=0; j < nWave; j++) {

    // Function that ensures damping of the spectrum near the Kolmogorov scale: Eq.9
    su2double f_neta = exp(- pow(12.0 * WaveNumbers[j] / k_neta,2.0));
    f_neta_.push_back(f_neta);
    
    // (kcut) and function (fcut) that will dump the spectrum at wave
    // numbers larger than kcut: Eq.10 ---*/
    su2double f_cut = exp( - pow( 4.0 * max(WaveNumbers[j] - 0.9 * k_cut, 0.0) / k_cut, 3.0) );
    f_cut_.push_back(f_cut);
    
    // E_k is a prescribed spatial spectrum of the kinetic energy of turbulence represented
    // by a modified von Karman spectrum Eq.7
    su2double E_k_aux = (pow(WaveNumbers[j] / k_e, 4.0) / pow(1.0+2.4 * pow(WaveNumbers[j] / k_e,2.0),17.0/6.0)) * f_cut * f_neta;
    E_k.push_back(E_k_aux);
    E_k_sum += E_k_aux * DeltaWave[j];
    
  }
  
  // Loop again over the wave number and calculate the velocity fluctuations with the
  // contributions of all normalized amplitudes.
  su2double NormalizedAmplitude_sum = 0.0;
  for(int j=0; j < nWave; j++) {
  
    su2double ConvectiveTerm[3] = {0.0, 0.0, 0.0};
    
    // Convective term
    ConvectiveTerm[0] = (2.0 * PI_NUMBER * (val_coords[0] - max_velocity * CurrentTime)) / (WaveNumbers[j] * max_lengthEnergetic);
    ConvectiveTerm[1] = val_coords[1];
    ConvectiveTerm[2] = val_coords[2];
    
    // Dot product between the random vector and convective term.
    su2double dot_prod = 0.0;
    for (unsigned short iDim=0; iDim<nDim; iDim++) dot_prod += RandUnitVec[j*nDim+iDim]*ConvectiveTerm[iDim];
    dot_prod = WaveNumbers[j] * dot_prod;
    
    //  Normalized Amplitude of mode n
    su2double NormalizedAmplitude_n = (E_k[j] * DeltaWave[j]) / max(E_k_sum, 1e-10);
    NormalizedAmplitude_n_.push_back(NormalizedAmplitude_n);
    NormalizedAmplitude_sum += NormalizedAmplitude_n;
        
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      VelAuxTurb[iDim] += sqrt(NormalizedAmplitude_n) * RandUnitNormal[j*nDim+iDim] * cos(dot_prod + PhaseMode[j]);
  }
  
  
  // Calculate the auxiliary vector of velocity fluctuations.
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    VelAuxTurb[iDim] =  2.0 * sqrt(1.5) * VelAuxTurb[iDim];
  
  // Now multiply it by the Cholesky decomposition.
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      VelTurb[iDim] += a_ij[iDim][jDim] * VelAuxTurb[jDim];
    
  // TODO: Check why VelTurb is nan at the wall
  if (std::isnan(VelTurb[0]) || std::isnan(VelTurb[1]) || std::isnan(VelTurb[2])){
    VelTurb[0] = 0.0; VelTurb[1] = 0.0; VelTurb[2] = 0.0;
  }
  
  E_k.clear();
  f_cut_.clear();
  f_neta_.clear();
  NormalizedAmplitude_n_.clear();
  
  /* Compute the conservative variables. Note that both 2D and 3D
   cases are treated correctly. */
  rho /= rhoRef;
  u = (u + Constant_Turb*VelTurb[0])/uRef;
  v = (v + Constant_Turb*VelTurb[1])/uRef;
  w = (w + Constant_Turb*VelTurb[2])/uRef;
  
  val_solution[0]      = rho;
  val_solution[1]      = rho * u;
  val_solution[2]      = rho * v;
  val_solution[3]      = rho * w;
  val_solution[nVar-1] = ovGm1 + 0.5*rho*(u*u + v*v + w*w);
}


//void CTurbChannelSolution::GetSolution(const su2double *val_coords,
//                               const su2double val_t,
//                               su2double       *val_solution) {
//
//  /* The initial conditions are set for the turbulent channel flow test case.*/
//
//  su2double val_coordsZ      = 0.0;
//  if (nDim == 3) val_coordsZ = val_coords[2];
//
//  // Determine the y-coordinate:- Remember that the channel grid is from 0 -> 2*delta
//  //                            - RANS equations are only solved for the upper half of the domain;
//
//  su2double y = fabs(val_coords[1]/halfChan - halfChan);
//
//  // Search for the y-coordinate in yRANS.
//  unsigned int i;
//  for(i=1; i<yRANS.size(); ++i)
//    if((y >= yRANS[i-1]) && (y <= yRANS[i])) break;
//
//  // Determine the interpolation weights.
//  const su2double wi   = (y - yRANS[i-1])/(yRANS[i] - yRANS[i-1]);
//  const su2double wim1 = 1.0 - wi;
//
//  /* Compute the primitive variables. */
//  su2double rho = (wim1*rhoRANS[i-1] + wi*rhoRANS[i])/rhoRef;
//  su2double u   = (wim1*uRANS[i-1]   + wi*uRANS[i])/uRef;
//
//  // Determine the seed of the random generator based on the coordinates.
//  // This means that the random generator will always produce the same
//  // number for the same coordinates.
//  const int seed = (int) (36*val_coords[0] + 821*val_coords[1] + 18955*val_coords[2]);
//  std::default_random_engine rand_gen(seed);
//  std::uniform_real_distribution<su2double> u01(-0.1,0.1);
//
//  /* Compute the conservative variables. Note that both 2D and 3D
//   cases are treated correctly. */
//
//  val_solution[0]      = rho;
//  val_solution[1]      = rho*(u + u*u01(rand_gen));
//  val_solution[2]      = rho*(u*u01(rand_gen));
//  val_solution[3]      = rho*(u*u01(rand_gen));
//  val_solution[nVar-1] = ovGm1 + 0.5*rho*(u*u);
//}

bool CTurbChannelSolution::ExactSolutionKnown(void) {return false;}

// Function, which computes the point distribution of a 1D grid using the
// TANH stretching algorithm of Vinokur.
void CTurbChannelSolution::PointDistributionVinokur(const su2double len,
                              const su2double        spacingEnds,
                               vector<su2double> &yPoints)
{
  //----------------------------------------------------------------------------------
  // The TANH distribution of Vinokur is used, see Vinokur, Marcel, On One-Dimensional
  // Stretching Functions for Finite-Difference Calculations, JCP, 50, 215, 1983.
  // See also http://www.cfd-online.com/Wiki/Structured_mesh_generation.
  //----------------------------------------------------------------------------------

  // Easier storage of the number of points and the inverse of the number
  // of elements.
  const int nPoints        = (int) yPoints.size();
  const su2double nElemInv = 1.0/(nPoints-1);

  // Determine the scaled spacing at the beginning and end of the line. For this
  // particular case this spacing is the same, but the algorithm is capable of
  // handling different values of the begin and end spacing.
  const su2double sBeg = spacingEnds/len;
  const su2double sEnd = sBeg;

  // Compute the parameters A and B. The former is used in the actual distribution
  // function, the latter is used to control the stretching factor.
  const su2double A = sqrt(sEnd/sBeg);
  const su2double B = nElemInv/sqrt(sBeg*sEnd);
  su2double delta;

  // Solve the transcendental equation sinh(delta)/delta = B to obtain
  // the stretching factor delta. Several regions for the value of B must
  // be distinguished.
  if(B < ((su2double) 0.9999999))
  {
    // The solution of the transcendental equation is purely imaginary.
    // Substitution of delta = I*delta transforms the transcendental equation
    // into sin(delta)/delta = B. The Newton algorithm converges quickly.
    delta = 1.0;
    for(;;)
    {
      const su2double tmp = sin(delta)/delta;
      const su2double f   = tmp - B;
      const su2double df  = (cos(delta) - tmp)/delta;

      const su2double ddelta = -f/df;
      delta += ddelta;
      if(fabs(ddelta) < convergenceThreshold) break;
    }

    // Compute the scaled point distribution.
    const su2double f = 1.0/tan(0.5*delta);

    for(int i=0; i<nPoints; ++i)
    {
      const su2double xi = i*nElemInv - 0.5;
      const su2double u  = 0.5*(1.0 + f*tan(delta*xi));
      yPoints[i] = u/(A + (1.0-A)*u);
    }
  }
  else if(B > ((su2double) 1.0000001))
  {
    // The solution of the transcendental equation is real.
    // The Newton algorithm converges quickly.
    delta = 5.0;
    for(;;)
    {
      const su2double tmp = sinh(delta)/delta;
      const su2double f   = tmp - B;
      const su2double df  = (cosh(delta) - tmp)/delta;

      const su2double ddelta = -f/df;
      delta += ddelta;
      if(fabs(ddelta) < convergenceThreshold) break;
    }

    // Compute the scaled point distribution.
    const su2double f = 1.0/tanh(0.5*delta);

    for(int i=0; i<nPoints; ++i)
    {
      const su2double xi = i*nElemInv - 0.5;
      const su2double u  = 0.5*(1.0 + f*tanh(delta*xi));
      yPoints[i] = u/(A + (1.0-A)*u);
    }
  }
  else
  {
    // Special situation. B == 1 and the solution of the transcendental
    // equation is delta = 1. This is a limit situation in which the
    // parameter u is equal to xi.
    for(int i=0; i<nPoints; ++i)
    {
      const su2double xi = i*nElemInv;
      yPoints[i] = xi/(A + (1.0-A)*xi);
    }
  }

  // Scale yPoints, such that the range corresponds to -1 to 1.
  for(int i=0; i<nPoints; ++i)
    yPoints[i] = 2.0*yPoints[i] - 1.0;

  // Make sure that the point distribution is exactly symmetric.
  int ii = nPoints-1;
  for(int i=0; i<=ii; ++i, --ii)
  {
    const su2double val = 0.5*(yPoints[ii] - yPoints[i]);
    yPoints[i]  = -val;
    yPoints[ii] =  val;
  }

  // Make sure that the first point is -1 and the last is 1.0.
  yPoints[0]         = -1.0;
  yPoints[nPoints-1] =  1.0;

  // Multiply the point distribution by 0.5*len to obtain the
  // correct dimensional values.
  for(int i=0; i<nPoints; ++i)
    yPoints[i] *= 0.5*len;
}

// Function, which computes the numerical solution of the fully developed RANS
// equations for a channel flow.
void CTurbChannelSolution::ComputeFullyDevelopedRANS(CConfig*  config,
                                vector<su2double> &yRANS,
                                vector<su2double> &rhoRANS,
                                vector<su2double> &uRANS,
                                vector<su2double> &kRANS,
                                vector<su2double> &omegaRANS,
                                vector<su2double> &mutRANS)
{
  
  
  su2double muLam = config->GetViscosity_FreeStream();
  su2double mu    = config->GetViscosity_FreeStreamND();
  su2double pDim  = config->GetPressure_FreeStream();
  su2double fBody = fBodyX;
  
  su2double Prandtl_Lam = config->GetPrandtl_Lam();
  su2double Prandtl_Turb = config->GetPrandtl_Turb();
  su2double Gamma = config->GetGamma();
  su2double Cv    = RGas/(Gamma - 1.0);
  su2double kCondLam = muLam*Gamma*Cv/Prandtl_Lam;
  
  // Allocate the memory for the vectors.
  yRANS.resize(nGridPoints);
  rhoRANS.resize(nGridPoints);
  uRANS.resize(nGridPoints);
  kRANS.resize(nGridPoints);
  omegaRANS.resize(nGridPoints);
  mutRANS.resize(nGridPoints);
  dudyRANS.resize(nGridPoints);
  distRANS.resize(nGridPoints);

  // Only the master rank computes the numerical solution.
  if(rank == MASTER_NODE)
  {
    // Allocate the memory for TRANS.
     vector<su2double> TRANS(nGridPoints);

    // Determine the required spacing at the wall.
    su2double dyWall = yPlusWall*lTau;

    // Determine the distribution of the grid points. Use the existing function
    // to compute the point distribution over the full channel.
    vector<su2double> yFullChannel(2*nGridPoints-1);
    PointDistributionVinokur(config->GetLength_Reynolds(), dyWall, yFullChannel);

    // Copy the relevant part of the point distribution in yRANS.
    for(int i=0; i<nGridPoints; ++i)
      yRANS[i] = yFullChannel[i+nGridPoints-1];

    // Compute the true wall spacing.
    dyWall = yRANS[nGridPoints-1] - yRANS[nGridPoints-2];

    // Set the wall boundary conditions for the last point.
    rhoRANS[nGridPoints-1]   = pDim/(RGas*TWall);
    uRANS[nGridPoints-1]     = 0.0;
    TRANS[nGridPoints-1]     = TWall;
    kRANS[nGridPoints-1]     = 0.0;
    mutRANS[nGridPoints-1]   = 0.0;
    omegaRANS[nGridPoints-1] = 60.0*mu/(rhoRANS[nGridPoints-1]*beta1*dyWall*dyWall);

    // Initial guess of the density, velocity and turbulent kinetic energy.
    for(int i=0; i<(nGridPoints-1); ++i)
    {
      // Determine the non-dimensional y-coordinate.
      const su2double y = yRANS[i]/halfChan;

      // Compute the velocity.
      uRANS[i] = a0*(1.0 - pow(y,alpha));

      // Compute the dimensional temperature and the density.
      TRANS[i]   = TWall + (TMiddle-TWall)*(1.0 - pow(y,(alpha+2.0)));
      rhoRANS[i] = pDim/(RGas*TRANS[i]);

      // Compute the y+ of this point.
      const su2double dist  = yRANS[nGridPoints-1] - yRANS[i];
      const su2double yPlus = dist/lTau;

      // Crude estimate for the turbulent kinetic energy.
      su2double kPlus = 1.0;
      if(yPlus < ((su2double) 25.0))
        kPlus = ((su2double) 0.16)*yPlus;
      else if(yPlus < ((su2double) 125.0))
        kPlus = 4.0 - ((su2double) 0.03)*(yPlus - ((su2double) 25.0));

      kRANS[i] = kPlus*uTau*uTau;
    }

    // Initial guess of omega.
    for(int i=0; i<(nGridPoints-1); ++i)
    {
      // Compute dudy for this grid point.
      su2double dudy = 0.0;
      if( i ) dudy = (uRANS[i+1] - uRANS[i-1])/(yRANS[i+1] - yRANS[i-1]);

      // Estimate omega from the Wray Agarwal model.
      su2double omega = 0.5*dudy*dudy/sqrt((su2double) 0.09);
      omega =  max(omega, (su2double) 1.e-10);

      // Compute the eddy viscosity. Make sure it is not more than
      // 1000 times the laminar viscosity.
      su2double muT = rhoRANS[i]*kRANS[i]/omega;
      muT =  min(muT, 1000*muLam);

      // Compute the final value of omega.
      omegaRANS[i] = rhoRANS[i]*kRANS[i]/muT;
    }

    // Define the values for the minimum and maximum value of omega.
    const su2double omegaMax = omegaRANS[nGridPoints-1];
    const su2double omegaMin = omegaRANS[0]/100;

    // Define the vectors to store the right hand side and the Jacobian matrix.
     vector<su2double> RHS(2*(nGridPoints-1));
     vector<su2double> aJac(4*(nGridPoints-1)), bJac(4*(nGridPoints-1)),
                           cJac(4*(nGridPoints-1));

    // Start of the iterative loop to compute the numerical solution.
    su2double duMaxOld = 100, ratioMuTurbMuLamMax;
    for(int iter=0;; ++iter)
    {
      //------------------------------------------------------------------------
      // Part 1: Iteration on the turbulence equations.
      //------------------------------------------------------------------------

      // Loop over the grid points to compute the source terms and the
      // contribution of the source terms to the central Jacobian.
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        // Set the pointers for residual and central Jacobian.
        su2double *res = RHS.data()  + 2*i;
        su2double *b   = bJac.data() + 4*i;

        // Compute the derivatives the u, k and ln(omega) in this grid point.
        // Note that for the first grid point the symmetry BC is applied.
        su2double dyCell = yRANS[i+1] - yRANS[ max(i-1,0)];
        if( i ) dyCell *= 0.5;

        const su2double dyInv = 1.0/(yRANS[i+1] - yRANS[ max(i-1,0)]);

        const su2double dudy       = i ? dyInv*(uRANS[i+1] - uRANS[i-1]) : 0.0;
        const su2double dkdy       = i ? dyInv*(kRANS[i+1] - kRANS[i-1]) : 0.0;
        const su2double dlnomegady = i ? dyInv*(log(omegaRANS[i+1]) - log(omegaRANS[i-1])) : 0.0;

        // Determine the distance to the wall.
        const su2double dist = yRANS[nGridPoints-1] - yRANS[i];

        // Compute the value of the cross diffusion term.
        const su2double CDOmega = 2.0*rhoRANS[i]*sigmaOm2*dkdy*dlnomegady;
        su2double CD =  max(CDOmega, ((su2double) 1.e-20));

        // Determine the three expressions, which determine the arguments
        // for the blending functions.
        const su2double tt1 = sqrt(kRANS[i])/(betaStar*omegaRANS[i]*dist);
        const su2double tt2 = 500*muLam/(rhoRANS[i]*omegaRANS[i]*dist*dist);
        const su2double tt3 = 4.0*rhoRANS[i]*sigmaOm2*kRANS[i]/(CD*dist*dist);

        // Compute the arguments for the blending functions.
        const su2double arg1 =  min( max(tt1,tt2),tt3);
        const su2double arg2 =  max(2.0*tt1,tt2);

        // Compute the blending functions F1 and F2.
        const su2double F1 = tanh(arg1*arg1*arg1*arg1);
        const su2double F2 = tanh(arg2*arg2);

        // Compute the blended values of beta and gam.
        const su2double beta = F1*beta1 + (1.0-F1)*beta2;
        const su2double gam  = F1*gam1  + (1.0-F1)*gam2;

        // Compute the eddy viscosity.
        const su2double vortMag = fabs(dudy);
        const su2double muTurb  = a1*rhoRANS[i]*kRANS[i]/ max(a1*omegaRANS[i],vortMag*F2);

        // Compute the destruction and production term of the k-equation.
        // The latter should be limited.
        const su2double destrK = betaStar*rhoRANS[i]*omegaRANS[i]*kRANS[i];
        const su2double prodK  =  min(muTurb*dudy*dudy, 20*destrK);

        // Compute the destruction and production term of the omega equation.
        const su2double destrOmega = beta*rhoRANS[i]*omegaRANS[i]*omegaRANS[i];
        const su2double prodOmega  = gam*rhoRANS[i]*dudy*dudy;

        // Compute the source terms.
        res[0] = dyCell*(prodK - destrK);
        res[1] = dyCell*((1.0-F1)*CDOmega + prodOmega - destrOmega);

        // Compute the central Jacobian. Note that only the destruction
        // terms are linearized to improve stability.
        b[0] = -dyCell*betaStar*rhoRANS[i]*omegaRANS[i];
        b[1] =  0.0;
        b[2] = -dyCell*betaStar*rhoRANS[i]*kRANS[i];
        b[3] = -dyCell*2.0*beta*rhoRANS[i]*omegaRANS[i];
      }

      // Loop over the faces to compute the diffusive fluxes.
      // These are scattered to the adjacent points to compute the residual.
      // The corresponding Jacobians are updated accordingly.
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        // This is the face between points i and i+1. Compute the average
        // quantities and the gradients, which are needed to compute the
        // diffusive flux of the flow equations.
        const su2double T   = 0.5*(TRANS[i] + TRANS[i+1]);
        const su2double rho = pDim/(RGas*T);

        const su2double k     = 0.5*(kRANS[i]     + kRANS[i+1]);
        const su2double omega = 0.5*(omegaRANS[i] + omegaRANS[i+1]);

        const su2double dyInv      = 1.0/(yRANS[i+1] - yRANS[i]);
        const su2double dudy       = dyInv*(uRANS[i+1] - uRANS[i]);
        const su2double dkdy       = dyInv*(kRANS[i+1] - kRANS[i]);
        const su2double domegady   = dyInv*(omegaRANS[i+1] - omegaRANS[i]);
        const su2double dlnomegady = dyInv*(log(omegaRANS[i+1]) - log(omegaRANS[i]));

        // Determine the distance to the wall.
        const su2double dist = yRANS[nGridPoints-1] - 0.5*(yRANS[i]+yRANS[i+1]);
        
        distRANS[i] = dist;
        dudyRANS[i] = dudy;

        // Compute the value of the cross diffusion term.
        const su2double CDOmega = 2.0*rho*sigmaOm2*dkdy*dlnomegady;
        su2double CD =  max(CDOmega, ((su2double) 1.e-20));

        // Determine the three expressions, which determine the arguments
        // for the blending functions.
        const su2double tt1 = sqrt(k)/(betaStar*omega*dist);
        const su2double tt2 = 500*muLam/(rho*omega*dist*dist);
        const su2double tt3 = 4.0*rho*sigmaOm2*k/(CD*dist*dist);

        // Compute the arguments for the blending functions.
        const su2double arg1 =  min( max(tt1,tt2),tt3);
        const su2double arg2 =  max(2.0*tt1,tt2);

        // Compute the blending functions F1 and F2.
        const su2double F1 = tanh(arg1*arg1*arg1*arg1);
        const su2double F2 = tanh(arg2*arg2);

        // Compute the blended values of sigmaK and sigmaOm
        const su2double sigmaK  = F1*sigmaK1  + (1.0-F1)*sigmaK2;
        const su2double sigmaOm = F1*sigmaOm1 + (1.0-F1)*sigmaOm2;

        // Compute the eddy viscosity.
        const su2double vortMag = fabs(dudy);
        const su2double muTurb  = a1*rho*k/ max(a1*omega,vortMag*F2);

        // Compute the effective viscosities for the k- and omega-equation.
        const su2double muK  = muLam + sigmaK*muTurb;
        const su2double muOm = muLam + sigmaOm*muTurb;

        // Compute the fluxes at the face.
        const su2double fK  = muK*dkdy;
        const su2double fOm = muOm*domegady;

        // Compute the derivatives of fK and fOm w.r.t. the nodal values
        // of k and omega. The value of muTurb is frozen.
        const su2double dfKdkip1 =  dyInv*muK;
        const su2double dfKdki   = -dyInv*muK;

        const su2double dfOmdomip1 =  dyInv*muOm;
        const su2double dfOmdomi   = -dyInv*muOm;

        // Set the pointers for the residual and Jacobian matrices of point i.
        su2double *res = RHS.data()  + 2*i;
        su2double *b   = bJac.data() + 4*i;
        su2double *c   = cJac.data() + 4*i;

        // Update the residuals and Jacobians of point i.
        res[0] += fK;
        res[1] += fOm;

        b[0] += dfKdki;
        b[3] += dfOmdomi;

        c[0] = dfKdkip1;
        c[1] = 0.0;
        c[2] = 0.0;
        c[3] = dfOmdomip1;

        // Check for the first face.
        if(i == 0)
        {
          // The flux through the mirror face must be computed as well. When the
          // symmetry boundary condition is applied, the value of a quantity is
          // the same, while the y-derivative is the opposite. Hence the fluxes
          // are the opposite. As these fluxes are substracted the net effect
          // is adding the current fluxes once more.
          res[0] += fK;
          res[1] += fOm;

          b[0] += dfKdki;
          b[3] += dfOmdomi;

          c[0] += dfKdkip1;
          c[3] += dfOmdomip1;
        }

        // Check if this is not the last face.
        if(i < (nGridPoints-2))
        {
          // Set the pointers for the residual and Jacobian matrices of point i+1.
          res = RHS.data()  + 2*(i+1);
          b   = bJac.data() + 4*(i+1);
          su2double *a = aJac.data() + 4*(i+1);

          res[0] -= fK;
          res[1] -= fOm;

          a[0] = -dfKdki;
          a[1] =  0.0;
          a[2] =  0.0;
          a[3] = -dfOmdomi;

          b[0] -= dfKdkip1;
          b[3] -= dfOmdomip1;
        }
      }

      // Underrelaxation of the central Jacobian. The multiplication with 1.01
      // corresponds to a CFL of 100.
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        su2double *bb = bJac.data() + 4*i;
        bb[0] *= (su2double) 1.01;
        bb[3] *= (su2double) 1.01;
      }

      // Solve the block-tri-diagonal system to compute the updates of k and omega.
      BlockThomas(nGridPoints-1, aJac, bJac, cJac, RHS);

      // Update k and omega.
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        kRANS[i]     -= RHS[2*i];
        omegaRANS[i] -= RHS[2*i+1];

        kRANS[i]     =  max(kRANS[i],0.0);
        omegaRANS[i] =  min( max(omegaRANS[i],omegaMin),omegaMax);
      }

      //------------------------------------------------------------------------
      // Part 2: Iteration on the flow equations.
      //------------------------------------------------------------------------

      // Loop over the grid points to compute the source terms and the
      // contribution of the source terms to the central Jacobian.
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        // Set the pointers for residual and central Jacobian.
        su2double *res = RHS.data()  + 2*i;
        su2double *b   = bJac.data() + 4*i;

        // Compute the source terms.
        su2double dyCell = yRANS[i+1] - yRANS[ max(i-1,0)];
        if( i ) dyCell *= 0.5;

        res[0] = dyCell*fBody;
        res[1] = dyCell*fBody*uRANS[i];

        // Compute the contributions to the central Jacobian.
        b[0] = 0.0;
        b[1] = dyCell*fBody;
        b[2] = 0.0;
        b[3] = 0.0;
      }

      // Loop over the faces to compute the diffusive fluxes.
      // These are scattered to the adjacent points to compute the residual.
      // The corresponding Jacobians are updated accordingly.
      ratioMuTurbMuLamMax = 0.0;
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        // This is the face between points i and i+1. Compute the average
        // quantities and the gradients, which are needed to compute the
        // diffusive flux of the flow equations.
        const su2double u   = 0.5*(uRANS[i] + uRANS[i+1]);
        const su2double T   = 0.5*(TRANS[i] + TRANS[i+1]);
        const su2double rho = pDim/(RGas*T);

        const su2double k     = 0.5*(kRANS[i]     + kRANS[i+1]);
        const su2double omega = 0.5*(omegaRANS[i] + omegaRANS[i+1]);

        const su2double dyInv = 1.0/(yRANS[i+1] - yRANS[i]);
        const su2double dudy  = dyInv*(uRANS[i+1] - uRANS[i]);
        const su2double dTdy  = dyInv*(TRANS[i+1] - TRANS[i]);

        // Compute the arguments for the blending function F2.
        const su2double dist = yRANS[nGridPoints-1] - 0.5*(yRANS[i]+yRANS[i+1]);
        const su2double tt1  = sqrt(k)/(betaStar*omega*dist);
        const su2double tt2  = 500*muLam/(rho*omega*dist*dist);

        // Compute the eddy viscosity.
        const su2double arg2    =  max(2.0*tt1,tt2);
        const su2double F2      = tanh(arg2*arg2);
        const su2double vortMag = fabs(dudy);
        const su2double muTurb  = a1*rho*k/ max(a1*omega,vortMag*F2);
        mutRANS[i] = muTurb;

        // Update the ratio ratioMuTurbMuLamMax.
        const su2double ratioMuTurbMuLam = muTurb/muLam;
        ratioMuTurbMuLamMax =  max(ratioMuTurbMuLamMax, ratioMuTurbMuLam);

        // Compute the total viscosity and thermal conductivity.
        const su2double mu    = muLam + muTurb;
        const su2double kCond = kCondLam + muTurb*Gamma*Cv/Prandtl_Turb;

        // Compute the fluxes at the face.
        const su2double fMomX = mu*dudy;
        const su2double fEner = mu*u*dudy + kCond*dTdy;

        // Compute the derivatives of fMomX and fEner w.r.t. the nodal values
        // of u and T. The values of muTurb and kCondTurb are frozen.
        const su2double dfMomXduip1 =  dyInv*mu;
        const su2double dfMomXdui   = -dyInv*mu;

        const su2double dfEnerduip1 = 0.5*mu*dudy + u*dyInv*mu;
        su2double dfEnerdui   = 0.5*mu*dudy - u*dyInv*mu;

        const su2double dfEnerdTip1 =  dyInv*kCond;
        const su2double dfEnerdTi   = -dyInv*kCond;

        // Set the pointers for the residual and Jacobian matrices of point i.
        su2double *res = RHS.data()  + 2*i;
        su2double *b   = bJac.data() + 4*i;
        su2double *c   = cJac.data() + 4*i;

        // Update the residuals and Jacobians of point i.
        res[0] += fMomX;
        res[1] += fEner;

        b[0] += dfMomXdui;
        b[1] += dfEnerdui;
        b[3] += dfEnerdTi;

        c[0] = dfMomXduip1;
        c[1] = dfEnerduip1;
        c[2] = 0.0;
        c[3] = dfEnerdTip1;

        // Check for the first face.
        if(i == 0)
        {
          // The flux through the mirror face must be computed as well. When the
          // symmetry boundary condition is applied, the value of a quantity is
          // the same, while the y-derivative is the opposite. Hence the fluxes
          // are the opposite. As these fluxes are substracted the net effect
          // is adding the current fluxes once more.
          res[0] += fMomX;
          res[1] += fEner;

          b[0] += dfMomXdui;
          b[1] += dfEnerdui;
          b[3] += dfEnerdTi;

          c[0] += dfMomXduip1;
          c[1] += dfEnerduip1;
          c[3] += dfEnerdTip1;
        }

        // Check if this is not the last face.
        if(i < (nGridPoints-2))
        {
          // Set the pointers for the residual and Jacobian matrices of point i+1.
          res = RHS.data()  + 2*(i+1);
          b   = bJac.data() + 4*(i+1);
          su2double *a = aJac.data() + 4*(i+1);

          res[0] -= fMomX;
          res[1] -= fEner;

          a[0] = -dfMomXdui;
          a[1] = -dfEnerdui;
          a[2] =  0.0;
          a[3] = -dfEnerdTi;

          b[0] -= dfMomXduip1;
          b[1] -= dfEnerduip1;
          b[3] -= dfEnerdTip1;
        }
      }

      // Solve the block-tri-diagonal system to compute the updates of
      // the velocities and temperature.
      BlockThomas(nGridPoints-1, aJac, bJac, cJac, RHS);

      // Update the velocities and temperature and compute the new values of
      // the density. Also compute the maximum update of the velocity.
      su2double duMax = 0.0;
      for(int i=0; i<(nGridPoints-1); ++i)
      {
        uRANS[i]  -= RHS[2*i];
        TRANS[i]  -= RHS[2*i+1];
        rhoRANS[i] = pDim/(RGas*TRANS[i]);

        duMax =  max(duMax, fabs(RHS[2*i]));
      }

      // Exit criterion.
      if((duMax < 1.e-3) && ((duMax > duMaxOld) || (duMax < 1.e-6))) break;
      duMaxOld = duMax;
    }

     cout << "#" <<  endl
              << "# Fully developed RANS solution computed." <<  endl
              << "# Maximum value of muTurb/muLam: " << ratioMuTurbMuLamMax <<  endl
              << "#" <<  endl;
  }

  // In MPI mode the RANS solution must be broadcast to the other ranks.
#ifdef HAVE_MPI
  SU2_MPI::Bcast(yRANS.data(),     yRANS.size(),     MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(rhoRANS.data(),   rhoRANS.size(),   MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(uRANS.data(),     uRANS.size(),     MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(kRANS.data(),     kRANS.size(),     MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(omegaRANS.data(), omegaRANS.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(mutRANS.data(),   mutRANS.size(),     MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(distRANS.data(),  distRANS.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(dudyRANS.data(),  dudyRANS.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
#endif
}

// Function, which solves a block tri-diagonal matrix using the blocked
// version of the Thomas algorithm. The blocks are (2X2) matrices.
void CTurbChannelSolution::BlockThomas(const int              n,
                                       vector<su2double> &a,
                                       vector<su2double> &b,
                                       vector<su2double> &c,
                                       vector<su2double> &d)
{
  //----------------------------------------------------------------------------
  // Forward sweep. Create the upper diagonal block form.
  //----------------------------------------------------------------------------

  // Loop over the number of block rows.
  for(int i=1; i<n; ++i)
  {
    // Set the required pointers.
    su2double *ai   = a.data() + 4*i;
    su2double *bi   = b.data() + 4*i;
    su2double *bim1 = b.data() + 4*(i-1);
    su2double *cim1 = c.data() + 4*(i-1);
    su2double *di   = d.data() + 2*i;
    su2double *dim1 = d.data() + 2*(i-1);

    // Determine the inverse of bim1.
    const su2double detInv =  1.0/(bim1[0]*bim1[3] - bim1[1]*bim1[2]);
    const su2double b00Inv =  detInv*bim1[3];
    const su2double b10Inv = -detInv*bim1[1];
    const su2double b01Inv = -detInv*bim1[2];
    const su2double b11Inv =  detInv*bim1[0];

    // Store the inverse of bim1 for the backward sweep.
    bim1[0] = b00Inv; bim1[1] = b10Inv; bim1[2] = b01Inv; bim1[3] = b11Inv;

    // Determine the weight matrix, which is the multiplication matrix
    // of row (i-1), which is subtracted from the current row.
    const su2double w00 = ai[0]*b00Inv + ai[2]*b10Inv;
    const su2double w10 = ai[1]*b00Inv + ai[3]*b10Inv;
    const su2double w01 = ai[0]*b01Inv + ai[2]*b11Inv;
    const su2double w11 = ai[1]*b01Inv + ai[3]*b11Inv;

    // Update the diagonal matrix of the current row.
    bi[0] -= w00*cim1[0] + w01*cim1[1];
    bi[1] -= w10*cim1[0] + w11*cim1[1];
    bi[2] -= w00*cim1[2] + w01*cim1[3];
    bi[3] -= w10*cim1[2] + w11*cim1[3];

    // Update the rhigh hand side of the current row.
    di[0] -= w00*dim1[0] + w01*dim1[1];
    di[1] -= w10*dim1[0] + w11*dim1[1];
  }

  //----------------------------------------------------------------------------
  // Backward sweep. Compute the actual solution.
  //----------------------------------------------------------------------------

  // Compute the solution of the last row.
  su2double *bb = b.data() + 4*(n-1);
  su2double *dd = d.data() + 2*(n-1);

  su2double detInv = 1.0/(bb[0]*bb[3] - bb[1]*bb[2]);
  su2double dd0 = dd[0], dd1 = dd[1];
  dd[0] = detInv*(dd0*bb[3] - dd1*bb[2]);
  dd[1] = detInv*(dd1*bb[0] - dd0*bb[1]);

  // Backward loop over the block rows.
  for(int i=(n-2); i>=0; --i)
  {
    // Set the required pointers to the current and next row.
    su2double *bi   = b.data() + 4*i;
    su2double *ci   = c.data() + 4*i;
    su2double *di   = d.data() + 2*i;
    su2double *dip1 = d.data() + 2*(i+1);

    // Subtract the solution of the next row.
    di[0] -= ci[0]*dip1[0] + ci[2]*dip1[1];
    di[1] -= ci[1]*dip1[0] + ci[3]*dip1[1];

    // Determine the solution of the current row. Note that inverse of
    // bi is stored.
    dd0 = di[0]; dd1 = di[1];

    di[0] = bi[0]*dd0 + bi[2]*dd1;
    di[1] = bi[1]*dd0 + bi[3]*dd1;
  }
}

void CTurbChannelSolution::STG_Preprocessing(vector<su2double> &PhaseMode,
                                        vector<su2double> &RandUnitVec,
                                        vector<su2double> &RandUnitNormal)
{
  
   // Generate random numbers
    if (rank == MASTER_NODE){
      std::default_random_engine rand_gen;
      
      // a somewhat random seed
      rand_gen.seed((int)time(0));
      uniform_real_distribution<su2double> u02pi(0.,2.0*PI_NUMBER);
      uniform_real_distribution<su2double> u01(0.,1.0);
      
      su2double theta, phi, x, y, z;
      su2double thetan, phin, xn, yn, zn;
      for (unsigned long i = 0; i < NModes; ++i){
        PhaseMode.push_back(u02pi(rand_gen));
        
        // Generate random number over a sphere with radius 1.
        theta = 2. * PI_NUMBER * u01(rand_gen);
        phi   = acos(1. - 2. * u01(rand_gen));
        x = sin(phi) * cos(theta);
        y = sin(phi) * sin(theta);
        z = cos(phi);

        thetan = 2. * PI_NUMBER * u01(rand_gen);
        phin   = acos(1. - 2. * u01(rand_gen));
        xn = sin(phin) * cos(thetan);
        yn = sin(phin) * sin(thetan);
        zn = cos(phin);
        
        // Create a normal vector
        su2double norm = x * xn + y * yn + z * zn;
        xn = x - xn*norm;
        yn = y - yn*norm;
        zn = z - zn*norm;
        
        RandUnitVec.push_back(x);
        RandUnitVec.push_back(y);
        RandUnitVec.push_back(z);
        RandUnitNormal.push_back(xn);
        RandUnitNormal.push_back(yn);
        RandUnitNormal.push_back(zn);
        
      }
    }
    
    // Resize memory for vectors in other ranks
    if (rank != MASTER_NODE){
      PhaseMode.resize(NModes);
      RandUnitVec.resize(NModes*nDim);
      RandUnitNormal.resize(NModes*nDim);
    }

  #ifdef HAVE_MPI
    SU2_MPI::Bcast(PhaseMode.data(), PhaseMode.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(RandUnitVec.data(), RandUnitVec.size() , MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
    SU2_MPI::Bcast(RandUnitNormal.data(), RandUnitVec.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  #endif
  
  // Number of waves (modes)
  int nWave;
  
  // Calculate the most energetic scale, Nyquist wave length and maximum velocity.
  max_lengthEnergetic = 0.0; max_velocity = 0.0; min_lengthNyquist = 1e10;
  
  if (rank == MASTER_NODE){

    for(int i=0; i<(nGridPoints-1); ++i){
      
      // Determine the distance to the wall.
      const su2double dist = yRANS[nGridPoints-1] - yRANS[i];
      
      // Determine the maximum length
      const su2double maxLength = max(max(dx,dz), yRANS[i+1]-yRANS[i]);
      
      // Determine K, Omega and Density
      const su2double turb_ke  = kRANS[i];
      const su2double omega    = omegaRANS[i];
      const su2double density  = rhoRANS[i];
      const su2double velocity = uRANS[i];
      
      //Determine the turbulent length scale for the SST model.
      const su2double lengthTurb = sqrt(turb_ke) / (0.09 * omega);
      su2double maxdelta = max(yRANS[i+1]-yRANS[i], dz);
      max_velocity = max(max_velocity, sqrt(pow(velocity,2.0)));
      max_lengthEnergetic = max(max_lengthEnergetic, min( 2.0 * dist, 3.0 * lengthTurb ));
      min_lengthNyquist   = min(2.0 * min( max(maxdelta, 0.3*maxLength) + 0.1 * dist, maxLength ), min_lengthNyquist);
      
    }
    
    // Determine the Wave Numbers;
    vector<su2double> WaveNumbersFace; // Aux Vector.
    
    // Calculate the local most energetic wave number
    su2double k_min = 2.0 * PI_NUMBER / max_lengthEnergetic;
    
    // Calculate the local Nyquist wave number
    su2double k_cut = 2.0 * PI_NUMBER / min_lengthNyquist;
    
    cout << "# Most Energetic wave number: " << k_min << endl;
    cout << "# Nyquist wave number: " << k_cut << endl;
    cout << "# Maximum Velocity at the interface: " <<  max_velocity << endl;
    
    /*--- Wave numbers at faces---*/
    su2double alpha = 0.01;
    unsigned short iMode = -1;
    while (k_min < 1.5*k_cut){
      iMode += 1;
      k_min = k_min * pow(alpha + 1.0, iMode);
      WaveNumbersFace.push_back(k_min);
    }
    
    // Wave numbers at the center---*/
    for(int j = 0; j < ((int)WaveNumbersFace.size() - 1); j++){
      WaveNumbers.push_back(0.5 * (WaveNumbersFace[j] + WaveNumbersFace[j+1]));
      DeltaWave.push_back(WaveNumbersFace[j+1] - WaveNumbersFace[j]);
    }
    
    nWave = (int)WaveNumbers.size();
    cout << "# Number of modes: " << nWave << endl;
    cout << "#" << endl;
        
  }

#ifdef HAVE_MPI
  
  SU2_MPI::Bcast(&max_lengthEnergetic, 1, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(&max_velocity, 1, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(&min_lengthNyquist, 1, MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(&nWave, 1, MPI_INT,  MASTER_NODE, MPI_COMM_WORLD);

  if (rank!=MASTER_NODE){
    WaveNumbers.resize(nWave);
    DeltaWave.resize(nWave);
  }

  SU2_MPI::Bcast(WaveNumbers.data(), WaveNumbers.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
  SU2_MPI::Bcast(DeltaWave.data(), DeltaWave.size(), MPI_DOUBLE,  MASTER_NODE, MPI_COMM_WORLD);
#endif
  
}
