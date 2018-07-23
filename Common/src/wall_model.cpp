/*!
 * \file wall_model.cpp
 * \brief File, which contains the implementation for the wall model functions
 *        for large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/wall_model.hpp"

/* MKL or LAPACK include files, if supported. */
#ifdef HAVE_MKL
#include "mkl.h"
#elif HAVE_LAPACK
#include "lapacke.h"
#endif

void CWallModel1DEQ::Initialize(const unsigned short *intInfo,
                                const su2double      *doubleInfo){

  /* Copy the data from the arguments into the member variables. */
  numPoints      = intInfo[0];
  thickness      = doubleInfo[0];
  expansionRatio = doubleInfo[1];

  /* Allocate the memory for the coordinates of the grid points used
     in the 1D equilibrium wall model. */
  coorGridPoints.resize(numPoints);

  /* Determine the scaled version of the normal coordinates, where the
     first normal coordinate is simply 1.0. */
  su2double currentHeight = 1.0;
  coorGridPoints[0] = 0.0;

  for(unsigned short i=1; i<numPoints; ++i) {
    coorGridPoints[i] = coorGridPoints[i-1] + currentHeight;
    currentHeight    *= expansionRatio;
  }

  /* Determine the scaling factor of the normal coordinates and
     apply the scaling to obtain the correct coordinates. */
  const su2double scaleFact = thickness/coorGridPoints[numPoints-1];

  for(unsigned short i=0; i<numPoints; ++i)
    coorGridPoints[i] *= scaleFact;
}

void CWallModel1DEQ::WallShearStressAndHeatFlux(const su2double rhoExchange,
                                                const su2double velExchange,
                                                const su2double pExchange,
                                                const su2double Wall_HeatFlux,
                                                const bool      HeatFlux_Prescribed,
                                                const su2double Wall_Temperature,
                                                const bool      Temperature_Prescribed,
                                                      su2double &tauWall,
                                                      su2double &qWall,
                                                      su2double &ViscosityWall,
                                                      su2double &kOverCvWall) {

  /*--- Set up vectors ---*/
  vector<su2double> y(numPoints,0.0);
  vector<su2double> u(numPoints,0.0);
  vector<su2double> T(numPoints,0.0);
  vector<su2double> rho(numPoints,0.0);
  vector<su2double> mu(numPoints,0.0);
  vector<su2double> muTurb(numPoints,0.0);
  vector<su2double> muTurb_new(numPoints,0.0);
  vector<su2double> nu(numPoints,0.0);

  vector<su2double> lower(numPoints-1,0.0);
  vector<su2double> upper(numPoints-1,0.0);
  vector<su2double> diagonal(numPoints,0.0);
  vector<su2double> rhs(numPoints,0.0);

  // Set some constants, assuming air at standard conditions
  // TO DO: Get these values from solver or config classes
  su2double C_1 = 1.458e-6;
  su2double S = 110.4;
  su2double R = 287.058;
  su2double kappa = 0.41;
  su2double A = 17;
  su2double gamma = 1.4;
  su2double Pr_lam = 0.7;
  su2double Pr_turb = 0.9;
  su2double c_p = (gamma*R)/(gamma-1);
  su2double c_v = R/(gamma-1);

  // Set tau wall to zero
  tauWall = 0.0;

  // Set exchange temperature based on pressure and density
  su2double tExchange = pExchange/(rhoExchange*R);

  // Set booleans for run control
  bool converged = false;
  bool initCondExists = false;

  // Set the exchange height equal to the top (last) wall model grid point
  su2double exchangeHeight = coorGridPoints[numPoints-1];

  unsigned short j = 0;
  while( converged == false ){

    // Make sure the tri-diagonal system data is zero'd
    lower.assign(numPoints-1,0.0);
    upper.assign(numPoints-1,0.0);
    diagonal.assign(numPoints,0.0);
    rhs.assign(numPoints,0.0);

    /*--- Set initial condition if it doesn't already exist ---*/
    if( initCondExists == false ){
      for(unsigned short i = 0; i<numPoints; i++){
        // Set y coordinates
        y[i] = coorGridPoints[i];

        // Set linear velocity profile
        u[i] = y[i] * velExchange/exchangeHeight;

        // Set constant temperature profile
        T[i] = tExchange;

        // Set density profile
        rho[i] = pExchange / (R*T[i]);

        // Set the viscosity profile, based on Sutherland's law
        mu[i] = C_1 * pow(T[i],1.5) / (T[i] + S);
        nu[i] = mu[i]/rho[i];
      }
      // Set the initial friction length based on wall shear with linear
      // velocity profile
      tauWall = mu[0] * (u[1]-u[0])/y[1];
      su2double u_tau = sqrt(tauWall/rho[0]);
      su2double l_tau = nu[0] / u_tau;
      for(unsigned short i = 0; i<numPoints; i++){
        /*--- Set the turbulent viscosity ---*/
        su2double y_plus = y[i]/l_tau;
        su2double D = pow(1-exp((-y_plus)/A),2.0);
        muTurb[i] = kappa * rhoExchange * y[i] * u_tau * D;
        //muTurb[i] = 0.0;
      }
      initCondExists = true;
    }

//    // Debugging output
//    for(unsigned short i=0; i<numPoints; i++){
//      cout << u[i] << ", ";
//    }
//    cout << endl;

    /*--- Solve the differential equation
     * d/dy[ (mu + mu_turb) * du/dy) = 0
     * Re-write as
     * d^2u/dy^2 + g(y)*du/dy = 0
     * where g(y) = f'(y) / f(y)
     * f(y) = (mu + mu_turb) ---*/

    for(unsigned short i=1; i<numPoints-1; i++)
    {
      // For equally spaced points, these are the same
      su2double dy_minus = (y[i] - y[i-1]);
      su2double dy_plus = (y[i+1] - y[i]);

      su2double f_minus = mu[i-1] + muTurb[i-1];
      su2double f = mu[i] + muTurb[i];
      su2double f_plus = mu[i+1] + muTurb[i+1];

      su2double f_prime = (f_plus - f_minus) / (dy_minus + dy_plus);

      su2double g = f_prime/f;

      lower[i-1] = (1/(dy_minus*dy_plus) - g/(dy_minus+dy_plus));
      upper[i] = (1/(dy_plus*dy_plus) + g/(dy_minus+dy_plus));
      diagonal[i] = -2.0/(dy_minus*dy_plus);
      rhs[i] = 0.0;
    }

    /*--- Set boundary conditions in matrix and rhs ---*/
    diagonal[0] = 1.0;
    rhs[0] = 0.0;
    diagonal[numPoints-1] = 1.0;
    rhs[numPoints-1] = velExchange;

    // Solve the matrix problem to get the velocity field
    //********LAPACK CALL*******
#if defined (HAVE_LAPACK) || defined(HAVE_MKL)
    int info = 0;
    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,numPoints,1,lower.data(),diagonal.data(),upper.data(),rhs.data(),numPoints);
#else
    SU2_MPI::Error("Not compiled with LAPACK support", CURRENT_FUNCTION);
#endif

    u = rhs;


    lower.assign(numPoints-1,0.0);
    upper.assign(numPoints-1,0.0);
    diagonal.assign(numPoints,0.0);
    rhs.assign(numPoints,0.0);

    // Set up energy equation matrix and rhs
    for(unsigned short i=1; i<numPoints-1; i++)
    {
      su2double dy_minus = (y[i] - y[i-1]);
      su2double dy_plus = (y[i+1] - y[i]);

      su2double f_minus = c_p * (mu[i-1]/Pr_lam + muTurb[i-1]/Pr_turb);
      su2double f = c_p * (mu[i]/Pr_lam + muTurb[i]/Pr_turb);
      su2double f_plus = c_p * (mu[i+1]/Pr_lam + muTurb[i+1]/Pr_turb);

      su2double f_prime = (f_plus - f_minus) / (dy_minus + dy_plus);

      su2double g = f_prime/f;

      lower[i-1] = (1/(dy_minus*dy_plus) - g/(dy_minus+dy_plus));
      upper[i] = (1/(dy_plus*dy_plus) + g/(dy_minus+dy_plus));
      diagonal[i] = -2.0/(dy_minus*dy_plus);

      su2double u_prime = (u[i+1] - u[i-1])/(dy_minus+dy_plus);
      su2double u_prime_prime = (u[i-1] - 2*u[i] + u[i+1])/(dy_minus*dy_plus);

      su2double rhs_1 = -f_prime * u[i] * u_prime;
      su2double rhs_2 = -f * u_prime * u_prime;
      su2double rhs_3 = -f * u[i] * u_prime_prime;
      rhs[i] = rhs_1 + rhs_2 + rhs_3;
    }

    // Set up energy boundary conditions
    if(Temperature_Prescribed == true)
    {
      // Temperature specified at wall
      diagonal[0] = 1.0;
      rhs[0] = Wall_Temperature;

      // Temperature specified by exchange
      diagonal[numPoints-1] = 1.0;
      rhs[numPoints-1] = tExchange;
    }
    else if(HeatFlux_Prescribed == true)
    {
      // Heat flux
      su2double f_zero = c_p * (mu[0]/Pr_lam + muTurb[0]/Pr_turb);
      diagonal[0] = -f_zero/(y[1]-y[0]);
      su2double f_one = c_p * (mu[1]/Pr_lam + muTurb[1]/Pr_turb);
      upper[0] = f_one/(y[1]-y[0]);
      rhs[0] = Wall_HeatFlux;

      // Temperature specified by exchange
      diagonal[numPoints-1] = 1.0;
      rhs[numPoints-1] = tExchange;
    }

    // Solve the matrix problem to get the temperature field
    // *******LAPACK CALL********
#if defined (HAVE_LAPACK) || defined(HAVE_MKL)
    info = 0;
    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,numPoints,1,lower.data(),diagonal.data(),upper.data(),rhs.data(),numPoints);
#else
    SU2_MPI::Error("Not compiled with LAPACK support", CURRENT_FUNCTION);
#endif

    T = rhs;

    // Update solution with new u and T profiles
    for(unsigned short i=0; i<numPoints; i++){
      rho[i] = pExchange / (R*T[i]);
      mu[i] = C_1 * pow(T[i],1.5)/(T[i] + S);
      nu[i] = mu[i]/rho[i];
    }

    tauWall = mu[0] * ( (u[1] - u[0])/(y[1] - y[0]) );
    su2double u_tau = sqrt(tauWall/rho[0]);
    su2double l_tau = nu[0]/u_tau;

    for(unsigned short i=0; i<numPoints; i++){
      su2double y_plus = y[i]/l_tau;
      su2double D = pow((1 - exp(-y_plus/A)),2.0);
      muTurb_new[i] = kappa * u_tau * rho[i] * y[i] * D;
    }

    // Calculate norms and residual
    su2double norm_diff_mu_turb = 0.0;
    su2double norm_mu_turb = 0.0;
    for(unsigned short i=0; i<numPoints; i++){
      norm_diff_mu_turb = norm_diff_mu_turb + abs(muTurb_new[i] - muTurb[i]);
      norm_mu_turb = norm_mu_turb + abs(muTurb[i]);
    }

    su2double residual = norm_diff_mu_turb/norm_mu_turb;

    muTurb = muTurb_new;

    if(residual < 0.000001){
      converged = true;
      tauWall = (mu[0] + muTurb[0]) * ( (u[1] - u[0])/(y[1] - y[0]) );
      qWall = c_p * (mu[0]/Pr_lam + muTurb[0]/Pr_turb) * ( (T[1] - T[0]) / (y[1] - y[0]));
      ViscosityWall = mu[0] + muTurb[0];
      kOverCvWall = c_p / c_v * (mu[0]/Pr_lam + muTurb[0]/Pr_turb);
    }
    else if(j == 50){
      cout << "CWallModel1DEQ::WallShearStressAndHeatFlux: Wall Model did not converge" << endl;
      //      // Debugging output
      //      cout << "tauWall = " << tauWall << endl;
      //      cout << "qWall = " << qWall << endl;
      //      cout << "ViscosityWall = " << ViscosityWall << endl;
      //      cout << "kOverCvWall = " << kOverCvWall << endl;
      //      cout << "y, u, T, mu, muTurb, rho" << endl;
      //      for(unsigned short i=0; i<numPoints; i++){
      //        cout << y[i] << ", ";
      //        cout << u[i] << ", ";
      //        cout << T[i] << ", ";
      //        cout << mu[i] << ", ";
      //        cout << muTurb[i] << ", ";
      //        cout << rho[i] << ", ";
      //        cout << endl;
      //	  }
      SU2_MPI::Error("Did not converge", CURRENT_FUNCTION);
    }

    j++;
  }
}
