/*!
 * \file wall_model.cpp
 * \brief File, which contains the implementation for the wall model functions
 *        for large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 6.2.0 "Falcon"
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

/* Prototypes for Lapack functions, if MKL or LAPACK is used. */
#if defined (HAVE_MKL) || defined(HAVE_LAPACK)
extern "C" void dgtsv_(int*, int*, passivedouble*, passivedouble*,
                       passivedouble*, passivedouble*, int*, int*);
#endif

void CWallModel1DEQ::Initialize(const unsigned short *intInfo,
                                const su2double      *doubleInfo){

  /* Copy the data from the arguments into the member variables. */
  numPoints      = intInfo[0];
  h_wm           = doubleInfo[0];
  expansionRatio = doubleInfo[1];

  unsigned short nfa = numPoints + 1;
  su2double currentHeight = 1.0;
  
  /* Allocate the memory for the coordinates of the grid points used
     in the 1D equilibrium wall model. */
  y_cv.resize(numPoints);
  y_fa.resize(numPoints+1);

  /* Determine the scaled version of the normal coordinates, where the
   first normal coordinate is simply 1.0. */
  
  y_fa[0] = 0.0;
  for(unsigned short i=1; i<nfa; ++i) {
    y_fa[i] = y_fa[i-1] + currentHeight;
    y_cv[i-1] = 0.5 * (y_fa[i] + y_fa[i-1]);
    currentHeight = currentHeight*expansionRatio;
  }
  
  su2double y_max = y_cv[numPoints-1];
  for(unsigned short i=0; i<nfa; ++i) {
    y_fa[i] = y_fa[i]/y_max * h_wm;
  }
  for(unsigned short i=0; i<numPoints; ++i) {
    y_cv[i] = y_cv[i]/y_max * h_wm;
  }
}

void CWallModel1DEQ::WallShearStressAndHeatFlux(const su2double tExchange,
                                                const su2double velExchange,
                                                const su2double muExchange,
                                                const su2double pExchange,
                                                const su2double Wall_HeatFlux,
                                                const bool      HeatFlux_Prescribed,
                                                const su2double TWall,
                                                const bool      Temperature_Prescribed,
                                                      su2double &tauWall,
                                                      su2double &qWall,
                                                      su2double &ViscosityWall,
                                                      su2double &kOverCvWall) {

  
  // Set tau wall to initial guess
  tauWall = 0.5;
  qWall = 0.0;
  ViscosityWall = 0.0;
  kOverCvWall = 0.0;
  
  // Set some constants, assuming air at standard conditions
  // TO DO: Get these values from solver or config classes
  //su2double C_1 = 2.03929e-04;
  su2double C_1 = 1.716e-5;
  su2double S = 110.4;
  su2double T_ref = 273.15;
  su2double R = 287.058;
  su2double kappa = 0.41;
  su2double A = 17;
  su2double gamma = 1.4;
  su2double Pr_lam = 0.7;
  su2double Pr_turb = 0.9;
  su2double c_p = (gamma*R)/(gamma-1);
  su2double c_v = R/(gamma-1);
  su2double h_wall = c_p * TWall;
  su2double h_bc   = c_p * tExchange;
  su2double tauWall_lam = muExchange * velExchange / h_wm;
  unsigned short nfa = numPoints + 1;

  /*--- Set up vectors ---*/
  
  vector<su2double> T(nfa, tExchange);
  vector<su2double> mu_fa(nfa, muExchange);
  vector<su2double> tmp(nfa, 0.0);
  vector<su2double> u(numPoints, 0.0);
  vector<su2double> lower(numPoints-1,0.0);
  vector<su2double> upper(numPoints-1,0.0);
  vector<su2double> diagonal(numPoints,0.0);
  vector<su2double> rhs(numPoints,0.0);
  
  // Set parameters for control
  bool converged = false;
  unsigned short iter = 0, max_iter = 25;
  su2double tauWall_prev = 0.0, tol = 1e-3,  aux_rhs=0.0;
  //su2double qWall_prev=0.0;
  su2double mut, nu, mu_lam, rho, utau, y_plus, D;
  
  while (converged == false){
    
    iter += 1;
    if (iter == max_iter) converged = true;
    
    tauWall_prev = tauWall;
    //qWall_prev = qWall;
    
    // total viscosity
    // note: rho and mu_lam will be a function of temperature when solving an energy equation
    for(unsigned short i=0; i < nfa; ++i){
      mu_lam = C_1 * pow(T[i]/T_ref, 1.5) * ((T_ref + S)/ (T[i] + S));
      mu_fa[i] = mu_lam;
    }
    
    for(unsigned short i=1; i < nfa; ++i){
      rho = pExchange / (R*T[i]);
      nu = mu_fa[i] / rho;
      utau = sqrt(tauWall / rho);
      y_plus = y_fa[i] * utau / nu;
      D = pow(1.0 - exp((-y_plus)/A),2.0);
      mut = rho * kappa * y_fa[i] * utau * D;
      mu_fa[i] += mut;
    }
    
    // Momentum matrix
    // solution vector is u at y_cv
    lower.assign(numPoints-1,0.0);
    upper.assign(numPoints-1,0.0);
    diagonal.assign(numPoints,0.0);
    rhs.assign(numPoints,0.0);
    
    // top bc
    diagonal[numPoints - 1] = 1.0;
    rhs[numPoints - 1] = velExchange;
    
    // internal cvs
    for (unsigned short i=1; i < (numPoints - 1); ++i){
      upper[i]  =  mu_fa[i + 1] / (y_cv[i + 1] -  y_cv[i] );
      lower[i-1]  = mu_fa[i] / (y_cv[i] -  y_cv[i - 1] );
      diagonal[i] = -1.0 * (upper[i] + lower[i - 1]);
    }
    
    // wall bc
    upper[0] = mu_fa[1]/(y_cv[1] - y_cv[0]);
    diagonal[0] = -1.0 * (upper[0] + mu_fa[0]/(y_cv[0]-y_fa[0]) );
    rhs[0] = 0.0;
    
    // Solve the matrix problem to get the velocity field
    // rhs returned the solution
    //********LAPACK CALL*******
#if (defined(HAVE_MKL) || defined(HAVE_LAPACK)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
    int info, nrhs = 1;

    dgtsv_(&numPoints,&nrhs,lower.data(),diagonal.data(),upper.data(),rhs.data(),&numPoints, &info);
    if (info != 0)
      SU2_MPI::Error("Unsuccessful call to dgtsv_", CURRENT_FUNCTION);
#else
    SU2_MPI::Error("Not compiled with MKL or LAPACK support", CURRENT_FUNCTION);
#endif

    u = rhs;
    
    // update total viscosity
    // note: rho and mu_lam will be a function of temperature when solving an energy equation
    for(unsigned short i=0; i < nfa; ++i){
      mu_lam = C_1 * pow(T[i]/T_ref, 1.5) * ((T_ref + S)/ (T[i] + S));
      mu_fa[i] = mu_lam;
      tmp[i] = mu_lam/Pr_lam;
    }
    // update tauWall
    tauWall = mu_fa[0] * (u[0] - 0.0)/(y_cv[0]-y_fa[0]);
    for(unsigned short i=1; i < nfa; ++i){
      rho = pExchange / (R*T[i]);
      nu = mu_fa[i] / rho;
      utau = sqrt(tauWall / rho);
      y_plus = y_fa[i] * utau / nu;
      D = pow(1.0 - exp((-y_plus)/A),2.0);
      mut = rho * kappa * y_fa[i] * utau * D;
      mu_fa[i] += mut;
      tmp[i] += mut/Pr_turb;
      tmp[i] *= c_p;
    }
    tmp[0] = tmp[0]*c_p;
    
    // Energy matrix
    // solution vector is T at y_cv
    lower.assign(numPoints-1,0.0);
    upper.assign(numPoints-1,0.0);
    diagonal.assign(numPoints,0.0);
    rhs.assign(numPoints,0.0);
    
    // internal cvs
    for (unsigned short i=1; i < (numPoints - 1); ++i){
      upper[i]  =  tmp[i + 1] / (y_cv[i + 1] -  y_cv[i] );
      lower[i-1]  = tmp[i] / (y_cv[i] -  y_cv[i - 1] );
      diagonal[i] = -1.0 * (upper[i] + lower[i - 1]);
    }
    
    // top bc
    diagonal[numPoints - 1] = 1.0;
    
    // wall bc
    upper[0] = tmp[1]/(y_cv[1] - y_cv[0]);
    diagonal[0] = -1.0 * (upper[0] + tmp[0]/(y_cv[0]-y_fa[0]) );
    aux_rhs = tmp[0]/(y_cv[0]-y_fa[0]);
    
    // RHS Energy
    /* compute flux -- (mu + mu_t) * u * du/dy -- */
    
    /* zero flux at the wall */
    tmp[0] = 0. ;
    for (unsigned short i=1; i < numPoints; ++i){
      tmp[i] = 0.5* (mu_fa[i]) * (u[i] + u[i-1]) * (u[i] -u[i-1])/(y_cv[i] -  y_cv[i - 1] )  ;
    }
    for (unsigned short i=0; i < (numPoints - 1); ++i){
      rhs[i] = -tmp[i+1] + tmp[i];
    }
    // END RHS Energy
    
    if (HeatFlux_Prescribed == true){
      /* dT/dy = 0 -> Twall = T[1] */
      h_wall = c_p * T[1];
    }
    
    rhs[0] -= aux_rhs * h_wall;
    rhs[numPoints-1] = h_bc;
    
    // Solve the matrix problem to get the temperature field
    // *******LAPACK CALL********
#if (defined(HAVE_MKL) || defined(HAVE_LAPACK)) && !(defined(CODI_REVERSE_TYPE) || defined(CODI_FORWARD_TYPE))
    dgtsv_(&numPoints,&nrhs,lower.data(),diagonal.data(),upper.data(),rhs.data(),&numPoints, &info);
    if (info != 0)
      SU2_MPI::Error("Unsuccessful call to dgtsv_", CURRENT_FUNCTION);
#else
    SU2_MPI::Error("Not compiled with MKL or LAPACK support", CURRENT_FUNCTION);
#endif

    // Get Temperature from enthalpy
    // Temperature will be at cv or face?
    T[0] = h_wall/c_p;
    T[numPoints] = h_bc/c_p;
    for (unsigned short i=0; i < numPoints-1; i++){
      T[i+1] = 0.5 * (rhs[i] + rhs[i+1])/c_p;
    }
    
    // final update tauWall
    mu_lam = C_1 * pow(T[0]/T_ref, 1.5) * ((T_ref + S)/ (T[0] + S));
    
    // These quantities will be returned.
    tauWall = mu_lam * (u[0] - 0.0)/(y_cv[0]-y_fa[0]);
    qWall = mu_lam * (c_p / Pr_lam) * -(T[1] - T[0]) / (y_cv[0]-y_fa[0]);
    ViscosityWall = mu_lam;
    //kOverCvWall = c_p / c_v * (mu[0]/Pr_lam + muTurb[0]/Pr_turb);
    kOverCvWall = c_p / c_v * (mu_lam/Pr_lam);
    
    // define a norm
    //if (abs( (tauWall - tauWall_prev)/tauWall_lam ) < tol && abs(qWall - qWall_prev) < tol){
    if (abs( (tauWall - tauWall_prev)/tauWall_lam ) < tol){
      converged = true;
//      mu_lam = C_1 * pow(T[0]/T_ref, 1.5) * ((T_ref + S)/ (T[0] + S));
//      rho = pExchange / (R * T[0]);
//      cout << "# " << tExchange << " " << velExchange << " " << muExchange << " " << pExchange << endl;
//      cout << "#Number of iterations: " << iter << endl;
//      cout << "#u tau:  " << sqrt(tauWall/rho) << endl;
//      cout << "#tau wall: " << tauWall << endl;
//      cout << "#qWall: " << qWall << endl;
//      cout << "#Re_tau: " << sqrt(tauWall/rho)*1.0/(mu_lam/rho) << endl;
//      cout << "#Y+  u+ u  T T+" << endl;
//      for(unsigned short i=0; i < numPoints; ++i){
//        cout<< y_cv[i] * sqrt(tauWall/rho) / (mu_lam/rho) << " " << u[i]/sqrt(tauWall/rho)  << " " << u[i] << " " << T[i+1] << " " << (TWall - T[i+1]) * (pExchange / (R * T[i+1])) *  sqrt(tauWall/rho) * c_p / qWall<< endl;
//      }
    }
  }
}

void CWallModelLogLaw::Initialize(const unsigned short *intInfo,
                                const su2double      *doubleInfo){
  
  /* Copy the data from the arguments into the member variables. */
  numPoints      = intInfo[0];
  h_wm           = doubleInfo[0];
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
  const su2double scaleFact = h_wm/coorGridPoints[numPoints-1];
  
  for(unsigned short i=0; i<numPoints; ++i)
    coorGridPoints[i] *= scaleFact;
}

void CWallModelLogLaw::WallShearStressAndHeatFlux(const su2double rhoExchange,
                                                const su2double velExchange,
                                                const su2double muExchange,
                                                const su2double pExchange,
                                                const su2double Wall_HeatFlux,
                                                const bool      HeatFlux_Prescribed,
                                                const su2double Wall_Temperature,
                                                const bool      Temperature_Prescribed,
                                                su2double &tauWall,
                                                su2double &qWall,
                                                su2double &ViscosityWall,
                                                su2double &kOverCvWall) {
  

  // Note: I replaced pExchange by muExchange
  su2double fval, fprime, newton_step, u_tau = 0.0;
  su2double k = 0.38;
  su2double C = 4.1;
  su2double u_tau0 = 0.001;
  su2double nuExchange = muExchange / rhoExchange;
  
  //su2double R = 287.058;
  su2double gamma = 1.4;
  su2double Pr_lam = 0.7;
  //su2double Pr_turb = 0.9;
  //su2double c_p = (gamma*R)/(gamma-1);
  //su2double c_v = R/(gamma-1);
  const su2double factHeatFlux_Lam  = HeatFlux_Prescribed ? su2double(0.0): gamma/Pr_lam;
  unsigned short maxIter = 50, counter = -1;
  
  SU2_MPI::Error("Implementation not finished yet.\nPlease comment this line if you want to test", CURRENT_FUNCTION);
  
  for (unsigned short i=0; i < maxIter; i++){
    counter += 1;
    //y_plus = h_wm * u_tau0 / nuExchange;
    //fval = (velExchange /u_tau0) - ((1.0/k) * log(y_plus) + C);
    //fprime = -velExchange/pow(u_tau0,2) - 1.0/(k*u_tau0);

    fval = velExchange/u_tau0 - (C - 1.0*log(k)/k)*(1.0 - exp(-0.0909090909090909*u_tau*h_wm/nuExchange) - 0.0909090909090909*u_tau0*h_wm*exp(-0.333333333333333*u_tau0*h_wm/nuExchange)/nuExchange) - 1.0*log(k*u_tau*h_wm/nuExchange + 1.0)/k;
    
    //-u_pw/u_tau**2 - (C - 1.0*log(k)/k)*(-0.0909090909090909*y*exp(-0.333333333333333*u_tau*y/nu)/nu + 0.0909090909090909*y*exp(-0.0909090909090909*u_tau*y/nu)/nu + 0.0303030303030303*u_tau*y**2*exp(-0.333333333333333*u_tau*y/nu)/nu**2) - 1.0*y/(nu*(k*u_tau*y/nu + 1.0))
    
    fprime = -velExchange/pow(u_tau0,2.0) - (C - 1.0*log(k)/k)*(-0.0909090909090909*h_wm*exp(-0.333333333333333*u_tau0*h_wm/nuExchange)/nuExchange + 0.0909090909090909*h_wm*exp(-0.0909090909090909*u_tau0*h_wm/nuExchange)/nuExchange + 0.0303030303030303*u_tau0*pow(h_wm,2.0)*exp(-0.333333333333333*u_tau0*h_wm/nuExchange)/pow(nuExchange, 2.0)) - 1.0*h_wm/(nuExchange*(k*u_tau0*h_wm/nuExchange + 1.0));
    
    newton_step = fval/fprime;
    
    u_tau = u_tau0 - newton_step;
    
    if (abs(u_tau - u_tau0) < 1e-6){
      break;
    }
    if (counter == 49) {
      cout << counter << " " <<  u_tau0 << " " << u_tau << endl;
      cout << fprime << " " << fval << " " << newton_step << endl;
      cout << rhoExchange << " " << velExchange << " " << muExchange << endl;
    }
    u_tau0 = u_tau;
  }
  //cout << counter << " " <<  u_tau0 << " " << u_tau << endl;
  //cout << thickness * u_tau0 / nuExchange << " " << velExchange /u_tau << endl;
  tauWall = rhoExchange * pow(u_tau,2.0);
  qWall = 0.0;
  ViscosityWall = muExchange;
  kOverCvWall = ViscosityWall * factHeatFlux_Lam;
  //kOverCvWall = c_p / c_v * (muExchange/Pr_lam);
  
  if(Temperature_Prescribed == true)
  {
    
    
  }
  else if(HeatFlux_Prescribed == true){
  }
  
}
