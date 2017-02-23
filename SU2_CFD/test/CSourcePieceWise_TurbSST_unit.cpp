/*!
 * \file CSourcePieceWise_TurbSST_unit
 * \brief Unit testing for CSourcePieceWise_TurbSST class
 * \author T. A. Oliver
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
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#define TINY 10e-12

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include "config_structure.hpp"
#include "numerics_structure.hpp"

int main() {
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  // Instantiate default configuration
  int return_flag = 0;
  int nDim = 2;
  CConfig testConfig;
  
  su2double paramsSST[10];
  // TODO: Set params to nominal values

  //  sigma_k1 = paramsSST[0];
  //  sigma_k2 = paramsSST[1];
  //  beta_star     = paramsSST[6];
  //  sigma_omega_1 = paramsSST[2];
  //  sigma_omega_2 = paramsSST[3];
  //  beta_1        = paramsSST[4];
  //  beta_2        = paramsSST[5];
  //  alfa_1        = paramsSST[8];
  //  alfa_2        = paramsSST[9];
  //  a1            = paramsSST[7];

  paramsSST[0] = 0.85;
  paramsSST[1] = 1.0;
  paramsSST[2] = 0.5;
  paramsSST[3] = 0.856;
  paramsSST[4] = 0.075;
  paramsSST[5] = 0.0828;
  paramsSST[6] = 0.09;
  paramsSST[7] = 1.0;
  paramsSST[8] = 0.5555555555555555555;
  paramsSST[9] = 0.44;

  // Data to feed to source term calc
  // TODO: Set values here to meaningful values
  su2double StrainMag_i = 1.0;
  // 0 -- k  1 -- omega
  su2double TurbVar_i[2] = {1.0, 1.0};

  su2double** TurbVar_Grad_i = new su2double* [2];
  TurbVar_Grad_i[0] = new su2double[2];
  TurbVar_Grad_i[1] = new su2double[2];

  su2double Volume = 1.0;
  su2double dist_i = 1.0;   // should be greater than 10e-10

  su2double F1_i   = 0.5;
  su2double F2_i   = 1.0;
  su2double CDkw_i = 0.6;



  su2double** PrimVar_Grad_i = new su2double*[4]; // size?
  PrimVar_Grad_i[0] = new su2double[2];
  PrimVar_Grad_i[1] = new su2double[2];
  PrimVar_Grad_i[2] = new su2double[2];
  PrimVar_Grad_i[3] = new su2double[2];

  PrimVar_Grad_i[0][0] = 0.30;
  PrimVar_Grad_i[0][1] = 0.020;
  PrimVar_Grad_i[0][2] = 0.180;
  PrimVar_Grad_i[1][0] = 0.540;
  PrimVar_Grad_i[1][1] = 0.340;
  PrimVar_Grad_i[1][2] = 0.230;
  PrimVar_Grad_i[2][0] = 0.110;
  PrimVar_Grad_i[2][1] = 0.30;
  PrimVar_Grad_i[2][2] = 0.80;
  PrimVar_Grad_i[3][0] = 0.50;
  PrimVar_Grad_i[3][1] = 0.90;
  PrimVar_Grad_i[3][2] = 0.01;


  su2double V_i[9];
  //  comressible flow
  //  0-2 u v w  ndim+2 density  ndim+5 laminar viscosity ndim+6 eddy viscosity
  V_i[0] = 0.230; 
  V_i[1] = 0.20; 
  V_i[2] = 1.0; 
  V_i[3] = 0.340; 
  V_i[4] = 3.0; 
  V_i[5] = 0.0; 
  V_i[6] = 1.20; 
  V_i[7] = 0.230; 
  V_i[8] = 0.0; 
  V_i[9] = 0.770; 

  // Instantiate class to compute SST source terms
  CSourcePieceWise_TurbSST testSST(nDim, 6, paramsSST, &testConfig);

  // Set CNumerics data members used in source term evaluation
  // NB: This mimics what happens in void CTurbSSTSolver::Source_Residual

  // set Primitive variables
  testSST.SetPrimitive(V_i, NULL);

  // set Primitive varaible gradients
  testSST.SetPrimVarGradient(PrimVar_Grad_i, NULL);

  // set Turbulent variables
  testSST.SetTurbVar(TurbVar_i, NULL);

  // set Turbulent variable gradients
  testSST.SetTurbVarGradient(TurbVar_Grad_i, NULL);

  // set Volume
  testSST.SetVolume(Volume);
    
  // set Distance
  testSST.SetDistance(dist_i, 0.0);

  // set blending functions
  testSST.SetF1blending(F1_i,0.0);
  testSST.SetF2blending(F2_i,0.0);

  //testSST.SetVorticity(NULL, NULL);
  
  // set Strain magnitude
  testSST.SetStrainMag(StrainMag_i, 0.0);
  testSST.SetCrossDiff(CDkw_i, 0.0);


  // Evaluate residual and jacobian
  su2double target[2];
  target[0] = target[1] = 0.0;
  su2double res[2];
  su2double** jac_i;
  jac_i    = new su2double* [2];
  jac_i[0] = new su2double[2];
  jac_i[1] = new su2double[2];
 

  // Compute targets
  // step 0 : compute F1_i, F2_i
  // step 1 : compute diverg
  // step 2 : compute pk
  // step 3 : compute pw
  // step 4 : assemble

  // step 0 F1_i, F2_i
  su2double alfa_blended = (1.0 - F1_i) * paramsSST[8] + F1_i * paramsSST[9];
  su2double beta_blended = (1.0 - F1_i) * paramsSST[4] + F1_i * paramsSST[5];

  // step 1  diverg
  su2double diverg = 0.0;
  for (int i=0;i<nDim; i++)
      diverg += PrimVar_Grad_i[i+1][i];

  // step 2 pk
  su2double pk = 0.0;
  pk = V_i[9] * StrainMag_i * StrainMag_i - 2.0/3.0*V_i[nDim+2]*TurbVar_i[0] * diverg;
  pk = min(pk, 20.0 * paramsSST[6] * V_i[nDim+2] * TurbVar_i[0]*TurbVar_i[1]);
  pk = (pk > 0.0) ? pk : 0.0;

  // step 3 pw
  su2double pw = 0.0;
  su2double zeta = 0.0;
  zeta = max(TurbVar_i[1], StrainMag_i * F2_i/paramsSST[7]);
  pw = StrainMag_i * StrainMag_i - 2.0/3.0 * zeta * diverg;
  pw = pw > 0.0 ? pw:0.0;

  // step 4 assemble
  target[0] = (pk - paramsSST[6] * V_i[nDim+2] * TurbVar_i[1] * TurbVar_i[0]) * Volume;
  target[1] = (alfa_blended * V_i[nDim+2] * pw - beta_blended*V_i[nDim+2]*TurbVar_i[1]*TurbVar_i[1] + (1.0 - F1_i)*CDkw_i) * Volume;
  // RES_omega = gamma * omega * P/k - beta * rho * omega^2 + 
  testSST.ComputeResidual(res, jac_i, NULL, &testConfig);
  // TODO: Check that residual and Jacobians are as expected.

  if (fabs(res[0] - target[0]) > TINY || fabs(res[1]-target[1]) > TINY){
      return_flag = 1;
      cout<<"residual 1 = "<<res[0]<<"\t target 1 = "<<target[0]<<endl;
      cout<<"residual 2 = "<<res[1]<<"\t target 2 = "<<target[1]<<endl;
  }
  // Clean up
  delete [] jac_i[1]; delete [] jac_i[0]; delete [] jac_i;
  delete [] PrimVar_Grad_i[3]; 
  delete [] PrimVar_Grad_i[2]; delete [] PrimVar_Grad_i[1];
  delete [] PrimVar_Grad_i[0]; delete [] PrimVar_Grad_i;
  delete [] TurbVar_Grad_i[1]; delete [] TurbVar_Grad_i[0];
  delete [] TurbVar_Grad_i;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return return_flag;

}
