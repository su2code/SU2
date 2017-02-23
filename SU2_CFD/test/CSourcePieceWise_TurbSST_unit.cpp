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
  CConfig testConfig;

  su2double paramsSST[10];
  // TODO: Set params to nominal values


  // Data to feed to source term calc
  // TODO: Set values here to meaningful values
  su2double StrainMag_i = 0.0;
  su2double TurbVar_i[2] = {0.0, 0.0};

  su2double** TurbVar_Grad_i = new su2double* [2];
  TurbVar_Grad_i[0] = new su2double[2];
  TurbVar_Grad_i[1] = new su2double[2];

  su2double Volume = 0.0;
  su2double dist_i = 0.0;

  su2double F1_i = 0.0;
  su2double F2_i = 0.0;
  su2double CDkw_i = 0.0;

  su2double** PrimVar_Grad_i = new su2double* [3]; // size?
  PrimVar_Grad_i[0] = new su2double[2];
  PrimVar_Grad_i[1] = new su2double[2];
  PrimVar_Grad_i[2] = new su2double[2];

  su2double V_i[9];

  // Instantiate class to compute SST source terms
  CSourcePieceWise_TurbSST testSST(2, 6, paramsSST, &testConfig);

  // Set CNumerics data members used in source term evaluation
  // NB: This mimics what happens in void CTurbSSTSolver::Source_Residual

  testSST.SetPrimitive(V_i, NULL);
  testSST.SetPrimVarGradient(PrimVar_Grad_i, NULL);
  testSST.SetTurbVar(TurbVar_i, NULL);
  testSST.SetTurbVarGradient(TurbVar_Grad_i, NULL);
  testSST.SetVolume(Volume);
  testSST.SetDistance(dist_i, 0.0);
  testSST.SetF1blending(F1_i,0.0);
  testSST.SetF2blending(F2_i,0.0);
  //testSST.SetVorticity(NULL, NULL);
  testSST.SetStrainMag(StrainMag_i, 0.0);
  testSST.SetCrossDiff(CDkw_i,0.0);


  // Evaluate residual and jacobian
  su2double res[2];
  su2double** jac_i;
  jac_i = new su2double* [2];
  jac_i[0] = new su2double[2];
  jac_i[1] = new su2double[2];

  testSST.ComputeResidual(res, jac_i, NULL, &testConfig);

  // TODO: Check that residual and Jacobians are as expected.


  // Clean up
  delete [] jac_i[1]; delete [] jac_i[0]; delete [] jac_i;
  delete [] PrimVar_Grad_i[2]; delete [] PrimVar_Grad_i[1];
  delete [] PrimVar_Grad_i[0]; delete [] PrimVar_Grad_i;
  delete [] TurbVar_Grad_i[1]; delete [] TurbVar_Grad_i[0];
  delete [] TurbVar_Grad_i;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0;
}
