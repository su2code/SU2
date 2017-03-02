/*!
 * \file resolution_integration_test.cpp
 * \brief This test checks whether the resolution tensor is correctly set for a grid
 * of quadrilateral cells.
 * \author C. Pederson
 * \version 4.3.0 "Cardinal"
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

#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

inline void CHECK() {
  std::cout << "Got to: " << __FILE__ << " : " << __LINE__ << std::endl;
};

const int nDim = 3;
const int nVar = 1;

class TestNumerics : public CNumerics {
 public:
  TestNumerics(CConfig* config) : CNumerics(nDim, nVar, config) {
  };
};

int main() {

  //---------------------------------------------------------------------------
  // Setup
  //---------------------------------------------------------------------------
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  CConfig* test_config = new CConfig();

  TestNumerics numerics(test_config);
  CHECK();
  su2double** gradprimvar = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    gradprimvar[iDim] = new su2double[nDim+1];
    for (int jDim = 0; jDim < nDim+1; jDim++)
      gradprimvar[iDim][jDim] = 0.0;
  }
  CHECK();
  gradprimvar[1][0] =  1; // dU/dy
  gradprimvar[2][1] = -1; // dV/dx
  CHECK();
  su2double normal[nDim] = {1.0, 0.0, 0.0};
  su2double laminar_viscosity = 0.005;
  CHECK();
  su2double eddy_viscosity = 0.015;
  su2double** aniso_eddy_viscosity = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    aniso_eddy_viscosity[iDim] = new su2double[nDim+1];
    for (int jDim = 0; jDim < nDim+1; jDim++)
      aniso_eddy_viscosity[iDim][jDim] = 0.0;
  }

  //---------------------------------------------------------------------------
  // Test
  //---------------------------------------------------------------------------
  CHECK();
  numerics.GetViscousArtCompProjFlux(gradprimvar, normal, laminar_viscosity,
                                     eddy_viscosity);


  //---------------------------------------------------------------------------
  // Teardown
  //---------------------------------------------------------------------------
  delete test_config;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





