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

#define BOOST_TEST_MODULE ViscousModelSplit
#include "MPI_global_fixture.hpp"


#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../../Common/include/config_structure.hpp"
#include "../include/numerics_structure.hpp"

const unsigned short nDim = 3;
const unsigned short nVar = nDim+2;
// We don't need all the primitive variables
const unsigned short nPrimVar = nDim+3;

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

#ifdef BUILD_TESTS

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_AUTO_TEST_CASE(RANSStressTensorMatchesModelSplit) {

  /*---
   * SETUP
   * ---*/

  CConfig* test_config = new CConfig();

  CAvgGrad_Base numerics(nDim, nVar, nVar+3, false, test_config);

  su2double primvar[nPrimVar];
  primvar[1] = 1.0;  // u
  primvar[2] = 2.0;  // v
  primvar[3] = 3.0;  // w
  primvar[nDim+2] = 5.0; // Density
  su2double** gradprimvar = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    gradprimvar[iVar] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      gradprimvar[iVar][jDim] = 0.0;
  }

  gradprimvar[1][0] = 1; // dU/dx
  gradprimvar[2][1] = 2; // dV/dy
  gradprimvar[2][1] = 1; // dV/dy
  gradprimvar[3][2] = 1; // dW/dz
  gradprimvar[3][1] = 4; // dW/dy
  const su2double turb_ke = 3.0;
  const su2double laminar_viscosity = 7.0;
  const su2double eddy_viscosity = 11.0;
  su2double **rans_tau = new su2double*[nDim];
  su2double **model_split_tau = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    rans_tau[iDim] = new su2double[nDim];
    model_split_tau[iDim] = new su2double[nDim];
  }

  /*---
   * TEST
   * ---*/

  numerics.GetStressTensor(primvar, gradprimvar, turb_ke,
                           laminar_viscosity, eddy_viscosity);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      rans_tau[iDim][jDim] = numerics.GetStressTensor(iDim, jDim);
    }
  }

  numerics.GetLaminarStressTensor(gradprimvar, laminar_viscosity);
  numerics.AddTauSGS(primvar, gradprimvar, turb_ke, eddy_viscosity);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      model_split_tau[iDim][jDim] = numerics.GetStressTensor(iDim, jDim);
    }
  }

  // Compare tau
  const su2double tolerance = 2*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_CLOSE_FRACTION(rans_tau[iDim][jDim],
                                 model_split_tau[iDim][jDim], tolerance);
    }
  }


  /*---
   * TEARDOWN
   * ---*/

  delete test_config;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] gradprimvar[iVar];
  }
  delete [] gradprimvar;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] model_split_tau[iDim];
    delete [] rans_tau[iDim];
  }
  delete [] model_split_tau;
  delete [] rans_tau;

}

BOOST_AUTO_TEST_CASE(RansStressMatchesIsotropicEddyViscosityStress) {

  /*---
   * SETUP
   * ---*/

  CConfig* test_config = new CConfig();

  CAvgGrad_Base numerics(nDim, nVar, nVar+3, false, test_config);

  su2double primvar[nPrimVar];
  primvar[1] = 1.0;  // u
  primvar[2] = 2.0;  // v
  primvar[3] = 3.0;  // w
  primvar[nDim+2] = 5.0; // Density
  su2double** gradprimvar = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    gradprimvar[iVar] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      gradprimvar[iVar][jDim] = 0.0;
  }

  gradprimvar[1][0] = 1; // dU/dx
  gradprimvar[2][1] = 2; // dV/dy
  gradprimvar[2][1] = 1; // dV/dy
  gradprimvar[3][2] = 1; // dW/dz
  gradprimvar[3][1] = 4; // dW/dy
  const su2double turb_ke = 3.0;
  const su2double laminar_viscosity = 7.0;
  const su2double eddy_viscosity = 11.0;
  su2double **rans_tau = new su2double*[nDim];
  su2double **model_split_tau = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    rans_tau[iDim] = new su2double[nDim];
    model_split_tau[iDim] = new su2double[nDim];
  }

  su2double **aniso_eddy_viscosity = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    aniso_eddy_viscosity[iDim] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      // FIXME: Does this need to be rescaled by a factor of sqrt(3) or something?
      if (iDim == jDim) gradprimvar[iDim][jDim] = eddy_viscosity;
      else gradprimvar[iDim][jDim] = 0.0;
    }
  }

  /*---
   * TEST
   * ---*/

  numerics.GetStressTensor(primvar, gradprimvar, turb_ke,
                           laminar_viscosity, eddy_viscosity);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      rans_tau[iDim][jDim] = numerics.GetStressTensor(iDim, jDim);
    }
  }

  numerics.GetLaminarStressTensor(gradprimvar, laminar_viscosity);
  numerics.AddTauSGET(gradprimvar, aniso_eddy_viscosity);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      model_split_tau[iDim][jDim] = numerics.GetStressTensor(iDim, jDim);
    }
  }

  // Compare tau
  const su2double tolerance = 2*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_CLOSE_FRACTION(rans_tau[iDim][jDim],
                                 model_split_tau[iDim][jDim], tolerance);
    }
  }


  /*---
   * TEARDOWN
   * ---*/

  delete test_config;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] gradprimvar[iVar];
  }
  delete [] gradprimvar;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] model_split_tau[iDim];
    delete [] rans_tau[iDim];
    delete [] aniso_eddy_viscosity[iDim];
  }
  delete [] model_split_tau;
  delete [] rans_tau;
  delete [] aniso_eddy_viscosity;

}

#endif
