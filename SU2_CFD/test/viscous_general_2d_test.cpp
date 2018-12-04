/*!
 * \file SU2_CFD.cpp
 * \brief Main file of the SU2 Computational Fluid Dynamics code
 * \author F. Palacios, T. Economon
 * \version 6.1.0 "Falcon"
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
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
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

#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <string.h>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "../include/numerics_structure.hpp"

unsigned short error_count = 0;

void BOOST_CHECK_CLOSE_FRACTION(const su2double& u, const su2double& v,
                                const su2double& tol) {
  if (u != 0 && v != 0) {
    const su2double rel_error_u = std::abs((u - v)/u);
    const su2double rel_error_v = std::abs((u - v)/v);
    if (rel_error_u > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Relative difference exceeds tolerance ";
      std::cout << "[" << rel_error_u << " > " << tol << "]\n";
      error_count += 1;
    } else if (rel_error_v > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Relative difference exceeds tolerance ";
      std::cout << "[" << rel_error_v << " > " << tol << "]\n";
      error_count += 1;
    }
  } else {
    const su2double abs_error = std::abs(u - v);
    if (abs_error > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Absolute difference exceeds tolerance ";
      std::cout << "[" << abs_error << " > " << tol << "]\n";
      error_count += 1;
    }
  }
}

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
void WriteCfgFile() {

  std::ofstream cfg_file;

  cfg_file.open("test.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "TIME_DISCRE_FLOW= EULER_IMPLICIT" << std::endl;

  cfg_file.close();

}

int main(int argc, char *argv[]) {

  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
  const unsigned int BUFSIZE = 3000000;        /*!< \brief MPI buffer. */
  int  buffsize;
  char *buffptr;
  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
#endif

  error_count = 0;

  /*---
   * Test Variables
   * ---*/

  WriteCfgFile();

  const unsigned short nDim = 2;
  const unsigned short nVar = nDim+2;
  const unsigned short nPrimVar = nDim+9;
  const unsigned short nSecVar = 4;

  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig("test.cfg", SU2_CFD, iZone, nZone, 2, VERB_NONE);
  config->SetGas_ConstantND(287.058);

  const su2double distance = 1;
  const su2double area = 3;

  su2double coord_i[nDim], coord_j[nDim];
  su2double normal[nDim];
  su2double** primvar_grad_i, **primvar_grad_j;
  su2double primvar_i[nPrimVar];
  su2double primvar_j[nPrimVar];
  su2double secvar_i[nSecVar];
  su2double secvar_j[nSecVar];
  su2double** Jacobian_i, **Jacobian_j;
  su2double* residual_i;
  su2double** expected_jacobian_i, **expected_jacobian_j;

  /*--- Inputs ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    coord_i[iDim] = 0;
    coord_j[iDim] = 0;
  }
  coord_j[0] += distance;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    normal[iDim] = coord_j[iDim]/distance*area;
  }

  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_i[iVar] = 0.0; primvar_j[iVar] = 0.0;
  }
  for (unsigned short iVar = 0; iVar < nSecVar; iVar++) {
    secvar_i[iVar] = 0.0; secvar_j[iVar] = 0.0;
  }

  primvar_grad_i = new su2double*[nPrimVar];
  primvar_grad_j = new su2double*[nPrimVar];
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_grad_i[iVar] = new su2double[nDim];
    primvar_grad_j[iVar] = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_i[iVar][iDim] = 0.0;
    }
  }

  /*--- Outputs ---*/

  Jacobian_i = new su2double*[nVar];
  Jacobian_j = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double[nVar];
    Jacobian_j[iVar] = new su2double[nVar];
  }

  residual_i = new su2double[nVar];

  expected_jacobian_i = new su2double*[nVar];
  expected_jacobian_j = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    expected_jacobian_i[iVar] = new su2double[nVar];
    expected_jacobian_j[iVar] = new su2double[nVar];
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      expected_jacobian_i[iVar][jVar] = 0;
      expected_jacobian_j[iVar][jVar] = 0;
    }
  }

  /*---
   * SETUP
   * ---*/

  CNumerics* numerics = new CGeneralAvgGrad_Flow(nDim, nVar, false, config);

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 2.0; // density
  primvar_i[nDim+5] = 1.0; // laminar viscosity
  primvar_i[nDim+6] = 5.0; // turbulent viscosity
  for (unsigned short iVar = 1; iVar < nDim+1; iVar++) {
    primvar_i[iVar] = iVar*iVar; // Velocities
  }
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  secvar_i[2] = 2.0;
  secvar_j[3] = 5.0;
  for (unsigned short iVar = 0; iVar < nSecVar; iVar++) {
    secvar_j[iVar] = secvar_i[iVar];
  }

  primvar_grad_i[1][0] =  1.0; // du/dx
  primvar_grad_i[2][1] =  2.0; // dv/dy
  primvar_grad_i[1][1] =  1.0; // du/dy
  primvar_grad_i[2][0] =  1.0; // dv/dx
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_j[iVar][iDim] = primvar_grad_i[iVar][iDim];
    }
  }

  const su2double tke = 3; // 3 cancels out 3 in denominator

  /*---
   * TEST
   * ---*/

  numerics->SetCoord(coord_i, coord_j);
  numerics->SetNormal(normal);
  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetSecondary(secvar_i, secvar_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  numerics->SetTauWall(0, 0);
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar];
  expected_residual[0] = 0;
  expected_residual[1] = -12;
  expected_residual[2] = 36;
  expected_residual[3] = 132;
  expected_jacobian_i[0][0] = 0;
  expected_jacobian_i[0][1] = 0;
  expected_jacobian_i[0][2] = 0;
  expected_jacobian_i[0][3] = 0;
  expected_jacobian_i[1][0] = 12;
  expected_jacobian_i[1][1] = -12;
  expected_jacobian_i[1][2] = -0;
  expected_jacobian_i[1][3] = 0;
  expected_jacobian_i[2][0] = 36;
  expected_jacobian_i[2][1] = -0;
  expected_jacobian_i[2][2] = -9;
  expected_jacobian_i[2][3] = 0;
  expected_jacobian_i[3][0] = 123;
  expected_jacobian_i[3][1] = -15;
  expected_jacobian_i[3][2] = -27;
  expected_jacobian_i[3][3] = -0;
  expected_jacobian_j[0][0] = -0;
  expected_jacobian_j[0][1] = -0;
  expected_jacobian_j[0][2] = -0;
  expected_jacobian_j[0][3] = -0;
  expected_jacobian_j[1][0] = -12;
  expected_jacobian_j[1][1] = 12;
  expected_jacobian_j[1][2] = 0;
  expected_jacobian_j[1][3] = -0;
  expected_jacobian_j[2][0] = -36;
  expected_jacobian_j[2][1] = 0;
  expected_jacobian_j[2][2] = 9;
  expected_jacobian_j[2][3] = -0;
  expected_jacobian_j[3][0] = -189;
  expected_jacobian_j[3][1] = 9;
  expected_jacobian_j[3][2] = 45;
  expected_jacobian_j[3][3] = 0;

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    BOOST_CHECK_CLOSE_FRACTION(expected_residual[iVar], residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_CLOSE_FRACTION(expected_jacobian_i[iVar][jVar], Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_CLOSE_FRACTION(expected_jacobian_j[iVar][jVar], Jacobian_j[iVar][jVar], tolerance);
    }
  }

  if (error_count > 0) {
    cout << "\n\n";
    cout << "If the calculated values are actually the correct values, then update the test with these calculated values:\n\n";
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      cout << "expected_residual[" << iVar << "] = ";
      cout << residual_i[iVar] << ";\n";
    }
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        cout << "expected_jacobian_i[" << iVar << "][" << jVar << "] = ";
        cout << Jacobian_i[iVar][jVar] << ";\n";
      }
    }
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        cout << "expected_jacobian_j[" << iVar << "][" << jVar << "] = ";
        cout << Jacobian_j[iVar][jVar] << ";\n";
      }
    }
  }

  /*---
   * TEARDOWN
   * ---*/

  delete numerics;

  delete config;
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    delete[] primvar_grad_i[iVar];
    delete[] primvar_grad_j[iVar];
  }
  delete[] primvar_grad_i;
  delete[] primvar_grad_j;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete[] Jacobian_i[iVar];
    delete[] Jacobian_j[iVar];
  }
  delete[] Jacobian_i;
  delete[] Jacobian_j;
  delete[] residual_i;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete[] expected_jacobian_i[iVar];
    delete[] expected_jacobian_j[iVar];
  }
  delete[] expected_jacobian_i;
  delete[] expected_jacobian_j;

  /*--- Finalize MPI parallelization ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif

  if (error_count == 0)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;

}
