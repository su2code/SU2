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

#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <string.h>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "../include/numerics_structure.hpp"

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

  int opt;
  unsigned long iterations = 1E6;
  for (int i = 1; i < argc; i++) {
    if (not(strcmp(argv[i], "-i"))) {
      iterations = atoi(argv[i+1]);
      i++;
    } else {
      std::cout << "Unrecognized option: " << argv[i] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "Number of iterations: " << iterations << std::endl;

  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
  const unsigned int BUFSIZE = 3000000;        /*!< \brief MPI buffer. */
  int  buffsize;
  char *buffptr;
  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
#endif


  /*---
   * Test Variables
   * ---*/

  WriteCfgFile();

  const unsigned short nDim = 3;
  const unsigned short nVar = nDim+2;
  const unsigned short nPrimVar = nDim+9;

  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig("test.cfg", SU2_CFD, iZone, nZone, 2, VERB_NONE);

  const su2double distance = 1;
  const su2double area = 3;

  su2double coord_i[nDim], coord_j[nDim];
  su2double normal[nDim];
  su2double** primvar_grad_i, **primvar_grad_j;
  su2double primvar_i[nPrimVar];
  su2double primvar_j[nPrimVar];
  su2double** Jacobian_i, **Jacobian_j;
  su2double* residual_i;

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

  /*---
   * SETUP
   * ---*/

  CNumerics* numerics = new CAvgGrad_Flow(nDim, nVar, false, config);

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 1.0; // density
  primvar_i[nDim+5] = 1.0; // laminar viscosity
  primvar_i[nDim+6] = 1.0; // turbulent viscosity
  for (unsigned short iVar = 1; iVar < nDim+1; iVar++) {
    primvar_i[iVar] = iVar; // Velocities
  }
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  primvar_grad_i[1][0] =  1.0; // du/dx
  primvar_grad_i[2][1] =  2.0; // dv/dy
  primvar_grad_i[3][2] =  3.0; // dw/dz
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
  numerics->SetSecondary(NULL, NULL);
  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  clock_t begin = clock();
  for (unsigned long i = 0; i < iterations; i++)
    numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  std::cout << "Elapsed time: " << elapsed_secs << std::endl;

  /*---
   * TEARDOWN
   * ---*/

  delete numerics;

  delete config;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
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

  /*--- Finalize MPI parallelization ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif

  return EXIT_SUCCESS;

}
