/*!
 * \file zeta_transform_test.cpp
 * \brief Checks the zeta transformation for the hybrid filtering
 * \author C. Pederson
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

#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/hybrid_RANS_LES_model.hpp"

void print_matrix(vector<vector<su2double> > M) {
  std::cout << "  [[";
  for (unsigned int i = 0; i < M.size(); i++) {
    for (unsigned int j = 0; j < M[i].size(); j++) {
      std::cout << M[i][j];
      if (j<M[i].size()-1) std::cout << ", ";
    }
    if (i < M.size()-1) std::cout << "]" << std::endl << "   [";
  }
  std::cout << "]]" << std::endl;
}

int main() {

  /**-------------------------------------------------------------------------
   * SETUP
   *
   *--------------------------------------------------------------------------
   */
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  int nDim = 3;
  int nConstants = 15;

  CConfig* test_config = new CConfig();

  su2double eigvalues_M[nDim];
  eigvalues_M[0] = 2.0;
  eigvalues_M[1] = 2.0;
  eigvalues_M[2] = 2.0;
  su2double** eigvecs_M = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    eigvecs_M[iDim] = new su2double[nDim];
    for (int jDim = 0; jDim< nDim; jDim++) {
      eigvecs_M[iDim][jDim] = (iDim == jDim);
    }
  }

  /**-------------------------------------------------------------------------
   * TEST
   *
   *--------------------------------------------------------------------------
   */
  CHybrid_Mediator mediator(nDim, test_config);

  vector<vector<su2double> > zeta = mediator.BuildZeta(eigvalues_M, eigvecs_M);

  su2double tol = 1e-8;

  // Check the eigenvectors
  for (int iDim = 0; iDim < nDim; iDim++) {
    for (int jDim = 0; jDim < nDim; jDim++) {
      if (abs(zeta[iDim][jDim] - (iDim == jDim)) > tol) {
        cout.precision(16);
        cout << "ERROR: Zeta was not calculated correctly." << endl;
        cout << "    Calculated:" << endl;
        print_matrix(zeta);
        return_flag = 1;
        break;
      }
    }
    if (return_flag != 0) break;
  }


  /**-------------------------------------------------------------------------
   * TEARDOWN
   *
   *--------------------------------------------------------------------------
   */
  for (int iDim = 0; iDim < nDim; iDim++)
    delete[] eigvecs_M[iDim];
  delete[] eigvecs_M;

  delete test_config;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}
