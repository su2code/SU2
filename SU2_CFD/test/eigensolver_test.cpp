/*!
 * \file eigensolver_test.cpp
 * \brief 
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

  CConfig* test_config = new CConfig();

  // Allocate and create input matrix
  su2double** matrix = new su2double*[nDim];
  for (int iDim=0; iDim<nDim; iDim++) {
    matrix[iDim] = new su2double[nDim];
    for (int jDim=0; jDim<nDim; jDim++) {
      matrix[iDim][jDim] = su2double(iDim == jDim)*(iDim+1);
    }
  }

  // Construct the test and correct values
  vector<su2double> eigvalues(nDim), expected_values(nDim);
  vector<vector<su2double> > eigvectors(nDim, vector<su2double>(nDim)),
                             expected_vectors(nDim, vector<su2double>(nDim));

  // Fill in the expected results
  for (int iDim=0; iDim<nDim; iDim++) {
    expected_values[iDim] = (iDim+1);
    for (int jDim=0; jDim<nDim; jDim++) {
      expected_vectors[iDim][jDim] = su2double(iDim == jDim);
    }
  }

  /**-------------------------------------------------------------------------
   * TEST
   *
   *--------------------------------------------------------------------------
   */
  CHybrid_Mediator mediator(nDim, test_config);
  mediator.SolveEigen(matrix, eigvalues, eigvectors);

  su2double tol = 1e-8;

  // Check the eigenvectors
  for (int iDim = 0; iDim < nDim; iDim++) {
    for (int jDim = 0; jDim < nDim; jDim++) {
      if (abs(expected_vectors[iDim][jDim] - eigvectors[iDim][jDim]) > tol) {
        cout << "ERROR: Eigenvectors were not calculated correctly." << endl;
        cout << "    Calculated:" << endl;
        print_matrix(eigvectors);
        cout << "    Expected:" << endl;
        print_matrix(expected_vectors);
        return_flag = 1;
        break;
      }
    }
    if (return_flag != 0) break;
  }

  // Check the eigenvalues
  for (int iDim = 0; iDim < nDim; iDim++) {
      if (abs(expected_values[iDim] - eigvalues[iDim]) > tol) {
        cout << "ERROR: Eigenvalues were not calculated correctly." << endl;
        cout << "    Calculated:" << endl;
        cout << "    [";
        cout << eigvalues[0] << ", ";
        cout << eigvalues[1] << ", ";
        cout << eigvalues[2] << "]" << endl;
        cout << "    Expected:" << endl;
        cout << "    [";
        cout << expected_values[0] << ", ";
        cout << expected_values[1] << ", ";
        cout << expected_values[2] << "]" << endl;
        return_flag = 1;
        break;
    }
    if (return_flag != 0) break;
  }

  /**-------------------------------------------------------------------------
   * TEARDOWN
   *
   *--------------------------------------------------------------------------
   */
   delete test_config;
    


#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}
