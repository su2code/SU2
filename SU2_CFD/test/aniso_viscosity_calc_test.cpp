/*!
 * \file aniso_viscosity_calc_test.cpp
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

#include "../include/numerics_structure.hpp"

const unsigned short nDim = 3;
const unsigned short nVar = 6;

template<std::size_t M, std::size_t N>
void print_matrix(su2double** A) {
  std::cout << "  [[";
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < N; j++) {
      std::cout << A[i][j];
      if (j<N-1) std::cout << ", ";
    }
    if (i < M-1) std::cout << "]" << std::endl << "   [";
  }
  std::cout << "]]" << std::endl;
}

class TestNumerics : public CAvgGrad_Flow {
 public:
  TestNumerics(CConfig* config) : CAvgGrad_Flow(3, 6, config) {
  }
};

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

  CConfig* test_config = new CConfig();

  TestNumerics numerics(test_config);

  /*XXX: These will have to be allocated and initialized when the
   * ComputeAnisoEddy function actually calls for them.
   */
  su2double** primvar_grad;
  su2double** resolution;

  su2double scalar_eddy_viscosity = 7.0;
  su2double** eddy_viscosity = new su2double*[nDim];
  su2double** correct_vals = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    eddy_viscosity[iDim] = new su2double[nDim];
    correct_vals[iDim] = new su2double[nDim];
  }

  /**-------------------------------------------------------------------------
   * TEST
   *
   *--------------------------------------------------------------------------
   */
  numerics.ComputeAnisoEddyViscosity(primvar_grad,
                                     resolution,
                                     scalar_eddy_viscosity,
                                     eddy_viscosity);

  su2double tol = 10*std::numeric_limits<su2double>::epsilon();
  // This was left general (not in a loop) for easy modification later
  correct_vals[0][0] = 7;
  correct_vals[0][1] = 0;
  correct_vals[0][2] = 0;
  correct_vals[1][0] = 0;
  correct_vals[1][1] = 7;
  correct_vals[1][2] = 0;
  correct_vals[2][0] = 0;
  correct_vals[2][1] = 0;
  correct_vals[2][2] = 7;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      if (std::abs(correct_vals[iDim][jDim] -
                   eddy_viscosity[iDim][jDim]) > tol)
        return_flag = 1;
    }
  }

  if (return_flag == 1) {
    std::cout << "ERROR: The anisotropic viscosity was computed incorrectly!" << std::endl;
    std::cout << "  Expected:" << std::endl;
    print_matrix<3,3>(correct_vals);
    std::cout << "  Found:" << std::endl;
    print_matrix<3,3>(eddy_viscosity);
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
