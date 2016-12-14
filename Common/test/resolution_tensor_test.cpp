/*!
 * \file resolution_tensor_test.cpp
 * \brief 
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

#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"

double dot_prod(double v[3], double w[3]) {
  double output = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

void ReduceVectors(double (&vecs_in)[6][3], double (&vecs_out)[3][3]) {
  // Use the first vector
  for (int i = 0; i<3; ++i) vecs_out[0][i] = vecs_in[0][i];

  // Find the vectors that are the closest to being orthogonal
  double norm_0 = 0.0;
  norm_0 = dot_prod(vecs_in[0], vecs_in[0]);
  double min_dp = norm_0;
  double second_dp = norm_0;
  int min_dp_index, second_dp_index;
  double dot_product;
  for (int i = 1; i<6; ++i) {
    dot_product = dot_prod(vecs_in[0], vecs_in[i]);
    if (std::abs(dot_product) < min_dp) {
      min_dp = dot_product;
      min_dp_index = i;
    } else if (std::abs(dot_product) < second_dp) {
      second_dp = dot_product;
      second_dp_index = i;
    }
  }

  // Construct the second vector
  double factor0 = min_dp/norm_0;
  std::cout << factor0 << std::endl;
  for (int i=0; i<3; ++i) vecs_out[1][i] = vecs_in[min_dp_index][i] - factor0*vecs_in[0][i];
  double norm_1 = dot_prod(vecs_out[1], vecs_out[1]);

  // Construct the third vector
  factor0 = second_dp/norm_0;
  double factor1 = dot_prod(vecs_in[second_dp_index], vecs_out[1])/norm_1;
  for (int i = 0; i < 3; ++i) {
    vecs_out[2][i] = vecs_in[second_dp_index][i];
    vecs_out[2][i] -= factor0*vecs_out[0][i];
    vecs_out[2][i] -= factor1*vecs_out[1][i];
  }

}

int main() {

  //---------------------------------------------------------------------------
  // Setup
  //---------------------------------------------------------------------------
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  CConfig* test_config;
  test_config = new CConfig();

  //---------------------------------------------------------------------------
  // Test
  //---------------------------------------------------------------------------

  double vecs_in[6][3] = {{ 1, 1, -1},
      { 1,  0,  2},
      { 2, -2,  3},
      { 1, -1, -1},
      { 1,  0, -2},
      {-2,  -2, 3}};


  double vecs_out[3][3] = {{0, 0, 0}, {0,0,0}, {0,0,0}};

  ReduceVectors(vecs_in, vecs_out);

  const double eps = std::numeric_limits<double>::epsilon();
  su2double tol = 10*eps;

  // Check that the return vectors are nonzero.
  su2double sum;
  for (int i=0; i<3; ++i) {
    sum = vecs_out[i][0] + vecs_out[i][1] + vecs_out[i][2];
    if (not(std::abs(sum) > tol))  {
      std::cout << "ERROR: Calculation of vectors for resolution tensor failed." << std::endl;
      std::cout << "Norm of computed vector #" << i+1 << " was less than "
          << eps << std::endl;
      std::cout << "Vector: [" << vecs_out[i][0] << " , "
          << vecs_out[i][1] << " , " << vecs_out[i][2] << "]" << std::endl;
      return_flag = 1;
      break;
    }
  }

  // Check that the return values are orthogonal.
  su2double dot_product;
  for (int i=1; i<3; ++i) {
    dot_product = dot_prod(vecs_out[0], vecs_out[i]);
    if (std::abs(dot_product) > eps) {
      std::cout << "ERROR: Calculation of vectors for resolution tensor failed." << std::endl;
      std::cout << "Computed vector were non-orthogonal." << std::endl;
      std::cout << "Vec 1:  [" << vecs_out[0][0] << " , " << vecs_out[0][1] << " , "
          << vecs_out[0][2] << "]" << std::endl;
      std::cout << "Vec 2: [" << vecs_out[1][0] << " , " << vecs_out[1][1] << " , "
          << vecs_out[1][2] << "]" << std::endl;
      std::cout << "Vec 3: [" << vecs_out[2][0] << " , " << vecs_out[2][1] << " , "
          << vecs_out[2][2] << "]" << std::endl;
      return_flag = 1;
      break;
    }
  }

  // Check that the first vector points in the same direction as the original
  dot_product = dot_prod(vecs_out[0], vecs_in[0]);
  double norm = dot_prod(vecs_in[0], vecs_in[0]);
  if (std::abs(dot_product - norm) > eps) {
    std::cout << "ERROR: Calculation of vectors for resolution tensor failed." << std::endl;
    std::cout << "First vector (the reference vector) changed directions." << std::endl;
    std::cout << "In:  [" << vecs_in[0][0] << " , " << vecs_in[0][1] << " , "
        << vecs_in[0][2] << "]" << std::endl;
    std::cout << "Out: [" << vecs_out[0][0] << " , " << vecs_out[0][1] << " , "
        << vecs_out[0][2] << "]" << std::endl;
    return_flag = 1;
  }

  //---------------------------------------------------------------------------
  // Teardown
  //---------------------------------------------------------------------------
  delete test_config;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





