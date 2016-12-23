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

class TestGeometry : public CPhysicalGeometry {
 protected:
 public:
  TestGeometry() {
    nElem = 4;
    for (int iElem=0; iElem<nElem; iElem++) {
      elem[iElem] = new CQuadrilateral(0,0,0,0,2);
    }
  }
  void SetResolutionTensor(void) {
  };
};

inline su2double dot_prod(su2double v[3], su2double w[3]) {
  su2double output = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

void ReduceVectors(su2double (&vecs_in)[6][3], su2double (&vecs_out)[3][3]) {
  // Use the first vector
  for (int i = 0; i<3; ++i) vecs_out[0][i] = vecs_in[0][i];

  // Find the vectors that are the closest to being orthogonal
  // FIXME: The second vector must be close to orthogonal with the third
  su2double norm_0 = 0.0;
  norm_0 = dot_prod(vecs_in[0], vecs_in[0]);
  su2double min_dp = norm_0;
  su2double second_dp = norm_0;
  int min_dp_index, second_dp_index;
  su2double dot_product;
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
  su2double factor0 = min_dp/norm_0;
  std::cout << factor0 << std::endl;
  for (int i=0; i<3; ++i) vecs_out[1][i] = vecs_in[min_dp_index][i] - factor0*vecs_in[0][i];
  su2double norm_1 = dot_prod(vecs_out[1], vecs_out[1]);

  // Construct the third vector
  factor0 = second_dp/norm_0;
  su2double factor1 = dot_prod(vecs_in[second_dp_index], vecs_out[1])/norm_1;
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
  int nDim = 2;
  CQuadrilateral* elem = new CQuadrilateral(1,1,1,1,nDim);

  const su2double eps = std::numeric_limits<su2double>::epsilon();
  su2double tol = 10*eps;

  // FIXME: The quadirlateral object does not yet have face CGs
  elem->SetResolutionTensor();
  vector<vector<su2double> > Mij = elem->GetResolutionTensor();

  // Check that the resolution tensor is nonzero
  su2double sum = 0.0;
  for (int i=0; i<nDim; ++i) {
    for (int j=0; j<nDim; ++j) {
      sum += std::abs(Mij[i][j]);
    }
  }
  if (sum < eps) {
    std::cout << "ERROR: The resolution tensor for a quadrilateral was not correct."
        << std::endl;
    std::cout << "The sum was within machine precision of zero." << std::endl;
    std::cout << "Sum: " << sum << std::endl;
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





