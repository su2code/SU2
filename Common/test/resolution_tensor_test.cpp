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

inline su2double dot_prod(su2double v[2], su2double w[2]) {
  return v[0]*w[0] + v[1]*w[1];
}

inline su2double dot_prod(vector<su2double> v, vector<su2double> w) {
  return v[0]*w[0] + v[1]*w[1];
}

inline su2double magnitude(su2double v[2]) {
  return std::sqrt(dot_prod(v,v));
}

inline void print_matrix(su2double v[2][2]) {
  std::cout << "[[" << v[0][0] << "," << v[0][1] << "],[";
  std::cout << v[1][0] << "," << v[1][1] << "]]" << std::endl;
}

class Test2DElem : public CPrimalGrid {
 protected:
  su2double **Mij;
  static const unsigned short nFaces = 4;
 public:
  Test2DElem() : CPrimalGrid() {
    unsigned short iDim, iFace;
    nDim = 2;
    /*--- Allocate CG coordinates ---*/
    Coord_CG = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++)
      Coord_CG[iDim] = 0.0;
    Coord_FaceElems_CG = new su2double* [nFaces];
    for (iFace = 0; iFace < nFaces; iFace++) {
      Coord_FaceElems_CG[iFace] = new su2double [nDim];
    }
    Coord_FaceElems_CG[0][0] =  8.0;
    Coord_FaceElems_CG[0][1] =  0.0;
    Coord_FaceElems_CG[1][0] = -8.0;
    Coord_FaceElems_CG[1][1] =  0.0;
    Coord_FaceElems_CG[2][0] =  0.0;
    Coord_FaceElems_CG[2][1] =  2.0;
    Coord_FaceElems_CG[3][0] =  0.0;
    Coord_FaceElems_CG[3][1] = -2.0;
  };
  ~Test2DElem(void) {};

  //---------------------------------------------------------------------------
  void SetResolutionTensor(void) {
    unsigned short iDim, jDim, kDim, lDim;

    // Allocate Mij
    Mij = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Mij[iDim] = new su2double [nDim];
      for (jDim = 0; jDim < nDim; ++jDim) {
        Mij[iDim][jDim] = 0.0;
      }
    }

    su2double vecs[nDim][nDim];

    // First vector
    for (iDim = 0; iDim < nDim; ++iDim) {
      vecs[0][iDim] = Coord_FaceElems_CG[0][iDim] - Coord_FaceElems_CG[1][iDim];
    }

    // Second vector
    for (iDim = 0; iDim < nDim; ++iDim) {
      vecs[1][iDim] = Coord_FaceElems_CG[2][iDim] - Coord_FaceElems_CG[3][iDim];
    }

    // Normalized vectors
    su2double eigvalues[nDim][nDim];
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        eigvalues[iDim][jDim] = 0.0;
      }
      eigvalues[iDim][iDim] = magnitude(vecs[iDim]);
      for (jDim = 0; jDim < nDim; ++jDim) {
        vecs[iDim][jDim] /= eigvalues[iDim][iDim];
      }
    }

    // TODO: Gram-Schmidt Orthogonalization

    // Perform matrix multiplication
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        for (kDim = 0; kDim < nDim; ++kDim) {
          for (lDim = 0; lDim < nDim; ++lDim) {
            Mij[iDim][jDim] += vecs[kDim][iDim]*
                std::sqrt(eigvalues[kDim][lDim])*vecs[lDim][jDim];
          }
        }
      }
    }

  };

  //---------------------------------------------------------------------------
  vector<vector<su2double> > GetResolutionTensor(void) {
    vector<vector<su2double> > output(nDim, vector<su2double>(nDim));
    for (unsigned short iDim = 0; iDim < nDim; ++iDim) {
      for (unsigned short jDim = 0; jDim < nDim; ++jDim) {
        output[iDim][jDim] = Mij[iDim][jDim];
      }
    }
    return output;
  }

  //---------------------------------------------------------------------------
  // Unused functions; Virtual in CPrimalGrid, so they must be implemented
  void SetDomainElement(unsigned long val_domainelement) {};
  unsigned long GetDomainElement(void) {return 0;};
  void Change_Orientation(void) {};
  unsigned short GetVTK_Type(void) {return 0;};
  unsigned short GetRotation_Type(void) {return 0;};
  void SetRotation_Type(unsigned short val_rotation_type) {};
  unsigned short GetnNeighbor_Nodes(unsigned short val_node) {return 0;};
  unsigned short GetnNeighbor_Elements(void) {return 0;};
  unsigned short GetnNodes(void)  {return 0;};
  unsigned short GetnFaces(void)  {return 0;};
  unsigned short GetnNodesFace(unsigned short val_face)  {return 0;};
  unsigned short GetMaxNodesFace(void)  {return 0;};
  unsigned long GetNode(unsigned short val_node)  {return 0;};
  void SetNode(unsigned short val_node, unsigned long val_point) {};
  unsigned short GetFaces(unsigned short val_face, unsigned short val_index)  {return 0;};
  unsigned short GetNeighbor_Nodes(unsigned short val_node, unsigned short val_index)  {return 0;};

};

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
  // Tests
  //---------------------------------------------------------------------------
  // Error tolerances
  const su2double eps = std::numeric_limits<su2double>::epsilon();
  su2double tol = 10*eps;

  Test2DElem* elem = new Test2DElem();
  static int nDim = 2;

  elem->SetResolutionTensor();
  vector<vector<su2double> > Mij = elem->GetResolutionTensor();

  // ---------------------------------------------------------------------------
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

  // ---------------------------------------------------------------------------
  // Check that all columns of resolution tensor are orthogonal
  su2double dp = 0;
  dp += dot_prod(Mij[0],Mij[1]);

  if (std::abs(dp) > tol) {
    std::cout << "ERROR: The resolution tensor for a quadrilateral was not correct."
        << std::endl;
    std::cout << "The column vectors are not orthogonal." << std::endl;
    std::cout << "Sum of dot products: " << dp << std::endl;
    return_flag = 1;
  }

  // ---------------------------------------------------------------------------
  // Check that the values of Mij are correct
  bool entries_correct = true;
  if (Mij[0][0] != 4.0 ) entries_correct = false;
  if (Mij[0][1] != 0.0 ) entries_correct = false;
  if (Mij[1][0] != 0.0 ) entries_correct = false;
  if (Mij[1][1] != 2.0 ) entries_correct = false;

  if (not(entries_correct)) {
    std::cout << "ERROR: The resolution tensor for a quadrilateral was not correct."
        << std::endl;
    std::cout << "The elements of the array were incorrect." << std::endl;
    std::cout << "Array elements: [[" << Mij[0][0] << "," << Mij[0][1] << "],[";
    std::cout << Mij[1][0] << "," << Mij[1][1] << "]]" << std::endl;
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





