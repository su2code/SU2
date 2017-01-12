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

template<std::size_t nDim>
inline su2double dot_prod(su2double v[nDim], su2double w[nDim]) {
  su2double dot_product = 0.0;
  for (unsigned int iDim = 0; iDim < nDim; ++iDim) {
    dot_product += v[iDim]*w[iDim];
  }
  return dot_product;
}

template<std::size_t nDim>
inline su2double dot_prod(vector<su2double> v, vector<su2double> w) {
  su2double dot_product = 0.0;
  for (unsigned int iDim = 0; iDim < nDim; ++iDim) {
    dot_product += v[iDim]*w[iDim];
  }
  return dot_product;
}

template<std::size_t nDim>
inline su2double magnitude(su2double v[nDim]) {
  return std::sqrt(dot_prod<nDim>(v,v));
}

template<std::size_t nDim>
inline void print_matrix(su2double v[nDim][nDim]) {
  std::cout << "[[";
  for (unsigned int jDim = 0; jDim < nDim; ++jDim) {
    for (unsigned int iDim = 0; iDim < nDim; ++iDim) {
      std::cout << v[iDim][jDim];
      if (iDim != nDim-1) std::cout << ",";
    }
    if (jDim != nDim-1)  {
      std::cout << "],[";
    } else {
      std::cout << "]]" << std::endl;
    }
  }
}

template<std::size_t nDim>
void GramSchmidt(su2double w[nDim][nDim], su2double v[nDim][nDim]) {
  unsigned short iDim, jDim;

  // Set the first basis vector to the first input vector
  for (iDim = 0; iDim < nDim; ++iDim) {
    v[0][iDim] = w[0][iDim]/magnitude<nDim>(w[0]);
  }

  // Compute the next orthogonal vector
  for (iDim = 0; iDim < nDim; ++iDim) {
    v[1][iDim] = w[1][iDim] - dot_prod<nDim>(w[1],v[0])*v[0][iDim];
  }
  // Normalize the 2nd vector
  su2double mag = magnitude<nDim>(v[1]);
  for (iDim = 0; iDim < nDim; ++iDim) {
    v[1][iDim] /= mag;
  }

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
    Coord_FaceElems_CG[0][0] =  2.0;
    Coord_FaceElems_CG[0][1] =  0.0;
    Coord_FaceElems_CG[2][0] = -2.0;
    Coord_FaceElems_CG[2][1] =  0.0;
    Coord_FaceElems_CG[1][0] =  0.0;
    Coord_FaceElems_CG[1][1] = -1.0;
    Coord_FaceElems_CG[3][0] =  0.0;
    Coord_FaceElems_CG[3][1] =  1.0;
  };
  ~Test2DElem(void) {};

  //---------------------------------------------------------------------------
  void SetResolutionTensor(void) {
    unsigned short iDim, jDim, kDim, lDim;
    unsigned short iFace;
    unsigned short* paired_faces;
    // paired_faces is used to sort the faces into matching pairs.
    // The code will look for pairs of faces that are mostly opposite, then
    // sort them so that the face indices in paired_faces[0] and paired_faces[1]
    // match, then paired_faces[2] and paired_faces[3] match, etc.


    // Allocate Mij
    Mij = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++) {
      Mij[iDim] = new su2double [nDim];
      for (jDim = 0; jDim < nDim; ++jDim) {
        Mij[iDim][jDim] = 0.0;
      }
    }

    paired_faces = new unsigned short [nFaces];

    su2double eigvecs[nDim][nDim];

    // Create cell center to face vectors
    su2double center2face[nFaces][nDim];
    for (iFace = 0; iFace < nFaces; ++iFace) {
      for (iDim = 0; iDim < nDim; ++iDim) {
        center2face[iFace][iDim] = Coord_FaceElems_CG[iFace][iDim] - Coord_CG[iDim];
      }
    }

    // First vector
    paired_faces[0] = 0; // Choose vector 1 as our first vector to pair up
    // Find vector mostly parallel to first
    su2double min_dp = 1.0;
    su2double current_dp;
    for (iFace = 1; iFace < nFaces; ++iFace) {
      // NOTE: You cannot use nDim in the template functions for the array size;
      // nDim is not const
      current_dp = dot_prod<2>(center2face[0],center2face[iFace])
          /(magnitude<2>(center2face[0])*magnitude<2>(center2face[iFace]));
      if (current_dp < min_dp) {
        min_dp = current_dp;
        paired_faces[1] = iFace;
      }
    }

    for (iDim = 0; iDim < nDim; ++iDim) {
      eigvecs[0][iDim] = Coord_FaceElems_CG[0][iDim] -
          Coord_FaceElems_CG[paired_faces[1]][iDim];
    }

    // Second vector
    paired_faces[2] = 0;
    paired_faces[3] = 0;
    for (iFace = 1; iFace < nFaces; ++iFace) {
      if (iFace != paired_faces[1]) {
        if (paired_faces[2] == 0) {
          paired_faces[2] = iFace;
        } else {
          paired_faces[3] = iFace;
        }
      }
    }

    for (iDim = 0; iDim < nDim; ++iDim) {
      eigvecs[1][iDim] = Coord_FaceElems_CG[paired_faces[2]][iDim] -
          Coord_FaceElems_CG[paired_faces[3]][iDim];
    }

    // Get magnitudes
    su2double eigvalues[nDim][nDim];
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        eigvalues[iDim][jDim] = 0.0;
      }
      eigvalues[iDim][iDim] = magnitude<2>(eigvecs[iDim]);
      for (jDim = 0; jDim < nDim; ++jDim) {
        eigvecs[iDim][jDim] /= eigvalues[iDim][iDim];
      }
    }
    std::cout << eigvalues[0][0] << "," << eigvalues[1][1] << std::endl;


    // Gram-Schmidt Process to make the vectors orthogonal
    su2double temp_eigvecs[nDim][nDim];
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        temp_eigvecs[iDim][jDim] = eigvecs[iDim][jDim];
      }
    }
    GramSchmidt<2>(temp_eigvecs, eigvecs);

    // Perform matrix multiplication
    for (iDim = 0; iDim < nDim; ++iDim) {
      for (jDim = 0; jDim < nDim; ++jDim) {
        for (kDim = 0; kDim < nDim; ++kDim) {
          for (lDim = 0; lDim < nDim; ++lDim) {
            Mij[iDim][jDim] += eigvecs[kDim][iDim]*
                eigvalues[kDim][lDim]*eigvecs[lDim][jDim];
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
  static const int nDim = 2;

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





