/*!
 * \file aniso_viscous_jacobian_compare.cpp
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
const unsigned short nVar = 7;

template<std::size_t M, std::size_t N>
void PrintMatrix(su2double** A) {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "  [[";
  for (unsigned int i = 0; i < M; i++) {
    for (unsigned int j = 0; j < N; j++) {
      std::cout << A[i][j];
      if (j<N-1) std::cout << ",\t";
    }
    if (i < M-1) std::cout << "]" << std::endl << "   [";
  }
  std::cout << "]]" << std::endl;
}

class TestNumerics : public CAvgGrad_Flow {
 public:
  TestNumerics(CConfig* config, bool hasAnisoViscosity)
    : CAvgGrad_Flow(3, 6, config, hasAnisoViscosity) {
    implicit = true;
  };
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

  TestNumerics isotropic_numerics  (test_config, false);
  TestNumerics anisotropic_numerics(test_config, true);

  /*--- Inputs ---*/

  su2double primvar_i[nVar+5];
  su2double primvar_j[nVar+5];
  for (int iVar = 0; iVar < nVar+5; iVar++) {
    primvar_i[iVar] = 1.0; primvar_j[iVar] = 1.0;
  }

  su2double** primvar_grad_i = new su2double*[nVar];
  su2double** primvar_grad_j = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    primvar_grad_i[iVar] = new su2double[nDim];
    primvar_grad_j[iVar] = new su2double[nDim];
    for (int iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_i[iVar][iDim] = 0.0;
    }
  }
  primvar_grad_i[1][1] = 1.0;
  primvar_grad_j[1][1] = 1.0;
  primvar_grad_i[2][2] = 0.3;
  primvar_grad_j[2][2] = 0.3;


  su2double** resolution_i = new su2double*[nDim];
  su2double** resolution_j = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    resolution_i[iDim] = new su2double[nDim];
    resolution_j[iDim] = new su2double[nDim];
    resolution_i[iDim][iDim] = 1.0;
    resolution_j[iDim][iDim] = 1.0;
  }

  su2double* hybrid_param_i = new su2double[1];
  su2double* hybrid_param_j = new su2double[1];
  hybrid_param_i[0] = 1.0;
  hybrid_param_j[0] = 1.0;

  su2double** anisotropy_i = new su2double*[nDim];
  su2double** anisotropy_j = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    anisotropy_i[iDim] = new su2double[nDim];
    anisotropy_j[iDim] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++) {
      anisotropy_i[iDim][jDim] = (iDim == jDim);
      anisotropy_j[iDim][jDim] = (iDim == jDim);
    }
  }

  su2double normal[nDim] = {0.577350269, 0.577350269, 0.577350269};

  su2double coord_i[nDim] = {0.0, 0.0, 0.0};
  su2double coord_j[nDim] = {0.577350269, 0.577350269, 0.577350269};

  /*--- Outputs ---*/

  su2double** Jacobian_i_a = new su2double*[nVar];
  su2double** Jacobian_j_a = new su2double*[nVar];
  su2double** Jacobian_i_b = new su2double*[nVar];
  su2double** Jacobian_j_b = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i_a[iVar] = new su2double[nVar];
    Jacobian_j_a[iVar] = new su2double[nVar];
    Jacobian_i_b[iVar] = new su2double[nVar];
    Jacobian_j_b[iVar] = new su2double[nVar];
  }

  su2double* residual = new su2double[nVar];

  /**-------------------------------------------------------------------------
   * TEST
   *
   *--------------------------------------------------------------------------
   */
   su2double tol = 10*std::numeric_limits<su2double>::epsilon();
   isotropic_numerics.SetNormal(normal);
   isotropic_numerics.SetCoord(coord_i, coord_j);
   isotropic_numerics.SetResolutionTensor(resolution_i, resolution_j);
   isotropic_numerics.SetPrimitive(primvar_i, primvar_j);
   isotropic_numerics.SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
   isotropic_numerics.ComputeResidual(residual, Jacobian_i_a, Jacobian_j_a, test_config);

   anisotropic_numerics.SetHybridParameter(hybrid_param_i, hybrid_param_j);
   anisotropic_numerics.SetNormal(normal);
   anisotropic_numerics.SetCoord(coord_i, coord_j);
   anisotropic_numerics.SetResolutionTensor(resolution_i, resolution_j);
   anisotropic_numerics.SetPrimitive(primvar_i, primvar_j);
   anisotropic_numerics.SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
   anisotropic_numerics.SetEddyViscAnisotropy(anisotropy_i, anisotropy_j);
   anisotropic_numerics.ComputeResidual(residual, Jacobian_i_b, Jacobian_j_b, test_config);

   for (int iVar = 0; iVar < nVar; iVar++) {
     for (int jVar = 0; jVar < nVar; jVar++) {
       if (fabs(Jacobian_i_a[iVar][jVar] - Jacobian_i_b[iVar][jVar]) > 1e-7)
        return_flag = 1;
     }
   }
   if (return_flag == 1) {
     std::cout << "ERROR: The anisotropic jacobian did not match an isotropic jacobian." << std::endl;
     std::cout << "  Isotropic:" << std::endl;
     PrintMatrix<nVar,nVar>(Jacobian_i_a);
     std::cout << endl;
     std::cout << "  Found:" << std::endl;
     PrintMatrix<nVar,nVar>(Jacobian_i_b);
   }

  /**-------------------------------------------------------------------------
   * TEARDOWN
   *
   *--------------------------------------------------------------------------
   */
   delete test_config;
   for (int iDim = 0; iDim < nDim; iDim++) {
     delete[] resolution_i[iDim];
     delete[] resolution_j[iDim];
   }
   delete[] resolution_i;
   delete[] resolution_j;
   for (int iVar = 0; iVar < nVar; iVar++) {
     delete[] Jacobian_i_a[iVar];
     delete[] Jacobian_j_a[iVar];
     delete[] Jacobian_i_b[iVar];
     delete[] Jacobian_j_b[iVar];
     delete[] primvar_grad_i[iVar];
     delete[] primvar_grad_j[iVar];
   }
   delete[] Jacobian_i_a;
   delete[] Jacobian_j_a;
   delete[] Jacobian_i_b;
   delete[] Jacobian_j_b;
   delete[] primvar_grad_i;
   delete[] primvar_grad_j;


#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return return_flag;
}
