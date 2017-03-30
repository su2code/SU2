/*!
 * \file aniso_viscous_flux.cpp
 * \brief checks whether the viscous projected flux is computed correctly for
 *        anisotropic eddy viscosities in compressible flow.
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
#include <typeinfo>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

const unsigned short nDim = 3;
const unsigned short nVar = 6;

class TestNumerics : public CNumerics {
 public:
  TestNumerics(CConfig* config) : CNumerics(3, 6, config) {
  };
};

int main() {

  //---------------------------------------------------------------------------
  // Setup
  //---------------------------------------------------------------------------
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
#endif

  int return_flag=0;

  CConfig* test_config = new CConfig();

  TestNumerics numerics(test_config);

  su2double normal[nDim] = {1.0, 0.0, 0.0};
  su2double prim_var[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  const su2double dist_ij = 1.0;
  const su2double dS = 3.0;
  const su2double laminar_viscosity = 0.000;

  /*-- Build the eddy viscosity --*/
  su2double** eddy_viscosity = new su2double*[nDim];
  for (int iDim = 0; iDim < nDim; iDim++) {
    eddy_viscosity[iDim] = new su2double[nDim];
  }
  eddy_viscosity[0][0] = 1.0;
  eddy_viscosity[0][1] = 2.0;
  eddy_viscosity[0][2] = 3.0;
  eddy_viscosity[1][0] = 2.0;
  eddy_viscosity[1][1] = 4.0;
  eddy_viscosity[1][2] = 5.0;
  eddy_viscosity[2][0] = 3.0;
  eddy_viscosity[2][1] = 5.0;
  eddy_viscosity[2][2] = 6.0;

  /*-- Build the results vectors --*/
  su2double* Proj_Visc_Flux = new su2double[nVar];
  su2double** Proj_Jac_Tensor_i = new su2double*[nVar];
  su2double** Proj_Jac_Tensor_j = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    Proj_Jac_Tensor_i[iVar] = new su2double[nVar];
    Proj_Jac_Tensor_j[iVar] = new su2double[nVar];
  }

  /**---------------------------------------------------------------------------
   * Test
   * 
   *  We set up:
   *  dS = 3.0
   *  dist = 1.0
   *  \mu = [[1.0, 2.0, 3.0],[2.0, 4.0, 5.0],[3.0, 5.0, 6.0]]
   *  n = [1, 0, 0]
   *
   *  So the end Jacobian should be:
   *  jac = [[-20, 12, 18], [-2, -3, 0], [-3, 0, -3]]
   *
   *---------------------------------------------------------------------------
   */
  const su2double tol = 5*std::numeric_limits<su2double>::epsilon();
  numerics.GetViscousProjJacs(prim_var, laminar_viscosity, eddy_viscosity,
                              dist_ij, normal, dS, Proj_Visc_Flux,
                              Proj_Jac_Tensor_i, Proj_Jac_Tensor_j);
  if (std::abs(Proj_Jac_Tensor_i[1][1] - (-20)) > tol ||
      std::abs(Proj_Jac_Tensor_i[1][2] - ( 12)) > tol ||
      std::abs(Proj_Jac_Tensor_i[1][3] - ( 18)) > tol ||
      std::abs(Proj_Jac_Tensor_i[2][1] - (- 2)) > tol ||
      std::abs(Proj_Jac_Tensor_i[2][2] - (- 3)) > tol ||
      std::abs(Proj_Jac_Tensor_i[2][3] - (  0)) > tol ||
      std::abs(Proj_Jac_Tensor_i[3][1] - (- 3)) > tol ||
      std::abs(Proj_Jac_Tensor_i[3][2] - (  0)) > tol ||
      std::abs(Proj_Jac_Tensor_i[3][3] - (- 3)) > tol ) return_flag = 1;
  if (return_flag == 1) {
    std::cout << "The projected flux jacobian for an anisotropic eddy";
    std::cout << " viscosity was incorrect" << std::endl;
    std::cout << "  The test case was: compressible flow." << std::endl;
    std::cout << "  Expected:" << std::endl;
    std::cout << "     [[-34, -12, -18]," << std::endl;
    std::cout << "      [ -6,  -3,   0]," << std::endl;
    std::cout << "      [ -9,   0,  -3]]" << std::endl;
    std::cout << "  Found:" << std::endl;
    std::cout << "    [[";
    std::cout << Proj_Jac_Tensor_i[1][1] << ", ";
    std::cout << Proj_Jac_Tensor_i[1][2] << ", ";
    std::cout << Proj_Jac_Tensor_j[1][3] << "]" << std::endl;
    std::cout << "     [";
    std::cout << Proj_Jac_Tensor_i[2][1] << ", ";
    std::cout << Proj_Jac_Tensor_i[2][2] << ", ";
    std::cout << Proj_Jac_Tensor_j[2][3] << "]" << std::endl;
    std::cout << "     [";
    std::cout << Proj_Jac_Tensor_i[3][1] << ", ";
    std::cout << Proj_Jac_Tensor_i[3][2] << ", ";
    std::cout << Proj_Jac_Tensor_j[3][3] << "]]" << std::endl;
  }


  //---------------------------------------------------------------------------
  // Teardown
  //---------------------------------------------------------------------------
  delete test_config;
//  for (int iVar = 0; iVar < nVar; iVar++) {
//    delete [] Proj_Jac_Tensor_i[iVar];
//    delete [] Proj_Jac_Tensor_j[iVar];
//  }
//  delete [] Proj_Jac_Tensor_i;
//  delete [] Proj_Jac_Tensor_j;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





