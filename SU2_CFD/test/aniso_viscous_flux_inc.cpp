/*!
 * \file aniso_viscous_flux.cpp
 * \brief checks whether the viscous projected flux is computed correctly for
 *        anisotropic eddy viscosities in incompressible flow.
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
const unsigned short nVar = 4;

class TestNumerics : public CNumerics {
 public:
  TestNumerics(CConfig* config) : CNumerics(3, 4, config) {
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

  su2double** gradprimvar = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    gradprimvar[iVar] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++)
      gradprimvar[iVar][jDim] = 0.0;
  }


  gradprimvar[1][1] =  1; // dU/dy
  gradprimvar[2][0] = -1; // dV/dx
  su2double normal[nDim] = {1.0, 0.0, 0.0};
  su2double laminar_viscosity = 0.000;
  su2double eddy_viscosity = 0.015;
  su2double** Aniso_Eddy_Viscosity = new su2double*[nDim];
  int counter = 1;
  for (int iDim = 0; iDim < nDim; iDim++) {
    Aniso_Eddy_Viscosity[iDim] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++) {
      Aniso_Eddy_Viscosity[iDim][jDim] = counter;
      counter++;
    }
  }
  //---------------------------------------------------------------------------
  // Test
  // We set up \nu_{ij} = [[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0]]
  //       and \tau_{ij} = [[0, 1, 0],[-1, 0, 0],[0, 0, 0]]
  // So that \tau_{ij} = \nu_{ik}*\pderiv{u_j}{x_k} + \nu_{jk}*\pderiv{u_i}{x_k}
  //         \tau_{11} = 2*\nu_{12}*\pderiv{u_1}{x_2} =  4.0
  //         \tau_{12} =   \nu_{11}*\pderiv{u_2}{x_1}
  //                     + \nu_{22}*\pderiv{u_1}{x_2} =  4.0
  //         \tau_{13} =   \nu_{32}*\pderiv{u_1}{x_2} =  8.0
  //         \tau_{22} = 2*\nu_{21}*\pderiv{u_2}{x_1} = -8.0
  //         \tau_{23} =   \nu_{31}*\pderiv{u_2}{x_1} = -7.0
  //         \tau_{33} =                              =  0.0
  // So the normal in the {1.0, 0.0, 0.0} direction is:
  //         \tau_{ij}*e_1 = [tau_{11}, tau_{12}, tau_{13}]
  //                       = [4.0, 4.0, 8.0]
  //---------------------------------------------------------------------------
  numerics.GetViscousArtCompProjFlux(gradprimvar, normal, laminar_viscosity,
                                     Aniso_Eddy_Viscosity);
  su2double* output = numerics.Proj_Flux_Tensor;
  su2double correct_output[nVar] = {0.0, 4.0, 4.0, 8.0};
  for (int iVar = 0; iVar < nVar; iVar++) {
    if (output[iVar] != correct_output[iVar]) {
      std::cout << "The projected flux tensor for an anisotropic eddy";
      std::cout << " viscosity was incorrect" << std::endl;
      std::cout << "    The test case was: incompressible flow." << std::endl;
      std::cout << "  Expected:" << std::endl;
      std::cout << "    [ " << correct_output[0] << ", ";
      std::cout << correct_output[1] << ", " << correct_output[2] << ", ";
      std::cout << correct_output[3] << "]" << std::endl;
      std::cout << "  Found:" << std::endl;
      std::cout << "    [ " << output[0] << ", " << output[1] << ", ";
      std::cout << output[2] << ", " << output[3] << "]" << std::endl;
      return_flag = 1;
      break;
    }
  }


  //---------------------------------------------------------------------------
  // Teardown
  //---------------------------------------------------------------------------
  delete test_config;
  for (int iVar = 0; iVar < nVar; iVar++) {
    delete [] gradprimvar[iVar];
  }
  delete [] gradprimvar;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





