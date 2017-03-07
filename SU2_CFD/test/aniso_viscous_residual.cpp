/*!
 * \file aniso_viscous_residual.cpp
 * \brief This test checks whether the correct viscous residual is calculated
 *        for anisotropic viscosities on incompressible flows.
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

#include <iomanip>
#include <iostream>
#include <typeinfo>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

const unsigned short nDim = 3;
const unsigned short nVar = 4;

class TestNumerics : public CAvgGradArtComp_Flow {
 public:
  bool useAnisoEddyViscosity;

  TestNumerics(CConfig* config) : CAvgGradArtComp_Flow(3, 4, config) {
    V_i = new su2double[nDim+4];
    V_j = new su2double[nDim+4];

    implicit = false;
  };

  void ComputeMyResidual(su2double *val_residual, su2double **val_Jacobian_i,
                         su2double **val_Jacobian_j, CConfig *config) {
    unsigned short jDim;
    /*--- Normalized normal vector ---*/
    Area = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      Area += Normal[iDim]*Normal[iDim];
    Area = sqrt(Area);

    for (iDim = 0; iDim < nDim; iDim++)
      UnitNormal[iDim] = Normal[iDim]/Area;

    /*--- Laminar and Eddy viscosity ---*/

    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];

    /*--- Mean Viscosities ---*/

    Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
    if (useAnisoEddyViscosity) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Mean_Aniso_Eddy_Viscosity[iDim][jDim] =
              0.5*(Aniso_Eddy_Viscosity_i[iDim][jDim] +
                   Aniso_Eddy_Viscosity_j[iDim][jDim]);
    } else {
      Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
    }

    /*--- Mean gradient approximation ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);

    /*--- Get projected flux tensor ---*/

    if (useAnisoEddyViscosity) {
      GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal,
                                Mean_Laminar_Viscosity,
                                Mean_Aniso_Eddy_Viscosity);
    } else {

      GetViscousArtCompProjFlux(Mean_GradPrimVar, Normal,
                                Mean_Laminar_Viscosity, Mean_Eddy_Viscosity);
    }

    /*--- Update viscous residual ---*/

    for (iVar = 0; iVar < nVar; iVar++)
      val_residual[iVar] = Proj_Flux_Tensor[iVar];

    /*--- Implicit part ---*/

    if (implicit) {

      dist_ij = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        dist_ij += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);
      dist_ij = sqrt(dist_ij);

      if (dist_ij == 0.0) {
        for (iVar = 0; iVar < nVar; iVar++) {
          for (jVar = 0; jVar < nVar; jVar++) {
            val_Jacobian_i[iVar][jVar] = 0.0;
            val_Jacobian_j[iVar][jVar] = 0.0;
          }
        }
      }
      else {
        GetViscousArtCompProjJacs(Mean_Laminar_Viscosity, Mean_Eddy_Viscosity, dist_ij, UnitNormal,
                                  Area, val_Jacobian_i, val_Jacobian_j);
      }

    }

  }
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

  su2double* residual = new su2double[nVar];
  su2double** Jacobian_i;
  su2double** Jacobian_j;

  su2double** eddy_viscosity_i = new su2double*[nDim];
  int counter = 2;
  for (int iDim = 0; iDim < nDim; iDim++) {
    eddy_viscosity_i[iDim] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity_i[iDim][jDim] = counter;
      counter++;
    }
  }

  su2double** eddy_viscosity_j = new su2double*[nDim];
  counter = 0;
  for (int iDim = 0; iDim < nDim; iDim++) {
    eddy_viscosity_j[iDim] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity_j[iDim][jDim] = counter;
      counter++;
    }
  }

  su2double** gradprimvar = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    gradprimvar[iVar] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++)
      gradprimvar[iVar][jDim] = 0.0;
  }
  gradprimvar[1][1] =  1; // dU/dy
  gradprimvar[2][0] = -1; // dV/dx

  su2double normal[3] = {1.0, 0.0, 0.0};

  numerics.useAnisoEddyViscosity = true;
  numerics.SetNormal(normal);
  numerics.SetPrimVarGradient(gradprimvar, gradprimvar);
  numerics.SetEddyViscosity(0.0, 0.0);
  numerics.SetLaminarViscosity(0.0, 0.0);
  numerics.SetAnisoEddyViscosity(eddy_viscosity_i, eddy_viscosity_j);

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
  numerics.ComputeMyResidual(residual, Jacobian_i, Jacobian_j, test_config);
  su2double* output = residual;
  su2double correct_output[nVar] = {0.0, 4.0, 4.0, 8.0};
  for (int iVar = 0; iVar < nVar; iVar++) {
    if (output[iVar] != correct_output[iVar]) {
      std::cout << "The computed viscous residual for an anisotropic eddy";
      std::cout << " viscosity was incorrect" << std::endl;
      std::cout << "    The test case was: incompressible flow." << std::endl;
      std::cout << "  Expected:" << std::endl;
      std::cout << "    [ " << correct_output[0] << ", ";
      std::cout << correct_output[1] << ", " << correct_output[2] << ", ";
      std::cout << correct_output[3] << "]" << std::endl;
      std::cout << "  Found:" << std::endl;
      std::cout << "    [ " << output[0] << ", " << output[1] << ", ";
      std:;cout << output[2] << ", " << output[3] << "]" << std::endl;
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

  for (int iDim = 0; iDim < nDim; iDim++) {
    delete [] eddy_viscosity_i[iDim];
    delete [] eddy_viscosity_j[iDim];
  }
  delete [] eddy_viscosity_i;
  delete [] eddy_viscosity_j;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return return_flag;
}





