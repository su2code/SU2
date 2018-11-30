/*!
 * \file viscous_ideal_vs_general.cpp
 * \brief This test checks whether the viscous residuals and jacobians match
 *        for ideal gases, regardless of which version of the numerics used.
 * \author C. Pederson
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <cmath>
#include <string.h>

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "../include/numerics_direct_mean.hpp"

unsigned short error_count;

void BOOST_CHECK_CLOSE_FRACTION(const su2double& u, const su2double& v,
                                const su2double& tol) {
  if (u != 0 && v != 0) {
    const su2double rel_error_u = std::abs((u - v)/u);
    const su2double rel_error_v = std::abs((u - v)/v);
    if (rel_error_u > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Relative difference exceeds tolerance ";
      std::cout << "[" << rel_error_u << " > " << tol << "]\n";
      error_count += 1;
    } else if (rel_error_v > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Relative difference exceeds tolerance ";
      std::cout << "[" << rel_error_v << " > " << tol << "]\n";
      error_count += 1;
    }
  } else {
    const su2double abs_error = std::abs(u - v);
    if (abs_error > tol) {
      std::cout << "error in test: [" << u << " != " << v << "]. ";
      std::cout << "Absolute difference exceeds tolerance ";
      std::cout << "[" << abs_error << " > " << tol << "]\n";
      error_count += 1;
    }
  }
}

const unsigned short nDim = 3;
const unsigned short nVar = nDim+2;
const unsigned short nPrimVar = nDim+9;
const unsigned short nSecVar = 4;

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
void WriteCfgFile() {

  std::ofstream cfg_file;

  cfg_file.open("test.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "TIME_DISCRE_FLOW= EULER_IMPLICIT" << std::endl;

  cfg_file.close();

}

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

/**
 * Sets up the test and runs it.
 *
 * This class is used to ensure that the same conditions are applied to
 * both types of numerics.  It encapsulates a lot of the test settings and
 * provides a method to run the test on arbitrary numerics.
 */
class TestRunner {
 public:

  /**
   * Construct the test problem and setup necessary variables
   * @param config - The problem defintion
   */
  TestRunner(CConfig* config);

  ~TestRunner();

  /**
   * Compute residuals and Jacobians given a set of numerics.
   *
   * The test runner is designed to be used for viscous residual calculations
   * only.
   *
   * \param[in] numerics - A viscous residual class.
   * \param[out] residual_i - The calculated residual
   * \param[out] Jacobian_i - The calculated jacobian at point i
   * \param[out] Jacobian_j - The calculated jacobian at point j
   * \param[in] config - The problem definition
   */
  void Run(CNumerics* numerics, su2double* residual_i,
           su2double** Jacobian_i, su2double** Jacobian_j, CConfig* config);

  // Geometry definition
  static const su2double distance = 1;
  static const su2double area = 1;

  // Gas properties
  static const su2double density = 1.0;
  static const su2double temperature = 300;

  // Extra variables to be defined
  static const su2double velocity_squared = 1.0;
  static const su2double turb_ke = 0.0;
  static const su2double eddy_viscosity = 0.0;
  static const su2double laminar_viscosity = 1.0;
  static const su2double laminar_prandtl = 1.0;

  // Models
  CIdealGas ideal_gas;
  CConstantPrandtl conductivity_model;

  // Velocity
  su2double velocity[nDim];

  // Variables set in the numerics class
  su2double coord_i[nDim], coord_j[nDim];
  su2double normal[nDim];
  su2double** primvar_grad;
  su2double primvar[nPrimVar];
  su2double secvar[nSecVar];
};

TestRunner::TestRunner(CConfig* config)
  : ideal_gas(config->GetGamma(), config->GetGas_ConstantND()),
    conductivity_model(config->GetPrandtl_Lam()) {

  ideal_gas.SetTDState_rhoT(density, temperature);
  conductivity_model.SetConductivity(temperature, density, laminar_viscosity,
                                     eddy_viscosity, ideal_gas.GetCp());

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    velocity[iDim] = sqrt(velocity_squared/3.0);
  }

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    coord_i[iDim] = 0; coord_j[iDim] = 0;
  }
  coord_j[0] += distance;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    normal[iDim] = coord_j[iDim]/distance*area;
  }

  primvar_grad = new su2double*[nPrimVar];
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_grad[iVar] = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad[iVar][iDim] = 0.0;
    }
  }

  primvar_grad[1][0] =  1.0; // du/dx
  primvar_grad[2][1] =  2.0; // dv/dy
  primvar_grad[1][1] =  1.0; // du/dy
  primvar_grad[2][0] =  1.0; // dv/dx
  primvar_grad[3][0] =  1.0; // dw/dx
  primvar_grad[1][2] =  1.0; // du/dz
  primvar_grad[0][0] =  1.0; // dT/dx
  primvar_grad[0][1] =  1.0; // dT/dy
  primvar_grad[0][2] =  1.0; // dT/dz

  const su2double static_energy = ideal_gas.GetStaticEnergy();
  const su2double total_energy = static_energy + 0.5*velocity_squared + 0.5*turb_ke;
  const su2double pressure = ideal_gas.GetPressure();
  const su2double enthalpy = total_energy + pressure / density;
  primvar[0] = ideal_gas.GetTemperature();
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    primvar[iDim+1] = velocity[iDim];
  }
  primvar[nDim+1] = ideal_gas.GetPressure();
  primvar[nDim+2] = ideal_gas.GetDensity();
  primvar[nDim+3] = enthalpy;
  primvar[nDim+5] = laminar_viscosity;
  primvar[nDim+6] = eddy_viscosity;
  primvar[nDim+7] = conductivity_model.GetConductivity();
  primvar[nDim+8] = ideal_gas.GetCp();

  secvar[2] = ideal_gas.GetdTdrho_e();
  secvar[3] = ideal_gas.GetdTde_rho();
}

TestRunner::~TestRunner() {

  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    delete [] primvar_grad[iVar];
  }
  delete [] primvar_grad;
}

void TestRunner::Run(CNumerics* numerics, su2double* residual_i,
                     su2double** Jacobian_i, su2double** Jacobian_j,
                     CConfig* config) {

  // Setup the variables needed for a viscous residual calculation

  numerics->SetCoord(coord_i, coord_j);
  numerics->SetNormal(normal);
  numerics->SetSecondary(secvar, secvar);
  numerics->SetPrimitive(primvar, primvar);
  numerics->SetPrimVarGradient(primvar_grad, primvar_grad);
  numerics->SetTurbKineticEnergy(turb_ke, turb_ke);

  // Compute the residual and jacobians

  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

}

/**
 * Compare the viscous numerics for the ideal gas and the generalized gas
 * law variants.  They **should** be identical, if an ideal gas law is used
 * for both versions.
 */
int main(int argc, char *argv[]) {

  /*--- MPI initialization, and buffer setting ---*/

#ifdef HAVE_MPI
  const unsigned int BUFSIZE = 3000000;        /*!< \brief MPI buffer. */
  int  buffsize;
  char *buffptr;
  SU2_MPI::Init(&argc, &argv);
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
#else
  SU2_Comm MPICommunicator(0);
#endif

  const int rank = SU2_MPI::GetRank();
  const int size = SU2_MPI::GetSize();

  error_count = 0;


  /*--- Setup ---*/

  WriteCfgFile();
  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig("test.cfg", SU2_CFD, iZone, nZone, 2, VERB_NONE);
  config->SetGas_ConstantND(287.058);

  CNumerics* ideal_numerics = new CAvgGrad_Flow(3, 5, config);
  CNumerics* general_numerics = new CGeneralAvgGrad_Flow(3, 5, config);

  su2double ideal_residual_i[nVar] = {0, 0, 0, 0, 0};
  su2double general_residual_i[nVar] = {0, 0, 0, 0, 0};

  su2double** ideal_jacobian_i = new su2double*[nVar];
  su2double** ideal_jacobian_j = new su2double*[nVar];
  su2double** general_jacobian_i = new su2double*[nVar];
  su2double** general_jacobian_j = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    ideal_jacobian_i[iVar] = new su2double[nVar];
    ideal_jacobian_j[iVar] = new su2double[nVar];
    general_jacobian_i[iVar] = new su2double[nVar];
    general_jacobian_j[iVar] = new su2double[nVar];
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      ideal_jacobian_i[iVar][jVar] = 0;
      ideal_jacobian_j[iVar][jVar] = 0;
      general_jacobian_i[iVar][jVar] = 0;
      general_jacobian_j[iVar][jVar] = 0;
    }
  }

  /*--- Test ---*/

  TestRunner test_runner(config);

  test_runner.Run(ideal_numerics, ideal_residual_i, ideal_jacobian_i,
                  ideal_jacobian_j, config);
  test_runner.Run(general_numerics, general_residual_i, general_jacobian_i,
                  general_jacobian_j, config);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar=0; iVar < nVar; iVar++) {
    BOOST_CHECK_CLOSE_FRACTION(ideal_residual_i[iVar],
                               general_residual_i[iVar], tolerance);
  }
  for (unsigned short iVar=0; iVar < nVar; iVar++) {
    for (unsigned short jVar=0; jVar < nVar; jVar++) {
      std::cout << "jacobian[" << iVar << "][" << jVar << "]\n";
      BOOST_CHECK_CLOSE_FRACTION(ideal_jacobian_i[iVar][jVar],
                                 general_jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_CLOSE_FRACTION(ideal_jacobian_j[iVar][jVar],
                                 general_jacobian_j[iVar][jVar], tolerance);
    }
  }

  /*--- Teardown ---*/

  delete config;
  delete ideal_numerics;
  delete general_numerics;

  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    delete [] ideal_jacobian_i[iVar];
    delete [] ideal_jacobian_j[iVar];
    delete [] general_jacobian_i[iVar];
    delete [] general_jacobian_j[iVar];
  }
  delete [] ideal_jacobian_i;
  delete [] ideal_jacobian_j;
  delete [] general_jacobian_i;
  delete [] general_jacobian_j;

  /*--- Finalize MPI parallelization ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif

  if (error_count == 0)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
