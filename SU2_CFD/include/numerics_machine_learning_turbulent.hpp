/*!
 * \file numerics_machine_learning_direct_turbulent.hpp
 * \brief Header for caller functions of the turbulence models.
 * \author B. Tracey
 * \version 3.2.9 "eagle"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <cmath>
#include <iostream>

using namespace std;

class SpalartAllmarasOtherOutputs{
public:
  SpalartAllmarasOtherOutputs();
  ~SpalartAllmarasOtherOutputs();
  double fw;
  double mul_production; // multiplier of OmegaNu
  double mul_destruction; // multiplier of (Nu/ d)^2
  double mul_crossproduction; // multiplier of dnuhat_i / dx_i
  double Omega; // Sqrt vorticity
};

class SpalartAllmarasConstants{
public:
  SpalartAllmarasConstants();
  ~SpalartAllmarasConstants();
  double  cv1_3,
          k2,
          cb1,
          cw2,
          cw3_6,
          sigma,
          cb2,
          cb2_sigma,
          cw1;
};

class SpalartAllmarasInputs{
private:
  int  nDim;
  double limiter; // How close to the wall should the turbulence model be turned off
  double**  DUiDXj;  // Mean flow derivative
  double* DTurb_Kin_Visc_DXj; // D NuTilde D X
  void init(int nDim, double limiter);
public:
  SpalartAllmarasInputs(int nDim);
  SpalartAllmarasInputs(int nDim, double limiter);
  ~SpalartAllmarasInputs();
  void Set(double** DUiDXj, double* DTurb_Kin_Visc_DXj, bool rotating_frame, bool transition, double dist, double Laminar_Viscosity, double Density, double Turbulent_Kinematic_Viscosity, double intermittency);
  int GetNumDim();
  double GetLimiter();
  double** GetMeanFlowGradient();
  double* GetTurbKinViscGradient();
  bool rotating_frame;
  bool transition;
  double Omega;
  double dist; // Wall distance
  double Laminar_Viscosity;
  double Density;
  double Turbulent_Kinematic_Viscosity;
  double intermittency; // Used for transition
};

/* \brief computes spalart allmaras source term. See
  numerics_machine_learning_direct_turbulent.cpp */
void SpalartAllmarasSourceTerm(SpalartAllmarasInputs* inputs, SpalartAllmarasConstants* constants, double* output_residual, double* output_jacobian, SpalartAllmarasOtherOutputs* otherOutput);

/* \brief Computes the vorticity from the velocity gradient
 tensor */
double ComputeVorticity(int nDim, double** DUiDXj);
