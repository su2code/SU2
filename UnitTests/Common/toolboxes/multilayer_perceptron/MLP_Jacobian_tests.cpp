/*!
 * \file MLP_Jacobian_tests.cpp
 * \brief Check if MLP is differentiable
 * \author E.C.Bunschoten
 * \version 8.5.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

#include "catch.hpp"
#include "../../../../Common/include/CConfig.hpp"
#if defined(HAVE_MLPCPP)
#define MLP_CUSTOM_TYPE su2double
#include "../../../../subprojects/MLPCpp/include/CLookUp_ANN.hpp"
#define USE_MLPCPP
#endif
#include <vector>

#ifdef USE_MLPCPP
TEST_CASE("MLP Jacobian test", "[LookUpANN]") {
  /* Create network with random weights */
  std::vector<size_t> NN = {3, 10, 10, 1};
  MLPToolbox::CNeuralNetwork mlp = MLPToolbox::CNeuralNetwork(NN);
  mlp.SetActivationFunction("gelu");
  mlp.RandomWeights();

  /* Enable Jacobian and Hessian computation */
  mlp.CalcJacobian(true);
  mlp.CalcHessian(true);

  /* Network input values */
  su2double val_in_1, val_in_2, val_in_3;
  val_in_1 = 0.2;
  val_in_2 = 0.3;
  val_in_3 = 0.9;

  /* Calculate Jacobians with AD */
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(val_in_1);
  AD::RegisterInput(val_in_2);
  AD::RegisterInput(val_in_3);

  mlp.SetInput(0, val_in_1);
  mlp.SetInput(1, val_in_2);
  mlp.SetInput(2, val_in_3);
  mlp.Predict();

  su2double val_out = mlp.GetOutput(0);

  AD::RegisterOutput(val_out);
  AD::StopRecording();
  SU2_TYPE::SetDerivative(val_out, 1.0);
  AD::ComputeAdjoint();
  const su2double jacobian_AD_1 = SU2_TYPE::GetDerivative(val_in_1);
  const su2double jacobian_AD_2 = SU2_TYPE::GetDerivative(val_in_2);
  const su2double jacobian_AD_3 = SU2_TYPE::GetDerivative(val_in_3);

  /* Retrieve analytical Jacobians */
  const su2double jacobian_ref_1 = mlp.GetJacobian(0, 0);
  const su2double jacobian_ref_2 = mlp.GetJacobian(0, 1);
  const su2double jacobian_ref_3 = mlp.GetJacobian(0, 2);

  /* Unit tests for Jacobians */
  CHECK(SU2_TYPE::GetValue(jacobian_AD_1) == Approx(SU2_TYPE::GetValue(jacobian_ref_1)));
  CHECK(SU2_TYPE::GetValue(jacobian_AD_2) == Approx(SU2_TYPE::GetValue(jacobian_ref_2)));
  CHECK(SU2_TYPE::GetValue(jacobian_AD_3) == Approx(SU2_TYPE::GetValue(jacobian_ref_3)));

  /* Calculate Hessians with AD */
  AD::Reset();
  AD::StartRecording();
  AD::RegisterInput(val_in_1);
  AD::RegisterInput(val_in_2);
  AD::RegisterInput(val_in_3);

  mlp.SetInput(0, val_in_1);
  mlp.SetInput(1, val_in_2);
  mlp.SetInput(2, val_in_3);
  mlp.Predict();

  su2double jacobian_analytical = mlp.GetJacobian(0, 0);

  AD::RegisterOutput(jacobian_analytical);
  AD::StopRecording();
  SU2_TYPE::SetDerivative(jacobian_analytical, 1.0);
  AD::ComputeAdjoint();
  const su2double hessian_AD_1 = SU2_TYPE::GetDerivative(val_in_1);
  const su2double hessian_AD_2 = SU2_TYPE::GetDerivative(val_in_2);
  const su2double hessian_AD_3 = SU2_TYPE::GetDerivative(val_in_3);

  /* Retrieve analytical Hessians */
  const su2double hessian_ref_1 = mlp.GetHessian(0, 0, 0);
  const su2double hessian_ref_2 = mlp.GetHessian(0, 0, 1);
  const su2double hessian_ref_3 = mlp.GetHessian(0, 0, 2);

  /* Unit tests on Hessians */
  CHECK(SU2_TYPE::GetValue(hessian_AD_1) == Approx(SU2_TYPE::GetValue(hessian_ref_1)));
  CHECK(SU2_TYPE::GetValue(hessian_AD_2) == Approx(SU2_TYPE::GetValue(hessian_ref_2)));
  CHECK(SU2_TYPE::GetValue(hessian_AD_3) == Approx(SU2_TYPE::GetValue(hessian_ref_3)));
}
#endif
