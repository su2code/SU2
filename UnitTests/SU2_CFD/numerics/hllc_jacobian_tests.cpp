/*!
 * \file hllc_jacobian_tests.cpp
 * \brief Finite-difference verification of the HLLC Jacobians (ideal and general gas, 2D and 3D), on fixed and moving
 *        faces, exact (USE_ACCURATE_FLUX_JACOBIANS) or default.
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
#include <array>
#include <cmath>
#include <memory>
#include <sstream>
#include "../../../SU2_CFD/include/numerics/flow/convection/hllc.hpp"

namespace {

constexpr su2double gamma = 1.4;
constexpr unsigned short maxVar = 5;

using State = std::array<su2double, maxVar>;  // conservative variables (rho, rho u, rho v, [rho w], rho E)

/*--- Primitive (T, u, ..., p, rho, h, c) and secondary (dp/drho_e, dp/de_rho) variables of an ideal gas whose total
 energy contains the turbulent kinetic energy k. ---*/
void Primitives(unsigned short nDim, const State& U, su2double k, su2double* V, su2double* S) {
  const su2double rho = U[0];
  su2double q2 = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    V[iDim + 1] = U[iDim + 1] / rho;
    q2 += V[iDim + 1] * V[iDim + 1];
  }
  const su2double e = U[nDim + 1] / rho - 0.5 * q2 - k;
  const su2double p = (gamma - 1) * rho * e;
  V[0] = p / (rho * 287.058);
  V[nDim + 1] = p;
  V[nDim + 2] = rho;
  V[nDim + 3] = (U[nDim + 1] + p) / rho;
  V[nDim + 4] = sqrt(gamma * p / rho);
  S[0] = (gamma - 1) * e;
  S[1] = (gamma - 1) * rho;
}

/*--- Conservative state from density, velocity, pressure and k. ---*/
State Conservative(unsigned short nDim, su2double rho, const su2double* vel, su2double p, su2double k) {
  State U{};
  su2double q2 = 0.0;
  U[0] = rho;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    U[iDim + 1] = rho * vel[iDim];
    q2 += vel[iDim] * vel[iDim];
  }
  U[nDim + 1] = p / (gamma - 1) + rho * (0.5 * q2 + k);
  return U;
}

struct Result {
  State flux;
  su2double jac_i[maxVar][maxVar], jac_j[maxVar][maxVar];
};

Result Evaluate(unsigned short nDim, CNumerics& numerics, const CConfig* config, const State& Ui, const State& Uj,
                su2double k_i, su2double k_j) {
  su2double Vi[20] = {0.0}, Vj[20] = {0.0}, Si[2], Sj[2];
  Primitives(nDim, Ui, k_i, Vi, Si);
  Primitives(nDim, Uj, k_j, Vj, Sj);
  numerics.SetPrimitive(Vi, Vj);
  numerics.SetSecondary(Si, Sj);
  numerics.SetTurbKineticEnergy(k_i, k_j);
  const auto res = numerics.ComputeResidual(config);
  Result r;
  for (unsigned short iVar = 0; iVar < nDim + 2; iVar++) {
    r.flux[iVar] = res.residual[iVar];
    for (unsigned short jVar = 0; jVar < nDim + 2; jVar++) {
      r.jac_i[iVar][jVar] = res.jacobian_i[iVar][jVar];
      r.jac_j[iVar][jVar] = res.jacobian_j[iVar][jVar];
    }
  }
  return r;
}

/*--- Largest difference between the Jacobians and central differences of the flux, relative to the largest entry. ---*/
su2double JacobianError(unsigned short nDim, CNumerics& numerics, const CConfig* config, const State& Ui,
                        const State& Uj, su2double k_i, su2double k_j) {
  const unsigned short nVar = nDim + 2;
  const auto ref = Evaluate(nDim, numerics, config, Ui, Uj, k_i, k_j);
  su2double maxErr = 0.0, maxJac = 0.0;
  for (int side = 0; side < 2; side++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      State Up = side ? Uj : Ui, Um = Up;
      const su2double h = 1e-6 * std::max(std::abs(Up[jVar]), su2double(1.0));
      Up[jVar] += h;
      Um[jVar] -= h;
      const auto fp = side ? Evaluate(nDim, numerics, config, Ui, Up, k_i, k_j)
                           : Evaluate(nDim, numerics, config, Up, Uj, k_i, k_j);
      const auto fm = side ? Evaluate(nDim, numerics, config, Ui, Um, k_i, k_j)
                           : Evaluate(nDim, numerics, config, Um, Uj, k_i, k_j);
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        const su2double fd = (fp.flux[iVar] - fm.flux[iVar]) / (2 * h);
        const su2double jac = side ? ref.jac_j[iVar][jVar] : ref.jac_i[iVar][jVar];
        maxErr = std::max(maxErr, std::abs(jac - fd));
        maxJac = std::max(maxJac, std::abs(jac));
      }
    }
  }
  return maxErr / maxJac;
}

/*--- Contact speed sM of the scheme, in the frame of the face (w = grid velocity along the unit normal n). ---*/
su2double ContactSpeed(unsigned short nDim, const State& Ui, const State& Uj, su2double k_i, su2double k_j,
                       const su2double* n, su2double w, su2double sL, su2double sR) {
  su2double Vi[20] = {0.0}, Vj[20] = {0.0}, S[2];
  Primitives(nDim, Ui, k_i, Vi, S);
  Primitives(nDim, Uj, k_j, Vj, S);
  su2double u_i = -w, u_j = -w;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    u_i += Vi[iDim + 1] * n[iDim];
    u_j += Vj[iDim + 1] * n[iDim];
  }
  const su2double p_i = Vi[nDim + 1], p_j = Vj[nDim + 1], rho_i = Vi[nDim + 2], rho_j = Vj[nDim + 2];
  return (p_i - p_j - rho_i * u_i * (sL - u_i) + rho_j * u_j * (sR - u_j)) / (rho_j * (sR - u_j) - rho_i * (sL - u_i));
}

std::unique_ptr<CConfig> MakeConfig(bool moving, bool accurate) {
  std::stringstream options;
  options << "SOLVER= EULER\nCONV_NUM_METHOD_FLOW= HLLC\nTIME_DISCRE_FLOW= EULER_IMPLICIT\n";
  options << "USE_ACCURATE_FLUX_JACOBIANS= " << (accurate ? "YES" : "NO") << "\n";
  if (moving) options << "GRID_MOVEMENT= STEADY_TRANSLATION\nTRANSLATION_RATE= 0.0 0.0 0.0\n";
  return std::make_unique<CConfig>(options, SU2_COMPONENT::SU2_CFD, false);
}

/*--- With USE_ACCURATE_FLUX_JACOBIANS= YES the Jacobians are exact for fixed wave speeds in all branches. By default
 the derivative of p* is approximated (star states), the supersonic branches remain exact. ---*/
template <class HLLC>
void CheckAllBranches(unsigned short nDim, bool moving, su2double k, bool accurate) {
  auto config = MakeConfig(moving, accurate);
  HLLC numerics(nDim, nDim + 2, config.get());

  /*--- Oblique unit normal (face area 1.3) and, on moving faces, a grid velocity with a normal component. ---*/
  const su2double n2[3] = {0.8, 0.6, 0.0}, n3[3] = {2.0 / 7, 3.0 / 7, 6.0 / 7};
  const su2double* n = (nDim == 2) ? n2 : n3;
  su2double normal[3] = {0.0}, gridVel[3] = {0.0};
  for (unsigned short iDim = 0; iDim < nDim; iDim++) normal[iDim] = 1.3 * n[iDim];
  if (moving) {
    gridVel[0] = 0.3;
    gridVel[1] = -0.1;
    gridVel[2] = 0.2;
  }
  numerics.SetNormal(normal);
  numerics.SetGridVel(gridVel, gridVel);
  su2double w = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) w += gridVel[iDim] * n[iDim];

  /*--- Two state pairs with different density, pressure and velocity (also tangential): flow along +n (sM > 0,
   left branches) and along -n (sM < 0, right branches). ---*/
  const su2double tangent[3] = {0.1, -0.2, 0.15};
  for (const su2double dir : {1.0, -1.0}) {
    su2double vel_i[3], vel_j[3];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      vel_i[iDim] = dir * 0.5 * n[iDim] + tangent[iDim];
      vel_j[iDim] = dir * 0.25 * n[iDim] - tangent[iDim];
    }
    const State Ui = Conservative(nDim, 1.2, vel_i, 1.1, k);
    const State Uj = Conservative(nDim, 0.9, vel_j, 0.8, 0.8 * k);

    /*--- Wave speeds relative to the face, held fixed (the Jacobians treat them as constants): the supersonic and
     the star branch on the side of the flow. ---*/
    const su2double supersonicSpeeds[2] = {dir > 0 ? 0.02 : -2.0, dir > 0 ? 2.0 : -0.02};
    const su2double starSpeeds[2] = {-2.0, 2.0};
    for (const auto* s : {supersonicSpeeds, starSpeeds}) {
      numerics.SetFixedWaveSpeeds(true, s[0], s[1]);
      const bool supersonic = (s[0] > 0.0) || (s[1] < 0.0);
      const su2double sM = ContactSpeed(nDim, Ui, Uj, k, 0.8 * k, n, w, s[0], s[1]);
      CAPTURE(nDim, moving, k, accurate, dir, s[0], s[1], sM);

      /*--- Make sure the intended branch is the one selected. ---*/
      if (supersonic)
        REQUIRE(((dir > 0) ? s[0] : -s[1]) > 0.0);
      else
        REQUIRE(sM * dir > 0.0);

      const su2double error = JacobianError(nDim, numerics, config.get(), Ui, Uj, k, 0.8 * k);
      if (accurate || supersonic)
        CHECK(error < 1e-6);
      else
        CHECK(error > 1e-3);  // the approximation is active
    }
  }
}

}  // namespace

TEST_CASE("HLLC Jacobians match finite differences", "[HLLC]") {
  for (const unsigned short nDim : {2, 3}) {
    for (const bool moving : {false, true}) {
      for (const bool accurate : {false, true}) {
        CheckAllBranches<CUpwHLLC_Flow>(nDim, moving, 0.0, accurate);
        CheckAllBranches<CUpwGeneralHLLC_Flow>(nDim, moving, 0.0, accurate);
      }
    }
  }
}
