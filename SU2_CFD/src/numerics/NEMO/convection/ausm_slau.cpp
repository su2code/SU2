/*!
 * \file ausm_slau.cpp
 * \brief Implementations of the AUSM-family of schemes in NEMO.
 * \author F. Palacios, S.R. Copeland, W. Maier, C. Garbacz
 * \version 8.0.0 "Harrier"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

#include "../../../../include/numerics/NEMO/convection/ausm_slau.hpp"

#include "../../../../../Common/include/toolboxes/geometry_toolbox.hpp"

CUpwAUSM_SLAU_Base_NEMO::CUpwAUSM_SLAU_Base_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                                                 unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                                                 const CConfig* config)
    : CNEMONumerics(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {
  if (config->GetDynamic_Grid() && (SU2_MPI::GetRank() == MASTER_NODE))
    cout << "WARNING: Grid velocities are NOT yet considered in AUSM-type schemes." << endl;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);

  rhos_i = new su2double[nSpecies];
  rhos_j = new su2double[nSpecies];

  Fc_L = new su2double[nVar];
  Fc_R = new su2double[nVar];
  dM_LP = new su2double[nVar];
  dM_RM = new su2double[nVar];
  dP_LP = new su2double[nVar];
  dP_RM = new su2double[nVar];
  da_L = new su2double[nVar];
  da_R = new su2double[nVar];

  Flux = new su2double[nVar];
  Jacobian_i = new su2double*[nVar];
  Jacobian_j = new su2double*[nVar];
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double[nVar];
    Jacobian_j[iVar] = new su2double[nVar];
  }
}

CUpwAUSM_SLAU_Base_NEMO::~CUpwAUSM_SLAU_Base_NEMO() {
  delete[] rhos_i;
  delete[] rhos_j;

  delete[] Fc_L;
  delete[] Fc_R;
  delete[] dM_LP;
  delete[] dM_RM;
  delete[] dP_LP;
  delete[] dP_RM;
  delete[] da_L;
  delete[] da_R;

  delete[] Flux;
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    delete[] Jacobian_i[iVar];
    delete[] Jacobian_j[iVar];
  }
  delete[] Jacobian_i;
  delete[] Jacobian_j;
}

void CUpwAUSM_SLAU_Base_NEMO::ComputeInterfaceQuantities(const CConfig* config, su2double* pressure,
                                                         su2double& interface_mach, su2double* interface_soundspeed) {
  /*--- For schemes that fit in the general form of Advection Upstream Splitting Method (AUSM) schemes you can inherit
  from this class and implement only the specifics, which should be the face pressure flux(es), interface Mach number
  and the interface soundspeed(s). For implicit solution methods this class will use the analytic AUSM Jacobians, until
  more variant Jacobians have been added.
  ---*/
}

void CUpwAUSM_SLAU_Base_NEMO::ComputeJacobian(su2double** val_Jacobian_i, su2double** val_Jacobian_j) {
  const auto& Ms = fluidmodel->GetSpeciesMolarMass();
  const auto& Cvtr = fluidmodel->GetSpeciesCvTraRot();
  const su2double Ru = 1000.0 * UNIVERSAL_GAS_CONSTANT;
  rhoCvtr_i = V_i[RHOCVTR_INDEX];
  rhoCvtr_j = V_j[RHOCVTR_INDEX];

  /*--- Initialize the Jacobians ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] = 0.0;
      Jacobian_j[iVar][jVar] = 0.0;
    }
  }

  /*--- Determine proper flux and soundspeed ---*/
  if (M_F >= 0.0)
    Fc_LR = Fc_L;
  else
    Fc_LR = Fc_R;

  su2double A_LR = 0;
  if (M_F >= 0.0)
    A_LR = A_F[0];
  else
    A_LR = A_F[1];

  /*--- Sound speed derivatives: Species density ---*/

  // Electrons only
  for (auto iSpecies = 0ul; iSpecies < nEl; iSpecies++) {
    da_L[iSpecies] = 1.0 / (2.0 * SoundSpeed_i * Density_i) * (1 + dPdU_i[nSpecies + nDim]) *
                     (dPdU_i[iSpecies] - Pressure_i / Density_i);
    da_R[iSpecies] = 1.0 / (2.0 * SoundSpeed_j * Density_j) * (1 + dPdU_j[nSpecies + nDim]) *
                     (dPdU_j[iSpecies] - Pressure_j / Density_j);
  }

  // Heavy species
  for (auto iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
    da_L[iSpecies] =
        1.0 / (2.0 * SoundSpeed_i) *
        (1 / rhoCvtr_i * (Ru / Ms[iSpecies] - Cvtr[iSpecies] * dPdU_i[nSpecies + nDim]) * Pressure_i / Density_i +
         1.0 / Density_i * (1.0 + dPdU_i[nSpecies + nDim]) * (dPdU_i[iSpecies] - Pressure_i / Density_i));
    da_R[iSpecies] =
        1.0 / (2.0 * SoundSpeed_j) *
        (1 / rhoCvtr_j * (Ru / Ms[iSpecies] - Cvtr[iSpecies] * dPdU_j[nSpecies + nDim]) * Pressure_j / Density_j +
         1.0 / Density_j * (1.0 + dPdU_j[nSpecies + nDim]) * (dPdU_j[iSpecies] - Pressure_j / Density_j));
  }

  /*--- Sound speed derivatives: Momentum ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    da_L[nSpecies + iDim] = -1.0 / (2.0 * Density_i * SoundSpeed_i) *
                            ((1.0 + dPdU_i[nSpecies + nDim]) * dPdU_i[nSpecies + nDim]) * Velocity_i[iDim];
    da_R[nSpecies + iDim] = -1.0 / (2.0 * Density_j * SoundSpeed_j) *
                            ((1.0 + dPdU_j[nSpecies + nDim]) * dPdU_j[nSpecies + nDim]) * Velocity_j[iDim];
  }

  /*--- Sound speed derivatives: Energy ---*/
  da_L[nSpecies + nDim] =
      1.0 / (2.0 * Density_i * SoundSpeed_i) * ((1.0 + dPdU_i[nSpecies + nDim]) * dPdU_i[nSpecies + nDim]);
  da_R[nSpecies + nDim] =
      1.0 / (2.0 * Density_j * SoundSpeed_j) * ((1.0 + dPdU_j[nSpecies + nDim]) * dPdU_j[nSpecies + nDim]);

  /*--- Sound speed derivatives: Vib-el energy ---*/
  da_L[nSpecies + nDim + 1] =
      1.0 / (2.0 * Density_i * SoundSpeed_i) * ((1.0 + dPdU_i[nSpecies + nDim]) * dPdU_i[nSpecies + nDim + 1]);
  da_R[nSpecies + nDim + 1] =
      1.0 / (2.0 * Density_j * SoundSpeed_j) * ((1.0 + dPdU_j[nSpecies + nDim]) * dPdU_j[nSpecies + nDim + 1]);

  /*--- Left state Jacobian ---*/
  if (M_F >= 0) {
    /*--- Jacobian contribution: dFc terms ---*/
    for (auto iVar = 0ul; iVar < nSpecies + nDim; iVar++) {
      for (auto jVar = 0ul; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] += M_F * Fc_L[iVar] * da_L[jVar];
      }
      Jacobian_i[iVar][iVar] += M_F * SoundSpeed_i;
    }
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      Jacobian_i[nSpecies + nDim][iSpecies] +=
          M_F * (dPdU_i[iSpecies] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[iSpecies]);
    }
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      Jacobian_i[nSpecies + nDim][nSpecies + iDim] +=
          M_F *
          (-dPdU_i[nSpecies + nDim] * Velocity_i[iDim] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + iDim]);
    }
    Jacobian_i[nSpecies + nDim][nSpecies + nDim] +=
        M_F * ((1.0 + dPdU_i[nSpecies + nDim]) * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + nDim]);
    Jacobian_i[nSpecies + nDim][nSpecies + nDim + 1] +=
        M_F * (dPdU_i[nSpecies + nDim + 1] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + nDim + 1]);
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Jacobian_i[nSpecies + nDim + 1][jVar] += M_F * Fc_L[nSpecies + nDim + 1] * da_L[jVar];
    }
    Jacobian_i[nSpecies + nDim + 1][nSpecies + nDim + 1] += M_F * SoundSpeed_i;
  }

  /*--- Calculate derivatives of the split pressure flux ---*/
  if ((M_F >= 0) || ((M_F < 0) && (fabs(M_F) <= 1.0))) {
    if (fabs(M_L) <= 1.0) {
      /*--- Mach number ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dM_LP[iSpecies] =
            0.5 * (M_L + 1.0) *
            (-ProjVelocity_i / (Density_i * SoundSpeed_i) - ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2)));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dM_LP[nSpecies + iDim] = 0.5 * (M_L + 1.0) *
                                 (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + iDim] +
                                  UnitNormal[iDim] / (Density_i * SoundSpeed_i));
      dM_LP[nSpecies + nDim] = 0.5 * (M_L + 1.0) * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim]);
      dM_LP[nSpecies + nDim + 1] =
          0.5 * (M_L + 1.0) * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim + 1]);

      /*--- Pressure ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dP_LP[iSpecies] = 0.25 * (M_L + 1.0) *
                          (dPdU_i[iSpecies] * (M_L + 1.0) * (2.0 - M_L) +
                           Pressure_i *
                               (-ProjVelocity_i / (Density_i * SoundSpeed_i) -
                                ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2))) *
                               (3.0 - 3.0 * M_L));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dP_LP[nSpecies + iDim] = 0.25 * (M_L + 1.0) *
                                 (-Velocity_i[iDim] * dPdU_i[nSpecies + nDim] * (M_L + 1.0) * (2.0 - M_L) +
                                  Pressure_i *
                                      (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + iDim] +
                                       UnitNormal[iDim] / (Density_i * SoundSpeed_i)) *
                                      (3.0 - 3.0 * M_L));
      dP_LP[nSpecies + nDim] =
          0.25 * (M_L + 1.0) *
          (dPdU_i[nSpecies + nDim] * (M_L + 1.0) * (2.0 - M_L) +
           Pressure_i * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim]) * (3.0 - 3.0 * M_L));
      dP_LP[nSpecies + nDim + 1] =
          0.25 * (M_L + 1.0) *
          (dPdU_i[nSpecies + nDim + 1] * (M_L + 1.0) * (2.0 - M_L) +
           Pressure_i * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim + 1]) * (3.0 - 3.0 * M_L));
    } else {
      /*--- Mach number ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dM_LP[iSpecies] =
            -ProjVelocity_i / (Density_i * SoundSpeed_i) - ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dM_LP[nSpecies + iDim] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + iDim] +
                                 UnitNormal[iDim] / (Density_i * SoundSpeed_i);
      dM_LP[nSpecies + nDim] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim];
      dM_LP[nSpecies + nDim + 1] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim + 1];

      /*--- Pressure ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) dP_LP[iSpecies] = dPdU_i[iSpecies];
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dP_LP[nSpecies + iDim] = (-Velocity_i[iDim] * dPdU_i[nSpecies + nDim]);
      dP_LP[nSpecies + nDim] = dPdU_i[nSpecies + nDim];
      dP_LP[nSpecies + nDim + 1] = dPdU_i[nSpecies + nDim + 1];
    }

    /*--- dM contribution ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
      for (auto jVar = 0ul; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] += dM_LP[jVar] * Fc_LR[iVar] * A_LR;
      }
    }

    /*--- Jacobian contribution: dP terms ---*/
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      for (auto iVar = 0ul; iVar < nVar; iVar++) {
        Jacobian_i[nSpecies + iDim][iVar] += dP_LP[iVar] * UnitNormal[iDim];
      }
    }
  }

  /*--- Right state Jacobian ---*/
  if (M_F < 0) {
    /*--- Jacobian contribution: dFc terms ---*/
    for (auto iVar = 0ul; iVar < nSpecies + nDim; iVar++) {
      for (auto jVar = 0ul; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] += M_F * Fc_R[iVar] * da_R[jVar];
      }
      Jacobian_j[iVar][iVar] += M_F * SoundSpeed_j;
    }
    for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
      Jacobian_j[nSpecies + nDim][iSpecies] +=
          M_F * (dPdU_j[iSpecies] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[iSpecies]);
    }
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      Jacobian_j[nSpecies + nDim][nSpecies + iDim] +=
          M_F *
          (-dPdU_j[nSpecies + nDim] * Velocity_j[iDim] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + iDim]);
    }
    Jacobian_j[nSpecies + nDim][nSpecies + nDim] +=
        M_F * ((1.0 + dPdU_j[nSpecies + nDim]) * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + nDim]);
    Jacobian_j[nSpecies + nDim][nSpecies + nDim + 1] +=
        M_F * (dPdU_j[nSpecies + nDim + 1] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + nDim + 1]);
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Jacobian_j[nSpecies + nDim + 1][jVar] += M_F * Fc_R[nSpecies + nDim + 1] * da_R[jVar];
    }
    Jacobian_j[nSpecies + nDim + 1][nSpecies + nDim + 1] += M_F * SoundSpeed_j;
  }

  /*--- Calculate derivatives of the split pressure flux ---*/
  if ((M_F < 0) || ((M_F >= 0) && (fabs(M_F) <= 1.0))) {
    if (fabs(M_R) <= 1.0) {
      /*--- Mach ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dM_RM[iSpecies] =
            -0.5 * (M_R - 1.0) *
            (-ProjVelocity_j / (Density_j * SoundSpeed_j) - ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2)));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dM_RM[nSpecies + iDim] = -0.5 * (M_R - 1.0) *
                                 (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + iDim] +
                                  UnitNormal[iDim] / (Density_j * SoundSpeed_j));
      dM_RM[nSpecies + nDim] = -0.5 * (M_R - 1.0) * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim]);
      dM_RM[nSpecies + nDim + 1] =
          -0.5 * (M_R - 1.0) * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim + 1]);

      /*--- Pressure ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dP_RM[iSpecies] = 0.25 * (M_R - 1.0) *
                          (dPdU_j[iSpecies] * (M_R - 1.0) * (2.0 + M_R) +
                           Pressure_j *
                               (-ProjVelocity_j / (Density_j * SoundSpeed_j) -
                                ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2))) *
                               (3.0 + 3.0 * M_R));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dP_RM[nSpecies + iDim] = 0.25 * (M_R - 1.0) *
                                 ((-Velocity_j[iDim] * dPdU_j[nSpecies + nDim]) * (M_R - 1.0) * (2.0 + M_R) +
                                  Pressure_j *
                                      (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + iDim] +
                                       UnitNormal[iDim] / (Density_j * SoundSpeed_j)) *
                                      (3.0 + 3.0 * M_R));
      dP_RM[nSpecies + nDim] =
          0.25 * (M_R - 1.0) *
          (dPdU_j[nSpecies + nDim] * (M_R - 1.0) * (2.0 + M_R) +
           Pressure_j * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim]) * (3.0 + 3.0 * M_R));
      dP_RM[nSpecies + nDim + 1] =
          0.25 * (M_R - 1.0) *
          (dPdU_j[nSpecies + nDim + 1] * (M_R - 1.0) * (2.0 + M_R) +
           Pressure_j * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim + 1]) * (3.0 + 3.0 * M_R));

    } else {
      /*--- Mach ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++)
        dM_RM[iSpecies] =
            -ProjVelocity_j / (Density_j * SoundSpeed_j) - ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2));
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dM_RM[nSpecies + iDim] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + iDim] +
                                 UnitNormal[iDim] / (Density_j * SoundSpeed_j);
      dM_RM[nSpecies + nDim] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim];
      dM_RM[nSpecies + nDim + 1] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim + 1];

      /*--- Pressure ---*/
      for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) dP_RM[iSpecies] = dPdU_j[iSpecies];
      for (auto iDim = 0ul; iDim < nDim; iDim++)
        dP_RM[nSpecies + iDim] = -Velocity_j[iDim] * dPdU_j[nSpecies + nDim];
      dP_RM[nSpecies + nDim] = dPdU_j[nSpecies + nDim];
      dP_RM[nSpecies + nDim + 1] = dPdU_j[nSpecies + nDim + 1];
    }

    /*--- Jacobian contribution: dM terms ---*/
    for (auto iVar = 0ul; iVar < nVar; iVar++) {
      for (auto jVar = 0ul; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] += dM_RM[jVar] * Fc_LR[iVar] * A_LR;
      }
    }

    /*--- Jacobian contribution: dP terms ---*/
    for (auto iDim = 0ul; iDim < nDim; iDim++) {
      for (auto iVar = 0ul; iVar < nVar; iVar++) {
        Jacobian_j[nSpecies + iDim][iVar] += dP_RM[iVar] * UnitNormal[iDim];
      }
    }
  }

  /*--- Integrate over dual-face area ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++) {
    for (auto jVar = 0ul; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] *= Area;
      Jacobian_j[iVar][jVar] *= Area;
    }
  }
}

CNumerics::ResidualType<> CUpwAUSM_SLAU_Base_NEMO::ComputeResidual(const CConfig* config) {
  /*--- Compute geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (auto iDim = 0ul; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim] / Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX + iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX + iSpecies];
  }
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[VEL_INDEX + iDim];
    Velocity_j[iDim] = V_j[VEL_INDEX + iDim];
  }

  Pressure_i = V_i[P_INDEX];
  Enthalpy_i = V_i[H_INDEX];
  Density_i = V_i[RHO_INDEX];
  SoundSpeed_i = V_i[A_INDEX];

  Pressure_j = V_j[P_INDEX];
  Enthalpy_j = V_j[H_INDEX];
  Density_j = V_j[RHO_INDEX];
  SoundSpeed_j = V_j[A_INDEX];

  e_ve_i = e_ve_j = 0;
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    e_ve_i += (V_i[RHOS_INDEX + iSpecies] * eve_i[iSpecies]) / Density_i;
    e_ve_j += (V_j[RHOS_INDEX + iSpecies] * eve_j[iSpecies]) / Density_j;
  }

  /*--- Projected velocities ---*/
  ProjVelocity_i = GeometryToolbox::DotProduct(nDim, Velocity_i, UnitNormal);
  ProjVelocity_j = GeometryToolbox::DotProduct(nDim, Velocity_j, UnitNormal);

  /*--- Compute mass and pressure fluxes of specific scheme ---*/
  ComputeInterfaceQuantities(config, PressureFlux, M_F, A_F);

  const su2double MassFlux_i = M_F * A_F[0];
  const su2double MassFlux_j = M_F * A_F[1];

  const su2double DissFlux_i = fabs(MassFlux_i);
  const su2double DissFlux_j = fabs(MassFlux_j);

  /*--- Assign left & right convective flux vectors ---*/
  for (auto iSpecies = 0ul; iSpecies < nSpecies; iSpecies++) {
    Fc_L[iSpecies] = rhos_i[iSpecies];
    Fc_R[iSpecies] = rhos_j[iSpecies];
  }
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    Fc_L[nSpecies + iDim] = Density_i * Velocity_i[iDim];
    Fc_R[nSpecies + iDim] = Density_j * Velocity_j[iDim];
  }
  Fc_L[nSpecies + nDim] = Density_i * Enthalpy_i;
  Fc_R[nSpecies + nDim] = Density_j * Enthalpy_j;
  Fc_L[nSpecies + nDim + 1] = Density_i * e_ve_i;
  Fc_R[nSpecies + nDim + 1] = Density_j * e_ve_j;

  /*--- Compute numerical flux ---*/
  for (auto iVar = 0ul; iVar < nVar; iVar++)
    Flux[iVar] = 0.5 * ((MassFlux_i + DissFlux_i) * Fc_L[iVar] + (MassFlux_j - DissFlux_j) * Fc_R[iVar]) * Area;

  for (auto iDim = 0ul; iDim < nDim; iDim++) Flux[nSpecies + iDim] += PressureFlux[iDim] * Area;

  /*--- If required, compute Jacobians (approximated using AUSM) ---*/
  if (implicit) ComputeJacobian(Jacobian_i, Jacobian_j);

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

CUpwAUSM_NEMO::CUpwAUSM_NEMO(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                             unsigned short val_nPrimVarGrad, const CConfig* config)
    : CUpwAUSM_SLAU_Base_NEMO(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {}

void CUpwAUSM_NEMO::ComputeInterfaceQuantities(const CConfig* config, su2double* pressure, su2double& interface_mach,
                                               su2double* interface_soundspeed) {
  /*--- Calculate L/R Mach numbers ---*/
  interface_soundspeed[0] = SoundSpeed_i;
  interface_soundspeed[1] = SoundSpeed_j;
  M_L = ProjVelocity_i / interface_soundspeed[0];
  M_R = ProjVelocity_j / interface_soundspeed[1];

  /*--- Calculate split numerical fluxes ---*/
  su2double M_LP, M_RM, P_LP, P_RM;

  if (fabs(M_L) <= 1.0) {
    M_LP = 0.25 * (M_L + 1.0) * (M_L + 1.0);
    P_LP = 0.25 * Pressure_i * (M_L + 1.0) * (M_L + 1.0) * (2.0 - M_L);
  } else {
    M_LP = 0.5 * (M_L + fabs(M_L));
    P_LP = 0.5 * Pressure_i * (M_L + fabs(M_L)) / M_L;
  }

  if (fabs(M_R) <= 1.0) {
    M_RM = -0.25 * (M_R - 1.0) * (M_R - 1.0);
    P_RM = 0.25 * Pressure_j * (M_R - 1.0) * (M_R - 1.0) * (2.0 + M_R);
  } else {
    M_RM = 0.5 * (M_R - fabs(M_R));
    P_RM = 0.5 * Pressure_j * (M_R - fabs(M_R)) / M_R;
  }

  // M_Lnferface = M(1/2) = (M_LP + M_RM)
  interface_mach = (M_LP + M_RM);

  // Split pressure P(1/2)
  // Note: only a single pressure flux in AUSM
  for (auto iDim = 0ul; iDim < nDim; iDim++) pressure[iDim] = (P_LP + P_RM) * UnitNormal[iDim];
}

CUpwAUSMPLUSUP2_NEMO::CUpwAUSMPLUSUP2_NEMO(unsigned short val_nDim, unsigned short val_nVar,
                                           unsigned short val_nPrimVar, unsigned short val_nPrimVarGrad,
                                           const CConfig* config)
    : CUpwAUSM_SLAU_Base_NEMO(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {
  Minf = config->GetMach();
  Kp = 0.25;
  sigma = 1.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+-Up2 requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

void CUpwAUSMPLUSUP2_NEMO::ComputeInterfaceQuantities(const CConfig* config, su2double* pressure,
                                                      su2double& interface_mach, su2double* interface_soundspeed) {
  const su2double sq_veli = GeometryToolbox::SquaredNorm(nDim, Velocity_i);
  const su2double sq_velj = GeometryToolbox::SquaredNorm(nDim, Velocity_j);

  /*--- Compute C*  ---*/
  const su2double CstarL = sqrt(2.0 * (Gamma_i - 1.0) / (Gamma_i + 1.0) * Enthalpy_i);
  const su2double CstarR = sqrt(2.0 * (Gamma_j - 1.0) / (Gamma_j + 1.0) * Enthalpy_j);

  /*--- Compute C^ ---*/
  const su2double ChatL = CstarL * CstarL / max(CstarL, ProjVelocity_i);
  const su2double ChatR = CstarR * CstarR / max(CstarR, -ProjVelocity_j);

  /*--- Interface speed of sound ---*/
  const su2double aF = min(ChatL, ChatR);
  interface_soundspeed[0] = interface_soundspeed[1] = aF;

  const su2double M_L = ProjVelocity_i / aF;
  const su2double M_R = ProjVelocity_j / aF;

  const su2double rhoF = 0.5 * (Density_i + Density_j);
  const su2double MFsq = 0.5 * (M_L * M_L + M_R * M_R);

  const su2double param1 = max(MFsq, Minf * Minf);
  const su2double Mrefsq = (min(1.0, param1));
  const su2double fa = 2.0 * sqrt(Mrefsq) - Mrefsq;

  const su2double alpha = 3.0 / 16.0 * (-4.0 + 5.0 * fa * fa);
  const su2double beta = 1.0 / 8.0;

  /*--- Pressure diffusion term ---*/
  const su2double Mp = -(Kp / fa) * max((1.0 - sigma * MFsq), 0.0) * (Pressure_j - Pressure_i) / (rhoF * aF * aF);

  su2double M_LP, P_LP, M_RM, P_RM;

  if (fabs(M_L) <= 1.0) {
    M_LP = 0.25 * (M_L + 1.0) * (M_L + 1.0) + beta * (M_L * M_L - 1.0) * (M_L * M_L - 1.0);
    P_LP = (0.25 * (M_L + 1.0) * (M_L + 1.0) * (2.0 - M_L) + alpha * M_L * (M_L * M_L - 1.0) * (M_L * M_L - 1.0));
  } else {
    M_LP = 0.5 * (M_L + fabs(M_L));
    P_LP = 0.5 * (M_L + fabs(M_L)) / M_L;
  }

  if (fabs(M_R) <= 1.0) {
    M_RM = -0.25 * (M_R - 1.0) * (M_R - 1.0) - beta * (M_R * M_R - 1.0) * (M_R * M_R - 1.0);
    P_RM = (0.25 * (M_R - 1.0) * (M_R - 1.0) * (2.0 + M_R) - alpha * M_R * (M_R * M_R - 1.0) * (M_R * M_R - 1.0));
  } else {
    M_RM = 0.5 * (M_R - fabs(M_R));
    P_RM = 0.5 * (M_R - fabs(M_R)) / M_R;
  }

  /*--- Interface Mach number ---*/
  interface_mach = (M_LP + M_RM + Mp);

  /*--- Modified pressure flux ---*/
  const su2double pFi = sqrt(0.5 * (sq_veli + sq_velj)) * (P_LP + P_RM - 1.0) * 0.5 * (Density_j + Density_i) * aF;

  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    pressure[iDim] =
        (0.5 * (Pressure_j + Pressure_i) + 0.5 * (P_LP - P_RM) * (Pressure_i - Pressure_j) + pFi) * UnitNormal[iDim];
  }
}

CUpwAUSMPLUSM_NEMO::CUpwAUSMPLUSM_NEMO(unsigned short val_nDim, unsigned short val_nVar, unsigned short val_nPrimVar,
                                       unsigned short val_nPrimVarGrad, const CConfig* config)
    : CUpwAUSM_SLAU_Base_NEMO(val_nDim, val_nVar, val_nPrimVar, val_nPrimVarGrad, config) {
  Minf = config->GetMach();
  beta = 1.0 / 8.0;

  if (Minf < EPS)
    SU2_MPI::Error("AUSM+M requires a reference Mach number (\"MACH_NUMBER\") greater than 0.", CURRENT_FUNCTION);
}

void CUpwAUSMPLUSM_NEMO::ComputeInterfaceQuantities(const CConfig* config, su2double* pressure,
                                                    su2double& interface_mach, su2double* interface_soundspeed) {
  const su2double sq_veli = GeometryToolbox::SquaredNorm(nDim, Velocity_i);
  const su2double sq_velj = GeometryToolbox::SquaredNorm(nDim, Velocity_j);

  /*--- Calculate interface numerical gammas and speed of sound ---*/
  const su2double Enthalpy_norm = 0.5 * (Enthalpy_i + Enthalpy_j);
  const su2double Gamma_ij = 0.5 * (Gamma_i + Gamma_j);
  const su2double A_s = sqrt(2.0 * Enthalpy_norm * (Gamma_ij - 1.0) / (Gamma_ij + 1.0));

  su2double A_F;
  if (0.5 * (ProjVelocity_i + ProjVelocity_j) >= 0.0)
    A_F = A_s * A_s / max(fabs(ProjVelocity_i), A_s);
  else
    A_F = A_s * A_s / max(fabs(ProjVelocity_j), A_s);

  /*--- Compute L/R Mach numbers ---*/
  interface_soundspeed[0] = interface_soundspeed[1] = A_F;
  M_L = ProjVelocity_i / interface_soundspeed[0];
  M_R = ProjVelocity_j / interface_soundspeed[1];

  /*--- Interface mach number w/ pressure diffusion term (M_p)---*/
  const su2double Density_F = 0.5 * (Density_i + Density_j);
  const su2double MF_sq = 0.5 * (sq_veli + sq_velj) / (A_F * A_F);

  const su2double param1 = max(MF_sq, Minf * Minf);
  const su2double Mrefsq = (min(1.0, param1));
  const su2double fa = 2.0 * sqrt(Mrefsq) - Mrefsq;

  const su2double alpha = 3.0 / 16.0 * (-4.0 + 5.0 * fa * fa);
  const su2double f = 0.5 * (1 - cos(PI_NUMBER * min(1.0, max(abs(M_L), abs(M_R)))));

  /*--- Pressure sensor terms ---*/
  const su2double h = min(Sensor_i, Sensor_j);
  const su2double g = 0.5 * (1 + cos(PI_NUMBER * h));
  const su2double f0 = min(1.0, max(f, Minf * Minf));

  /*--- Pressure diffusion term ---*/
  const su2double M_p = -0.5 * (1.0 - f) * (Pressure_j - Pressure_i) / (Density_F * A_F * A_F) * (1.0 - g);

  /*--- Compute base split mach and pressure fluxes ---*/
  su2double M_LP, M_RM, P_LP, P_RM;

  if (fabs(M_L) <= 1.0) {
    M_LP = 0.25 * (M_L + 1.0) * (M_L + 1.0) + beta * (M_L * M_L - 1.0) * (M_L * M_L - 1.0);
    P_LP = (0.25 * (M_L + 1.0) * (M_L + 1.0) * (2.0 - M_L) + alpha * M_L * (M_L * M_L - 1.0) * (M_L * M_L - 1.0));
  } else {
    M_LP = 0.5 * (M_L + fabs(M_L));
    P_LP = 0.5 * (M_L + fabs(M_L)) / M_L;
  }

  if (fabs(M_R) <= 1.0) {
    M_RM = -0.25 * (M_R - 1.0) * (M_R - 1.0) - beta * (M_R * M_R - 1.0) * (M_R * M_R - 1.0);
    P_RM = (0.25 * (M_R - 1.0) * (M_R - 1.0) * (2.0 + M_R) - alpha * M_R * (M_R * M_R - 1.0) * (M_R * M_R - 1.0));
  } else {
    M_RM = 0.5 * (M_R - fabs(M_R));
    P_RM = 0.5 * (M_R - fabs(M_R)) / M_R;
  }

  interface_mach = M_LP + M_RM + M_p;

  /*--- Compute and add pressure sensor term to pressure flux ---*/
  const su2double pFi = f0 * (P_LP + P_RM - 1.0) * 0.5 * (Pressure_i + Pressure_j);

  /*--- Velocity diffusion term---*/
  su2double P_un[MAXNDIM] = {0.0};

  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    su2double Vel_L = (Velocity_i[iDim] + ProjVelocity_i * UnitNormal[iDim]);
    su2double Vel_R = (Velocity_j[iDim] - ProjVelocity_j * UnitNormal[iDim]);

    P_un[iDim] = -g * Gamma_ij * 0.5 * (Pressure_i + Pressure_j) / A_F * P_LP * P_RM * (Vel_R - Vel_L);
  }

  /*--- Pressure flux ---*/
  for (auto iDim = 0ul; iDim < nDim; iDim++) {
    pressure[iDim] =
        (0.5 * (Pressure_j + Pressure_i) + 0.5 * (P_LP - P_RM) * (Pressure_i - Pressure_j) + pFi) * UnitNormal[iDim] +
        P_un[iDim];
  }
}
