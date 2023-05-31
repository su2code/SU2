/*!
 * \file ausm_slau.cpp
 * \brief Implementations of the AUSM-family of schemes in NEMO.
 * \author F. Palacios, S.R. Copeland, W. Maier, C. Garbacz
 * \version 7.5.1 "Blackbird"
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
#include <iomanip>

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
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Jacobian_i[iVar] = new su2double[nVar];
    Jacobian_j[iVar] = new su2double[nVar];
  }
}

CUpwAUSM_SLAU_Base_NEMO::~CUpwAUSM_SLAU_Base_NEMO(void) {
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
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
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
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
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
  for (unsigned short iSpecies = 0; iSpecies < nEl; iSpecies++) {
    da_L[iSpecies] = 1.0 / (2.0 * SoundSpeed_i * Density_i) * (1 + dPdU_i[nSpecies + nDim]) *
                     (dPdU_i[iSpecies] - Pressure_i / Density_i);
    da_R[iSpecies] = 1.0 / (2.0 * SoundSpeed_j * Density_j) * (1 + dPdU_j[nSpecies + nDim]) *
                     (dPdU_j[iSpecies] - Pressure_j / Density_j);
  }

  // Heavy species
  for (unsigned short iSpecies = nEl; iSpecies < nSpecies; iSpecies++) {
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
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
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
    for (auto iVar = 0u; iVar < nSpecies + nDim; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] += M_F * Fc_L[iVar] * da_L[jVar];
      }
      Jacobian_i[iVar][iVar] += M_F * SoundSpeed_i;
    }
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Jacobian_i[nSpecies + nDim][iSpecies] +=
          M_F * (dPdU_i[iSpecies] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[iSpecies]);
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Jacobian_i[nSpecies + nDim][nSpecies + iDim] +=
          M_F *
          (-dPdU_i[nSpecies + nDim] * Velocity_i[iDim] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + iDim]);
    }
    Jacobian_i[nSpecies + nDim][nSpecies + nDim] +=
        M_F * ((1.0 + dPdU_i[nSpecies + nDim]) * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + nDim]);
    Jacobian_i[nSpecies + nDim][nSpecies + nDim + 1] +=
        M_F * (dPdU_i[nSpecies + nDim + 1] * SoundSpeed_i + Density_i * Enthalpy_i * da_L[nSpecies + nDim + 1]);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[nSpecies + nDim + 1][jVar] += M_F * Fc_L[nSpecies + nDim + 1] * da_L[jVar];
    }
    Jacobian_i[nSpecies + nDim + 1][nSpecies + nDim + 1] += M_F * SoundSpeed_i;
  }

  /*--- Calculate derivatives of the split pressure flux ---*/
  if ((M_F >= 0) || ((M_F < 0) && (fabs(M_F) <= 1.0))) {
    if (fabs(M_L) <= 1.0) {
      /*--- Mach number ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dM_LP[iSpecies] =
            0.5 * (M_L + 1.0) *
            (-ProjVelocity_i / (Density_i * SoundSpeed_i) - ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2)));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dM_LP[nSpecies + iDim] = 0.5 * (M_L + 1.0) *
                                 (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + iDim] +
                                  UnitNormal[iDim] / (Density_i * SoundSpeed_i));
      dM_LP[nSpecies + nDim] = 0.5 * (M_L + 1.0) * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim]);
      dM_LP[nSpecies + nDim + 1] =
          0.5 * (M_L + 1.0) * (-ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim + 1]);

      /*--- Pressure ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dP_LP[iSpecies] = 0.25 * (M_L + 1.0) *
                          (dPdU_i[iSpecies] * (M_L + 1.0) * (2.0 - M_L) +
                           Pressure_i *
                               (-ProjVelocity_i / (Density_i * SoundSpeed_i) -
                                ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2))) *
                               (3.0 - 3.0 * M_L));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
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
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dM_LP[iSpecies] =
            -ProjVelocity_i / (Density_i * SoundSpeed_i) - ProjVelocity_i * da_L[iSpecies] / (pow(SoundSpeed_i, 2));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dM_LP[nSpecies + iDim] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + iDim] +
                                 UnitNormal[iDim] / (Density_i * SoundSpeed_i);
      dM_LP[nSpecies + nDim] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim];
      dM_LP[nSpecies + nDim + 1] = -ProjVelocity_i / (pow(SoundSpeed_i, 2)) * da_L[nSpecies + nDim + 1];

      /*--- Pressure ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) dP_LP[iSpecies] = dPdU_i[iSpecies];
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dP_LP[nSpecies + iDim] = (-Velocity_i[iDim] * dPdU_i[nSpecies + nDim]);
      dP_LP[nSpecies + nDim] = dPdU_i[nSpecies + nDim];
      dP_LP[nSpecies + nDim + 1] = dPdU_i[nSpecies + nDim + 1];
    }

    /*--- dM contribution ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        Jacobian_i[iVar][jVar] += dM_LP[jVar] * Fc_LR[iVar] * A_LR;
      }
    }

    /*--- Jacobian contribution: dP terms ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Jacobian_i[nSpecies + iDim][iVar] += dP_LP[iVar] * UnitNormal[iDim];
      }
    }
  }

  /*--- Right state Jacobian ---*/
  if (M_F < 0) {
    /*--- Jacobian contribution: dFc terms ---*/
    for (auto iVar = 0u; iVar < nSpecies + nDim; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] += M_F * Fc_R[iVar] * da_R[jVar];
      }
      Jacobian_j[iVar][iVar] += M_F * SoundSpeed_j;
    }
    for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
      Jacobian_j[nSpecies + nDim][iSpecies] +=
          M_F * (dPdU_j[iSpecies] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[iSpecies]);
    }
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Jacobian_j[nSpecies + nDim][nSpecies + iDim] +=
          M_F *
          (-dPdU_j[nSpecies + nDim] * Velocity_j[iDim] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + iDim]);
    }
    Jacobian_j[nSpecies + nDim][nSpecies + nDim] +=
        M_F * ((1.0 + dPdU_j[nSpecies + nDim]) * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + nDim]);
    Jacobian_j[nSpecies + nDim][nSpecies + nDim + 1] +=
        M_F * (dPdU_j[nSpecies + nDim + 1] * SoundSpeed_j + Density_j * Enthalpy_j * da_R[nSpecies + nDim + 1]);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      Jacobian_j[nSpecies + nDim + 1][jVar] += M_F * Fc_R[nSpecies + nDim + 1] * da_R[jVar];
    }
    Jacobian_j[nSpecies + nDim + 1][nSpecies + nDim + 1] += M_F * SoundSpeed_j;
  }

  /*--- Calculate derivatives of the split pressure flux ---*/
  if ((M_F < 0) || ((M_F >= 0) && (fabs(M_F) <= 1.0))) {
    if (fabs(M_R) <= 1.0) {
      /*--- Mach ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dM_RM[iSpecies] =
            -0.5 * (M_R - 1.0) *
            (-ProjVelocity_j / (Density_j * SoundSpeed_j) - ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2)));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dM_RM[nSpecies + iDim] = -0.5 * (M_R - 1.0) *
                                 (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + iDim] +
                                  UnitNormal[iDim] / (Density_j * SoundSpeed_j));
      dM_RM[nSpecies + nDim] = -0.5 * (M_R - 1.0) * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim]);
      dM_RM[nSpecies + nDim + 1] =
          -0.5 * (M_R - 1.0) * (-ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim + 1]);

      /*--- Pressure ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dP_RM[iSpecies] = 0.25 * (M_R - 1.0) *
                          (dPdU_j[iSpecies] * (M_R - 1.0) * (2.0 + M_R) +
                           Pressure_j *
                               (-ProjVelocity_j / (Density_j * SoundSpeed_j) -
                                ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2))) *
                               (3.0 + 3.0 * M_R));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
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
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++)
        dM_RM[iSpecies] =
            -ProjVelocity_j / (Density_j * SoundSpeed_j) - ProjVelocity_j * da_R[iSpecies] / (pow(SoundSpeed_j, 2));
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dM_RM[nSpecies + iDim] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + iDim] +
                                 UnitNormal[iDim] / (Density_j * SoundSpeed_j);
      dM_RM[nSpecies + nDim] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim];
      dM_RM[nSpecies + nDim + 1] = -ProjVelocity_j / (pow(SoundSpeed_j, 2)) * da_R[nSpecies + nDim + 1];

      /*--- Pressure ---*/
      for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) dP_RM[iSpecies] = dPdU_j[iSpecies];
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        dP_RM[nSpecies + iDim] = -Velocity_j[iDim] * dPdU_j[nSpecies + nDim];
      dP_RM[nSpecies + nDim] = dPdU_j[nSpecies + nDim];
      dP_RM[nSpecies + nDim + 1] = dPdU_j[nSpecies + nDim + 1];
    }

    /*--- Jacobian contribution: dM terms ---*/
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        Jacobian_j[iVar][jVar] += dM_RM[jVar] * Fc_LR[iVar] * A_LR;
      }
    }

    /*--- Jacobian contribution: dP terms ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Jacobian_j[nSpecies + iDim][iVar] += dP_RM[iVar] * UnitNormal[iDim];
      }
    }
  }

  /*--- Integrate over dual-face area ---*/
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      Jacobian_i[iVar][jVar] *= Area;
      Jacobian_j[iVar][jVar] *= Area;
    }
  }
}

CNumerics::ResidualType<> CUpwAUSM_SLAU_Base_NEMO::ComputeResidual(const CConfig* config) {
  /*--- Compute geometric quantities ---*/
  Area = GeometryToolbox::Norm(nDim, Normal);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) UnitNormal[iDim] = Normal[iDim] / Area;

  /*--- Pull stored primitive variables ---*/
  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    rhos_i[iSpecies] = V_i[RHOS_INDEX + iSpecies];
    rhos_j[iSpecies] = V_j[RHOS_INDEX + iSpecies];
  }
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
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

  e_ve_i = 0;
  e_ve_j = 0;
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    //std::cout << "V_i[RHOS_INDEX + iSpecies]=" << V_i[RHOS_INDEX + iSpecies] << std::endl;
    //std::cout << "V_j[RHOS_INDEX + iSpecies]=" << V_j[RHOS_INDEX + iSpecies] << std::endl;
    //std::cout << "eve_i[iSpecies]=" << eve_i[iSpecies] << std::endl;
    //std::cout << "eve_j[iSpecies]=" << eve_j[iSpecies] << std::endl;
    //std::cout << "Density_i=" << Density_i << std::endl;
    //std::cout << "Density_j=" << Density_j << std::endl;
    e_ve_i += (V_i[RHOS_INDEX + iSpecies] * eve_i[iSpecies]) / Density_i;
    e_ve_j += (V_j[RHOS_INDEX + iSpecies] * eve_j[iSpecies]) / Density_j;
  }

  //for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
  //  std::cout << std::setprecision(20) << "rhos_i[iSpecies]=" << rhos_i[iSpecies] << std::endl;
  //  std::cout << std::setprecision(20) << "rhos_j[iSpecies]=" << rhos_j[iSpecies] << std::endl;
  //}
  //for (unsigned short iDim = 0; iDim < nDim; iDim++) {
  //  std::cout << std::setprecision(20) << "u_i[iDim]=" << Velocity_i[iDim] << std::endl;
  //  std::cout << std::setprecision(20) << "u_j[iDim]=" << Velocity_j[iDim] << std::endl;
  //}
  //std::cout << std::setprecision(20) << "P_i=" << Pressure_i << std::endl;
  //std::cout << std::setprecision(20) << "P_j=" << Pressure_j << std::endl;
  //std::cout << std::setprecision(20) << "h_i=" << Enthalpy_i << std::endl;
  //std::cout << std::setprecision(20) << "h_j=" << Enthalpy_j << std::endl;
  //std::cout << std::setprecision(20) << "a_i=" << SoundSpeed_i << std::endl;
  //std::cout << std::setprecision(20) << "a_j=" << SoundSpeed_j << std::endl;
  //std::cout << std::setprecision(20) << "eve_i=" << e_ve_i << std::endl;
  //std::cout << std::setprecision(20) << "eve_j=" << e_ve_j << std::endl;  

  /*--- Projected velocities ---*/
  ProjVelocity_i = GeometryToolbox::DotProduct(nDim, Velocity_i, UnitNormal);
  ProjVelocity_j = GeometryToolbox::DotProduct(nDim, Velocity_j, UnitNormal);

  /*--- Compute mass and pressure fluxes of specific scheme ---*/
  ComputeInterfaceQuantities(config, PressureFlux, M_F, A_F);

  su2double MassFlux_i = M_F * A_F[0];
  su2double MassFlux_j = M_F * A_F[1];

  su2double DissFlux_i = fabs(MassFlux_i);
  su2double DissFlux_j = fabs(MassFlux_j);

  /*--- Assign left & right convective flux vectors ---*/
  for (unsigned short iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
    Fc_L[iSpecies] = rhos_i[iSpecies];
    Fc_R[iSpecies] = rhos_j[iSpecies];
  }
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Fc_L[nSpecies + iDim] = Density_i* Velocity_i[iDim];
    Fc_R[nSpecies + iDim] = Density_j* Velocity_j[iDim];
  }
  Fc_L[nSpecies + nDim] =     Density_i* Enthalpy_i;
  Fc_R[nSpecies + nDim] =     Density_j* Enthalpy_j;
  Fc_L[nSpecies + nDim + 1] = Density_i* e_ve_i;
  Fc_R[nSpecies + nDim + 1] = Density_j* e_ve_j;

  /*--- Compute numerical flux ---*/
  for (unsigned short iVar = 0; iVar < nVar; iVar++){
  //  std::cout << std::setprecision(20) << "MassFlux_i=" << MassFlux_i << std::endl;
  //  std::cout << std::setprecision(20) << "DissFlux_i=" << DissFlux_i << std::endl;
  //  std::cout << std::setprecision(20) << "Fc_L[" << iVar <<"]=" << Fc_L[iVar] << std::endl;
  //  std::cout << std::setprecision(20) << "MassFlux_j=" << MassFlux_j << std::endl;
  //  std::cout << std::setprecision(20) << "DissFlux_j=" << DissFlux_j << std::endl;
  //  std::cout << std::setprecision(20) << "Fc_R[" << iVar <<"]=" << Fc_R[iVar] << std::endl;
  //  std::cout << std::setprecision(20) << "Area=" << Area << std::endl;
    Flux[iVar] = 0.5 * ((MassFlux_i + DissFlux_i) * Fc_L[iVar] + (MassFlux_j - DissFlux_j) * Fc_R[iVar]) * Area;
  //  std::cout << std::setprecision(20) << "Flux[" << iVar << "]=" << Flux[iVar] << std::endl;
  }
  
  //exit(0);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) Flux[nSpecies + iDim] += PressureFlux[iDim] * Area;

  /*--- If required, compute Jacobians (approximated using AUSM) ---*/
  if (implicit) ComputeJacobian(Jacobian_i, Jacobian_j);

  return ResidualType<>(Flux, Jacobian_i, Jacobian_j);
}

//CNumerics::ResidualType<> CUpwAUSM_SLAU_Base_NEMO::ComputeResidual(const CConfig *config) {
//
//  unsigned short iDim, iVar, iSpecies;
//  su2double Density_i, Density_j, 
//  e_ve_i, e_ve_j, mL, mR, mLP, mRM, mF, pLP, pRM, pF, Phi;
//
//  /*--- Compute geometric quantities ---*/
//  Area = GeometryToolbox::Norm(nDim, Normal);
//
//  for (iDim = 0; iDim < nDim; iDim++)
//    UnitNormal[iDim] = Normal[iDim]/Area;
//
//  /*--- Pull stored primitive variables ---*/
//  // Primitives: [rho1,...,rhoNs, T, Tve, u, v, w, P, rho, h, a, c]
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    rhos_i[iSpecies] = V_i[RHOS_INDEX+iSpecies];
//    rhos_j[iSpecies] = V_j[RHOS_INDEX+iSpecies];
//  }
//  for (iDim = 0; iDim < nDim; iDim++) {
//    Velocity_i[iDim] = V_i[VEL_INDEX+iDim];
//    Velocity_j[iDim] = V_j[VEL_INDEX+iDim];
//  }
//
//  Pressure_i   = V_i[P_INDEX];   Pressure_j   = V_j[P_INDEX];
//  Enthalpy_i   = V_i[H_INDEX];   Enthalpy_j   = V_j[H_INDEX];
//  SoundSpeed_i   = V_i[A_INDEX];   SoundSpeed_j   = V_j[A_INDEX];
//  Density_i = V_i[RHO_INDEX]; Density_j = V_j[RHO_INDEX];
//  
//  e_ve_i  = 0; e_ve_j  = 0;
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    //std::cout << "eve_i[iSpecies]=" << eve_i[iSpecies] << std::endl;
//    e_ve_i += (V_i[RHOS_INDEX+iSpecies]*eve_i[iSpecies])/Density_i;
//    e_ve_j += (V_j[RHOS_INDEX+iSpecies]*eve_j[iSpecies])/Density_j;
//  }
//
//  /*--- Projected velocities ---*/
//  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0;
//  for (iDim = 0; iDim < nDim; iDim++) {
//    ProjVelocity_i += Velocity_i[iDim]*UnitNormal[iDim];
//    ProjVelocity_j += Velocity_j[iDim]*UnitNormal[iDim];
//  }
//
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    std::cout << std::setprecision(20)<< "rhos_i[iSpecies]=" << rhos_i[iSpecies] << std::endl;
//    std::cout << std::setprecision(20)<< "rhos_j[iSpecies]=" << rhos_j[iSpecies] << std::endl;
//  }
//  for (iDim = 0; iDim < nDim; iDim++) {
//    std::cout << std::setprecision(20)<< "u_i[iDim]=" << Velocity_i[iDim] << std::endl;
//    std::cout << std::setprecision(20)<< "u_j[iDim]=" << Velocity_j[iDim] << std::endl;
//  }
//  std::cout << std::setprecision(20) << "P_i=" << Pressure_i << std::endl;
//  std::cout << std::setprecision(20) << "P_j=" << Pressure_j << std::endl;
//  std::cout << std::setprecision(20) << "h_i=" << Enthalpy_i << std::endl;
//  std::cout << std::setprecision(20) << "h_j=" << Enthalpy_j << std::endl;
//  std::cout << std::setprecision(20) << "a_i=" << SoundSpeed_i << std::endl;
//  std::cout << std::setprecision(20) << "a_j=" << SoundSpeed_j << std::endl;
//  std::cout << std::setprecision(20) << "eve_i=" << e_ve_i << std::endl;
//  std::cout << std::setprecision(20) << "eve_j=" << e_ve_j << std::endl;
//
//  /*--- Calculate L/R Mach numbers ---*/
//  mL = ProjVelocity_i/SoundSpeed_i;
//  mR = ProjVelocity_j/SoundSpeed_j;
//
//  /*--- Calculate split numerical fluxes ---*/
//  if (fabs(mL) <= 1.0) mLP = 0.25*(mL+1.0)*(mL+1.0);
//  else                 mLP = 0.5*(mL+fabs(mL));
//
//  if (fabs(mR) <= 1.0) mRM = -0.25*(mR-1.0)*(mR-1.0);
//  else                 mRM = 0.5*(mR-fabs(mR));
//
//  mF = mLP + mRM;
//
//  if (fabs(mL) <= 1.0) pLP = 0.25*Pressure_i*(mL+1.0)*(mL+1.0)*(2.0-mL);
//  else                 pLP = 0.5*Pressure_i*(mL+fabs(mL))/mL;
//
//  if (fabs(mR) <= 1.0) pRM = 0.25*Pressure_j*(mR-1.0)*(mR-1.0)*(2.0+mR);
//  else                 pRM = 0.5*Pressure_j*(mR-fabs(mR))/mR;
//
//  pF = pLP + pRM;
//  Phi = fabs(mF);
//
//  /*--- Assign left & right convective vectors ---*/
//  for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
//    Fc_L[iSpecies] = rhos_i[iSpecies]*SoundSpeed_i;
//    Fc_R[iSpecies] = rhos_j[iSpecies]*SoundSpeed_j;
//  }
//  for (iDim = 0; iDim < nDim; iDim++) {
//    Fc_L[nSpecies+iDim] = Density_i*SoundSpeed_i*Velocity_i[iDim];
//    Fc_R[nSpecies+iDim] = Density_j*SoundSpeed_j*Velocity_j[iDim];
//  }
//  Fc_L[nSpecies+nDim]   = Density_i*SoundSpeed_i*Enthalpy_i;
//  Fc_R[nSpecies+nDim]   = Density_j*SoundSpeed_j*Enthalpy_j;
//  Fc_L[nSpecies+nDim+1] = Density_i*SoundSpeed_i*e_ve_i;
//  Fc_R[nSpecies+nDim+1] = Density_j*SoundSpeed_j*e_ve_j;
//
//  /*--- Compute numerical flux ---*/
//  for (iVar = 0; iVar < nVar; iVar++){
//    Flux[iVar] = 0.5*((mF+Phi)*Fc_L[iVar]+(mF-Phi)*Fc_R[iVar])*Area;
//    std::cout << std::setprecision(20)<< "Flux[" << iVar << "]=" << Flux[iVar] << std::endl;
//  }
//  
// // exit(0);
//
//
//  for (iDim = 0; iDim < nDim; iDim++)
//    Flux[nSpecies+iDim] += pF*UnitNormal[iDim]*Area;
//
////  if (implicit)
//
////    /*--- Initialize the Jacobians ---*/
////    for (iVar = 0; iVar < nVar; iVar++) {
////      for (jVar = 0; jVar < nVar; jVar++) {
////        val_Jacobian_i[iVar][jVar] = 0.0;
////        val_Jacobian_j[iVar][jVar] = 0.0;
////      }
////    }
////
////    if (mF >= 0.0) FcLR = FcL;
////    else           FcLR = FcR;
////
////    /*--- Sound speed derivatives: Species density ---*/
////    for (iSpecies = 0; iSpecies < nHeavy; iSpecies++) {
////      Cvtrs = (3.0/2.0+xi[iSpecies]/2.0)*Ru/Ms[iSpecies];
////      daL[iSpecies] = 1.0/(2.0*a_i) * (1/rhoCvtr_i*(Ru/Ms[iSpecies] - Cvtrs*dPdU_i[nSpecies+nDim])*P_i/rho_i
////          + 1.0/rho_i*(1.0+dPdU_i[nSpecies+nDim])*(dPdU_i[iSpecies] - P_i/rho_i));
////      daR[iSpecies] = 1.0/(2.0*a_j) * (1/rhoCvtr_j*(Ru/Ms[iSpecies] - Cvtrs*dPdU_j[nSpecies+nDim])*P_j/rho_j
////          + 1.0/rho_j*(1.0+dPdU_j[nSpecies+nDim])*(dPdU_j[iSpecies] - P_j/rho_j));
////    }
////    for (iSpecies = 0; iSpecies < nEl; iSpecies++) {
////      daL[nSpecies-1] = 1.0/(2.0*a_i*rho_i) * (1+dPdU_i[nSpecies+nDim])*(dPdU_i[nSpecies-1] - P_i/rho_i);
////      daR[nSpecies-1] = 1.0/(2.0*a_j*rho_j) * (1+dPdU_j[nSpecies+nDim])*(dPdU_j[nSpecies-1] - P_j/rho_j);
////    }
////
////    /*--- Sound speed derivatives: Momentum ---*/
////    for (iDim = 0; iDim < nDim; iDim++) {
////      daL[nSpecies+iDim] = -1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim])*u_i[iDim];
////      daR[nSpecies+iDim] = -1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim])*u_j[iDim];
////    }
////
////    /*--- Sound speed derivatives: Energy ---*/
////    daL[nSpecies+nDim]   = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim]);
////    daR[nSpecies+nDim]   = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim]);
////
////    /*--- Sound speed derivatives: Vib-el energy ---*/
////    daL[nSpecies+nDim+1] = 1.0/(2.0*rho_i*a_i) * ((1.0+dPdU_i[nSpecies+nDim])*dPdU_i[nSpecies+nDim+1]);
////    daR[nSpecies+nDim+1] = 1.0/(2.0*rho_j*a_j) * ((1.0+dPdU_j[nSpecies+nDim])*dPdU_j[nSpecies+nDim+1]);
////
////    /*--- Left state Jacobian ---*/
////    if (mF >= 0) {
////
////      /*--- Jacobian contribution: dFc terms ---*/
////      for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
////        for (jVar = 0; jVar < nVar; jVar++) {
////          val_Jacobian_i[iVar][jVar] += mF * FcL[iVar]/a_i * daL[jVar];
////        }
////        val_Jacobian_i[iVar][iVar] += mF * a_i;
////      }
////      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
////        val_Jacobian_i[nSpecies+nDim][iSpecies] += mF * (dPdU_i[iSpecies]*a_i + rho_i*h_i*daL[iSpecies]);
////      }
////      for (iDim = 0; iDim < nDim; iDim++) {
////        val_Jacobian_i[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_i[nSpecies+nDim]*u_i[iDim]*a_i + rho_i*h_i*daL[nSpecies+iDim]);
////      }
////      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_i[nSpecies+nDim])*a_i + rho_i*h_i*daL[nSpecies+nDim]);
////      val_Jacobian_i[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_i[nSpecies+nDim+1]*a_i + rho_i*h_i*daL[nSpecies+nDim+1]);
////      for (jVar = 0; jVar < nVar; jVar++) {
////        val_Jacobian_i[nSpecies+nDim+1][jVar] +=  mF * FcL[nSpecies+nDim+1]/a_i * daL[jVar];
////      }
////      val_Jacobian_i[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_i;
////    }
////
////    /*--- Calculate derivatives of the split pressure flux ---*/
////    if ( (mF >= 0) || ((mF < 0)&&(fabs(mF) <= 1.0)) ) {
////      if (fabs(mL) <= 1.0) {
////
////        /*--- Mach number ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dmLP[iSpecies] = 0.5*(mL+1.0) * (-ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i));
////        for (iDim = 0; iDim < nDim; iDim++)
////          dmLP[nSpecies+iDim] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i));
////        dmLP[nSpecies+nDim]   = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim]);
////        dmLP[nSpecies+nDim+1] = 0.5*(mL+1.0) * (-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1]);
////
////        /*--- Pressure ---*/
////        for(iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dpLP[iSpecies] = 0.25*(mL+1.0) * (dPdU_i[iSpecies]*(mL+1.0)*(2.0-mL)
////                                            + P_i*(-ProjVel_i/(rho_i*a_i)
////                                                   -ProjVel_i*daL[iSpecies]/(a_i*a_i))*(3.0-3.0*mL));
////        for (iDim = 0; iDim < nDim; iDim++)
////          dpLP[nSpecies+iDim] = 0.25*(mL+1.0) * (-u_i[iDim]*dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
////              + P_i*( -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim]
////              + UnitNormal[iDim]/(rho_i*a_i))*(3.0-3.0*mL));
////        dpLP[nSpecies+nDim]   = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim]*(mL+1.0)*(2.0-mL)
////            + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim])*(3.0-3.0*mL));
////        dpLP[nSpecies+nDim+1] = 0.25*(mL+1.0) * (dPdU_i[nSpecies+nDim+1]*(mL+1.0)*(2.0-mL)
////            + P_i*(-ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1])*(3.0-3.0*mL));
////      } else {
////
////        /*--- Mach number ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dmLP[iSpecies]      = -ProjVel_i/(rho_i*a_i) - ProjVel_i*daL[iSpecies]/(a_i*a_i);
////        for (iDim = 0; iDim < nDim; iDim++)
////          dmLP[nSpecies+iDim] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+iDim] + UnitNormal[iDim]/(rho_i*a_i);
////        dmLP[nSpecies+nDim]   = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim];
////        dmLP[nSpecies+nDim+1] = -ProjVel_i/(a_i*a_i) * daL[nSpecies+nDim+1];
////
////        /*--- Pressure ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dpLP[iSpecies] = dPdU_i[iSpecies];
////        for (iDim = 0; iDim < nDim; iDim++)
////          dpLP[nSpecies+iDim] = (-u_i[iDim]*dPdU_i[nSpecies+nDim]);
////        dpLP[nSpecies+nDim]   = dPdU_i[nSpecies+nDim];
////        dpLP[nSpecies+nDim+1] = dPdU_i[nSpecies+nDim+1];
////      }
////
////      /*--- dM contribution ---*/
////      for (iVar = 0; iVar < nVar; iVar++) {
////        for (jVar = 0; jVar < nVar; jVar++) {
////          val_Jacobian_i[iVar][jVar] += dmLP[jVar]*FcLR[iVar];
////        }
////      }
////
////      /*--- Jacobian contribution: dP terms ---*/
////      for (iDim = 0; iDim < nDim; iDim++) {
////        for (iVar = 0; iVar < nVar; iVar++) {
////          val_Jacobian_i[nSpecies+iDim][iVar] += dpLP[iVar]*UnitNormal[iDim];
////        }
////      }
////    }
////
////    /*--- Right state Jacobian ---*/
////    if (mF < 0) {
////
////      /*--- Jacobian contribution: dFc terms ---*/
////      for (iVar = 0; iVar < nSpecies+nDim; iVar++) {
////        for (jVar = 0; jVar < nVar; jVar++) {
////          val_Jacobian_j[iVar][jVar] += mF * FcR[iVar]/a_j * daR[jVar];
////        }
////        val_Jacobian_j[iVar][iVar] += mF * a_j;
////      }
////      for (iSpecies = 0; iSpecies < nSpecies; iSpecies++) {
////        val_Jacobian_j[nSpecies+nDim][iSpecies] += mF * (dPdU_j[iSpecies]*a_j + rho_j*h_j*daR[iSpecies]);
////      }
////      for (iDim = 0; iDim < nDim; iDim++) {
////        val_Jacobian_j[nSpecies+nDim][nSpecies+iDim] += mF * (-dPdU_j[nSpecies+nDim]*u_j[iDim]*a_j + rho_j*h_j*daR[nSpecies+iDim]);
////      }
////      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim]   += mF * ((1.0+dPdU_j[nSpecies+nDim])*a_j + rho_j*h_j*daR[nSpecies+nDim]);
////      val_Jacobian_j[nSpecies+nDim][nSpecies+nDim+1] += mF * (dPdU_j[nSpecies+nDim+1]*a_j + rho_j*h_j*daR[nSpecies+nDim+1]);
////      for (jVar = 0; jVar < nVar; jVar++) {
////        val_Jacobian_j[nSpecies+nDim+1][jVar] +=  mF * FcR[nSpecies+nDim+1]/a_j * daR[jVar];
////      }
////      val_Jacobian_j[nSpecies+nDim+1][nSpecies+nDim+1] += mF * a_j;
////    }
////
////    /*--- Calculate derivatives of the split pressure flux ---*/
////    if ( (mF < 0) || ((mF >= 0)&&(fabs(mF) <= 1.0)) ) {
////      if (fabs(mR) <= 1.0) {
////
////        /*--- Mach ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dmRM[iSpecies] = -0.5*(mR-1.0) * (-ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j));
////        for (iDim = 0; iDim < nDim; iDim++)
////          dmRM[nSpecies+iDim] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j));
////        dmRM[nSpecies+nDim]   = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim]);
////        dmRM[nSpecies+nDim+1] = -0.5*(mR-1.0) * (-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1]);
////
////        /*--- Pressure ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dpRM[iSpecies] = 0.25*(mR-1.0) * (dPdU_j[iSpecies]*(mR-1.0)*(2.0+mR)
////                                            + P_j*(-ProjVel_j/(rho_j*a_j)
////                                                   -ProjVel_j*daR[iSpecies]/(a_j*a_j))*(3.0+3.0*mR));
////        for (iDim = 0; iDim < nDim; iDim++)
////          dpRM[nSpecies+iDim] = 0.25*(mR-1.0) * ((-u_j[iDim]*dPdU_j[nSpecies+nDim])*(mR-1.0)*(2.0+mR)
////              + P_j*( -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim]
////              + UnitNormal[iDim]/(rho_j*a_j))*(3.0+3.0*mR));
////        dpRM[nSpecies+nDim]   = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim]*(mR-1.0)*(2.0+mR)
////            + P_j*(-ProjVel_j/(a_j*a_j)*daR[nSpecies+nDim])*(3.0+3.0*mR));
////        dpRM[nSpecies+nDim+1] = 0.25*(mR-1.0) * (dPdU_j[nSpecies+nDim+1]*(mR-1.0)*(2.0+mR)
////            + P_j*(-ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1])*(3.0+3.0*mR));
////
////      } else {
////
////        /*--- Mach ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dmRM[iSpecies]      = -ProjVel_j/(rho_j*a_j) - ProjVel_j*daR[iSpecies]/(a_j*a_j);
////        for (iDim = 0; iDim < nDim; iDim++)
////          dmRM[nSpecies+iDim] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+iDim] + UnitNormal[iDim]/(rho_j*a_j);
////        dmRM[nSpecies+nDim]   = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim];
////        dmRM[nSpecies+nDim+1] = -ProjVel_j/(a_j*a_j) * daR[nSpecies+nDim+1];
////
////        /*--- Pressure ---*/
////        for (iSpecies = 0; iSpecies < nSpecies; iSpecies++)
////          dpRM[iSpecies] = dPdU_j[iSpecies];
////        for (iDim = 0; iDim < nDim; iDim++)
////          dpRM[nSpecies+iDim] = -u_j[iDim]*dPdU_j[nSpecies+nDim];
////        dpRM[nSpecies+nDim]   = dPdU_j[nSpecies+nDim];
////        dpRM[nSpecies+nDim+1] = dPdU_j[nSpecies+nDim+1];
////      }
////
////      /*--- Jacobian contribution: dM terms ---*/
////      for (iVar = 0; iVar < nVar; iVar++) {
////        for (jVar = 0; jVar < nVar; jVar++) {
////          val_Jacobian_j[iVar][jVar] += dmRM[jVar] * FcLR[iVar];
////        }
////      }
////
////      /*--- Jacobian contribution: dP terms ---*/
////      for (iDim = 0; iDim < nDim; iDim++) {
////        for (iVar = 0; iVar < nVar; iVar++) {
////          val_Jacobian_j[nSpecies+iDim][iVar] += dpRM[iVar]*UnitNormal[iDim];
////        }
////      }
////    }
////
////    /*--- Integrate over dual-face area ---*/
////    for (iVar = 0; iVar < nVar; iVar++) {
////      for (jVar = 0; jVar < nVar; jVar++) {
////        val_Jacobian_i[iVar][jVar] *= Area;
////        val_Jacobian_j[iVar][jVar] *= Area;
////      }
////    }
////  }
//
//  return ResidualType<>(Flux, nullptr, nullptr);
//}


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
