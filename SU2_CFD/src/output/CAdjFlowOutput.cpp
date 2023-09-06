/*!
 * \file CAdjFlowOutput.cpp
 * \brief Main subroutines for flow discrete adjoint output
 * \author T. Kattmann
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

#include "../../include/output/CAdjFlowOutput.hpp"

#include "../../include/solvers/CSolver.hpp"

CAdjFlowOutput::CAdjFlowOutput(CConfig* config, unsigned short nDim)
    : COutput(config, nDim, false),
      turb_model(config->GetKind_Turb_Model()),
      cont_adj(config->GetContinuous_Adjoint()),
      frozen_visc((config->GetFrozen_Visc_Disc() && !cont_adj) || (config->GetFrozen_Visc_Cont() && cont_adj)) {}

// The "AddHistoryOutput(" must not be split over multiple lines to ensure proper python parsing
// clang-format off
void CAdjFlowOutput::AddHistoryOutputFields_AdjScalarRMS_RES(const CConfig* config) {
  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        /// DESCRIPTION: Root-mean square residual of the adjoint nu tilde.
        AddHistoryOutput("RMS_ADJ_NU_TILDE", "rms[A_nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::KW:
        /// DESCRIPTION: Root-mean square residual of the adjoint kinetic energy.
        AddHistoryOutput("RMS_ADJ_TKE", "rms[A_k]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
        /// DESCRIPTION: Root-mean square residual of the adjoint dissipation.
        AddHistoryOutput("RMS_ADJ_DISSIPATION", "rms[A_w]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddHistoryOutput("RMS_ADJ_SPECIES_" + std::to_string(iVar), "rms[A_rho*Y_" + std::to_string(iVar) + "]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint transported species.", HistoryFieldType::RESIDUAL);
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      AddHistoryOutput("RMS_ADJ_"+cv_name, "rms[" + cv_name + "]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint " + cv_name, HistoryFieldType::RESIDUAL);
    }
    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      AddHistoryOutput("RMS_ADJ_" + scalar_name, "rms[" + scalar_name + "]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of the adjoint of " + scalar_name + " .", HistoryFieldType::RESIDUAL);
    }
  }
}

void CAdjFlowOutput::AddHistoryOutputFields_AdjScalarMAX_RES(const CConfig* config) {
  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        /// DESCRIPTION: Maximum residual of the adjoint nu tilde.
        AddHistoryOutput("MAX_ADJ_NU_TILDE", "max[A_nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::KW:
        /// DESCRIPTION: Maximum residual of the adjoint kinetic energy.
        AddHistoryOutput("MAX_ADJ_TKE", "max[A_k]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
        /// DESCRIPTION: Maximum residual of the adjoint dissipation.
        AddHistoryOutput("MAX_ADJ_DISSIPATION", "max[A_w]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddHistoryOutput("MAX_ADJ_SPECIES_" + std::to_string(iVar), "max[A_rho*Y_" + std::to_string(iVar) + "]",ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint transported species.", HistoryFieldType::RESIDUAL);
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      AddHistoryOutput("MAX_ADJ_" + cv_name, "max["+cv_name +"]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint " + cv_name, HistoryFieldType::RESIDUAL);
    }
    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      AddHistoryOutput("MAX_ADJ_" + scalar_name, "max[scalar_" + scalar_name + "]",ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of the adjoint of " + scalar_name + " .", HistoryFieldType::RESIDUAL);
    }
  }
}

void CAdjFlowOutput::AddHistoryOutputFields_AdjScalarBGS_RES(const CConfig* config) {
  if (!multiZone) return;

  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        /// DESCRIPTION: BGS residual of the adjoint nu tilde.
        AddHistoryOutput("BGS_ADJ_NU_TILDE", "bgs[A_nu]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint nu tilde.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::KW:
        /// DESCRIPTION: BGS residual of the adjoint kinetic energy.
        AddHistoryOutput("BGS_ADJ_TKE", "bgs[A_k]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint kinetic energy.", HistoryFieldType::RESIDUAL);
        /// DESCRIPTION: BGS residual of the adjoint dissipation.
        AddHistoryOutput("BGS_ADJ_DISSIPATION", "bgs[A_w]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint dissipation.", HistoryFieldType::RESIDUAL);
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddHistoryOutput("BGS_ADJ_SPECIES_" + std::to_string(iVar), "bgs[A_rho*Y_" + std::to_string(iVar) + "]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint transported species.", HistoryFieldType::RESIDUAL);
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      AddHistoryOutput("BGS_ADJ_" + cv_name, "bgs[" + cv_name + "]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint " + cv_name, HistoryFieldType::RESIDUAL);
    }

    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      AddHistoryOutput("BGS_ADJ_" + scalar_name, "bgs[" + scalar_name + "]", ScreenOutputFormat::FIXED, "BGS_RES", "BGS residual of the adjoint of " + scalar_name + " .", HistoryFieldType::RESIDUAL);
    }
  }
}

void CAdjFlowOutput::AddHistoryOutputFieldsAdjScalarLinsol(const CConfig* config) {
  if ((turb_model != TURB_MODEL::NONE) && !frozen_visc) {
    AddHistoryOutput("LINSOL_ITER_TURB", "LinSolIterTurb", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for turbulence.");
    AddHistoryOutput("LINSOL_RESIDUAL_TURB", "LinSolResTurb", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for turbulence.");
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    AddHistoryOutput("LINSOL_ITER_SPECIES", "LinSolIterSpecies", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for species solver.");
    AddHistoryOutput("LINSOL_RESIDUAL_SPECIES", "LinSolResSpecies", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for species solver.");
  }
  
  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    AddHistoryOutput("LINSOL_ITER_FLAMELET", "LinSolIterScalar", ScreenOutputFormat::INTEGER, "LINSOL", "Number of iterations of the linear solver for scalar solver.");
    AddHistoryOutput("LINSOL_RESIDUAL_FLAMELET", "LinSolResScalar", ScreenOutputFormat::FIXED, "LINSOL", "Residual of the linear solver for scalar solver.");
  }
}
// clang-format on

void CAdjFlowOutput::LoadHistoryDataAdjScalar(const CConfig* config, const CSolver* const* solver) {
  const auto adjturb_solver = solver[ADJTURB_SOL];
  const auto adjspecies_solver = solver[ADJSPECIES_SOL];

  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        SetHistoryOutputValue("RMS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_RMS(0)));
        SetHistoryOutputValue("MAX_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_Max(0)));
        if (multiZone) {
          SetHistoryOutputValue("BGS_ADJ_NU_TILDE", log10(adjturb_solver->GetRes_BGS(0)));
        }
        break;
      case TURB_FAMILY::KW:
        SetHistoryOutputValue("RMS_ADJ_TKE", log10(adjturb_solver->GetRes_RMS(0)));
        SetHistoryOutputValue("RMS_ADJ_DISSIPATION", log10(adjturb_solver->GetRes_RMS(1)));
        SetHistoryOutputValue("MAX_ADJ_TKE", log10(adjturb_solver->GetRes_Max(0)));
        SetHistoryOutputValue("MAX_ADJ_DISSIPATION", log10(adjturb_solver->GetRes_Max(1)));
        if (multiZone) {
          SetHistoryOutputValue("BGS_ADJ_TKE", log10(adjturb_solver->GetRes_BGS(0)));
          SetHistoryOutputValue("BGS_ADJ_DISSIPATION", log10(adjturb_solver->GetRes_BGS(1)));
        }
        break;
      case TURB_FAMILY::NONE:
        break;
    }

    if (turb_model != TURB_MODEL::NONE) {
      SetHistoryOutputValue("LINSOL_ITER_TURB", adjturb_solver->GetIterLinSolver());
      SetHistoryOutputValue("LINSOL_RESIDUAL_TURB", log10(adjturb_solver->GetResLinSolver()));
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      SetHistoryOutputValue("RMS_ADJ_SPECIES_" + std::to_string(iVar), log10(adjspecies_solver->GetRes_RMS(iVar)));
      SetHistoryOutputValue("MAX_ADJ_SPECIES_" + std::to_string(iVar), log10(adjspecies_solver->GetRes_Max(iVar)));
      if (multiZone) {
        SetHistoryOutputValue("BGS_ADJ_SPECIES_" + std::to_string(iVar), log10(adjspecies_solver->GetRes_BGS(iVar)));
      }
    }

    SetHistoryOutputValue("LINSOL_ITER_SPECIES", adjspecies_solver->GetIterLinSolver());
    SetHistoryOutputValue("LINSOL_RESIDUAL_SPECIES", log10(adjspecies_solver->GetResLinSolver()));
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
      for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
        const auto& cv_name = config->GetControllingVariableName(iCV);
        SetHistoryOutputValue("RMS_ADJ_"+cv_name, log10(adjspecies_solver->GetRes_RMS(iCV)));
        SetHistoryOutputValue("MAX_ADJ_"+cv_name, log10(adjspecies_solver->GetRes_Max(iCV)));
        if (multiZone)
          SetHistoryOutputValue("BGS_ADJ_" + cv_name, log10(adjspecies_solver->GetRes_BGS(iCV)));
      }

    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      SetHistoryOutputValue("RMS_ADJ_" + scalar_name, log10(adjspecies_solver->GetRes_RMS(2 + i_scalar)));
      SetHistoryOutputValue("MAX_ADJ_" + scalar_name, log10(adjspecies_solver->GetRes_Max(2 + i_scalar)));
      if (multiZone) {
        SetHistoryOutputValue("BGS_ADJ_" + scalar_name, log10(adjspecies_solver->GetRes_BGS(2 + i_scalar)));
      }
    }

    SetHistoryOutputValue("LINSOL_ITER_FLAMELET", adjspecies_solver->GetIterLinSolver());
    SetHistoryOutputValue("LINSOL_RESIDUAL_FLAMELET", log10(adjspecies_solver->GetResLinSolver()));
  }

  ComputeSimpleCustomOutputs(config);
}

void CAdjFlowOutput::SetVolumeOutputFieldsAdjScalarSolution(const CConfig* config) {
  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        /// DESCRIPTION: Adjoint nu tilde.
        AddVolumeOutput("ADJ_NU_TILDE", "Adjoint_Nu_Tilde", "SOLUTION", "Adjoint Spalart-Allmaras variable");
        break;
      case TURB_FAMILY::KW:
        /// DESCRIPTION: Adjoint kinetic energy.
        AddVolumeOutput("ADJ_TKE", "Adjoint_TKE", "SOLUTION", "Adjoint turbulent kinetic energy");
        /// DESCRIPTION: Adjoint dissipation.
        AddVolumeOutput("ADJ_DISSIPATION", "Adjoint_Omega", "SOLUTION", "Adjoint rate of dissipation");
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddVolumeOutput("ADJ_SPECIES_" + std::to_string(iVar), "Adjoint_Species_" + std::to_string(iVar), "SOLUTION",
                      "Adjoint Species variable");
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      AddVolumeOutput("ADJ_" + cv_name, "Adjoint_" + cv_name, "SOLUTION", "Adjoint of the " + cv_name + " controlling variable.");
    }

    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      AddVolumeOutput("ADJ_" + scalar_name, "Adjoint_" + scalar_name, "SOLUTION", "Adjoint of " + scalar_name);
    }
  }
}

void CAdjFlowOutput::SetVolumeOutputFieldsAdjScalarResidual(const CConfig* config) {
  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        /// DESCRIPTION: Residual of the nu tilde.
        AddVolumeOutput("RES_ADJ_NU_TILDE", "Residual_Adjoint_Nu_Tilde", "RESIDUAL",
                        "Residual of the adjoint Spalart-Allmaras variable");
        break;
      case TURB_FAMILY::KW:
        /// DESCRIPTION: Residual of the adjoint kinetic energy.
        AddVolumeOutput("RES_ADJ_TKE", "Residual_Adjoint_TKE", "RESIDUAL",
                        "Residual of the adjoint turb. kinetic energy");
        /// DESCRIPTION: Residual of the adjoint dissipation.
        AddVolumeOutput("RES_ADJ_DISSIPATION", "Residual_Adjoint_Omega", "RESIDUAL",
                        "Residual of adjoint rate of dissipation");
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      AddVolumeOutput("RES_ADJ_SPECIES_" + std::to_string(iVar), "Residual_Adjoint_Species_" + std::to_string(iVar),
                      "RESIDUAL", "Residual of the adjoint Species variable");
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {
    for(auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      AddVolumeOutput("RES_ADJ_" + cv_name, "Residual_Adjoint_" + cv_name, "RESIDUAL", "Residual of the adjoint of " + cv_name);
    }
    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      AddVolumeOutput("RES_ADJ_" + scalar_name, "Residual_Adjoint_" + scalar_name, "RESIDUAL", "Residual of the adjoint of " + scalar_name);
    }
  }
}

void CAdjFlowOutput::LoadVolumeDataAdjScalar(const CConfig* config, const CSolver* const* solver,
                                              const unsigned long iPoint) {
  const auto Node_AdjTurb =
      ((turb_model != TURB_MODEL::NONE) && !frozen_visc) ? solver[ADJTURB_SOL]->GetNodes() : nullptr;
  const auto Node_AdjSpecies =
      (config->GetKind_Species_Model() != SPECIES_MODEL::NONE) ? solver[ADJSPECIES_SOL]->GetNodes() : nullptr;

  if (!frozen_visc) {
    switch (TurbModelFamily(turb_model)) {
      case TURB_FAMILY::SA:
        SetVolumeOutputValue("ADJ_NU_TILDE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
        SetVolumeOutputValue("RES_ADJ_NU_TILDE", iPoint,
                             Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
        break;
      case TURB_FAMILY::KW:
        SetVolumeOutputValue("ADJ_TKE", iPoint, Node_AdjTurb->GetSolution(iPoint, 0));
        SetVolumeOutputValue("ADJ_DISSIPATION", iPoint, Node_AdjTurb->GetSolution(iPoint, 1));
        SetVolumeOutputValue("RES_ADJ_TKE", iPoint,
                             Node_AdjTurb->GetSolution(iPoint, 0) - Node_AdjTurb->GetSolution_Old(iPoint, 0));
        SetVolumeOutputValue("RES_ADJ_DISSIPATION", iPoint,
                             Node_AdjTurb->GetSolution(iPoint, 1) - Node_AdjTurb->GetSolution_Old(iPoint, 1));
        break;
      case TURB_FAMILY::NONE:
        break;
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::SPECIES_TRANSPORT) {
    for (unsigned short iVar = 0; iVar < config->GetnSpecies(); iVar++) {
      SetVolumeOutputValue("ADJ_SPECIES_" + std::to_string(iVar), iPoint, Node_AdjSpecies->GetSolution(iPoint, iVar));
      SetVolumeOutputValue("RES_ADJ_SPECIES_" + std::to_string(iVar), iPoint,
                           Node_AdjSpecies->GetSolution(iPoint, iVar) - Node_AdjSpecies->GetSolution_Old(iPoint, iVar));
    }
  }

  if (config->GetKind_Species_Model() == SPECIES_MODEL::FLAMELET) {

    for (auto iCV=0u; iCV < config->GetNControlVars(); iCV++) {
      const auto& cv_name = config->GetControllingVariableName(iCV);
      SetVolumeOutputValue("ADJ_" + cv_name, iPoint, Node_AdjSpecies->GetSolution(iPoint, iCV));
      SetVolumeOutputValue("RES_ADJ_" + cv_name, iPoint, Node_AdjSpecies->GetSolution(iPoint, iCV) - Node_AdjSpecies->GetSolution_Old(iPoint, iCV));
    }
    for (unsigned short i_scalar = 0; i_scalar < config->GetNUserScalars(); i_scalar++) {
      const auto& scalar_name = config->GetUserScalarName(i_scalar);
      SetVolumeOutputValue("ADJ_" + scalar_name, iPoint, Node_AdjSpecies->GetSolution(iPoint, config->GetNControlVars() + i_scalar));
      SetVolumeOutputValue("RES_ADJ_" + scalar_name, iPoint, Node_AdjSpecies->GetSolution(iPoint, config->GetNControlVars() + i_scalar) - Node_AdjSpecies->GetSolution_Old(iPoint, config->GetNControlVars() + i_scalar));
    }
  }

}
