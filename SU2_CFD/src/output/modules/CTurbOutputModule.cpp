#include "../../../include/output/modules/CTurbOutputModule.hpp"
#include "../../../include/solvers/CSolver.hpp"


void CTurbOutputModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){
  switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
      historyFields.AddField("RMS_NU_TILDE", "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_NU_TILDE", "max[nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", FieldType::RESIDUAL);
      historyFields.AddField("BGS_NU_TILDE", "bgs[nu]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of nu tilde (SA model).", FieldType::RESIDUAL);
      break;
    case SST: case SST_SUST:
      /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
      historyFields.AddField("RMS_TKE", "rms[k]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_TKE", "max[k]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("BSG_TKE", "bgs[k]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of kinetic energy (SST model).", FieldType::RESIDUAL);

      /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
      historyFields.AddField("RMS_DISSIPATION", "rms[w]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_DISSIPATION", "max[w]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("BGS_DISSIPATION", "bgs[w]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of dissipation (SST model).", FieldType::RESIDUAL);
      break;
    default: break;
  }
}

void CTurbOutputModule::LoadHistoryData(CHistoryOutFieldManager &historyFields){
  CSolver* turb_solver = solverData.solver[TURB_SOL];

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    historyFields.SetFieldValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    historyFields.SetFieldValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    if (solverData.config->GetMultizone_Problem())
      historyFields.SetFieldValue("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
    break;
  case SST: case SST_SUST:
    historyFields.SetFieldValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    historyFields.SetFieldValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    if (solverData.config->GetMultizone_Problem())
      historyFields.SetFieldValue("BGS_TKE", log10(turb_solver->GetRes_BGS(0)));

    historyFields.SetFieldValue("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
    historyFields.SetFieldValue("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    if (solverData.config->GetMultizone_Problem())
      historyFields.SetFieldValue("BGS_DISSIPATION", log10(turb_solver->GetRes_BGS(1)));
    break;
  default: break;
  }
}

void CTurbOutputModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){
  // Turbulent Residuals
  switch(turb_model){
    case SST: case SST_SUST:
      volumeFields.AddField("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy", FieldType::DEFAULT);
      volumeFields.AddField("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation", FieldType::DEFAULT);
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      volumeFields.AddField("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable", FieldType::DEFAULT);
      break;
    case NONE:
      break;
  }
}

void CTurbOutputModule::LoadVolumeData(CVolumeOutFieldManager &volumeFields){

  const auto iPoint = solverData.iPoint;
  const auto *Node_Turb = solverData.solver[TURB_SOL]->GetNodes();

  switch(turb_model){
  case SST: case SST_SUST:
    volumeFields.SetFieldValue("TKE", Node_Turb->GetSolution(iPoint, 0));
    volumeFields.SetFieldValue("DISSIPATION",  Node_Turb->GetSolution(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    volumeFields.SetFieldValue("NU_TILDE", Node_Turb->GetSolution(iPoint, 0));
    break;
  case NONE:
    break;
  }
}
