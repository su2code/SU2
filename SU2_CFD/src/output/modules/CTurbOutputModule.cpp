#include "../../../include/output/modules/CTurbOutputModule.hpp"
#include "../../../include/solvers/CSolver.hpp"


void CTurbOutputModule::DefineHistoryFields(COutFieldCollection &fieldCollection){
  switch(turb_model){
    case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
      /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
      fieldCollection.AddItem("RMS_NU_TILDE", COutputField("rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of nu tilde (SA model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("MAX_NU_TILDE", COutputField("max[nu]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of nu tilde (SA model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("BGS_NU_TILDE", COutputField("bgs[nu]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of nu tilde (SA model).", FieldType::RESIDUAL));
      break;
    case SST: case SST_SUST:
      /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
      fieldCollection.AddItem("RMS_TKE", COutputField("rms[k]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of kinetic energy (SST model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("MAX_TKE", COutputField("max[k]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of kinetic energy (SST model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("BSG_TKE", COutputField("bgs[k]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of kinetic energy (SST model).", FieldType::RESIDUAL));

      /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
      fieldCollection.AddItem("RMS_DISSIPATION", COutputField("rms[w]", ScreenOutputFormat::FIXED, "RMS_RES", "Root-mean square residual of dissipation (SST model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("MAX_DISSIPATION", COutputField("max[w]", ScreenOutputFormat::FIXED, "MAX_RES", "Maximum residual of dissipation (SST model).", FieldType::RESIDUAL));
      fieldCollection.AddItem("BGS_DISSIPATION", COutputField("bgs[w]", ScreenOutputFormat::FIXED, "BGS_RES", "Block-Gauss-Seidel residual of dissipation (SST model).", FieldType::RESIDUAL));
      break;
    default: break;
  }
}

void CTurbOutputModule::LoadHistoryData(COutFieldCollection &fieldCollection){
  CSolver* turb_solver = solverData.solver[TURB_SOL];

  switch(turb_model){
  case SA: case SA_NEG: case SA_E: case SA_COMP: case SA_E_COMP:
    fieldCollection.SetValueByKey("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
    fieldCollection.SetValueByKey("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
    if (solverData.config->GetMultizone_Problem())
      fieldCollection.SetValueByKey("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
    break;
  case SST: case SST_SUST:
    fieldCollection.SetValueByKey("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
    fieldCollection.SetValueByKey("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
    if (solverData.config->GetMultizone_Problem())
      fieldCollection.SetValueByKey("BGS_TKE", log10(turb_solver->GetRes_BGS(0)));

    fieldCollection.SetValueByKey("RMS_DISSIPATION",    log10(turb_solver->GetRes_RMS(1)));
    fieldCollection.SetValueByKey("MAX_DISSIPATION",    log10(turb_solver->GetRes_Max(1)));
    if (solverData.config->GetMultizone_Problem())
      fieldCollection.SetValueByKey("BGS_DISSIPATION", log10(turb_solver->GetRes_BGS(1)));
    break;
  default: break;
  }
}

void CTurbOutputModule::DefineVolumeFields(COutFieldCollection &fieldCollection){
  // Turbulent Residuals
  switch(turb_model){
    case SST: case SST_SUST:
      fieldCollection.AddItem("TKE", COutputField("Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy", FieldType::DEFAULT));
      fieldCollection.AddItem("DISSIPATION", COutputField("Omega", "SOLUTION", "Rate of dissipation", FieldType::DEFAULT));
      break;
    case SA: case SA_COMP: case SA_E:
    case SA_E_COMP: case SA_NEG:
      fieldCollection.AddItem("NU_TILDE", COutputField("Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable", FieldType::DEFAULT));
      break;
    case NONE:
      break;
  }
}

void CTurbOutputModule::LoadVolumeData(COutFieldCollection &fieldCollection){

  unsigned long iPoint = solverData.iPoint;

  CVariable* Node_Turb = solverData.solver[TURB_SOL]->GetNodes();

  switch(turb_model){
  case SST: case SST_SUST:
    fieldCollection.SetValueByKey("TKE", Node_Turb->GetSolution(iPoint, 0));
    fieldCollection.SetValueByKey("DISSIPATION",  Node_Turb->GetSolution(iPoint, 1));
    break;
  case SA: case SA_COMP: case SA_E:
  case SA_E_COMP: case SA_NEG:
    fieldCollection.SetValueByKey("NU_TILDE", Node_Turb->GetSolution(iPoint, 0));
    break;
  case NONE:
    break;
  }
}
