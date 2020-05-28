#include "../../../include/output/modules/CTurbOutputModule.hpp"

#include "../../../include/solvers/CSolver.hpp"

void CTurbOutputModule::DefineHistoryFields(CHistoryOutFieldManager& historyFields) {
  switch (turb_model) {
    case SA:
    case SA_NEG:
    case SA_E:
    case SA_COMP:
    case SA_E_COMP:
      /// DESCRIPTION: Root-mean square residual of nu tilde (SA model).
      historyFields.AddField("RMS_NU_TILDE", "rms[nu]", ScreenOutputFormat::FIXED, "RMS_RES",
                             "Root-mean square residual of nu tilde (SA model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_NU_TILDE", "max[nu]", ScreenOutputFormat::FIXED, "MAX_RES",
                             "Maximum residual of nu tilde (SA model).", FieldType::RESIDUAL);
      historyFields.AddField("BGS_NU_TILDE", "bgs[nu]", ScreenOutputFormat::FIXED, "BGS_RES",
                             "Block-Gauss-Seidel residual of nu tilde (SA model).", FieldType::RESIDUAL);
      break;
    case SST:
    case SST_SUST:
      /// DESCRIPTION: Root-mean square residual of kinetic energy (SST model).
      historyFields.AddField("RMS_TKE", "rms[k]", ScreenOutputFormat::FIXED, "RMS_RES",
                             "Root-mean square residual of kinetic energy (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_TKE", "max[k]", ScreenOutputFormat::FIXED, "MAX_RES",
                             "Maximum residual of kinetic energy (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("BSG_TKE", "bgs[k]", ScreenOutputFormat::FIXED, "BGS_RES",
                             "Block-Gauss-Seidel residual of kinetic energy (SST model).", FieldType::RESIDUAL);

      /// DESCRIPTION: Root-mean square residual of the dissipation (SST model).
      historyFields.AddField("RMS_DISSIPATION", "rms[w]", ScreenOutputFormat::FIXED, "RMS_RES",
                             "Root-mean square residual of dissipation (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("MAX_DISSIPATION", "max[w]", ScreenOutputFormat::FIXED, "MAX_RES",
                             "Maximum residual of dissipation (SST model).", FieldType::RESIDUAL);
      historyFields.AddField("BGS_DISSIPATION", "bgs[w]", ScreenOutputFormat::FIXED, "BGS_RES",
                             "Block-Gauss-Seidel residual of dissipation (SST model).", FieldType::RESIDUAL);
      break;
    default:
      break;
  }
}

void CTurbOutputModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                        const IterationInfo& iterationInfo) {

  const auto* config      = solverData.config;
  const auto* turb_solver = solverData.solver[TURB_SOL];

  switch (turb_model) {
    case SA:
    case SA_NEG:
    case SA_E:
    case SA_COMP:
    case SA_E_COMP:
      historyFields.SetFieldValue("RMS_NU_TILDE", log10(turb_solver->GetRes_RMS(0)));
      historyFields.SetFieldValue("MAX_NU_TILDE", log10(turb_solver->GetRes_Max(0)));
      if (config->GetMultizone_Problem())
        historyFields.SetFieldValue("BGS_NU_TILDE", log10(turb_solver->GetRes_BGS(0)));
      break;
    case SST:
    case SST_SUST:
      historyFields.SetFieldValue("RMS_TKE", log10(turb_solver->GetRes_RMS(0)));
      historyFields.SetFieldValue("MAX_TKE", log10(turb_solver->GetRes_Max(0)));
      if (config->GetMultizone_Problem())
        historyFields.SetFieldValue("BGS_TKE", log10(turb_solver->GetRes_BGS(0)));

      historyFields.SetFieldValue("RMS_DISSIPATION", log10(turb_solver->GetRes_RMS(1)));
      historyFields.SetFieldValue("MAX_DISSIPATION", log10(turb_solver->GetRes_Max(1)));
      if (config->GetMultizone_Problem())
        historyFields.SetFieldValue("BGS_DISSIPATION", log10(turb_solver->GetRes_BGS(1)));
      break;
    default:
      break;
  }
}

void CTurbOutputModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {
  volumeFields.AddField("EDDY_VISCOSITY", "Eddy_Viscosity", "PRIMITIVE", "Turbulent eddy viscosity",
                        FieldType::DEFAULT);

  // Turbulent Residuals
  switch (turb_model) {
    case SST:
    case SST_SUST:
      volumeFields.AddField("TKE", "Turb_Kin_Energy", "SOLUTION", "Turbulent kinetic energy", FieldType::DEFAULT);
      volumeFields.AddField("DISSIPATION", "Omega", "SOLUTION", "Rate of dissipation", FieldType::DEFAULT);
      volumeFields.AddField("RES_TKE", "Residual_TKE", "RESIDUAL", "Residual of turbulent kinetic energy",
                            FieldType::DEFAULT);
      volumeFields.AddField("RES_DISS", "Residual_Omega", "RESIDUAL", "Residual of the rate of dissipation.",
                            FieldType::DEFAULT);
      volumeFields.AddField("LIM_TKE", "Limiter_TKE", "LIMITER", "Limiter value of turb. kinetic energy.",
                            FieldType::DEFAULT);
      volumeFields.AddField("LIM_DISS", "Limiter_Omega", "LIMITER", "Limiter value of dissipation rate.",
                            FieldType::DEFAULT);
      break;
    case SA:
    case SA_COMP:
    case SA_E:
    case SA_E_COMP:
    case SA_NEG:
      volumeFields.AddField("NU_TILDE", "Nu_Tilde", "SOLUTION", "Spalart-Allmaras variable", FieldType::DEFAULT);
      volumeFields.AddField("RES_NU_TILDE", "Residual_Nu_Tilde", "RESIDUAL",
                            "Residual of the Spalart–Allmaras variable", FieldType::DEFAULT);
      volumeFields.AddField("LIM_NU_TILDE", "Limiter_Nu_Tilde", "LIMITER",
                            "Limiter value of Spalart–Allmaras variable.", FieldType::DEFAULT);
      break;
    case NONE:
      break;
  }

  if (hybrid_RANSLES != NO_HYBRIDRANSLES) {
    volumeFields.AddField("LENGTH_SCALE", "DES_LengthScale", "DDES", "DES length scale value", FieldType::DEFAULT);
    volumeFields.AddField("WALL_DISTANCE", "Wall_Distance", "DDES", "Wall distance value", FieldType::DEFAULT);
  }

  if (transModel == BC) {
    volumeFields.AddField("INTERMITTENCY", "gamma_BC", "INTERMITTENCY", "Intermittency", FieldType::DEFAULT);
  }
}

void CTurbOutputModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) {

  const auto* geometry    = solverData.geometry;
  auto* turb_solver       = solverData.solver[TURB_SOL];
  auto* flow_solver       = solverData.solver[FLOW_SOL];
  const auto  iPoint      = pointInfo.iPoint;
  const auto* Node_Turb  = turb_solver->GetNodes();
  const auto* Node_Flow  = flow_solver->GetNodes();
  const auto* Node_Geo   = geometry->nodes;

  volumeFields.SetFieldValue("EDDY_VISCOSITY", Node_Flow->GetEddyViscosity(iPoint));

  switch (turb_model) {
    case SST:
    case SST_SUST:
      volumeFields.SetFieldValue("TKE", Node_Turb->GetSolution(iPoint, 0));
      volumeFields.SetFieldValue("DISSIPATION", Node_Turb->GetSolution(iPoint, 1));
      volumeFields.SetFieldValue("RES_TKE",  turb_solver->LinSysRes(iPoint, 0));
      volumeFields.SetFieldValue("RES_DISS", turb_solver->LinSysRes(iPoint, 1));
      volumeFields.SetFieldValue("LIM_TKE", Node_Turb->GetLimiter_Primitive(iPoint, 0));
      volumeFields.SetFieldValue("LIM_DISS", Node_Turb->GetLimiter_Primitive(iPoint, 1));
      break;
    case SA:
    case SA_COMP:
    case SA_E:
    case SA_E_COMP:
    case SA_NEG:
      volumeFields.SetFieldValue("NU_TILDE", Node_Turb->GetSolution(iPoint, 0));
      volumeFields.SetFieldValue("RES_NU_TILDE", turb_solver->LinSysRes(iPoint, 0));
      volumeFields.SetFieldValue("LIM_NU_TILDE", Node_Turb->GetLimiter_Primitive(iPoint, 0));
      break;
    case NONE:
      break;
  }
  if (hybrid_RANSLES != NO_HYBRIDRANSLES) {
    volumeFields.SetFieldValue("LENGTH_SCALE", Node_Flow->GetDES_LengthScale(iPoint));
    volumeFields.SetFieldValue("WALL_DISTANCE", Node_Geo->GetWall_Distance(iPoint));
  }

  if (transModel == BC) {
    volumeFields.SetFieldValue("INTERMITTENCY", Node_Turb->GetGammaBC(iPoint));
  }
}
