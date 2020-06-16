#include "../../../include/output/modules/CFixedCLModule.hpp"
#include "../../../include/solvers/CSolver.hpp"

void CFixedCLModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){
  /// DESCRIPTION: Difference between current and target CL
  historyFields.AddField("DELTA_CL","Delta_CL", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL",  "Difference between Target CL and current CL", FieldType::DEFAULT);
  /// DESCRIPTION: Angle of attack before the most recent update
  historyFields.AddField("PREV_AOA","Previous_AOA", ScreenOutputFormat::FIXED, "FIXED_CL", "Angle of Attack at the previous iteration of the Fixed CL driver", FieldType::DEFAULT);
  /// DESCRIPTION: Last change in angle of attack by the Fixed CL driver
  historyFields.AddField("CHANGE_IN_AOA", "Change_in_AOA", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "Last change in Angle of Attack by Fixed CL Driver", FieldType::DEFAULT);
  /// DESCRIPTION: AOA control command by the CL Driver
  historyFields.AddField("CL_DRIVER_COMMAND", "CL_Driver_Command", ScreenOutputFormat::SCIENTIFIC, "FIXED_CL", "CL Driver's control command", FieldType::DEFAULT);

  historyFields.AddField("AOA", "AoA", ScreenOutputFormat::FIXED, "FIXED_CL", "Angle of attack", FieldType::DEFAULT);
}

void CFixedCLModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                     const IterationInfo&){

  const auto* config      = solverData.config;
  const auto* flow_solver = solverData.solver[FLOW_SOL];

  historyFields.SetFieldValue("DELTA_CL", fabs(flow_solver->GetTotal_CL() - config->GetTarget_CL()));
  historyFields.SetFieldValue("PREV_AOA", flow_solver->GetPrevious_AoA());
  historyFields.SetFieldValue("CHANGE_IN_AOA", config->GetAoA()-flow_solver->GetPrevious_AoA());
  historyFields.SetFieldValue("CL_DRIVER_COMMAND", flow_solver->GetAoA_inc());
  historyFields.SetFieldValue("AOA", config->GetAoA());

}

