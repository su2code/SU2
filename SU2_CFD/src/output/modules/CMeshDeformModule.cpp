#include "../../../include/output/modules/CMeshDeformModule.hpp"
#include "../../../include/solvers/CMeshSolver.hpp"

void CMeshDeformModule::LoadHistoryData(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                        const IterationInfo&){

  const auto* mesh_solver = solverData.solver[MESH_SOL];

  historyFields.SetFieldValue("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
  historyFields.SetFieldValue("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
  historyFields.SetFieldValue("DEFORM_ITER", mesh_solver->GetIterLinSolver());
  historyFields.SetFieldValue("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
}

void CMeshDeformModule::DefineHistoryFields(CHistoryOutFieldManager &historyFields){

  historyFields.AddField("DEFORM_MIN_VOLUME", "MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh", FieldType::DEFAULT);
  historyFields.AddField("DEFORM_MAX_VOLUME", "MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh", FieldType::DEFAULT);
  historyFields.AddField("DEFORM_ITER",       "DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation", FieldType::DEFAULT);
  historyFields.AddField("DEFORM_RESIDUAL",   "DeformRes", ScreenOutputFormat::FIXED, "DEFORM",  "Residual of the linear solver for the mesh deformation", FieldType::DEFAULT);

}
