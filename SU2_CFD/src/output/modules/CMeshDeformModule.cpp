#include "../../../include/output/modules/CMeshDeformModule.hpp"
#include "../../../include/solvers/CMeshSolver.hpp"

void CMeshDeformModule::LoadHistoryData(COutFieldCollection& fieldCollection){


  CSolver* mesh_solver = solverData.solver[MESH_SOL];

  fieldCollection.SetValueByKey("DEFORM_MIN_VOLUME", mesh_solver->GetMinimum_Volume());
  fieldCollection.SetValueByKey("DEFORM_MAX_VOLUME", mesh_solver->GetMaximum_Volume());
  fieldCollection.SetValueByKey("DEFORM_ITER", mesh_solver->GetIterLinSolver());
  fieldCollection.SetValueByKey("DEFORM_RESIDUAL", log10(mesh_solver->GetResLinSolver()));
}

void CMeshDeformModule::DefineHistoryFields(COutFieldCollection &fieldCollection){

  fieldCollection.AddItem("DEFORM_MIN_VOLUME", COutputField("MinVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Minimum volume in the mesh", FieldType::DEFAULT));
  fieldCollection.AddItem("DEFORM_MAX_VOLUME", COutputField("MaxVolume", ScreenOutputFormat::SCIENTIFIC, "DEFORM", "Maximum volume in the mesh", FieldType::DEFAULT));
  fieldCollection.AddItem("DEFORM_ITER",       COutputField("DeformIter", ScreenOutputFormat::INTEGER, "DEFORM", "Linear solver iterations for the mesh deformation", FieldType::DEFAULT));
  fieldCollection.AddItem("DEFORM_RESIDUAL",   COutputField("DeformRes", ScreenOutputFormat::FIXED, "DEFORM",  "Residual of the linear solver for the mesh deformation", FieldType::DEFAULT));

}
