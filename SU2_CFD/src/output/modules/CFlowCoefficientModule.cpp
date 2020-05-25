#include "../../../include/output/modules/CFlowCoefficientModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CFlowCoefficientModule::CFlowCoefficientModule(CConfig *config, int nDim) : CSolverOutputModule(nDim){}

void CFlowCoefficientModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){



  volumeFields.AddField("MASSFLOW", "Mass flow", "FLOW_COEFFICIENT", "Mass flow rate", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("NORMAL_VELOCITY", "Normal velocity", "FLOW_COEFFICIENT", "Velocity component normal to the surface", FieldType::DEFAULT);
  volumeFields.AddField("TANG_VELOCITY", "Tangential velocity", "FLOW_COEFFICIENT", "Velocity component tangential to the surface", FieldType::DEFAULT);


}

void CFlowCoefficientModule::LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                            const IterationInfo& iterationInfo, const PointInfo& pointInfo){


}

void CFlowCoefficientModule::LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                                             const IterationInfo& iterationInfo, const PointInfo& pointInfo){
  const auto* config   = std::get<0>(solverData);
  const auto* geometry = std::get<1>(solverData);
  const auto* solver   = std::get<2>(solverData);
  const auto iPoint    = std::get<0>(pointInfo);
  const auto iVertex   = std::get<1>(pointInfo);
  const auto iMarker   = std::get<2>(pointInfo);
  const auto *Node_Flow = solver[FLOW_SOL]->GetNodes();

  const su2double Area = volumeFields.GetFieldValue("AREA");
  const su2double* Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
  su2double MassFlux = 0.0;
  su2double Vn = 0.0, Vt = 0.0;
  for (int iDim = 0; iDim < nDim; iDim++){
    MassFlux += Normal[iDim]*Node_Flow->GetDensity(iPoint)*
        Node_Flow->GetVelocity(iPoint, iDim)*config->GetDensity_Ref()*config->GetVelocity_Ref();
    Vn += Node_Flow->GetVelocity(iPoint, iDim) * Normal[iDim];
  }
  for (int iDim = 0; iDim < nDim; iDim++){
    Vt += Node_Flow->GetVelocity(iPoint, iDim) - Vn*Normal[iDim]/Area;
  }

  volumeFields.SetFieldValue("MASSFLOW", MassFlux);
  volumeFields.SetFieldValue("NORMAL_VELOCITY", Vn);
  volumeFields.SetFieldValue("TANG_VELOCITY", Vt);

}