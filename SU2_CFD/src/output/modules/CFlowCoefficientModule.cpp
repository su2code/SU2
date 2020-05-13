#include "../../../include/output/modules/CFlowCoefficientModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CFlowCoefficientModule::CFlowCoefficientModule(CConfig *config){}

void CFlowCoefficientModule::DefineVolumeFields(COutFieldCollection &fieldCollection){

  fieldCollection.AddItem("MASSFLOW", COutputField("Mass flow", "FLOW_COEFFICIENT", "Mass flow rate", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("AREA",     COutputField("Area", "FLOW_COEFFICIENT", "Area", FieldType::SURFACE_INTEGRATE));
  fieldCollection.AddItem("NORMAL_VELOCITY", COutputField("Normal velocity", "FLOW_COEFFICIENT", "Velocity component normal to the surface", FieldType::DEFAULT));
  fieldCollection.AddItem("TANGENTIAL_VELOCITY", COutputField("Tangential velocity", "FLOW_COEFFICIENT", "Velocity component tangential to the surface", FieldType::DEFAULT));

}


void CFlowCoefficientModule::LoadSurfaceData(COutFieldCollection &fieldCollection){

  su2double Area = 0.0;
  const su2double* Vector = solverData.vertex->GetNormal();
  su2double MassFlux = 0.0;
  su2double Vn = 0.0, Vt = 0.0;
  const unsigned long iPoint = solverData.iPoint;
  for (int iDim = 0; iDim < solverData.geometry->GetnDim(); iDim++){
    Area += Vector[iDim]*Vector[iDim];
    MassFlux += Vector[iDim]*solverData.solver[FLOW_SOL]->GetNodes()->GetDensity(iPoint)*
        solverData.solver[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim)*solverData.config->GetDensity_Ref()*solverData.config->GetVelocity_Ref();
    Vn += solverData.solver[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim) * Vector[iDim];
  }
  Area = sqrt(Area);
  for (int iDim = 0; iDim < solverData.geometry->GetnDim(); iDim++){
    Vt += solverData.solver[FLOW_SOL]->GetNodes()->GetVelocity(iPoint, iDim) - Vn*Vector[iDim]/Area;
  }

  fieldCollection.SetValueByKey("MASSFLOW", MassFlux);
  fieldCollection.SetValueByKey("AREA", Area);
  fieldCollection.SetValueByKey("NORMAL_VELOCITY", Vn);
  fieldCollection.SetValueByKey("TANGENTIAL_VELOCITY", Vt);

}