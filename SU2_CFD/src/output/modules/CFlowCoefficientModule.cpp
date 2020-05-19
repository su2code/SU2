#include "../../../include/output/modules/CFlowCoefficientModule.hpp"
#include "../../../../Common/include/CConfig.hpp"
#include "../../../../Common/include/geometry/CPhysicalGeometry.hpp"
#include "../../../include/solvers/CSolver.hpp"

CFlowCoefficientModule::CFlowCoefficientModule(CConfig *config, int nDim) : CSolverOutputModule(nDim){}

void CFlowCoefficientModule::DefineVolumeFields(CVolumeOutFieldManager &volumeFields){

  volumeFields.AddField("PRESSURE",    "Pressure",                "PRIMITIVE", "Pressure", FieldType::DEFAULT);
  volumeFields.AddField("TEMPERATURE", "Temperature",             "PRIMITIVE", "Temperature", FieldType::DEFAULT);
  volumeFields.AddField("MACH",        "Mach",                    "PRIMITIVE", "Mach number", FieldType::DEFAULT);
  volumeFields.AddField("PRESSURE_COEFF", "Pressure_Coefficient", "PRIMITIVE", "Pressure coefficient", FieldType::DEFAULT);



  volumeFields.AddField("MASSFLOW", "Mass flow", "FLOW_COEFFICIENT", "Mass flow rate", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("AREA",     "Area", "FLOW_COEFFICIENT", "Area", FieldType::SURFACE_INTEGRATE);
  volumeFields.AddField("NORMAL_VELOCITY", "Normal velocity", "FLOW_COEFFICIENT", "Velocity component normal to the surface", FieldType::DEFAULT);
  volumeFields.AddField("TANGENTIAL_VELOCITY", "Tangential velocity", "FLOW_COEFFICIENT", "Velocity component tangential to the surface", FieldType::DEFAULT);


}


void CFlowCoefficientModule::LoadSurfaceData(CVolumeOutFieldManager &volumeFields){

  const auto iPoint = solverData.iPoint;
  const auto *Node_Flow = solverData.solver[FLOW_SOL]->GetNodes();
  const auto *config = solverData.config;

  su2double Area = 0.0;
  const su2double* Vector = solverData.vertex->GetNormal();
  su2double MassFlux = 0.0;
  su2double Vn = 0.0, Vt = 0.0;
  for (int iDim = 0; iDim < solverData.geometry->GetnDim(); iDim++){
    Area += Vector[iDim]*Vector[iDim];
    MassFlux += Vector[iDim]*Node_Flow->GetDensity(iPoint)*
        Node_Flow->GetVelocity(iPoint, iDim)*solverData.config->GetDensity_Ref()*solverData.config->GetVelocity_Ref();
    Vn += Node_Flow->GetVelocity(iPoint, iDim) * Vector[iDim];
  }
  Area = sqrt(Area);
  for (int iDim = 0; iDim < solverData.geometry->GetnDim(); iDim++){
    Vt += Node_Flow->GetVelocity(iPoint, iDim) - Vn*Vector[iDim]/Area;
  }

  volumeFields.SetFieldValue("PRESSURE", Node_Flow->GetPressure(iPoint)*config->GetPressure_Ref());
  volumeFields.SetFieldValue("TEMPERATURE", Node_Flow->GetTemperature(iPoint)*config->GetTemperature_Ref());
  const su2double Mach = sqrt(Node_Flow->GetVelocity2(iPoint))/Node_Flow->GetSoundSpeed(iPoint);
  volumeFields.SetFieldValue("MACH", Mach);
  su2double VelMag = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++){
    VelMag += pow(solverData.solver[FLOW_SOL]->GetVelocity_Inf(iDim),2.0);
  }
  const su2double factor = 1.0/(0.5*solverData.solver[FLOW_SOL]->GetDensity_Inf()*VelMag);
  volumeFields.SetFieldValue("PRESSURE_COEFF",
                             (Node_Flow->GetPressure(iPoint) - solverData.solver[FLOW_SOL]->GetPressure_Inf())*factor);

  volumeFields.SetFieldValue("MASSFLOW", MassFlux);
  volumeFields.SetFieldValue("AREA", Area);
  volumeFields.SetFieldValue("NORMAL_VELOCITY", Vn);
  volumeFields.SetFieldValue("TANGENTIAL_VELOCITY", Vt);

}