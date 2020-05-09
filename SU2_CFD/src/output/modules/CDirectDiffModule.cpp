#include "../../../include/output/modules/CDirectDiffModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

constexpr static char prefix[] = "D_";

CDirectDiffModule::CDirectDiffModule(CConfig* config) : CSolverOutputModule(config->GetDirectDiff()){}

void CDirectDiffModule::DefineHistoryFieldModifier(COutFieldCollection &fieldCollection){

  modifiedFields.clear();
  baseFields = fieldCollection.GetFieldsByType({FieldType::COEFFICIENT, FieldType::PER_SURFACE_COEFFICIENT});

  for (const auto& field : baseFields){
    auto newField = fieldCollection.AddItem(prefix + field->first, COutputField("d[" + field->second.fieldName + "]",
                                                                      field->second.screenFormat,
                                                                      prefix + field->second.outputGroup,
                                                                      FieldType::AUTO_COEFFICIENT,
                                                                      "Derivative value (DIRECT_DIFF=YES)"));
    modifiedFields.push_back(newField);
  }

}

void CDirectDiffModule::LoadHistoryDataModifier(COutFieldCollection &fieldCollection){
  for (unsigned int iField = 0; iField < baseFields.size(); iField++){
    fieldCollection.SetValueByKey(modifiedFields[iField]->first, SU2_TYPE::GetDerivative(baseFields[iField]->second.value));
  }
}