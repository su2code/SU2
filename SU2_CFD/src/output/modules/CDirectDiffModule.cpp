#include "../../../include/output/modules/CDirectDiffModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

constexpr static char prefix[] = "D_";

CDirectDiffModule::CDirectDiffModule(CConfig* config, int nDim) : CSolverOutputModule(nDim, config->GetDirectDiff()){}

void CDirectDiffModule::DefineHistoryFieldModifier(CHistoryOutFieldManager &historyFields){

  modifiedFields.clear();
  baseFields = historyFields.GetCollection().GetFieldsByType({FieldType::COEFFICIENT, FieldType::PER_SURFACE_COEFFICIENT});

  for (const auto& field : baseFields){
    auto newField = historyFields.AddField(prefix + field->first, "d[" + field->second.fieldName + "]",
                                           field->second.screenFormat,
                                           prefix + field->second.outputGroup,
                                           "Derivative value (DIRECT_DIFF=YES)",
                                           FieldType::AUTO_COEFFICIENT );
    modifiedFields.push_back(newField);
  }

}

void CDirectDiffModule::LoadHistoryDataModifier(CHistoryOutFieldManager &historyFields){
  for (unsigned int iField = 0; iField < baseFields.size(); iField++){
    historyFields.SetFieldValue(modifiedFields[iField]->first, SU2_TYPE::GetDerivative(baseFields[iField]->second.value));
  }
}