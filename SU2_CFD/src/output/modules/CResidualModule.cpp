#include "../../../include/output/modules/CResidualModule.hpp"
#include "../../../../Common/include/CConfig.hpp"

CResidualModule::CResidualModule(CConfig *config){}


void CResidualModule::DefineHistoryFieldModifier(COutFieldCollection &fieldCollection){
  map<string, bool> Average;

  for (const auto& field : fieldCollection.GetFieldsByType({FieldType::RESIDUAL})){
    fieldCollection.AddItem("REL_" + field->first, COutputField("rel" + field->second.fieldName,
                                                                      field->second.screenFormat,
                                                                      "REL_" + field->second.outputGroup,
                                                                      FieldType::AUTO_RESIDUAL,
                                                                      "Relative residual."));
    Average[field->second.outputGroup] = true;
  }

  for (auto it = Average.begin(); it != Average.end(); it++){
    if (AverageGroupName.count(it->first) > 0) {
      fieldCollection.AddItem("AVG_" + it->first, COutputField("avg[" + AverageGroupName.at(it->first) + "]",
                                                                     ScreenOutputFormat::FIXED,
                                                                     "AVG_" + it->first ,
                                                                     FieldType::AUTO_RESIDUAL,
                                                                     "Average residual over all solution variables."));
    }
  }
}

void CResidualModule::LoadHistoryDataModifier(COutFieldCollection &fieldCollection){
  map<string, pair<su2double, int> > Average;

  bool InitResidual = solverData.config->GetTime_Domain() ? solverData.Iter == 0 : solverData.Iter < 2;

  for (const auto& field : fieldCollection.GetFieldsByType({FieldType::RESIDUAL})){
    if (InitResidual || (field->second.value > initialResiduals[field->first])) {
      initialResiduals[field->first] = field->second.value;
    }
    fieldCollection.SetValueByKey("REL_" + field->first, field->second.value - initialResiduals[field->first]);
    Average[field->second.outputGroup].first += field->second.value;
    Average[field->second.outputGroup].second++;
  }

  for (auto it = Average.begin(); it != Average.end(); it++){
    if (AverageGroupName.count(it->first) > 0) {
      const su2double& value = it->second.first;
      const int& count = it->second.second;
      const su2double average = value/count;
      fieldCollection.SetValueByKey("AVG_" + it->first, average);
    }
  }

}

