#include "../../../include/output/modules/CResidualModule.hpp"

#include "../../../../Common/include/CConfig.hpp"

CResidualModule::CResidualModule(CConfig* config, int nDim) : CSolverOutputModule(nDim) {}

void CResidualModule::DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) {
  map<string, bool> Average;

  for (const auto& field : historyFields.GetCollection().GetFieldsByType({FieldType::RESIDUAL})) {
    historyFields.AddField("REL_" + field->first, "rel" + field->second.fieldName, field->second.screenFormat,
                           "REL_" + field->second.outputGroup, "Relative residual.", FieldType::AUTO_RESIDUAL);
    Average[field->second.outputGroup] = true;
  }

  for (auto it = Average.begin(); it != Average.end(); it++) {
    if (AverageGroupName.count(it->first) > 0) {
      historyFields.AddField("AVG_" + it->first, "avg[" + AverageGroupName.at(it->first) + "]",
                             ScreenOutputFormat::FIXED, "AVG_" + it->first,
                             "Average residual over all solution variables.", FieldType::AUTO_RESIDUAL);
    }
  }

  for (const auto& field : historyFields.GetCollection().GetFieldsByType({FieldType::RESIDUAL})) {
    initialResiduals[field->first] = numeric_limits<su2double>::lowest();
  }
}

void CResidualModule::LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                                              const IterationInfo& iterationInfo) {

  const auto* config = solverData.config;
  const auto Iter    = iterationInfo.Iter;

  map<string, pair<su2double, int> > Average;

  bool InitResidual = config->GetTime_Domain() ? Iter == 0 : Iter < 2;

  for (const auto& field : historyFields.GetCollection().GetFieldsByType({FieldType::RESIDUAL})) {
    if (InitResidual || (field->second.value > initialResiduals[field->first])) {
      initialResiduals[field->first] = field->second.value;
    }
    historyFields.SetFieldValue("REL_" + field->first, field->second.value - initialResiduals[field->first]);
    Average[field->second.outputGroup].first += field->second.value;
    Average[field->second.outputGroup].second++;
  }

  for (auto it = Average.begin(); it != Average.end(); it++) {
    if (AverageGroupName.count(it->first) > 0) {
      const su2double& value = it->second.first;
      const int& count = it->second.second;
      const su2double average = value / count;
      historyFields.SetFieldValue("AVG_" + it->first, average);
    }
  }
}
