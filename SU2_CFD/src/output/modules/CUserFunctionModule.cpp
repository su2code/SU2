#include "../../../include/output/modules/CUserFunctionModule.hpp"

#include "../../../../Common/include/CConfig.hpp"

CUserFunctionModule::CUserFunctionModule(CConfig* config, int nDim) : CSolverOutputModule(nDim) {

  Parser::CFunctionParser funcParser;
  funcParser.ParseFunctionsFromFile("functions.su2x");

  historyUserFunctions = funcParser.GetExpressions({"historyfield"});
  volumeUserFunctions = funcParser.GetExpressions({"volumefield", "surfaceintegral"});

  auto addInlineUserFunction = [&](std::string field, std::vector<Parser::CFunctionParser::RawExpression>& userFuncsVector) {
    if (*field.begin() == '{' && *(field.end() - 1) == '}') {
      auto rawExpr = field;
      rawExpr.pop_back();
      rawExpr.erase(0, 1);
      userFuncsVector.push_back({"inline", field, rawExpr, {}});
    }
  };

  for (int iField = 0; iField < config->GetnHistoryOutput(); iField++) {
    addInlineUserFunction(config->GetHistoryOutput_Field(iField), historyUserFunctions);
  }
  for (int iField = 0; iField < config->GetnScreenOutput(); iField++) {
    addInlineUserFunction(config->GetScreenOutput_Field(iField), historyUserFunctions);
  }
  for (int iField = 0; iField < config->GetnVolumeOutput(); iField++) {
    addInlineUserFunction(config->GetVolumeOutput_Field(iField), volumeUserFunctions);
  }

  historyScope = std::unique_ptr<Parser::Scope>(new Parser::Scope());
  volumeScope = std::unique_ptr<Parser::Scope>(new Parser::Scope());

  for (int iMarker_CfgFile = 0; iMarker_CfgFile <  config->GetnMarker_CfgFile(); iMarker_CfgFile++){
    markerNames.push_back(config->GetMarker_CfgFile_TagBound(iMarker_CfgFile));
  }
}

void CUserFunctionModule::DefineVolumeFields(CVolumeOutFieldManager& volumeFields) {
  for (const auto& field : volumeFields.GetCollection().GetReferencesAll()) {
    volumeScope->addVariable(field->first, field->second.value);
  }

  for (uint iFunc = 0; iFunc < volumeUserFunctions.size(); iFunc++) {
    const std::string funcName = volumeUserFunctions[iFunc].name;
    FieldType fieldType;
    if (volumeUserFunctions[iFunc].type == "surfaceintegral") {
      fieldType = FieldType::SURFACE_INTEGRATE;
    } else {
      fieldType = FieldType::DEFAULT;
    }

    auto newField = volumeFields.AddField(funcName, funcName, "CUSTOM", "User-defined field", fieldType);

    volumeFunctionFields.push_back(newField);
    if (!volumeScope->addExpression(volumeUserFunctions[iFunc].expr)) {
      SU2_MPI::Error("In output expression " + volumeUserFunctions[iFunc].name + ": " + volumeScope->GetError(),
                     CURRENT_FUNCTION);
    }
  }
}

void CUserFunctionModule::LoadSurfaceData(CVolumeOutFieldManager&, const SolverData&, const IterationInfo&,
                                          const PointInfo&) {
  EvalUserFunctions(volumeScope.get(), volumeFunctionFields);
}

void CUserFunctionModule::LoadVolumeData(CVolumeOutFieldManager&, const SolverData&, const IterationInfo&,
                                         const PointInfo&) {
  EvalUserFunctions(volumeScope.get(), volumeFunctionFields);
}

void CUserFunctionModule::DefineHistoryFields(CHistoryOutFieldManager&) {}

void CUserFunctionModule::DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) {


  for (const auto& field : historyFields.GetCollection().GetReferencesAll()) {
    historyScope->addVariable(field->first, field->second.value);

    if (field->second.fieldType == FieldType::COEFFICIENT){
      historyScope->addCustomFunction(field->first + "@", Parser::CustomFunction::SURFACE_SUM);
    }
  }
  for (auto& markerName : markerNames){
    historyScope->addStringVar(markerName, markerName);
  }
  for (uint iFunc = 0; iFunc < historyUserFunctions.size(); iFunc++) {
    const auto funcName = historyUserFunctions[iFunc].name;
    auto newField = historyFields.AddField(funcName, funcName, ScreenOutputFormat::SCIENTIFIC, "CUSTOM",
                                           "User-defined field", FieldType::DEFAULT);

    historyFunctionFields.push_back(newField);
    if (!historyScope->addExpression(historyUserFunctions[iFunc].expr)) {
      SU2_MPI::Error("In output expression: " + historyUserFunctions[iFunc].name + ": " + historyScope->GetError(),
                     CURRENT_FUNCTION);
    }
  }
}

void CUserFunctionModule::LoadHistoryDataModifier(CHistoryOutFieldManager&, const SolverData&, const IterationInfo&) {
  EvalUserFunctions(historyScope.get(), historyFunctionFields);
}

void CUserFunctionModule::EvalUserFunctions(Parser::Scope* scope, std::vector<FieldRef>& functionFields) {
  auto values = scope->EvalExpressions();
  for (uint iField = 0; iField < functionFields.size(); iField++) {
    functionFields[iField]->second.value = values[iField];
  }
}
