#pragma once

#include "COutputModule.hpp"
#include "../../../../Common/include/toolboxes/parser/ExpressionParser.hpp"
#include "../../../../Common/include/toolboxes/parser/CFunctionParser.hpp"
#include <memory>

class CUserFunctionModule final : public CModifierModule {

  using FieldRef = typename COutFieldCollection::Map::iterator;
  std::vector<FieldRef> volumeFunctionFields;
  std::vector<FieldRef> historyFunctionFields;

  std::vector<Parser::CFunctionParser::RawFunction> historyUserFunctions;
  std::vector<Parser::CFunctionParser::RawFunction> volumeUserFunctions;

  std::unique_ptr<Parser::Scope> historyScope;
  std::unique_ptr<Parser::Scope> volumeScope;

  std::vector<std::string> markerNames;

public:
  CUserFunctionModule(CConfig *config, int nDim);

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields,
                               const IterationInfo& iterationInfo) override;

  void DefineVolumeFieldModifier(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceDataModifier(CVolumeOutFieldManager& volumeFields,
                       const IterationInfo& iterationInfo) override;

  void LoadVolumeDataModifier(CVolumeOutFieldManager& volumeFields, const IterationInfo& iterationInfo) override;

  static void EvalUserFunctions(Parser::Scope *scope, std::vector<FieldRef> &functionFields);

};


