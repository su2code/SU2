#pragma once

#include "COutputModule.hpp"
#include "../../../../Common/include/toolboxes/parser/ExpressionParser.hpp"
#include "../../../../Common/include/toolboxes/parser/CFunctionParser.hpp"
#include <memory>

class CUserFunctionModule final : public CSolverOutputModule {

  using FieldRef = typename COutFieldCollection::Map::iterator;
  std::vector<FieldRef> volumeFunctionFields;
  std::vector<FieldRef> historyFunctionFields;

//  using FieldRefFuncPair  = std::pair<FieldRef, Parser::Expression>;

  std::vector<Parser::CFunctionParser::RawFunction> historyUserFunctions;
  std::vector<Parser::CFunctionParser::RawFunction> volumeUserFunctions;

  std::unique_ptr<Parser::Scope> historyScope;
  std::unique_ptr<Parser::Scope> volumeScope;

  std::vector<std::string> markerNames;

public:
  explicit CUserFunctionModule(CConfig *config, int nDim);

  void DefineHistoryFieldModifier(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryDataModifier(CHistoryOutFieldManager& historyFields, const SolverData& solverData,
                               const IterationInfo& iterationInfo) override;

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void DefineVolumeFields(CVolumeOutFieldManager& volumeFields) override;

  void LoadSurfaceData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                       const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  void LoadVolumeData(CVolumeOutFieldManager& volumeFields, const SolverData& solverData,
                      const IterationInfo& iterationInfo, const PointInfo& pointInfo) override;

  static void EvalUserFunctions(Parser::Scope *scope, std::vector<FieldRef> &functionFields);

};


