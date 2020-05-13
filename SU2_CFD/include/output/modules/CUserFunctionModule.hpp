#pragma once

#include "COutputModule.hpp"

class CUserFunctionModule final : public CSolverOutputModule {

  using FieldRef = typename COutFieldCollection::Map::iterator;
  using FieldRefFuncPair  = std::pair<FieldRef, interpreter::UserFunction*>;
  using FieldRefTokenPair = std::pair<FieldRef, packToken*>;

  TokenMap historyScope;
  TokenMap volumeScope;

  std::vector<interpreter::UserFunction*> historyUserFunctions;
  std::vector<interpreter::UserFunction*> volumeUserFunctions;

  std::vector<FieldRefTokenPair> histFieldTokenRefs;
  std::vector<FieldRefTokenPair> volFieldTokenRefs;

  std::vector<FieldRefFuncPair> funcHistFieldRefs;
  std::vector<FieldRefFuncPair>  funcVolFieldRefs;

public:
  explicit CUserFunctionModule(CConfig *config);

  void DefineHistoryFieldModifier(COutFieldCollection& fieldCollection) override;

  void LoadHistoryDataModifier(COutFieldCollection& fieldCollection) override;

  void DefineHistoryFields(COutFieldCollection& fieldCollection) override;

  void DefineVolumeFields(COutFieldCollection& fieldCollection) override;

  void LoadSurfaceData(COutFieldCollection& fieldCollection) override;

  void LoadVolumeData(COutFieldCollection& fieldCollection) override;

  std::string GetName(const std::string& baseName);

  void EvalUserFunctions(const std::vector<FieldRefTokenPair> &fieldTokenRef,
                         const std::vector<FieldRefFuncPair> &funcFieldRef, TokenMap &scope);
};


