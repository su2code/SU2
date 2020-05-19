#ifndef CCOMMONMODULE_HPP
#define CCOMMONMODULE_HPP

#include "COutputModule.hpp"

class CCommonModule : public CSolverOutputModule{

public:
  CCommonModule(CConfig* config) {};

  void DefineHistoryFields(CHistoryOutFieldManager& historyFields) override;

  void LoadHistoryData(CHistoryOutFieldManager& historyFields) override;
};

#endif // CCOMMONMODULE_HPP
