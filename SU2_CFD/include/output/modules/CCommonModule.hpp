#ifndef CCOMMONMODULE_HPP
#define CCOMMONMODULE_HPP

#include "COutputModule.hpp"

class CCommonModule : public CSolverOutputModule{

public:
  CCommonModule(CConfig* config) {};

  void DefineHistoryFields(COutFieldCollection& fieldCollection) override;

  void LoadHistoryData(COutFieldCollection& fieldCollection) override;
};

#endif // CCOMMONMODULE_HPP
